//! Distributed generalized Poisson / linearized Poisson-Boltzmann (LPB) solver
//! built directly on BAANDI's `DistributedPairFft` with zero inner-loop heap
//! allocations and deterministic `Collective::tree_sum_f64` reductions.

use super::{
    IonicPseudoCharge, SolvationError, SolvationEvaluation, SolvationParams, SolveReport,
    SolverConfig, SolverKind, cavity,
};
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_numerics::real::{RealEighWorkspace, RealMatrix};
use dft_solvers::distributed_fft::DistributedPairFft;
use num_complex::Complex64;
use std::f64::consts::PI;

const FOUR_PI: f64 = 4.0 * PI;
const INV_EIGHT_PI: f64 = 1.0 / (8.0 * PI);
const DEFAULT_DIIS_CAP: usize = 8;

/// Core distributed FFT engine, reciprocal-space metadata, and plane-streaming
/// scratch buffers shared by `apply_operator`, `solve_pde`, and `evaluate`.
struct FftEngine {
    cell: Cell,
    fft: DistributedPairFft,
    dims: [usize; 3],
    total_points: usize,
    volume: f64,
    dv: f64,
    g_cart: Vec<[f64; 3]>,
    g2: Vec<f64>,
    g2_min: f64,
    zero_g2_indices: Vec<usize>,
    mirror_idx: Option<Vec<usize>>,
    batch_g: Vec<Complex64>,
    batch_out_g: Vec<Complex64>,
    slab_r0: Vec<f64>,
    slab_r1: Vec<f64>,
    slab_r2: Vec<f64>,
    slab_r3: Vec<f64>,
    rank_pair_counts: Vec<usize>,
    rank_pair_buf: Vec<f64>,
}

impl FftEngine {
    fn new<K: Collective>(
        cell: &Cell,
        mut fft: DistributedPairFft,
        comm: &K,
    ) -> Result<Self, SolvationError> {
        let dims = fft.shape();
        if dims.contains(&0) {
            return Err(SolvationError::InvalidParameter(
                "solvation grid dimensions must be positive",
            ));
        }
        let total_points = dims[0]
            .checked_mul(dims[1])
            .and_then(|n| n.checked_mul(dims[2]))
            .ok_or(SolvationError::Shape)?;
        let volume = cell.volume_bohr3();
        if !volume.is_finite() || volume <= 0.0 {
            return Err(SolvationError::InvalidParameter(
                "cell volume must be finite and positive",
            ));
        }

        let ll = fft.local_len();
        let pp = fft.plane_points();
        let out_of_range = fft.local_positions().iter().any(|&p| p >= total_points);
        let mut global_counts = [0.0_f64; 3];
        comm.tree_sum_f64(
            &[ll as f64, pp as f64, f64::from(out_of_range)],
            &mut global_counts,
        )?;
        if (global_counts[0] - total_points as f64).abs() > 0.5
            || (global_counts[1] - total_points as f64).abs() > 0.5
            || global_counts[2] != 0.0
        {
            return Err(SolvationError::Shape);
        }

        // Pre-warm FFT batch buffers once for np = 4 so greedy fallback (np = 3, 4)
        // and evaluate's np = 3 inverse batch never resize or mismatch send lengths.
        let batch_g = vec![Complex64::default(); 4 * ll];
        let batch_out_g = vec![Complex64::default(); 4 * ll];
        fft.begin_to_real_batch(comm, &batch_g, 4)?;

        let recip_cols = cell.reciprocal_columns_inverse_bohr();
        let local_pos = fft.local_positions();
        let mut g_cart = Vec::with_capacity(ll);
        let mut g2 = Vec::with_capacity(ll);
        let mut zero_g2_indices = Vec::new();
        let mut local_g2_min = f64::INFINITY;
        for (k, &g_idx) in local_pos.iter().enumerate() {
            let gc = cavity::skew_adjoint_g_cart(dims, g_idx, recip_cols);
            let norm2 = gc[0].mul_add(gc[0], gc[1].mul_add(gc[1], gc[2] * gc[2]));
            g_cart.push(gc);
            g2.push(norm2);
            if norm2 == 0.0 {
                zero_g2_indices.push(k);
            } else if norm2 < local_g2_min {
                local_g2_min = norm2;
            }
        }

        let ranks = comm.size();
        let rank_pair_counts = vec![2usize; ranks];
        let mut rank_pair_buf = vec![0.0_f64; 2 * ranks];
        comm.allgatherv_f64(&[local_g2_min, 0.0], &rank_pair_counts, &mut rank_pair_buf)?;
        let mut g2_min = f64::INFINITY;
        for r in 0..ranks {
            let v = rank_pair_buf[2 * r];
            if v.is_finite() && v > 0.0 && v < g2_min {
                g2_min = v;
            }
        }
        if !g2_min.is_finite() {
            g2_min = 1.0;
        }

        // Check whether every owned G-vector's mirror (-G mod N) is rank-local on ALL ranks.
        let [n0, n1, n2] = dims;
        let n12 = n1 * n2;
        let mut local_mirror = Vec::with_capacity(ll);
        let mut all_local = true;
        for &g_idx in local_pos {
            let ix = g_idx / n12;
            let iy = (g_idx / n2) % n1;
            let iz = g_idx % n2;
            let mx = (n0 - ix) % n0;
            let my = (n1 - iy) % n1;
            let mz = (n2 - iz) % n2;
            let mg = (mx * n1 + my) * n2 + mz;
            if let Ok(mk) = local_pos.binary_search(&mg) {
                local_mirror.push(mk);
            } else {
                all_local = false;
                break;
            }
        }
        let mut missing_mirrors = [0.0_f64];
        comm.tree_sum_f64(&[f64::from(!all_local)], &mut missing_mirrors)?;
        let mirror_idx = if missing_mirrors[0] == 0.0 {
            Some(local_mirror)
        } else {
            None
        };

        let slab_len = fft.rows() * dims[2];
        Ok(Self {
            cell: cell.clone(),
            fft,
            dims,
            total_points,
            volume,
            dv: volume / total_points as f64,
            g_cart,
            g2,
            g2_min,
            zero_g2_indices,
            mirror_idx,
            batch_g,
            batch_out_g,
            slab_r0: vec![0.0; slab_len],
            slab_r1: vec![0.0; slab_len],
            slab_r2: vec![0.0; slab_len],
            slab_r3: vec![0.0; slab_len],
            rank_pair_counts,
            rank_pair_buf,
        })
    }

    fn local_len(&self) -> usize {
        self.fft.local_len()
    }

    fn plane_points(&self) -> usize {
        self.fft.plane_points()
    }

    fn mirror_paired(&self) -> bool {
        self.mirror_idx.is_some()
    }

    fn global_min_max<K: Collective>(
        &mut self,
        comm: &K,
        values: &[f64],
    ) -> Result<(f64, f64), SolvationError> {
        let mut local_min = f64::INFINITY;
        let mut local_max = f64::NEG_INFINITY;
        for &v in values {
            if v < local_min {
                local_min = v;
            }
            if v > local_max {
                local_max = v;
            }
        }
        comm.allgatherv_f64(
            &[local_min, local_max],
            &self.rank_pair_counts,
            &mut self.rank_pair_buf,
        )?;
        let mut global_min = f64::INFINITY;
        let mut global_max = f64::NEG_INFINITY;
        for r in 0..self.rank_pair_counts.len() {
            let r_min = self.rank_pair_buf[2 * r];
            let r_max = self.rank_pair_buf[2 * r + 1];
            if r_min < global_min {
                global_min = r_min;
            }
            if r_max > global_max {
                global_max = r_max;
            }
        }
        Ok((global_min, global_max))
    }

    fn forward_real<K: Collective>(
        &mut self,
        comm: &K,
        real_planes: &[f64],
        out_g: &mut [Complex64],
    ) -> Result<(), SolvationError> {
        if real_planes.len() != self.plane_points() || out_g.len() != self.local_len() {
            return Err(SolvationError::Shape);
        }
        let slab_len = self.fft.rows() * self.dims[2];
        for i in 0..self.fft.plane_count() {
            let src = &real_planes[i * slab_len..(i + 1) * slab_len];
            let slab = self.fft.plane_mut();
            for (dst, &v) in slab.iter_mut().zip(src) {
                *dst = Complex64::new(v, 0.0);
            }
            self.fft.to_reciprocal_plane(i)?;
        }
        self.fft.finish_to_reciprocal(comm, out_g)?;
        Ok(())
    }

    fn inverse_to_real<K: Collective>(
        &mut self,
        comm: &K,
        coeff_g: &[Complex64],
        out_planes: &mut [f64],
    ) -> Result<(), SolvationError> {
        if coeff_g.len() != self.local_len() || out_planes.len() != self.plane_points() {
            return Err(SolvationError::Shape);
        }
        self.fft.begin_to_real(comm, coeff_g)?;
        let slab_len = self.fft.rows() * self.dims[2];
        for i in 0..self.fft.plane_count() {
            self.fft.to_real_plane(i)?;
            let slab = self.fft.plane();
            let dst = &mut out_planes[i * slab_len..(i + 1) * slab_len];
            for (d, s) in dst.iter_mut().zip(slab) {
                *d = s.re;
            }
        }
        Ok(())
    }

    /// Applies the exact discrete SPD operator `A x = sum_a D_a^T (eps D_a x) + kappa2 * x`
    /// and returns the number of directional 3D complex FFTs executed.
    fn apply_operator<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        has_kappa: bool,
        x_g: &[Complex64],
        out_g: &mut [Complex64],
    ) -> Result<usize, SolvationError> {
        let ll = self.local_len();
        let pp = self.plane_points();
        if eps_planes.len() != pp
            || kappa2_planes.len() != pp
            || x_g.len() != ll
            || out_g.len() != ll
        {
            return Err(SolvationError::Shape);
        }

        // G -> r leg (np = 2, valid on every stick partition):
        //   u0 = (i gx - gy) * x_g  =>  Re(u0) = D_x x, Im(u0) = D_y x
        //   u1 = i (gz + 1) * x_g   =>  Re(u1) = D_z x, Im(u1) = x
        let (b0, b1) = self.batch_g[..2 * ll].split_at_mut(ll);
        for (((d0, d1), &[gx, gy, gz]), &x) in
            b0.iter_mut().zip(b1.iter_mut()).zip(&self.g_cart).zip(x_g)
        {
            *d0 = Complex64::new(-gy, gx) * x;
            *d1 = Complex64::new(0.0, gz + 1.0) * x;
        }
        self.fft
            .begin_to_real_batch(comm, &self.batch_g[..2 * ll], 2)?;

        let slab_len = self.fft.rows() * self.dims[2];
        let plane_count = self.fft.plane_count();

        if let Some(mirror) = self.mirror_idx.as_deref() {
            // Mirror-stick r -> G leg (np = 2):
            //   v0 = eps * D_x x + i * eps * D_y x
            //   v1 = eps * D_z x + i * kappa2 * x
            for i in 0..plane_count {
                let off = i * slab_len;
                let eps = &eps_planes[off..off + slab_len];
                let kap = &kappa2_planes[off..off + slab_len];

                self.fft.to_real_plane_in(0, 2, i)?;
                {
                    let slab = self.fft.plane_mut();
                    for (s, &e) in slab.iter_mut().zip(eps) {
                        *s *= e;
                    }
                }
                self.fft.to_reciprocal_plane_in(0, 2, i)?;

                self.fft.to_real_plane_in(1, 2, i)?;
                {
                    let slab = self.fft.plane_mut();
                    for ((s, &e), &kp) in slab.iter_mut().zip(eps).zip(kap) {
                        *s = Complex64::new(e * s.re, kp * s.im);
                    }
                }
                self.fft.to_reciprocal_plane_in(1, 2, i)?;
            }
            self.fft
                .finish_to_reciprocal_batch(comm, &mut self.batch_out_g[..2 * ll], 2)?;

            let u0 = &self.batch_out_g[..ll];
            let u1 = &self.batch_out_g[ll..2 * ll];
            let neg_half_i = Complex64::new(0.0, -0.5);
            for (k, out) in out_g.iter_mut().enumerate() {
                let mk = mirror[k];
                let u0_k = u0[k];
                let u0_m = u0[mk].conj();
                let u1_k = u1[k];
                let u1_m = u1[mk].conj();
                let ax = 0.5 * (u0_k + u0_m);
                let ay = neg_half_i * (u0_k - u0_m);
                let az = 0.5 * (u1_k + u1_m);
                let kp = neg_half_i * (u1_k - u1_m);
                let [gx, gy, gz] = self.g_cart[k];
                *out = Complex64::new(0.0, -gx) * ax
                    + Complex64::new(0.0, -gy) * ay
                    + Complex64::new(0.0, -gz) * az
                    + kp;
            }
            Ok(4)
        } else {
            // Greedy fallback r -> G leg: 3 real fields (kappa = 0) or 4 real fields (kappa > 0).
            let np_out = if has_kappa { 4 } else { 3 };
            for i in 0..plane_count {
                let off = i * slab_len;
                let eps = &eps_planes[off..off + slab_len];
                let kap = &kappa2_planes[off..off + slab_len];

                self.fft.to_real_plane_in(0, 2, i)?;
                {
                    let slab = self.fft.plane();
                    for (((r0, r1), s), &e) in self
                        .slab_r0
                        .iter_mut()
                        .zip(&mut self.slab_r1)
                        .zip(slab)
                        .zip(eps)
                    {
                        *r0 = e * s.re;
                        *r1 = e * s.im;
                    }
                }
                self.fft.to_real_plane_in(1, 2, i)?;
                {
                    let slab = self.fft.plane();
                    for ((((r2, r3), s), &e), &kp) in self
                        .slab_r2
                        .iter_mut()
                        .zip(&mut self.slab_r3)
                        .zip(slab)
                        .zip(eps)
                        .zip(kap)
                    {
                        *r2 = e * s.re;
                        *r3 = kp * s.im;
                    }
                }

                for (dst, &v) in self.fft.plane_mut().iter_mut().zip(&self.slab_r0) {
                    *dst = Complex64::new(v, 0.0);
                }
                self.fft.to_reciprocal_plane_in(0, np_out, i)?;

                for (dst, &v) in self.fft.plane_mut().iter_mut().zip(&self.slab_r1) {
                    *dst = Complex64::new(v, 0.0);
                }
                self.fft.to_reciprocal_plane_in(1, np_out, i)?;

                for (dst, &v) in self.fft.plane_mut().iter_mut().zip(&self.slab_r2) {
                    *dst = Complex64::new(v, 0.0);
                }
                self.fft.to_reciprocal_plane_in(2, np_out, i)?;

                if has_kappa {
                    for (dst, &v) in self.fft.plane_mut().iter_mut().zip(&self.slab_r3) {
                        *dst = Complex64::new(v, 0.0);
                    }
                    self.fft.to_reciprocal_plane_in(3, np_out, i)?;
                }
            }
            self.fft.finish_to_reciprocal_batch(
                comm,
                &mut self.batch_out_g[..np_out * ll],
                np_out,
            )?;

            let ax = &self.batch_out_g[..ll];
            let ay = &self.batch_out_g[ll..2 * ll];
            let az = &self.batch_out_g[2 * ll..3 * ll];
            if has_kappa {
                let kp = &self.batch_out_g[3 * ll..4 * ll];
                for (k, out) in out_g.iter_mut().enumerate() {
                    let [gx, gy, gz] = self.g_cart[k];
                    *out = Complex64::new(0.0, -gx) * ax[k]
                        + Complex64::new(0.0, -gy) * ay[k]
                        + Complex64::new(0.0, -gz) * az[k]
                        + kp[k];
                }
            } else {
                for (k, out) in out_g.iter_mut().enumerate() {
                    let [gx, gy, gz] = self.g_cart[k];
                    *out = Complex64::new(0.0, -gx) * ax[k]
                        + Complex64::new(0.0, -gy) * ay[k]
                        + Complex64::new(0.0, -gz) * az[k];
                }
            }
            Ok(2 + np_out)
        }
    }

    /// Evaluates `phi_solv`, `|D phi_solv|^2`, and `phi_pol` on owned planes in one
    /// batched `np = 3` inverse FFT (1 `alltoall_f64`).
    fn reconstruct_potentials_and_grad_sq<K: Collective>(
        &mut self,
        comm: &K,
        phi_solv_g: &[Complex64],
        phi_pol_g: &[Complex64],
        phi_solv_planes: &mut [f64],
        grad_phi_solv_sq_planes: &mut [f64],
        phi_pol_planes: &mut [f64],
    ) -> Result<(), SolvationError> {
        let ll = self.local_len();
        let (b0, rest) = self.batch_g[..3 * ll].split_at_mut(ll);
        let (b1, b2) = rest.split_at_mut(ll);
        for ((((d0, d1), d2), &[gx, gy, gz]), (&ps, &pp)) in b0
            .iter_mut()
            .zip(b1.iter_mut())
            .zip(b2.iter_mut())
            .zip(&self.g_cart)
            .zip(phi_solv_g.iter().zip(phi_pol_g))
        {
            *d0 = Complex64::new(-gy, gx) * ps;
            *d1 = Complex64::new(0.0, gz + 1.0) * ps;
            *d2 = pp;
        }
        self.fft
            .begin_to_real_batch(comm, &self.batch_g[..3 * ll], 3)?;

        let slab_len = self.fft.rows() * self.dims[2];
        for i in 0..self.fft.plane_count() {
            let off = i * slab_len;
            self.fft.to_real_plane_in(0, 3, i)?;
            {
                let slab = self.fft.plane();
                let g2_dst = &mut grad_phi_solv_sq_planes[off..off + slab_len];
                for (g2, s) in g2_dst.iter_mut().zip(slab) {
                    *g2 = s.re.mul_add(s.re, s.im * s.im);
                }
            }

            self.fft.to_real_plane_in(1, 3, i)?;
            {
                let slab = self.fft.plane();
                let g2_dst = &mut grad_phi_solv_sq_planes[off..off + slab_len];
                let ps_dst = &mut phi_solv_planes[off..off + slab_len];
                for ((g2, ps), s) in g2_dst.iter_mut().zip(ps_dst).zip(slab) {
                    *g2 = s.re.mul_add(s.re, *g2);
                    *ps = s.im;
                }
            }

            self.fft.to_real_plane_in(2, 3, i)?;
            {
                let slab = self.fft.plane();
                let pp_dst = &mut phi_pol_planes[off..off + slab_len];
                for (pp, s) in pp_dst.iter_mut().zip(slab) {
                    *pp = s.re;
                }
            }
        }
        Ok(())
    }

    /// Computes discrete cavitation surface integral `sum_r (sqrt(|D S|^2 + delta_S^2) - delta_S)`
    /// (local unscaled sum) and `div_n_planes = sum_a D_a (D_a S / sqrt(|D S|^2 + delta_S^2))`
    /// using the exact same skew-adjoint `D_a`.
    fn cavitation_surface_grad_and_div<K: Collective>(
        &mut self,
        comm: &K,
        s_g: &[Complex64],
        delta_s: f64,
        div_g_scratch: &mut [Complex64],
        div_n_planes: &mut [f64],
    ) -> Result<f64, SolvationError> {
        let ll = self.local_len();
        let (b0, b1) = self.batch_g[..2 * ll].split_at_mut(ll);
        for (((d0, d1), &[gx, gy, gz]), &s) in
            b0.iter_mut().zip(b1.iter_mut()).zip(&self.g_cart).zip(s_g)
        {
            *d0 = Complex64::new(-gy, gx) * s;
            *d1 = Complex64::new(0.0, gz) * s;
        }
        self.fft
            .begin_to_real_batch(comm, &self.batch_g[..2 * ll], 2)?;

        let plane_count = self.fft.plane_count();
        let delta2 = delta_s * delta_s;
        let mut local_area_sum = 0.0_f64;

        if let Some(mirror) = self.mirror_idx.as_deref() {
            for i in 0..plane_count {
                self.fft.to_real_plane_in(0, 2, i)?;
                {
                    let slab = self.fft.plane();
                    for ((r0, r1), s) in self.slab_r0.iter_mut().zip(&mut self.slab_r1).zip(slab) {
                        *r0 = s.re;
                        *r1 = s.im;
                    }
                }
                self.fft.to_real_plane_in(1, 2, i)?;
                {
                    let slab = self.fft.plane();
                    for (r2, s) in self.slab_r2.iter_mut().zip(slab) {
                        *r2 = s.re;
                    }
                }

                // Normalize (D_x S, D_y S, D_z S) by mag = sqrt(|D S|^2 + delta_S^2).
                for ((r0, r1), r2) in self
                    .slab_r0
                    .iter_mut()
                    .zip(&mut self.slab_r1)
                    .zip(&mut self.slab_r2)
                {
                    let (dx, dy, dz) = (*r0, *r1, *r2);
                    let mag = (dx * dx + dy * dy + dz * dz + delta2).sqrt();
                    local_area_sum += mag - delta_s;
                    let inv_mag = 1.0 / mag;
                    *r0 = dx * inv_mag;
                    *r1 = dy * inv_mag;
                    *r2 = dz * inv_mag;
                }

                // Pack nx + i*ny in slot 0, nz + i*0 in slot 1.
                for ((dst, &nx), &ny) in self
                    .fft
                    .plane_mut()
                    .iter_mut()
                    .zip(&self.slab_r0)
                    .zip(&self.slab_r1)
                {
                    *dst = Complex64::new(nx, ny);
                }
                self.fft.to_reciprocal_plane_in(0, 2, i)?;

                for (dst, &nz) in self.fft.plane_mut().iter_mut().zip(&self.slab_r2) {
                    *dst = Complex64::new(nz, 0.0);
                }
                self.fft.to_reciprocal_plane_in(1, 2, i)?;
            }
            self.fft
                .finish_to_reciprocal_batch(comm, &mut self.batch_out_g[..2 * ll], 2)?;

            let u0 = &self.batch_out_g[..ll];
            let u1 = &self.batch_out_g[ll..2 * ll];
            let neg_half_i = Complex64::new(0.0, -0.5);
            for (k, dst) in div_g_scratch.iter_mut().enumerate() {
                let mk = mirror[k];
                let u0_k = u0[k];
                let u0_m = u0[mk].conj();
                let u1_k = u1[k];
                let u1_m = u1[mk].conj();
                let nx_g = 0.5 * (u0_k + u0_m);
                let ny_g = neg_half_i * (u0_k - u0_m);
                let nz_g = 0.5 * (u1_k + u1_m);
                let [gx, gy, gz] = self.g_cart[k];
                *dst = Complex64::new(0.0, gx) * nx_g
                    + Complex64::new(0.0, gy) * ny_g
                    + Complex64::new(0.0, gz) * nz_g;
            }
        } else {
            for i in 0..plane_count {
                self.fft.to_real_plane_in(0, 2, i)?;
                {
                    let slab = self.fft.plane();
                    for ((r0, r1), s) in self.slab_r0.iter_mut().zip(&mut self.slab_r1).zip(slab) {
                        *r0 = s.re;
                        *r1 = s.im;
                    }
                }
                self.fft.to_real_plane_in(1, 2, i)?;
                {
                    let slab = self.fft.plane();
                    for (r2, s) in self.slab_r2.iter_mut().zip(slab) {
                        *r2 = s.re;
                    }
                }

                for ((r0, r1), r2) in self
                    .slab_r0
                    .iter_mut()
                    .zip(&mut self.slab_r1)
                    .zip(&mut self.slab_r2)
                {
                    let (dx, dy, dz) = (*r0, *r1, *r2);
                    let mag = (dx * dx + dy * dy + dz * dz + delta2).sqrt();
                    local_area_sum += mag - delta_s;
                    let inv_mag = 1.0 / mag;
                    *r0 = dx * inv_mag;
                    *r1 = dy * inv_mag;
                    *r2 = dz * inv_mag;
                }

                for (dst, &nx) in self.fft.plane_mut().iter_mut().zip(&self.slab_r0) {
                    *dst = Complex64::new(nx, 0.0);
                }
                self.fft.to_reciprocal_plane_in(0, 3, i)?;

                for (dst, &ny) in self.fft.plane_mut().iter_mut().zip(&self.slab_r1) {
                    *dst = Complex64::new(ny, 0.0);
                }
                self.fft.to_reciprocal_plane_in(1, 3, i)?;

                for (dst, &nz) in self.fft.plane_mut().iter_mut().zip(&self.slab_r2) {
                    *dst = Complex64::new(nz, 0.0);
                }
                self.fft.to_reciprocal_plane_in(2, 3, i)?;
            }
            self.fft
                .finish_to_reciprocal_batch(comm, &mut self.batch_out_g[..3 * ll], 3)?;

            let nx_g = &self.batch_out_g[..ll];
            let ny_g = &self.batch_out_g[ll..2 * ll];
            let nz_g = &self.batch_out_g[2 * ll..3 * ll];
            for (k, dst) in div_g_scratch.iter_mut().enumerate() {
                let [gx, gy, gz] = self.g_cart[k];
                *dst = Complex64::new(0.0, gx) * nx_g[k]
                    + Complex64::new(0.0, gy) * ny_g[k]
                    + Complex64::new(0.0, gz) * nz_g[k];
            }
        }

        self.inverse_to_real(comm, div_g_scratch, div_n_planes)?;
        Ok(local_area_sum)
    }
}

/// Inner PDE solver state and preallocated scratch buffers (separated from `EvalScratch`
/// so `evaluate` can pass cavity/RHS slices directly into `solve_pde` with zero copies).
struct PdeSolver {
    engine: FftEngine,
    rhs_eff_g: Vec<Complex64>,
    r_g: Vec<Complex64>,
    p_g: Vec<Complex64>,
    z_g: Vec<Complex64>,
    w_g: Vec<Complex64>,
    aux_g: Vec<Complex64>,
    precond: Vec<f64>,
    prev_phi_g: Vec<Complex64>,
    prev_rhs_g: Vec<Complex64>,
    has_prev: bool,
    prev_has_kappa: bool,
    // Plane-local scratch used by FisicaroPcg, PolarizationDiis, and cavitation divergence.
    work_r0: Vec<f64>,
    work_r1: Vec<f64>,
    work_r2: Vec<f64>,
    work_r3: Vec<f64>,
    // Preallocated Anderson/DIIS history and symmetric eigensolver workspace.
    diis_capacity: usize,
    diis_prev_x: Vec<Complex64>,
    diis_prev_f: Vec<Complex64>,
    diis_dx: Vec<Vec<Complex64>>,
    diis_df: Vec<Vec<Complex64>>,
    diis_gram: RealMatrix,
    diis_q: RealMatrix,
    diis_lambda: Vec<f64>,
    diis_eigh: RealEighWorkspace,
    diis_forms_loc: Vec<f64>,
    diis_forms_glob: Vec<f64>,
    diis_gamma: Vec<f64>,
}

impl PdeSolver {
    fn new(engine: FftEngine) -> Result<Self, SolvationError> {
        let ll = engine.local_len();
        let pp = engine.plane_points();
        let cap = DEFAULT_DIIS_CAP;
        let max_forms = cap * (cap + 1) / 2 + cap;
        Ok(Self {
            engine,
            rhs_eff_g: vec![Complex64::default(); ll],
            r_g: vec![Complex64::default(); ll],
            p_g: vec![Complex64::default(); ll],
            z_g: vec![Complex64::default(); ll],
            w_g: vec![Complex64::default(); ll],
            aux_g: vec![Complex64::default(); ll],
            precond: vec![0.0; ll],
            prev_phi_g: vec![Complex64::default(); ll],
            prev_rhs_g: vec![Complex64::default(); ll],
            has_prev: false,
            prev_has_kappa: false,
            work_r0: vec![0.0; pp],
            work_r1: vec![0.0; pp],
            work_r2: vec![0.0; pp],
            work_r3: vec![0.0; pp],
            diis_capacity: cap,
            diis_prev_x: vec![Complex64::default(); ll],
            diis_prev_f: vec![Complex64::default(); ll],
            diis_dx: vec![vec![Complex64::default(); ll]; cap],
            diis_df: vec![vec![Complex64::default(); ll]; cap],
            diis_gram: RealMatrix::zeros(cap, cap)?,
            diis_q: RealMatrix::zeros(cap, cap)?,
            diis_lambda: vec![0.0; cap],
            diis_eigh: RealEighWorkspace::new(cap)?,
            diis_forms_loc: vec![0.0; max_forms],
            diis_forms_glob: vec![0.0; max_forms],
            diis_gamma: vec![0.0; cap],
        })
    }

    fn ensure_diis_capacity(&mut self, cap: usize) -> Result<(), SolvationError> {
        if cap <= self.diis_capacity {
            return Ok(());
        }
        let ll = self.engine.local_len();
        let max_forms = cap * (cap + 1) / 2 + cap;
        self.diis_dx
            .resize_with(cap, || vec![Complex64::default(); ll]);
        self.diis_df
            .resize_with(cap, || vec![Complex64::default(); ll]);
        self.diis_gram = RealMatrix::zeros(cap, cap)?;
        self.diis_q = RealMatrix::zeros(cap, cap)?;
        self.diis_lambda.resize(cap, 0.0);
        self.diis_eigh = RealEighWorkspace::new(cap)?;
        self.diis_forms_loc.resize(max_forms, 0.0);
        self.diis_forms_glob.resize(max_forms, 0.0);
        self.diis_gamma.resize(cap, 0.0);
        self.diis_capacity = cap;
        Ok(())
    }

    fn reset_warm_start(&mut self) {
        self.has_prev = false;
        self.prev_phi_g.fill(Complex64::default());
        self.prev_rhs_g.fill(Complex64::default());
    }

    fn save_warm_start(&mut self, phi_g: &[Complex64], has_kappa: bool) {
        self.prev_phi_g.copy_from_slice(phi_g);
        self.prev_rhs_g.copy_from_slice(&self.rhs_eff_g);
        self.prev_has_kappa = has_kappa;
        self.has_prev = true;
    }

    fn solve_pde<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        rhs_g: &[Complex64],
        config: &SolverConfig,
        out_phi_g: &mut [Complex64],
    ) -> Result<SolveReport, SolvationError> {
        config.validate()?;
        let ll = self.engine.local_len();
        let pp = self.engine.plane_points();

        // Collective-safe validation of shapes, finiteness, and physical bounds.
        let bad_shape = eps_planes.len() != pp
            || kappa2_planes.len() != pp
            || rhs_g.len() != ll
            || out_phi_g.len() != ll;
        let (nonfinite, bad_eps, bad_kappa, local_kappa2_sum, local_eps_sum) = if bad_shape {
            (false, false, false, 0.0, 0.0)
        } else {
            let nf_eps = eps_planes.iter().any(|v| !v.is_finite());
            let nf_kap = kappa2_planes.iter().any(|v| !v.is_finite());
            let nf_rhs = rhs_g.iter().any(|z| !z.re.is_finite() || !z.im.is_finite());
            let nonfinite = nf_eps || nf_kap || nf_rhs;
            let bad_eps = !nf_eps && eps_planes.iter().any(|&v| v < 1.0);
            let bad_kappa = !nf_kap && kappa2_planes.iter().any(|&v| v < 0.0);
            let kap_sum = if nf_kap {
                0.0
            } else {
                kappa2_planes.iter().copied().sum::<f64>()
            };
            let eps_sum = if nf_eps {
                0.0
            } else {
                eps_planes.iter().copied().sum::<f64>()
            };
            (nonfinite, bad_eps, bad_kappa, kap_sum, eps_sum)
        };
        let mut checks = [0.0_f64; 6];
        comm.tree_sum_f64(
            &[
                f64::from(bad_shape),
                f64::from(nonfinite),
                f64::from(bad_eps),
                f64::from(bad_kappa),
                local_kappa2_sum,
                local_eps_sum,
            ],
            &mut checks,
        )?;
        if checks[0] != 0.0 {
            return Err(SolvationError::Shape);
        }
        if checks[1] != 0.0 {
            return Err(SolvationError::NonFinite);
        }
        if checks[2] != 0.0 {
            return Err(SolvationError::InvalidParameter(
                "epsilon must be finite and >= 1.0",
            ));
        }
        if checks[3] != 0.0 {
            return Err(SolvationError::InvalidParameter(
                "kappa2 must be finite and >= 0.0",
            ));
        }

        let has_kappa = checks[4] > 0.0;
        let kappa2_bar = checks[4] / self.engine.total_points as f64;
        let eps_mean = checks[5] / self.engine.total_points as f64;
        let (eps_min, eps_max) = self.engine.global_min_max(comm, eps_planes)?;
        let eps_bar = (eps_min * eps_max).sqrt();

        // Prepare effective RHS (projected onto range(A) when kappa = 0).
        self.rhs_eff_g.copy_from_slice(rhs_g);
        if !has_kappa {
            for &k in &self.engine.zero_g2_indices {
                self.rhs_eff_g[k] = Complex64::default();
            }
        }
        let b_norm2_loc: f64 = self.rhs_eff_g.iter().map(Complex64::norm_sqr).sum();
        let mut b_norm2_glob = [0.0_f64];
        comm.tree_sum_f64(&[b_norm2_loc], &mut b_norm2_glob)?;
        let b_norm = (self.engine.volume * b_norm2_glob[0]).sqrt();
        if !b_norm.is_finite() {
            return Err(SolvationError::NonFinite);
        }
        if b_norm == 0.0 {
            out_phi_g.fill(Complex64::default());
            if matches!(
                config.kind,
                SolverKind::DirectPcg | SolverKind::FisicaroPcg | SolverKind::HybridFisicaroDirect
            ) {
                self.save_warm_start(out_phi_g, has_kappa);
            }
            return Ok(SolveReport {
                kind: config.kind,
                iterations: 0,
                converged: true,
                residual_norm: 0.0,
                relative_residual: 0.0,
                fft_count: 0,
                mirror_paired: self.engine.mirror_paired(),
                warm_started: false,
            });
        }

        let b_scale = b_norm.max(1.0);
        let tol_abs = config.tolerance * b_scale;

        match config.kind {
            SolverKind::DirectPcg => self.solve_direct_pcg(
                comm,
                eps_planes,
                kappa2_planes,
                has_kappa,
                eps_bar,
                eps_mean,
                kappa2_bar,
                b_scale,
                tol_abs,
                config,
                out_phi_g,
            ),
            SolverKind::FisicaroPcg => self.solve_fisicaro_pcg(
                comm,
                eps_planes,
                kappa2_planes,
                has_kappa,
                eps_max,
                b_scale,
                tol_abs,
                config,
                true,
                out_phi_g,
            ),
            SolverKind::PolarizationDiis => self.solve_polarization_diis(
                comm,
                eps_planes,
                kappa2_planes,
                has_kappa,
                eps_max,
                b_scale,
                tol_abs,
                config,
                out_phi_g,
            ),
            SolverKind::HybridFisicaroDirect => self.solve_hybrid_fisicaro_direct(
                comm,
                eps_planes,
                kappa2_planes,
                has_kappa,
                eps_bar,
                eps_max,
                kappa2_bar,
                b_scale,
                tol_abs,
                config,
                out_phi_g,
            ),
        }
    }

    fn fill_direct_precond(&mut self, eps_bar: f64, kappa2_bar: f64, has_kappa: bool) {
        // 0-FFT strictly-SPD G-space preconditioner with optimal geometric-mean `eps_bar`
        // (minimizing kappa_2(M^{-1} A) <= eps_max/eps_min) and exact `g2` (including the
        // 7 even-Nyquist corner modes where D_a = 0 so A acts solely via kappa2).
        for (p, &g2) in self.precond.iter_mut().zip(&self.engine.g2) {
            *p = if !has_kappa && g2 == 0.0 {
                0.0
            } else {
                1.0 / (eps_bar * g2 + kappa2_bar)
            };
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_direct_pcg<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        has_kappa: bool,
        eps_bar: f64,
        eps_mean: f64,
        kappa2_bar: f64,
        b_scale: f64,
        tol_abs: f64,
        config: &SolverConfig,
        out_phi_g: &mut [Complex64],
    ) -> Result<SolveReport, SolvationError> {
        self.fill_direct_precond(eps_bar, kappa2_bar, has_kappa);

        // Warm-start predictor uses `eps_mean` (matching bulk dielectric screening without
        // overshooting by sqrt(eps_b) in the solvent bulk).
        let warm_started = config.warm_start && self.has_prev && self.prev_has_kappa == has_kappa;
        if warm_started {
            for (((phi, &prev_phi), &g2), (&rhs, &prev_rhs)) in out_phi_g
                .iter_mut()
                .zip(&self.prev_phi_g)
                .zip(&self.engine.g2)
                .zip(self.rhs_eff_g.iter().zip(&self.prev_rhs_g))
            {
                let bulk_pre = if !has_kappa && g2 == 0.0 {
                    0.0
                } else {
                    1.0 / (eps_mean * g2 + kappa2_bar)
                };
                *phi = prev_phi + bulk_pre * (rhs - prev_rhs);
            }
        } else {
            for ((phi, &pre), &rhs) in out_phi_g.iter_mut().zip(&self.precond).zip(&self.rhs_eff_g)
            {
                *phi = pre * rhs;
            }
        }
        if !has_kappa {
            for &k in &self.engine.zero_g2_indices {
                out_phi_g[k] = Complex64::default();
            }
        }

        self.solve_direct_pcg_from_guess(
            comm,
            eps_planes,
            kappa2_planes,
            has_kappa,
            b_scale,
            tol_abs,
            config.max_iterations,
            SolverKind::DirectPcg,
            0,
            0,
            warm_started,
            out_phi_g,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_direct_pcg_from_guess<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        has_kappa: bool,
        b_scale: f64,
        tol_abs: f64,
        max_iterations: usize,
        report_kind: SolverKind,
        base_iterations: usize,
        base_ffts: usize,
        warm_started: bool,
        out_phi_g: &mut [Complex64],
    ) -> Result<SolveReport, SolvationError> {
        let volume = self.engine.volume;
        let mut fft_count = base_ffts
            + self.engine.apply_operator(
                comm,
                eps_planes,
                kappa2_planes,
                has_kappa,
                out_phi_g,
                &mut self.w_g,
            )?;
        if !has_kappa {
            for &k in &self.engine.zero_g2_indices {
                self.w_g[k] = Complex64::default();
            }
        }

        let mut r_norm2_loc = 0.0_f64;
        let mut rz_loc = 0.0_f64;
        for ((((r_dst, z_dst), p_dst), (&rhs, &w)), &pre) in self
            .r_g
            .iter_mut()
            .zip(&mut self.z_g)
            .zip(&mut self.p_g)
            .zip(self.rhs_eff_g.iter().zip(&self.w_g))
            .zip(&self.precond)
        {
            let r = rhs - w;
            let z = pre * r;
            *r_dst = r;
            *z_dst = z;
            *p_dst = z;
            r_norm2_loc += r.norm_sqr();
            rz_loc += r.re * z.re + r.im * z.im;
        }
        let mut sums = [0.0_f64; 2];
        comm.tree_sum_f64(&[r_norm2_loc, rz_loc], &mut sums)?;
        let mut r_norm = (volume * sums[0]).sqrt();
        let mut rz = volume * sums[1];
        if !r_norm.is_finite() || !rz.is_finite() {
            return Err(SolvationError::NonFinite);
        }
        if r_norm <= tol_abs {
            self.save_warm_start(out_phi_g, has_kappa);
            return Ok(SolveReport {
                kind: report_kind,
                iterations: base_iterations,
                converged: true,
                residual_norm: r_norm,
                relative_residual: r_norm / b_scale,
                fft_count,
                mirror_paired: self.engine.mirror_paired(),
                warm_started,
            });
        }

        let mut scalar_glob = [0.0_f64];
        for iter in 1..=max_iterations {
            fft_count += self.engine.apply_operator(
                comm,
                eps_planes,
                kappa2_planes,
                has_kappa,
                &self.p_g,
                &mut self.w_g,
            )?;

            let mut pap_loc = 0.0_f64;
            for (&p, &w) in self.p_g.iter().zip(&self.w_g) {
                pap_loc += p.re * w.re + p.im * w.im;
            }
            comm.tree_sum_f64(&[pap_loc], &mut scalar_glob)?;
            let pap = volume * scalar_glob[0];
            if !pap.is_finite() || pap <= 0.0 {
                return Err(SolvationError::Breakdown(
                    "direct PCG encountered non-positive p^T A p",
                ));
            }
            let alpha = rz / pap;
            if !alpha.is_finite() {
                return Err(SolvationError::Breakdown(
                    "direct PCG encountered non-finite step size",
                ));
            }

            let mut r_norm2_loc = 0.0_f64;
            let mut rz_new_loc = 0.0_f64;
            for ((((phi, r_dst), z_dst), (&p, &w)), &pre) in out_phi_g
                .iter_mut()
                .zip(&mut self.r_g)
                .zip(&mut self.z_g)
                .zip(self.p_g.iter().zip(&self.w_g))
                .zip(&self.precond)
            {
                *phi += alpha * p;
                let r = *r_dst - alpha * w;
                let z = pre * r;
                *r_dst = r;
                *z_dst = z;
                r_norm2_loc += r.norm_sqr();
                rz_new_loc += r.re * z.re + r.im * z.im;
            }
            comm.tree_sum_f64(&[r_norm2_loc, rz_new_loc], &mut sums)?;
            r_norm = (volume * sums[0]).sqrt();
            let rz_new = volume * sums[1];

            if !r_norm.is_finite() || !rz_new.is_finite() {
                return Err(SolvationError::NonFinite);
            }
            if r_norm <= tol_abs {
                if !has_kappa {
                    for &k in &self.engine.zero_g2_indices {
                        out_phi_g[k] = Complex64::default();
                    }
                }
                self.save_warm_start(out_phi_g, has_kappa);
                return Ok(SolveReport {
                    kind: report_kind,
                    iterations: base_iterations + iter,
                    converged: true,
                    residual_norm: r_norm,
                    relative_residual: r_norm / b_scale,
                    fft_count,
                    mirror_paired: self.engine.mirror_paired(),
                    warm_started,
                });
            }
            if rz <= 0.0 {
                return Err(SolvationError::Breakdown(
                    "direct PCG encountered non-positive r^T z",
                ));
            }
            let beta = rz_new / rz;
            for (p, &z) in self.p_g.iter_mut().zip(&self.z_g) {
                *p = z + beta * *p;
            }
            rz = rz_new;
        }

        Err(SolvationError::NotConverged {
            iterations: base_iterations + max_iterations,
            residual_norm: r_norm,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_hybrid_fisicaro_direct<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        has_kappa: bool,
        eps_bar: f64,
        eps_max: f64,
        kappa2_bar: f64,
        b_scale: f64,
        tol_abs: f64,
        config: &SolverConfig,
        out_phi_g: &mut [Complex64],
    ) -> Result<SolveReport, SolvationError> {
        let warm_started = config.warm_start && self.has_prev && self.prev_has_kappa == has_kappa;

        // Stage 1: fast 2-FFT/iter warm-started Fisicaro PCG predictor.
        let (stage1_iters, stage1_ffts) = match self.solve_fisicaro_pcg(
            comm,
            eps_planes,
            kappa2_planes,
            has_kappa,
            eps_max,
            b_scale,
            tol_abs,
            config,
            false,
            out_phi_g,
        ) {
            Ok(rep) => (rep.iterations, rep.fft_count),
            Err(SolvationError::Breakdown(_)) | Err(SolvationError::NotConverged { .. }) => {
                if warm_started {
                    out_phi_g.copy_from_slice(&self.prev_phi_g);
                } else {
                    out_phi_g.fill(Complex64::default());
                }
                if !has_kappa {
                    for &k in &self.engine.zero_g2_indices {
                        out_phi_g[k] = Complex64::default();
                    }
                }
                (0, 0)
            }
            Err(err) => return Err(err),
        };

        // Stage 2: exact discrete SPD DirectPcg verification & defect-correction polish.
        self.fill_direct_precond(eps_bar, kappa2_bar, has_kappa);
        self.solve_direct_pcg_from_guess(
            comm,
            eps_planes,
            kappa2_planes,
            has_kappa,
            b_scale,
            tol_abs,
            config.max_iterations,
            SolverKind::HybridFisicaroDirect,
            stage1_iters,
            stage1_ffts,
            warm_started,
            out_phi_g,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_fisicaro_pcg<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        has_kappa: bool,
        eps_max: f64,
        b_scale: f64,
        tol_abs: f64,
        config: &SolverConfig,
        save_warm: bool,
        out_phi_g: &mut [Complex64],
    ) -> Result<SolveReport, SolvationError> {
        let ll = self.engine.local_len();
        let volume = self.engine.volume;
        let eps_max_sqrt = eps_max.sqrt();
        let tol_tilde = tol_abs / eps_max_sqrt;

        // Setup (3 FFTs):
        //   work_r0 = sqrt(eps)
        //   aux_g   = F(sqrt(eps))
        //   work_r1 = q(r) = (-L_0 sqrt(eps) + kappa2 / sqrt(eps)) / sqrt(eps)
        //   work_r2 = b(r) / sqrt(eps)
        //   r_g     = b_tilde_g = F(b / sqrt(eps))
        let mut kap_over_eps_loc = 0.0_f64;
        for ((r0, &e), &k2) in self.work_r0.iter_mut().zip(eps_planes).zip(kappa2_planes) {
            *r0 = e.sqrt();
            kap_over_eps_loc += k2 / e;
        }
        self.engine
            .forward_real(comm, &self.work_r0, &mut self.aux_g)?;
        let mut fft_count = 1usize;

        // Paired single inverse FFT: Re = -L_0 sqrt(eps) = +nabla^2 sqrt(eps), Im = b(r).
        for (((dst, &g2), &aux), &rhs) in self.engine.batch_g[..ll]
            .iter_mut()
            .zip(&self.engine.g2)
            .zip(&self.aux_g)
            .zip(&self.rhs_eff_g)
        {
            *dst = -g2 * aux + Complex64::new(0.0, 1.0) * rhs;
        }
        self.engine
            .fft
            .begin_to_real(comm, &self.engine.batch_g[..ll])?;
        let slab_len = self.engine.fft.rows() * self.engine.dims[2];
        for i in 0..self.engine.fft.plane_count() {
            self.engine.fft.to_real_plane(i)?;
            let slab = self.engine.fft.plane();
            let off = i * slab_len;
            let range = off..off + slab_len;
            for ((((r1, r2), &s), &k2), &sl) in self.work_r1[range.clone()]
                .iter_mut()
                .zip(&mut self.work_r2[range.clone()])
                .zip(&self.work_r0[range.clone()])
                .zip(&kappa2_planes[range])
                .zip(slab)
            {
                *r1 = (sl.re + k2 / s) / s;
                *r2 = sl.im / s;
            }
        }
        fft_count += 1;

        self.engine
            .forward_real(comm, &self.work_r2, &mut self.r_g)?;
        fft_count += 1;

        // Global s_norm2 = ||sqrt(eps)||_2^2 and mean(kappa2 / eps).
        let s_norm2_loc: f64 = self.aux_g.iter().map(Complex64::norm_sqr).sum();
        let mut setup_glob = [0.0_f64; 2];
        comm.tree_sum_f64(&[s_norm2_loc, kap_over_eps_loc], &mut setup_glob)?;
        let s_norm2 = volume * setup_glob[0];
        let kappa2_over_eps_bar = setup_glob[1] / self.engine.total_points as f64;

        let g2_min = self.engine.g2_min;
        for (p, &g2) in self.precond.iter_mut().zip(&self.engine.g2) {
            *p = if has_kappa {
                1.0 / (g2 + kappa2_over_eps_bar)
            } else if g2 == 0.0 {
                1.0 / g2_min
            } else {
                1.0 / g2
            };
        }

        // Deflate b_tilde orthogonal to nullspace vector sqrt(eps) when kappa = 0.
        if !has_kappa && s_norm2 > 0.0 {
            let mut sb_loc = 0.0_f64;
            for (&s, &b) in self.aux_g.iter().zip(&self.r_g) {
                sb_loc += s.re * b.re + s.im * b.im;
            }
            let mut sb_glob = [0.0_f64];
            comm.tree_sum_f64(&[sb_loc], &mut sb_glob)?;
            let coeff = (volume * sb_glob[0]) / s_norm2;
            for (r, &s) in self.r_g.iter_mut().zip(&self.aux_g) {
                *r -= coeff * s;
            }
        }

        // Initial guess phi_tilde_0: warm-started F[sqrt(eps) * prev_phi] or M_0^{-1} b_tilde
        // (deflated orthogonal to sqrt(eps) when kappa = 0).
        let warm_started = config.warm_start && self.has_prev && self.prev_has_kappa == has_kappa;
        if warm_started {
            out_phi_g.copy_from_slice(&self.prev_phi_g);
            self.engine
                .fft
                .apply_local_batch_planes(comm, out_phi_g, &self.work_r0, 1)?;
            fft_count += 2;
        } else {
            for ((phi, &pre), &r) in out_phi_g.iter_mut().zip(&self.precond).zip(&self.r_g) {
                *phi = pre * r;
            }
        }
        if !has_kappa && s_norm2 > 0.0 {
            let mut sp_loc = 0.0_f64;
            for (&s, &p) in self.aux_g.iter().zip(out_phi_g.iter()) {
                sp_loc += s.re * p.re + s.im * p.im;
            }
            let mut sp_glob = [0.0_f64];
            comm.tree_sum_f64(&[sp_loc], &mut sp_glob)?;
            let coeff = (volume * sp_glob[0]) / s_norm2;
            for (phi, &s) in out_phi_g.iter_mut().zip(&self.aux_g) {
                *phi -= coeff * s;
            }
        }

        // Initial residual r_0 = b_tilde - (L_0 + q) phi_tilde_0.
        self.w_g.copy_from_slice(out_phi_g);
        self.engine
            .fft
            .apply_local_batch_planes(comm, &mut self.w_g, &self.work_r1, 1)?;
        fft_count += 2;
        for ((((w, r), z), (&g2, &pre)), &phi) in self
            .w_g
            .iter_mut()
            .zip(&mut self.r_g)
            .zip(&mut self.z_g)
            .zip(self.engine.g2.iter().zip(&self.precond))
            .zip(out_phi_g.iter())
        {
            *w += g2 * phi;
            *r -= *w;
            *z = pre * *r;
        }
        if !has_kappa && s_norm2 > 0.0 {
            let mut sr_loc = 0.0_f64;
            let mut sz_loc = 0.0_f64;
            for ((&s, &r), &z) in self.aux_g.iter().zip(&self.r_g).zip(&self.z_g) {
                sr_loc += s.re * r.re + s.im * r.im;
                sz_loc += s.re * z.re + s.im * z.im;
            }
            let mut proj = [0.0_f64; 2];
            comm.tree_sum_f64(&[sr_loc, sz_loc], &mut proj)?;
            let cr = (volume * proj[0]) / s_norm2;
            let cz = (volume * proj[1]) / s_norm2;
            for ((r, z), &s) in self.r_g.iter_mut().zip(&mut self.z_g).zip(&self.aux_g) {
                *r -= cr * s;
                *z -= cz * s;
            }
        }

        let mut r_norm2_loc = 0.0_f64;
        let mut rz_loc = 0.0_f64;
        for ((p, &r), &z) in self.p_g.iter_mut().zip(&self.r_g).zip(&self.z_g) {
            *p = z;
            r_norm2_loc += r.norm_sqr();
            rz_loc += r.re * z.re + r.im * z.im;
        }
        let mut sums = [0.0_f64; 2];
        comm.tree_sum_f64(&[r_norm2_loc, rz_loc], &mut sums)?;
        let mut r_tilde_norm = (volume * sums[0]).sqrt();
        let mut rz = volume * sums[1];
        if !r_tilde_norm.is_finite() || !rz.is_finite() {
            return Err(SolvationError::NonFinite);
        }

        let mut converged_iter = (r_tilde_norm <= tol_tilde).then_some(0usize);
        let mut pap_glob = [0.0_f64];

        if converged_iter.is_none() {
            for iter in 1..=config.max_iterations {
                self.w_g.copy_from_slice(&self.p_g);
                self.engine
                    .fft
                    .apply_local_batch_planes(comm, &mut self.w_g, &self.work_r1, 1)?;
                fft_count += 2;

                let mut pap_loc = 0.0_f64;
                for ((w_dst, &g2), &p) in self.w_g.iter_mut().zip(&self.engine.g2).zip(&self.p_g) {
                    let w = *w_dst + g2 * p;
                    *w_dst = w;
                    pap_loc += p.re * w.re + p.im * w.im;
                }
                comm.tree_sum_f64(&[pap_loc], &mut pap_glob)?;
                let pap = volume * pap_glob[0];
                if !pap.is_finite() || pap <= 0.0 {
                    return Err(SolvationError::Breakdown(
                        "Fisicaro operator became indefinite",
                    ));
                }
                let alpha = rz / pap;
                if !alpha.is_finite() {
                    return Err(SolvationError::Breakdown(
                        "Fisicaro operator became indefinite",
                    ));
                }

                for ((((phi, r_dst), z_dst), (&p, &w)), &pre) in out_phi_g
                    .iter_mut()
                    .zip(&mut self.r_g)
                    .zip(&mut self.z_g)
                    .zip(self.p_g.iter().zip(&self.w_g))
                    .zip(&self.precond)
                {
                    *phi += alpha * p;
                    let r = *r_dst - alpha * w;
                    *r_dst = r;
                    *z_dst = pre * r;
                }
                if !has_kappa && s_norm2 > 0.0 {
                    let mut sr_loc = 0.0_f64;
                    let mut sz_loc = 0.0_f64;
                    for ((&s, &r), &z) in self.aux_g.iter().zip(&self.r_g).zip(&self.z_g) {
                        sr_loc += s.re * r.re + s.im * r.im;
                        sz_loc += s.re * z.re + s.im * z.im;
                    }
                    let mut proj = [0.0_f64; 2];
                    comm.tree_sum_f64(&[sr_loc, sz_loc], &mut proj)?;
                    let cr = (volume * proj[0]) / s_norm2;
                    let cz = (volume * proj[1]) / s_norm2;
                    for ((r, z), &s) in self.r_g.iter_mut().zip(&mut self.z_g).zip(&self.aux_g) {
                        *r -= cr * s;
                        *z -= cz * s;
                    }
                }

                let mut r_norm2_loc = 0.0_f64;
                let mut rz_new_loc = 0.0_f64;
                for (&r, &z) in self.r_g.iter().zip(&self.z_g) {
                    r_norm2_loc += r.norm_sqr();
                    rz_new_loc += r.re * z.re + r.im * z.im;
                }
                comm.tree_sum_f64(&[r_norm2_loc, rz_new_loc], &mut sums)?;
                r_tilde_norm = (volume * sums[0]).sqrt();
                let rz_new = volume * sums[1];
                if !r_tilde_norm.is_finite() || !rz_new.is_finite() {
                    return Err(SolvationError::NonFinite);
                }
                if r_tilde_norm <= tol_tilde {
                    converged_iter = Some(iter);
                    break;
                }
                if rz <= 0.0 {
                    return Err(SolvationError::Breakdown(
                        "Fisicaro operator became indefinite",
                    ));
                }
                let beta = rz_new / rz;
                for (p, &z) in self.p_g.iter_mut().zip(&self.z_g) {
                    *p = z + beta * *p;
                }
                rz = rz_new;
            }
        }

        let Some(iterations) = converged_iter else {
            return Err(SolvationError::NotConverged {
                iterations: config.max_iterations,
                residual_norm: r_tilde_norm * eps_max_sqrt,
            });
        };

        // Reconstruct phi(r) = phi_tilde(r) / sqrt(eps(r)) and transform back to G-space (2 FFTs).
        self.engine
            .inverse_to_real(comm, out_phi_g, &mut self.work_r2)?;
        for (dst, &s) in self.work_r2.iter_mut().zip(&self.work_r0) {
            *dst /= s;
        }
        self.engine.forward_real(comm, &self.work_r2, out_phi_g)?;
        fft_count += 2;
        if !has_kappa {
            for &k in &self.engine.zero_g2_indices {
                out_phi_g[k] = Complex64::default();
            }
        }

        if save_warm {
            self.save_warm_start(out_phi_g, has_kappa);
        }

        let res_est = r_tilde_norm * eps_max_sqrt;
        Ok(SolveReport {
            kind: SolverKind::FisicaroPcg,
            iterations,
            converged: true,
            residual_norm: res_est,
            relative_residual: res_est / b_scale,
            fft_count,
            mirror_paired: self.engine.mirror_paired(),
            warm_started,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_polarization_diis<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        has_kappa: bool,
        eps_max: f64,
        b_scale: f64,
        tol_abs: f64,
        config: &SolverConfig,
        out_phi_g: &mut [Complex64],
    ) -> Result<SolveReport, SolvationError> {
        let m = config.diis_history;
        self.ensure_diis_capacity(m)?;
        let ll = self.engine.local_len();
        let volume = self.engine.volume;

        // Setup (4 FFTs):
        //   work_r0 = D_x(ln eps), work_r1 = D_y(ln eps), work_r2 = D_z(ln eps), work_r3 = b / eps
        let mut kap_over_eps_loc = 0.0_f64;
        for ((r0, &e), &k2) in self.work_r0.iter_mut().zip(eps_planes).zip(kappa2_planes) {
            *r0 = e.ln();
            kap_over_eps_loc += k2 / e;
        }
        let mut kap_glob = [0.0_f64];
        comm.tree_sum_f64(&[kap_over_eps_loc], &mut kap_glob)?;
        let kappa2_eff_bar = kap_glob[0] / self.engine.total_points as f64;

        self.engine
            .forward_real(comm, &self.work_r0, &mut self.aux_g)?;
        let mut fft_count = 1usize;

        // Paired np = 2 inverse FFT:
        //   slot 0: Re = D_x(ln eps), Im = D_y(ln eps)
        //   slot 1: Re = D_z(ln eps), Im = b(r)
        let (b0, b1) = self.engine.batch_g[..2 * ll].split_at_mut(ll);
        for ((((d0, d1), &[gx, gy, gz]), &le), &rhs) in b0
            .iter_mut()
            .zip(b1.iter_mut())
            .zip(&self.engine.g_cart)
            .zip(&self.aux_g)
            .zip(&self.rhs_eff_g)
        {
            *d0 = Complex64::new(-gy, gx) * le;
            *d1 = Complex64::new(0.0, gz) * le + Complex64::new(0.0, 1.0) * rhs;
        }
        self.engine
            .fft
            .begin_to_real_batch(comm, &self.engine.batch_g[..2 * ll], 2)?;
        let slab_len = self.engine.fft.rows() * self.engine.dims[2];
        let plane_count = self.engine.fft.plane_count();
        for i in 0..plane_count {
            let off = i * slab_len;
            let range = off..off + slab_len;
            self.engine.fft.to_real_plane_in(0, 2, i)?;
            {
                let slab = self.engine.fft.plane();
                for ((r0, r1), s) in self.work_r0[range.clone()]
                    .iter_mut()
                    .zip(&mut self.work_r1[range.clone()])
                    .zip(slab)
                {
                    *r0 = s.re;
                    *r1 = s.im;
                }
            }
            self.engine.fft.to_real_plane_in(1, 2, i)?;
            {
                let slab = self.engine.fft.plane();
                for (((r2, r3), &e), s) in self.work_r2[range.clone()]
                    .iter_mut()
                    .zip(&mut self.work_r3[range.clone()])
                    .zip(&eps_planes[range])
                    .zip(slab)
                {
                    *r2 = s.re;
                    *r3 = s.im / e;
                }
            }
        }
        fft_count += 2;

        // Forward FFT of b / eps for initial guess (L_0 + kappa2_eff_bar)^-1 (b / eps).
        self.engine
            .forward_real(comm, &self.work_r3, &mut self.aux_g)?;
        fft_count += 1;

        for ((phi, &g2), &aux) in out_phi_g.iter_mut().zip(&self.engine.g2).zip(&self.aux_g) {
            *phi = if !has_kappa && g2 == 0.0 {
                Complex64::default()
            } else {
                aux / (g2 + kappa2_eff_bar)
            };
        }

        let mut last_res = f64::INFINITY;
        let mut r_norm2_glob = [0.0_f64];

        for iter in 1..=config.max_iterations {
            // Fixed-point map (3 FFTs):
            //   (L_0 + kappa2_eff_bar) y = b/eps + (D ln eps).(D x)
            //                              - (kappa2/eps - kappa2_eff_bar) x
            let (b0, b1) = self.engine.batch_g[..2 * ll].split_at_mut(ll);
            for (((d0, d1), &[gx, gy, gz]), &x) in b0
                .iter_mut()
                .zip(b1.iter_mut())
                .zip(&self.engine.g_cart)
                .zip(out_phi_g.iter())
            {
                *d0 = Complex64::new(-gy, gx) * x;
                *d1 = Complex64::new(0.0, gz + 1.0) * x;
            }
            self.engine
                .fft
                .begin_to_real_batch(comm, &self.engine.batch_g[..2 * ll], 2)?;

            for i in 0..plane_count {
                let off = i * slab_len;
                let range = off..off + slab_len;
                self.engine.fft.to_real_plane_in(0, 2, i)?;
                {
                    let slab = self.engine.fft.plane();
                    for ((((s0, &r3), &r0), &r1), s) in self
                        .engine
                        .slab_r0
                        .iter_mut()
                        .zip(&self.work_r3[range.clone()])
                        .zip(&self.work_r0[range.clone()])
                        .zip(&self.work_r1[range.clone()])
                        .zip(slab)
                    {
                        *s0 = r3 + r0 * s.re + r1 * s.im;
                    }
                }
                self.engine.fft.to_real_plane_in(1, 2, i)?;
                {
                    let slab = self.engine.fft.plane_mut();
                    for ((((s, &s0), &r2), &k2), &e) in slab
                        .iter_mut()
                        .zip(&self.engine.slab_r0)
                        .zip(&self.work_r2[range.clone()])
                        .zip(&kappa2_planes[range.clone()])
                        .zip(&eps_planes[range])
                    {
                        let scr_shift = k2 / e - kappa2_eff_bar;
                        let rhs_tot = s0 + r2 * s.re - scr_shift * s.im;
                        *s = Complex64::new(rhs_tot, 0.0);
                    }
                }
                self.engine.fft.to_reciprocal_plane_in(0, 1, i)?;
            }
            self.engine.fft.finish_to_reciprocal(comm, &mut self.w_g)?;
            fft_count += 3;

            let mut r_norm2_loc = 0.0_f64;
            for (((r_dst, &g2), &w), &phi) in self
                .r_g
                .iter_mut()
                .zip(&self.engine.g2)
                .zip(&self.w_g)
                .zip(out_phi_g.iter())
            {
                if !has_kappa && g2 == 0.0 {
                    *r_dst = Complex64::default();
                    continue;
                }
                let denom = g2 + kappa2_eff_bar;
                let r_scaled = w - denom * phi;
                r_norm2_loc += r_scaled.norm_sqr();
                *r_dst = w / denom - phi;
            }
            comm.tree_sum_f64(&[r_norm2_loc], &mut r_norm2_glob)?;
            let res_est = eps_max * (volume * r_norm2_glob[0]).sqrt();
            if !res_est.is_finite() {
                return Err(SolvationError::NonFinite);
            }
            last_res = res_est;
            if res_est <= tol_abs {
                return Ok(SolveReport {
                    kind: SolverKind::PolarizationDiis,
                    iterations: iter,
                    converged: true,
                    residual_norm: res_est,
                    relative_residual: res_est / b_scale,
                    fft_count,
                    mirror_paired: self.engine.mirror_paired(),
                    warm_started: false,
                });
            }

            if iter == 1 {
                self.diis_prev_x.copy_from_slice(out_phi_g);
                self.diis_prev_f.copy_from_slice(&self.r_g);
                for (phi, &r) in out_phi_g.iter_mut().zip(&self.r_g) {
                    *phi += 0.5 * r;
                }
            } else {
                let slot = (iter - 2) % m;
                for ((((dx, df), &phi), &px), (&r, &pf)) in self.diis_dx[slot]
                    .iter_mut()
                    .zip(&mut self.diis_df[slot])
                    .zip(out_phi_g.iter())
                    .zip(&self.diis_prev_x)
                    .zip(self.r_g.iter().zip(&self.diis_prev_f))
                {
                    *dx = phi - px;
                    *df = r - pf;
                }
                self.diis_prev_x.copy_from_slice(out_phi_g);
                self.diis_prev_f.copy_from_slice(&self.r_g);

                let p = (iter - 1).min(m);
                let mut form_idx = 0usize;
                for i in 0..p {
                    for j in 0..=i {
                        let mut dot = 0.0_f64;
                        for (&a, &b) in self.diis_df[i].iter().zip(&self.diis_df[j]) {
                            dot += a.re * b.re + a.im * b.im;
                        }
                        self.diis_forms_loc[form_idx] = dot;
                        form_idx += 1;
                    }
                }
                for df_i in &self.diis_df[..p] {
                    let mut dot = 0.0_f64;
                    for (&a, &b) in df_i.iter().zip(&self.r_g) {
                        dot += a.re * b.re + a.im * b.im;
                    }
                    self.diis_forms_loc[form_idx] = dot;
                    form_idx += 1;
                }
                comm.tree_sum_f64(
                    &self.diis_forms_loc[..form_idx],
                    &mut self.diis_forms_glob[..form_idx],
                )?;

                self.diis_gram.reshape_within_storage(p, p)?;
                self.diis_q.reshape_within_storage(p, p)?;
                let mut idx = 0usize;
                for i in 0..p {
                    for j in 0..=i {
                        let v = volume * self.diis_forms_glob[idx];
                        self.diis_gram[(i, j)] = v;
                        self.diis_gram[(j, i)] = v;
                        idx += 1;
                    }
                }
                let rhs_start = idx;
                self.diis_eigh.solve_into(
                    &self.diis_gram,
                    &mut self.diis_lambda[..p],
                    &mut self.diis_q,
                )?;

                self.diis_gamma[..p].fill(0.0);
                let max_eig = self.diis_lambda[p - 1];
                if max_eig.is_finite() && max_eig > 0.0 {
                    let cutoff = 1.0e-12 * max_eig;
                    for j in 0..p {
                        let eig = self.diis_lambda[j];
                        if eig <= cutoff {
                            continue;
                        }
                        let mut proj = 0.0_f64;
                        for i in 0..p {
                            let rhs_i = volume * self.diis_forms_glob[rhs_start + i];
                            proj += self.diis_q[(i, j)] * rhs_i;
                        }
                        let scale = proj / eig;
                        for (i, gamma) in self.diis_gamma[..p].iter_mut().enumerate() {
                            *gamma += self.diis_q[(i, j)] * scale;
                        }
                    }
                }

                const MIX_BETA: f64 = 0.7;
                for (k, (phi, &r)) in out_phi_g.iter_mut().zip(&self.r_g).enumerate() {
                    let mut x_bar = *phi;
                    let mut f_bar = r;
                    for ((&g, dx_i), df_i) in self.diis_gamma[..p]
                        .iter()
                        .zip(&self.diis_dx[..p])
                        .zip(&self.diis_df[..p])
                    {
                        x_bar -= g * dx_i[k];
                        f_bar -= g * df_i[k];
                    }
                    *phi = x_bar + MIX_BETA * f_bar;
                }
            }

            if !has_kappa {
                for &k in &self.engine.zero_g2_indices {
                    out_phi_g[k] = Complex64::default();
                }
            }
        }

        Err(SolvationError::NotConverged {
            iterations: config.max_iterations,
            residual_norm: last_res,
        })
    }
}

/// Preallocated buffers for `SolvationWorkspace::evaluate`.
struct EvalScratch {
    s: Vec<f64>,
    ds_dn: Vec<f64>,
    eps: Vec<f64>,
    deps_dn: Vec<f64>,
    kappa2: Vec<f64>,
    dkappa2_dn: Vec<f64>,
    phi_solv: Vec<f64>,
    phi_pol: Vec<f64>,
    grad_phi_sq: Vec<f64>,
    rho_ion_g: Vec<Complex64>,
    n_g: Vec<Complex64>,
    rho_sol_g: Vec<Complex64>,
    phi_vac_g: Vec<Complex64>,
    b_pol_g: Vec<Complex64>,
    phi_pol_g: Vec<Complex64>,
    phi_solv_g: Vec<Complex64>,
}

impl EvalScratch {
    fn new(ll: usize, pp: usize) -> Self {
        Self {
            s: vec![0.0; pp],
            ds_dn: vec![0.0; pp],
            eps: vec![0.0; pp],
            deps_dn: vec![0.0; pp],
            kappa2: vec![0.0; pp],
            dkappa2_dn: vec![0.0; pp],
            phi_solv: vec![0.0; pp],
            phi_pol: vec![0.0; pp],
            grad_phi_sq: vec![0.0; pp],
            rho_ion_g: vec![Complex64::default(); ll],
            n_g: vec![Complex64::default(); ll],
            rho_sol_g: vec![Complex64::default(); ll],
            phi_vac_g: vec![Complex64::default(); ll],
            b_pol_g: vec![Complex64::default(); ll],
            phi_pol_g: vec![Complex64::default(); ll],
            phi_solv_g: vec![Complex64::default(); ll],
        }
    }
}

/// Reusable MPI-distributed implicit solvation workspace.
///
/// All FFT plans, transpose buffers, G-space vectors, DIIS history, and
/// plane-local scratch arrays are allocated once at construction.
pub struct SolvationWorkspace {
    solver: PdeSolver,
    eval: EvalScratch,
}

impl SolvationWorkspace {
    /// Constructs a full-box distributed FFT plan and solvation workspace.
    pub fn new<K: Collective>(
        cell: &Cell,
        dims: [usize; 3],
        comm: &K,
    ) -> Result<Self, SolvationError> {
        if dims.contains(&0) {
            return Err(SolvationError::InvalidParameter(
                "solvation grid dimensions must be positive",
            ));
        }
        let total = dims[0]
            .checked_mul(dims[1])
            .and_then(|n| n.checked_mul(dims[2]))
            .ok_or(SolvationError::Shape)?;
        let mut indices = Vec::with_capacity(total);
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    indices.push([i, j, k]);
                }
            }
        }
        let fft = DistributedPairFft::new_direct(dims, &indices, comm.rank(), comm.size())?;
        Self::new_with_fft(cell, fft, comm)
    }

    /// Constructs a solvation workspace from an existing full-box `DistributedPairFft`.
    pub fn new_with_fft<K: Collective>(
        cell: &Cell,
        fft: DistributedPairFft,
        comm: &K,
    ) -> Result<Self, SolvationError> {
        let engine = FftEngine::new(cell, fft, comm)?;
        let ll = engine.local_len();
        let pp = engine.plane_points();
        let solver = PdeSolver::new(engine)?;
        let eval = EvalScratch::new(ll, pp);
        Ok(Self { solver, eval })
    }

    pub fn dims(&self) -> [usize; 3] {
        self.solver.engine.dims
    }

    pub fn local_positions(&self) -> &[usize] {
        self.solver.engine.fft.local_positions()
    }

    pub fn plane_points(&self) -> usize {
        self.solver.engine.plane_points()
    }

    pub fn mirror_paired(&self) -> bool {
        self.solver.engine.mirror_paired()
    }

    /// Maximum real-space grid spacing `max_a |a_a| / N_a` in Bohr.
    ///
    /// Used to guard `SolverKind::FisicaroPcg`, whose product-rule aliasing error
    /// decays below `1e-5 Ha` once `max_grid_spacing_bohr() <= 0.25` Bohr.
    pub fn max_grid_spacing_bohr(&self) -> f64 {
        let cols = self.solver.engine.cell.lattice_columns_bohr();
        let dims = self.solver.engine.dims;
        (0..3)
            .map(|a| {
                let norm =
                    (cols[0][a] * cols[0][a] + cols[1][a] * cols[1][a] + cols[2][a] * cols[2][a])
                        .sqrt();
                norm / dims[a] as f64
            })
            .fold(0.0_f64, f64::max)
    }

    /// Global C-order grid indices of owned real-space points in exact
    /// `DistributedPairFft` plane-local order `[plane][row][z]`.
    pub fn owned_real_grid_indices(&self) -> Vec<usize> {
        let fft = &self.solver.engine.fft;
        let nz = self.solver.engine.dims[2];
        let rows = fft.rows();
        let mut out = Vec::with_capacity(fft.plane_points());
        for i in 0..fft.plane_count() {
            for y in 0..rows {
                let base = fft.grid_start(i, y);
                for z in 0..nz {
                    out.push(base + z);
                }
            }
        }
        out
    }

    /// Cartesian coordinates (Bohr) of owned real-space grid points in plane-local order.
    pub fn owned_real_coords_bohr(&self) -> Vec<[f64; 3]> {
        let [n0, n1, n2] = self.solver.engine.dims;
        let n12 = n1 * n2;
        let inv0 = 1.0 / n0 as f64;
        let inv1 = 1.0 / n1 as f64;
        let inv2 = 1.0 / n2 as f64;
        let cell = &self.solver.engine.cell;
        self.owned_real_grid_indices()
            .into_iter()
            .map(|g| {
                let ix = g / n12;
                let iy = (g / n2) % n1;
                let iz = g % n2;
                cell.fractional_to_cartesian_bohr([
                    ix as f64 * inv0,
                    iy as f64 * inv1,
                    iz as f64 * inv2,
                ])
            })
            .collect()
    }

    /// Forward normalized 3D FFT from owned real-space planes to owned G-space sticks.
    pub fn forward_real<K: Collective>(
        &mut self,
        comm: &K,
        real_planes: &[f64],
        out_g: &mut [Complex64],
    ) -> Result<(), SolvationError> {
        self.solver.engine.forward_real(comm, real_planes, out_g)
    }

    /// Inverse unscaled 3D FFT from owned G-space sticks to owned real-space planes.
    pub fn inverse_to_real<K: Collective>(
        &mut self,
        comm: &K,
        coeff_g: &[Complex64],
        out_planes: &mut [f64],
    ) -> Result<(), SolvationError> {
        self.solver
            .engine
            .inverse_to_real(comm, coeff_g, out_planes)
    }

    /// Clears retained warm-start state so the next `DirectPcg` solve cold-starts.
    pub fn reset_warm_start(&mut self) {
        self.solver.reset_warm_start();
    }

    /// Applies the exact discrete SPD operator `A x = sum_a D_a^T (eps D_a x) + kappa2 * x`.
    pub fn apply_operator<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        x_g: &[Complex64],
        out_g: &mut [Complex64],
    ) -> Result<(), SolvationError> {
        let ll = self.solver.engine.local_len();
        let pp = self.solver.engine.plane_points();
        let bad_shape = eps_planes.len() != pp
            || kappa2_planes.len() != pp
            || x_g.len() != ll
            || out_g.len() != ll;
        let local_has_kappa = !bad_shape && kappa2_planes.iter().any(|&v| v != 0.0);
        let mut flags = [0.0_f64; 2];
        comm.tree_sum_f64(
            &[f64::from(bad_shape), f64::from(local_has_kappa)],
            &mut flags,
        )?;
        if flags[0] != 0.0 {
            return Err(SolvationError::Shape);
        }
        let has_kappa = flags[1] != 0.0;
        let _ = self.solver.engine.apply_operator(
            comm,
            eps_planes,
            kappa2_planes,
            has_kappa,
            x_g,
            out_g,
        )?;
        Ok(())
    }

    /// Solves `A phi = rhs` in G-space using the solver selected in `config.kind`.
    pub fn solve_pde<K: Collective>(
        &mut self,
        comm: &K,
        eps_planes: &[f64],
        kappa2_planes: &[f64],
        rhs_g: &[Complex64],
        config: &SolverConfig,
        out_phi_g: &mut [Complex64],
    ) -> Result<SolveReport, SolvationError> {
        self.solver
            .solve_pde(comm, eps_planes, kappa2_planes, rhs_g, config, out_phi_g)
    }

    /// Evaluates the full variational implicit solvation functional (energies,
    /// Kohn-Sham potential components on owned planes, electrostatic potentials,
    /// and Gaussian ionic pseudo-charge forces).
    pub fn evaluate<K: Collective>(
        &mut self,
        comm: &K,
        density_planes: &[f64],
        ions: &[IonicPseudoCharge],
        params: &SolvationParams,
        config: &SolverConfig,
    ) -> Result<SolvationEvaluation, SolvationError> {
        self.evaluate_inner(comm, density_planes, None, ions, params, config)
    }

    /// Evaluates the full variational implicit solvation functional while reusing
    /// caller-supplied reciprocal-space electron density `density_g` (in owned
    /// `local_positions()` order), saving one full 3D forward FFT per SCF map.
    pub fn evaluate_with_density_g<K: Collective>(
        &mut self,
        comm: &K,
        density_planes: &[f64],
        density_g: &[Complex64],
        ions: &[IonicPseudoCharge],
        params: &SolvationParams,
        config: &SolverConfig,
    ) -> Result<SolvationEvaluation, SolvationError> {
        self.evaluate_inner(comm, density_planes, Some(density_g), ions, params, config)
    }

    fn evaluate_inner<K: Collective>(
        &mut self,
        comm: &K,
        density_planes: &[f64],
        density_g: Option<&[Complex64]>,
        ions: &[IonicPseudoCharge],
        params: &SolvationParams,
        config: &SolverConfig,
    ) -> Result<SolvationEvaluation, SolvationError> {
        params.validate()?;
        config.validate()?;

        let pp = self.plane_points();
        let ll = self.solver.engine.local_len();
        let bad_shape = density_planes.len() != pp || density_g.is_some_and(|dg| dg.len() != ll);
        let bad_nonfinite = !bad_shape
            && (density_planes.iter().any(|v| !v.is_finite())
                || density_g
                    .is_some_and(|dg| dg.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())));
        let local_n_sum = if bad_shape || bad_nonfinite {
            0.0
        } else {
            density_planes.iter().copied().sum::<f64>()
        };
        let mut check = [0.0_f64; 3];
        comm.tree_sum_f64(
            &[f64::from(bad_shape), f64::from(bad_nonfinite), local_n_sum],
            &mut check,
        )?;
        if check[0] != 0.0 {
            return Err(SolvationError::Shape);
        }
        if check[1] != 0.0 {
            return Err(SolvationError::NonFinite);
        }
        let n_el = self.solver.engine.dv * check[2];

        cavity::evaluate_cavity_fields(
            density_planes,
            params,
            &mut self.eval.s,
            &mut self.eval.ds_dn,
            &mut self.eval.eps,
            &mut self.eval.deps_dn,
            &mut self.eval.kappa2,
            &mut self.eval.dkappa2_dn,
        )?;

        let local_kappa2_sum: f64 = self.eval.kappa2.iter().copied().sum();
        let mut global_kappa2_sum = [0.0_f64];
        comm.tree_sum_f64(&[local_kappa2_sum], &mut global_kappa2_sum)?;
        let has_kappa = global_kappa2_sum[0] > 0.0;

        let sum_alpha_z_sm = cavity::build_ionic_pseudo_charge_g(
            &self.solver.engine.cell,
            self.solver.engine.dims,
            self.solver.engine.fft.local_positions(),
            ions,
            &mut self.eval.rho_ion_g,
        )?;
        if let Some(dg) = density_g {
            self.eval.n_g.copy_from_slice(dg);
        } else {
            self.solver
                .engine
                .forward_real(comm, density_planes, &mut self.eval.n_g)?;
        }
        for ((rho_sol, &rho_ion), &n) in self
            .eval
            .rho_sol_g
            .iter_mut()
            .zip(&self.eval.rho_ion_g)
            .zip(&self.eval.n_g)
        {
            *rho_sol = rho_ion - n;
        }

        let q_ion: f64 = ions.iter().map(|ion| ion.charge).sum();
        let background_correction_ha = (q_ion - n_el) * sum_alpha_z_sm;

        // Exact vacuum potential with the identical skew-adjoint L_0 operator.
        for ((phi_vac, &g2), &rho_sol) in self
            .eval
            .phi_vac_g
            .iter_mut()
            .zip(&self.solver.engine.g2)
            .zip(&self.eval.rho_sol_g)
        {
            *phi_vac = if g2 == 0.0 {
                Complex64::default()
            } else {
                (FOUR_PI / g2) * rho_sol
            };
        }

        // Fast exact vacuum path when all physical solvation couplings are off and neutral.
        if params.is_vacuum() && (q_ion - n_el).abs() < 1.0e-12 {
            self.eval.phi_pol_g.fill(Complex64::default());
            self.solver.engine.inverse_to_real(
                comm,
                &self.eval.phi_vac_g,
                &mut self.eval.phi_solv,
            )?;
            self.eval.phi_pol.fill(0.0);
            return Ok(SolvationEvaluation {
                electrostatic_ha: 0.0,
                background_correction_ha: 0.0,
                cavitation_ha: 0.0,
                total_ha: 0.0,
                v_es_planes: vec![0.0; pp],
                v_diel_planes: vec![0.0; pp],
                v_scr_planes: vec![0.0; pp],
                v_cav_planes: vec![0.0; pp],
                v_total_planes: vec![0.0; pp],
                phi_pol_planes: vec![0.0; pp],
                phi_solv_planes: self.eval.phi_solv.clone(),
                ionic_forces_ha_bohr: vec![[0.0; 3]; ions.len()],
                report: SolveReport {
                    kind: config.kind,
                    iterations: 0,
                    converged: true,
                    residual_norm: 0.0,
                    relative_residual: 0.0,
                    fft_count: 0,
                    mirror_paired: self.mirror_paired(),
                    warm_started: false,
                },
            });
        }

        // Polarization RHS:
        //   b_pol = 4*pi * rho_sol - A * phi_vac
        //         = 4*pi * rho_sol|_{ker(L_0)} - [sum_a D_a^T((eps - 1) D_a) + kappa2] phi_vac
        // Computing via (eps - 1) avoids O(Z / r_c) core-peak cancellation inside the cavity.
        for (dst, &eps) in self.solver.work_r0.iter_mut().zip(&self.eval.eps) {
            *dst = eps - 1.0;
        }
        let _ = self.solver.engine.apply_operator(
            comm,
            &self.solver.work_r0,
            &self.eval.kappa2,
            has_kappa,
            &self.eval.phi_vac_g,
            &mut self.eval.b_pol_g,
        )?;
        for v in &mut self.eval.b_pol_g {
            *v = -*v;
        }
        for &k in &self.solver.engine.zero_g2_indices {
            self.eval.b_pol_g[k] += FOUR_PI * self.eval.rho_sol_g[k];
        }

        let report = match self.solver.solve_pde(
            comm,
            &self.eval.eps,
            &self.eval.kappa2,
            &self.eval.b_pol_g,
            config,
            &mut self.eval.phi_pol_g,
        ) {
            Ok(rep) => rep,
            Err(SolvationError::Breakdown(_)) if config.kind == SolverKind::FisicaroPcg => {
                // Automatic guaranteed-SPD fallback when Fisicaro's product-rule
                // operator becomes indefinite on an under-resolved cavity wall.
                let fallback_cfg = SolverConfig {
                    kind: SolverKind::DirectPcg,
                    ..*config
                };
                self.solver.solve_pde(
                    comm,
                    &self.eval.eps,
                    &self.eval.kappa2,
                    &self.eval.b_pol_g,
                    &fallback_cfg,
                    &mut self.eval.phi_pol_g,
                )?
            }
            Err(err) => return Err(err),
        };

        for ((phi_solv, &phi_vac), &phi_pol) in self
            .eval
            .phi_solv_g
            .iter_mut()
            .zip(&self.eval.phi_vac_g)
            .zip(&self.eval.phi_pol_g)
        {
            *phi_solv = phi_vac + phi_pol;
        }

        // One np = 3 inverse batch yields phi_solv, |D phi_solv|^2, and phi_pol on owned planes.
        self.solver.engine.reconstruct_potentials_and_grad_sq(
            comm,
            &self.eval.phi_solv_g,
            &self.eval.phi_pol_g,
            &mut self.eval.phi_solv,
            &mut self.eval.grad_phi_sq,
            &mut self.eval.phi_pol,
        )?;

        // Dual electrostatic solvation energy: 0.5 * Omega * sum_G Re(conj(rho_sol) * phi_pol).
        let mut es_loc = 0.0_f64;
        for (&rho, &phi) in self.eval.rho_sol_g.iter().zip(&self.eval.phi_pol_g) {
            es_loc += rho.re * phi.re + rho.im * phi.im;
        }
        let mut es_glob = [0.0_f64];
        comm.tree_sum_f64(&[es_loc], &mut es_glob)?;
        let electrostatic_ha = 0.5 * self.solver.engine.volume * es_glob[0];

        // Unconstrained variational derivative of E_bg,corr = (Q_ion - N_e) * sum_I alpha_Z,I^sm
        // with respect to n(r) is -sum_I alpha_Z,I^sm for both kappa = 0 and kappa > 0.
        // At fixed N_e (kappa = 0), this uniform constant shifts KS eigenvalues without
        // altering orbitals or density while making dE_total/dn and Bregman remainders exact.
        let bg_pot_shift = sum_alpha_z_sm;
        let mut v_es_planes = vec![0.0; pp];
        let mut v_diel_planes = vec![0.0; pp];
        let mut v_scr_planes = vec![0.0; pp];
        for ((((v_es, v_diel), v_scr), (&phi_pol, &phi_solv)), ((&deps, &grad2), &dkap)) in
            v_es_planes
                .iter_mut()
                .zip(&mut v_diel_planes)
                .zip(&mut v_scr_planes)
                .zip(self.eval.phi_pol.iter().zip(&self.eval.phi_solv))
                .zip(
                    self.eval
                        .deps_dn
                        .iter()
                        .zip(&self.eval.grad_phi_sq)
                        .zip(&self.eval.dkappa2_dn),
                )
        {
            *v_es = dft_units::ElectrostaticPotentialHaPerE(phi_pol + bg_pot_shift)
                .to_electron_energy()
                .0;
            *v_diel = -INV_EIGHT_PI * deps * grad2;
            *v_scr = -INV_EIGHT_PI * dkap * phi_solv * phi_solv;
        }

        // Discrete cavitation (+ volume) energy and exact variational potential.
        let mut v_cav_planes = vec![0.0; pp];
        let vol_e_loc = if params.beta_ha_bohr3 != 0.0 {
            params.beta_ha_bohr3 * self.eval.s.iter().map(|&sv| 1.0 - sv).sum::<f64>()
        } else {
            0.0
        };
        let surf_e_loc = if params.tau_ha_bohr2 > 0.0 {
            self.solver
                .engine
                .forward_real(comm, &self.eval.s, &mut self.solver.aux_g)?;
            let area_loc = self.solver.engine.cavitation_surface_grad_and_div(
                comm,
                &self.solver.aux_g,
                params.cavity_grad_reg,
                &mut self.solver.w_g,
                &mut self.solver.work_r0,
            )?;
            for ((dst, &ds), &div_n) in v_cav_planes
                .iter_mut()
                .zip(&self.eval.ds_dn)
                .zip(&self.solver.work_r0)
            {
                *dst = -params.tau_ha_bohr2 * ds * div_n - params.beta_ha_bohr3 * ds;
            }
            params.tau_ha_bohr2 * area_loc
        } else {
            if params.beta_ha_bohr3 != 0.0 {
                for (dst, &ds) in v_cav_planes.iter_mut().zip(&self.eval.ds_dn) {
                    *dst = -params.beta_ha_bohr3 * ds;
                }
            }
            0.0
        };
        let mut cav_glob = [0.0_f64];
        comm.tree_sum_f64(&[surf_e_loc + vol_e_loc], &mut cav_glob)?;
        let cavitation_ha = self.solver.engine.dv * cav_glob[0];

        let mut v_total_planes = vec![0.0; pp];
        for ((((dst, &v_es), &v_diel), &v_scr), &v_cav) in v_total_planes
            .iter_mut()
            .zip(&v_es_planes)
            .zip(&v_diel_planes)
            .zip(&v_scr_planes)
            .zip(&v_cav_planes)
        {
            *dst = v_es + v_diel + v_scr + v_cav;
        }

        let ionic_forces_ha_bohr = cavity::ionic_solvation_forces(
            comm,
            &self.solver.engine.cell,
            self.solver.engine.dims,
            self.solver.engine.fft.local_positions(),
            ions,
            &self.eval.phi_pol_g,
        )?;

        let total_ha = electrostatic_ha + background_correction_ha + cavitation_ha;

        Ok(SolvationEvaluation {
            electrostatic_ha,
            background_correction_ha,
            cavitation_ha,
            total_ha,
            v_es_planes,
            v_diel_planes,
            v_scr_planes,
            v_cav_planes,
            v_total_planes,
            phi_pol_planes: self.eval.phi_pol.clone(),
            phi_solv_planes: self.eval.phi_solv.clone(),
            ionic_forces_ha_bohr,
            report,
        })
    }
}
