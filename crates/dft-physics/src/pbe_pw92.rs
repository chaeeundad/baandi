//! Unpolarized correlation from PW92 (PRB 45, 13244, Eq. 10/Table I)
//! and PBE96 (PRL 77, 3865, Eqs. 7–8). Hartree atomic units.
//! Only PW92 A is selectable; the exact PBE gamma is independent of it.
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Pw92Prefactor {
    #[default]
    Modified,
    Original,
}
impl Pw92Prefactor {
    pub fn name(self) -> &'static str {
        match self {
            Self::Modified => "modified",
            Self::Original => "original",
        }
    }
    pub fn value(self) -> f64 {
        match self {
            Self::Modified => 0.0310907,
            Self::Original => 0.031091,
        }
    }
}

/// Correlation energy per electron (Ha), d(n*exc)/dn (Ha), and
/// d(n*exc)/d sigma (Ha Bohr^5). The caller owns density screening.
#[derive(Clone, Copy, Debug)]
pub struct Correlation {
    pub exc: f64,
    pub vrho: f64,
    pub vsigma: f64,
}

/// Analytic chain rule at positive n and nonnegative sigma = |grad n|^2.
/// No finite differences or allocation are used. Negative sigma is allowed
/// locally for derivative oracles around zero, provided all logs stay defined.
pub fn correlation(n: f64, sigma: f64, prefactor: Pw92Prefactor) -> Correlation {
    let a = prefactor.value();
    let r = (3.0 / (4.0 * PI * n)).cbrt();
    let root = r.sqrt();
    let q = 2.0 * a * (7.5957 * root + 3.5876 * r + 1.6382 * r * root + 0.49294 * r * r);
    let dq = 2.0 * a * (7.5957 / (2.0 * root) + 3.5876 + 1.5 * 1.6382 * root + 2.0 * 0.49294 * r);
    let log = (1.0 / q).ln_1p();
    let ec = -2.0 * a * (1.0 + 0.21370 * r) * log;
    let dec_dr = -2.0 * a * (0.21370 * log - (1.0 + 0.21370 * r) * dq / (q * (q + 1.0)));
    let dec_dn = dec_dr * (-r / (3.0 * n));
    let kf = (3.0 * PI * PI * n).cbrt();
    let dt_dsigma = PI / (16.0 * kf * n * n);
    let t = sigma * dt_dsigma; // squared PBE reduced gradient
    let gamma = (1.0 - 2.0_f64.ln()) / (PI * PI);
    let beta_gamma = 0.06672455060314922 / gamma;
    let expm1 = (-ec / gamma).exp_m1();
    let aux = beta_gamma / expm1;
    let daux_dn = aux * (1.0 + 1.0 / expm1) * dec_dn / gamma;
    let u = aux * t;
    let denominator = 1.0 + u + u * u;
    let ratio = (1.0 + u) / denominator;
    let dratio_du = -u * (2.0 + u) / (denominator * denominator);
    let f = t * ratio;
    let df_dt = ratio + u * dratio_du;
    let df_daux = t * t * dratio_du;
    let h = gamma * (beta_gamma * f).ln_1p();
    let dh_df = gamma * beta_gamma / (1.0 + beta_gamma * f);
    let dh_dn = dh_df * (df_dt * (-7.0 * t / (3.0 * n)) + df_daux * daux_dn);
    Correlation {
        exc: ec + h,
        vrho: ec + h + n * (dec_dn + dh_dn),
        vsigma: n * dh_df * df_dt * dt_dsigma,
    }
}
