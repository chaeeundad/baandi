"""Experimental "A3" optimizer, re-expressed for the DFT harness.

A3 (exp-tabpfn/tools/icl/run_precon.py mode "A3", base metric from
exp-tabpfn/tools/icl/atomistic.py) is L-BFGS whose initial inverse Hessian is the
exact inverse of a FIXED bond-spring model Hessian built once from the starting
geometry, with Armijo backtracking. Units here are ASE's (eV, Angstrom).

Model Hessian B_p (free Cartesian coordinates only):
- bonds: pairs closer than 1.3 (r_cov,i + r_cov,j); 1.8 for a molecular element
  (H, C, N, O) bonded to a non-molecular one; minimum-image, periodic;
- stiffness k = 5 exp(-2.5 (d/r0 - 1)) eV/A^2, r0 = r_cov,i + r_cov,j;
- per bond the 3x3 block k (0.3 I + 0.7 n n^T) assembled as a graph Laplacian;
- plus 0.05 I eV/A^2 so the metric is SPD on the free subspace (removes the
  translation null space and isolated atoms).
L-BFGS: memory 10, two-loop recursion with h0 = B_p^{-1}; a non-descent
direction falls back to -B_p^{-1} g; curvature pairs with s.y <= 1e-12 |s||y| are
skipped. Line search: Armijo f(x + t p) <= f + 1e-4 t g.p, t halves, at most 25
trials, each trial a full DFT evaluation (counted). Step cap: the experiment
capped the 2-norm of the whole step (0.3 A on EMT); here the cap is the SAME
per-atom displacement cap that every other method of the benchmark receives
(deviation recorded in the results as ``cap_mode = "per_atom"``).
"""
import numpy as np
from ase.data import covalent_radii
from ase.neighborlist import neighbor_list

MOL = {1, 6, 7, 8}


def base_metric(atoms, free_mask):
    zs = sorted(set(atoms.numbers.tolist()))
    cut = {(p, q): (1.8 if ((p in MOL) != (q in MOL)) else 1.3) * (covalent_radii[p] + covalent_radii[q])
           for p in zs for q in zs}
    i, j, D = neighbor_list("ijD", atoms, cut)
    sel = i < j
    i, j, D = i[sel], j[sel], D[sel]
    d = np.linalg.norm(D, axis=1)
    z = atoms.numbers
    r0 = covalent_radii[z[i]] + covalent_radii[z[j]]
    k = 5.0 * np.exp(-2.5 * (d / r0 - 1.0))
    n = D / d[:, None]
    N = len(atoms)
    H = np.zeros((3 * N, 3 * N))
    for p, q, kk, nn in zip(i, j, k, n):
        blk = kk * (0.3 * np.eye(3) + 0.7 * np.outer(nn, nn))
        H[3 * p:3 * p + 3, 3 * p:3 * p + 3] += blk
        H[3 * q:3 * q + 3, 3 * q:3 * q + 3] += blk
        H[3 * p:3 * p + 3, 3 * q:3 * q + 3] -= blk
        H[3 * q:3 * q + 3, 3 * p:3 * p + 3] -= blk
    fm = np.repeat(free_mask, 3)
    H = H[np.ix_(fm, fm)] + 0.05 * np.eye(int(fm.sum()))
    return H, {"bonds": int(len(i)), "k_min_ev_a2": float(k.min()) if len(k) else None,
               "k_max_ev_a2": float(k.max()) if len(k) else None}


def run(atoms, free_mask, cap_a, accepted, memory=10, c1=1e-4, max_trials=25, max_iter=10 ** 6):
    """Minimise; ``accepted`` receives a copy of every accepted iterate (Angstrom).
    Stops by exceptions of the evaluator (Converged / BudgetExhausted) or returns a status."""
    fm = np.repeat(free_mask, 3)
    Bp, info = base_metric(atoms, free_mask)
    w, Q = np.linalg.eigh(Bp)
    with np.errstate(all="ignore"):  # numpy 2.0 + Accelerate emits spurious matmul warnings
        Binv = (Q / w) @ Q.T
    assert np.isfinite(Binv).all() and w.min() > 0
    info["metric_eig_min_ev_a2"], info["metric_eig_max_ev_a2"] = float(w.min()), float(w.max())

    def fg(x):
        pos = atoms.get_positions().ravel().copy()
        pos[fm] = x
        atoms.set_positions(pos.reshape(-1, 3))
        return atoms.get_potential_energy(), -atoms.get_forces().ravel()[fm]

    x = atoms.get_positions().ravel()[fm].copy()
    f, g = fg(x)
    accepted.append(atoms.get_positions().copy())
    S, Y = [], []
    for _ in range(max_iter):
        q = g.copy()
        alphas = []
        for s, y in zip(reversed(S), reversed(Y)):
            a = (s @ q) / (y @ s)
            alphas.append(a)
            q -= a * y
        z = Binv @ q
        for (s, y), a in zip(zip(S, Y), reversed(alphas)):
            b = (y @ z) / (y @ s)
            z += (a - b) * s
        p = -z
        if g @ p >= 0:
            p = -(Binv @ g)
        per_atom = np.zeros(len(atoms))
        per_atom[free_mask] = np.linalg.norm(p.reshape(-1, 3), axis=1)
        if per_atom.max() > cap_a:
            p *= cap_a / per_atom.max()
        t, ok = 1.0, False
        for _ in range(max_trials):
            fn, gn = fg(x + t * p)
            if fn <= f + c1 * t * (g @ p):
                ok = True
                break
            t *= 0.5
        if not ok:
            return "line_search_failed", info
        s, y = t * p, gn - g
        if s @ y > 1e-12 * np.linalg.norm(s) * np.linalg.norm(y):
            S.append(s)
            Y.append(y)
            if len(S) > memory:
                S.pop(0)
                Y.pop(0)
        x, f, g = x + s, fn, gn
        accepted.append(atoms.get_positions().copy())
    return "max_iter", info
