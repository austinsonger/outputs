#!/usr/bin/env python3
"""
Numerical inscribed-square finder for Jordan curves.
Toeplitz / Square Peg Problem experimental toolkit.

Method (diagonal formulation):
  A square inscribed in curve G is exactly a pair of parameters (t1, t2) such that,
  with p1 = G(t1), p2 = G(t2), center c = (p1+p2)/2, and w = rot90((p2-p1)/2),
  both q1 = c + w and q2 = c - w also lie on G.
  (p1, q1, p2, q2) then has equal perpendicular diagonals sharing a midpoint -> square.

Search:
  objective(t1,t2) = dist(q1, G)^2 + dist(q2, G)^2
  1. vectorized grid scan over the (t1,t2) torus (KD-tree distance queries)
  2. local minima below threshold -> Nelder-Mead refinement on the exact curve
  3. accept if corners verifiably on curve (segment-projected distance) and side
     above a degeneracy floor; dedupe by corner sets

Usage:
  python3 inscribed_squares.py            # gallery of 4 curves + roughness sweep
  python3 inscribed_squares.py gallery    # gallery only
  python3 inscribed_squares.py sweep      # roughness sweep only
"""

import sys
import time
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAU = 2.0 * np.pi
RNG_SEED = 7


# ----------------------------------------------------------------------------
# Curves: vectorized maps t in [0,1) -> points (n,2)
# ----------------------------------------------------------------------------

class Curve:
    def __init__(self, fn, name, n_dense=24000):
        self.fn = fn
        self.name = name
        ts = np.linspace(0.0, 1.0, n_dense, endpoint=False)
        self.dense = fn(ts)                      # (n,2) polyline sample
        self.n_dense = n_dense
        self.tree = cKDTree(self.dense)
        mins, maxs = self.dense.min(0), self.dense.max(0)
        self.diam = float(np.linalg.norm(maxs - mins))
        seg = np.linalg.norm(np.diff(np.vstack([self.dense, self.dense[:1]]), axis=0), axis=1)
        self.spacing = float(seg.mean())        # KD distance resolution floor

    def point(self, t):
        return self.fn(np.atleast_1d(np.asarray(t, dtype=float)) % 1.0)

    def dist_fast(self, pts):
        """KD-tree distance to dense sample cloud (fast, slight overestimate)."""
        d, _ = self.tree.query(np.atleast_2d(pts))
        return d

    def dist_exact(self, pts):
        """Point-to-polyline distance with segment projection near KD hits.
        Vectorized over a batch of points. Returns (m,) distances."""
        p = np.atleast_2d(np.asarray(pts, dtype=float))          # (m,2)
        _, idx = self.tree.query(p, k=6)                          # (m,6)
        idx = np.atleast_2d(idx)
        n = self.n_dense
        seg = np.concatenate([(idx - 1) % n, idx % n], axis=1)    # (m,12) seg starts
        a = self.dense[seg]                                       # (m,12,2)
        b = self.dense[(seg + 1) % n]
        ab = b - a
        L2 = np.einsum("mij,mij->mi", ab, ab)
        ap = p[:, None, :] - a
        s = np.clip(np.einsum("mij,mij->mi", ap, ab) / np.maximum(L2, 1e-300), 0.0, 1.0)
        proj = a + s[..., None] * ab
        d = np.linalg.norm(p[:, None, :] - proj, axis=2)          # (m,12)
        return d.min(axis=1)


def ellipse(a=1.6, b=1.0):
    def fn(t):
        th = TAU * t
        return np.column_stack([a * np.cos(th), b * np.sin(th)])
    return Curve(fn, f"Ellipse a={a}, b={b}")


def fourier_blob(seed=RNG_SEED, n_modes=6, amp=0.16):
    """Smooth generic Jordan curve: radial Fourier series with fast decay."""
    rng = np.random.default_rng(seed)
    ph = rng.uniform(0, TAU, n_modes)
    cf = amp * rng.uniform(0.5, 1.0, n_modes) / (np.arange(1, n_modes + 1) ** 1.3)

    def fn(t):
        th = TAU * t
        r = 1.0 + sum(cf[k] * np.cos((k + 1) * th + ph[k]) for k in range(n_modes))
        return np.column_stack([r * np.cos(th), r * np.sin(th)])
    return Curve(fn, "Smooth random blob")


def star_polygon(points=5, r_out=1.0, r_in=0.45):
    """Non-convex polygon (piecewise linear, corners everywhere it matters)."""
    k = 2 * points
    ang = TAU * np.arange(k) / k - TAU / 4
    rad = np.where(np.arange(k) % 2 == 0, r_out, r_in)
    vx = np.column_stack([rad * np.cos(ang), rad * np.sin(ang)])
    closed = np.vstack([vx, vx[:1]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = cum[-1]

    def fn(t):
        s = (np.asarray(t) % 1.0) * total
        i = np.clip(np.searchsorted(cum, s, side="right") - 1, 0, k - 1)
        frac = (s - cum[i]) / seg[i]
        return closed[i] + frac[:, None] * (closed[i + 1] - closed[i])
    return Curve(fn, f"{points}-point star polygon")


def weierstrass_curve(h=0.5, K=7, base=3, amp=0.5, seed=RNG_SEED, n_dense=60000):
    """
    Radial Weierstrass-type curve: r(th) = 1 + amp * norm * sum lam^k cos(base^k th + ph_k),
    lam = base^(-h). Truncated at K terms; as h -> 0 the limit curve loses all
    smoothness (nowhere differentiable at h < 1 in the K -> inf limit).
    Radial form guarantees a simple (Jordan) curve.
    """
    rng = np.random.default_rng(seed)
    ph = rng.uniform(0, TAU, K)
    lam = float(base) ** (-h)
    w = lam ** np.arange(K)
    norm = 1.0 / w.sum()

    def fn(t):
        th = TAU * np.asarray(t)
        r = 1.0 + amp * norm * sum(w[k] * np.cos((base ** k) * th + ph[k]) for k in range(K))
        return np.column_stack([r * np.cos(th), r * np.sin(th)])
    return Curve(fn, f"Weierstrass h={h}", n_dense=n_dense)


# ----------------------------------------------------------------------------
# Square finder
# ----------------------------------------------------------------------------

def square_corners(curve, t1, t2):
    """Corners (p1, q1, p2, q2) of the candidate square with diagonal p1-p2."""
    p1 = curve.point(t1)[0]
    p2 = curve.point(t2)[0]
    c = 0.5 * (p1 + p2)
    v = 0.5 * (p2 - p1)
    w = np.array([-v[1], v[0]])
    return np.array([p1, c + w, p2, c - w])


def find_squares(curve, n_grid=460, dt_min=0.004, side_floor_rel=0.01,
                 accept_tol_rel=1e-7, max_refine=400, verbose=False):
    """
    Returns list of dicts: corners (4,2), side, max corner distance to curve.
    side_floor_rel: reject squares smaller than this fraction of curve diameter
                    (the numerical degeneracy floor -- cf. the limiting-argument leak).
    """
    t0 = time.time()
    diam = curve.diam
    side_floor = side_floor_rel * diam
    accept_tol = accept_tol_rel * diam

    # --- vectorized grid scan over the (t1,t2) torus ---
    g = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    T1, T2 = np.meshgrid(g, g, indexing="ij")
    P1 = curve.point(T1.ravel())
    P2 = curve.point(T2.ravel())
    C = 0.5 * (P1 + P2)
    V = 0.5 * (P2 - P1)
    W = np.column_stack([-V[:, 1], V[:, 0]])
    d1 = curve.dist_fast(C + W)
    d2 = curve.dist_fast(C - W)
    OBJ = (d1 ** 2 + d2 ** 2).reshape(n_grid, n_grid)

    # mask the degenerate diagonal band t1 ~ t2 (p1=p2 gives a fake zero)
    dt = np.abs(T1 - T2)
    dt = np.minimum(dt, 1.0 - dt)
    OBJ[dt < dt_min] = np.inf
    # also mask geometrically tiny diagonals
    OBJ[(np.linalg.norm(2 * V, axis=1).reshape(n_grid, n_grid)) < side_floor] = np.inf

    # --- local minima of the grid (8-neighbor, torus wraparound) ---
    m = OBJ
    is_min = np.ones_like(m, dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            is_min &= m <= np.roll(np.roll(m, di, axis=0), dj, axis=1)
    thresh = (0.03 * diam) ** 2
    cand = np.argwhere(is_min & (m < thresh))
    cand = cand[np.argsort(m[cand[:, 0], cand[:, 1]])][:max_refine]

    # --- refine each candidate on the exact curve ---
    def objective(x):
        cs = square_corners(curve, x[0], x[1])
        return float(curve.dist_fast(cs[[1, 3]]) ** 2 @ np.ones(2))

    def objective_exact(x):
        c = square_corners(curve, x[0], x[1])
        d = curve.dist_exact(c[[1, 3]])
        return float(d @ d)

    def dedupe(items):
        uniq = []
        for f in items:
            dup = False
            for u in uniq:
                dmat = np.linalg.norm(f["corners"][:, None, :] - u["corners"][None, :, :],
                                      axis=2)
                if dmat.min(axis=1).max() < 0.02 * diam:
                    dup = True
                    if f["err"] < u["err"]:
                        u.update(f)
                    break
            if not dup:
                uniq.append(f)
        return uniq

    # --- stage 1: fast refinement of every grid candidate ---
    found = []
    for i, j in cand:
        res = minimize(objective, x0=[g[i], g[j]], method="Nelder-Mead",
                       options=dict(xatol=1e-10, fatol=1e-16, maxiter=400))
        cs = square_corners(curve, res.x[0], res.x[1])
        side = np.linalg.norm(cs[1] - cs[0])
        if side < side_floor:
            continue
        # pre-filter at the KD resolution floor; the exact polish decides for real
        worst = float(curve.dist_exact(cs[[1, 3]]).max())
        if worst > max(10.0 * curve.spacing, 4.0 * accept_tol):
            continue
        found.append(dict(corners=cs, side=float(side), err=worst, x=res.x))

    # --- stage 2: dedupe, then polish survivors against segment-exact distance
    # (KD point-cloud distance has a resolution floor ~ sample spacing) ---
    uniq = []
    for f in dedupe(found):
        res = minimize(objective_exact, x0=f["x"], method="Nelder-Mead",
                       options=dict(xatol=1e-12, fatol=1e-20, maxiter=250))
        cs = square_corners(curve, res.x[0], res.x[1])
        side = float(np.linalg.norm(cs[1] - cs[0]))
        worst = float(curve.dist_exact(cs[[1, 3]]).max())
        # scale-aware acceptance: small squares must be proportionally accurate
        if side < side_floor or worst > min(accept_tol, 0.01 * side):
            continue
        uniq.append(dict(corners=cs, side=side, err=worst, x=res.x.tolist()))
    uniq = dedupe(uniq)
    uniq.sort(key=lambda f: -f["side"])

    if verbose:
        print(f"  [{curve.name}] grid {n_grid}^2, {len(cand)} candidates -> "
              f"{len(uniq)} squares in {time.time()-t0:.1f}s")
    return uniq


# ----------------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------------

def draw_panel(ax, curve, squares, title_extra=""):
    ax.plot(*curve.dense.T, "k-", lw=1.0, zorder=1)
    cmap = plt.get_cmap("turbo")
    for k, sq in enumerate(squares):
        cs = np.vstack([sq["corners"], sq["corners"][:1]])
        col = cmap(0.1 + 0.8 * k / max(1, len(squares) - 1)) if len(squares) > 1 else cmap(0.25)
        ax.plot(*cs.T, "-", lw=1.6, color=col, zorder=2)
        ax.plot(*sq["corners"].T, "o", ms=3.5, color=col, zorder=3)
    ax.set_title(f"{curve.name}\n{len(squares)} inscribed square(s){title_extra}", fontsize=10)
    ax.set_aspect("equal")
    ax.axis("off")


def run_gallery(outdir="."):
    curves = [ellipse(), fourier_blob(), star_polygon(), weierstrass_curve(h=0.45)]
    fig, axes = plt.subplots(2, 2, figsize=(11, 11))
    results = {}
    for ax, curve in zip(axes.ravel(), curves):
        squares = find_squares(curve, verbose=True)
        results[curve.name] = squares
        draw_panel(ax, curve, squares)
    fig.suptitle("Inscribed squares found numerically (diagonal-pair search)",
                 fontsize=13, y=0.98)
    fig.tight_layout()
    path = f"{outdir}/squares_gallery.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"wrote {path}")

    # --- verification: ellipse analytic square ---
    a, b = 1.6, 1.0
    s = a * b / np.hypot(a, b)
    expected = np.array([[s, s], [-s, s], [-s, -s], [s, -s]])
    ell = results[f"Ellipse a={a}, b={b}"]
    ok = False
    for sq in ell:
        dmat = np.linalg.norm(sq["corners"][:, None, :] - expected[None, :, :], axis=2)
        if dmat.min(axis=1).max() < 1e-3:
            ok = True
    print(f"VERIFY ellipse: expected corners (+-{s:.6f}, +-{s:.6f}); "
          f"found {len(ell)} square(s); analytic match: {ok}")
    for name, sqs in results.items():
        sides = ", ".join(f"{q['side']:.3f}" for q in sqs)
        parity = "odd" if len(sqs) % 2 == 1 else "EVEN"
        print(f"  {name}: {len(sqs)} squares ({parity}) sides=[{sides}] "
              f"max_err={max((q['err'] for q in sqs), default=0):.2e}")
    return results


SWEEP_HS = [1.0, 0.8, 0.65, 0.5, 0.4, 0.3]


def sweep_part(hs, out_json):
    """Compute squares for a subset of roughness values, checkpoint to JSON."""
    import json
    rows = []
    for h in hs:
        curve = weierstrass_curve(h=h)
        squares = find_squares(curve, n_grid=520, side_floor_rel=0.008, verbose=True)
        rows.append(dict(h=h, diam=curve.diam,
                         squares=[dict(corners=s["corners"].tolist(),
                                       side=s["side"], err=s["err"]) for s in squares]))
    with open(out_json, "w") as f:
        json.dump(rows, f)
    print(f"wrote {out_json}")


def sweep_plot(part_jsons, outdir="."):
    import json
    rows = []
    for pj in part_jsons:
        with open(pj) as f:
            rows.extend(json.load(f))
    rows.sort(key=lambda r: -r["h"])

    fig = plt.figure(figsize=(13, 7.5))
    stats = []
    for k, r in enumerate(rows):
        curve = weierstrass_curve(h=r["h"])
        squares = [dict(corners=np.array(s["corners"]), side=s["side"], err=s["err"])
                   for s in r["squares"]]
        ax = fig.add_subplot(2, len(rows), k + 1)
        draw_panel(ax, curve, squares)
        ax.set_title(f"h={r['h']}   n={len(squares)}", fontsize=9)
        sides = np.array([s["side"] for s in squares])
        stats.append(dict(h=r["h"], count=len(squares),
                          min_side=sides.min() / r["diam"] if len(sides) else np.nan,
                          med_side=np.median(sides) / r["diam"] if len(sides) else np.nan,
                          max_side=sides.max() / r["diam"] if len(sides) else np.nan))

    ax2 = fig.add_subplot(2, 1, 2)
    H = [s["h"] for s in stats]
    ax2.semilogy(H, [s["min_side"] for s in stats], "o-",
                 label="smallest square / diameter")
    ax2.semilogy(H, [s["med_side"] for s in stats], "s--",
                 label="median square / diameter")
    ax2b = ax2.twinx()
    ax2b.plot(H, [s["count"] for s in stats], "^:", color="gray", label="# squares found")
    ax2b.set_ylabel("# squares found", color="gray")
    ax2.set_xlabel("roughness exponent h   (smooth  <--  1.0 ... 0.3  -->  rough)")
    ax2.set_ylabel("square side / curve diameter (log)")
    ax2.invert_xaxis()
    ax2.legend(loc="lower left", fontsize=9)
    ax2.set_title("Degeneration probe: smallest inscribed square shrinks as roughness grows",
                  fontsize=11)
    fig.tight_layout()
    path = f"{outdir}/roughness_sweep.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"wrote {path}")
    print("h, count, min_side/diam, med_side/diam, max_side/diam")
    for s in stats:
        print(f"  {s['h']:.2f}  {s['count']:3d}   {s['min_side']:.4f}   "
              f"{s['med_side']:.4f}   {s['max_side']:.4f}")
    return stats


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "gallery"
    np.set_printoptions(precision=5, suppress=True)
    if mode == "gallery":
        outdir = sys.argv[2] if len(sys.argv) > 2 else "."
        run_gallery(outdir)
    elif mode == "sweep_part":
        hs = [float(x) for x in sys.argv[2].split(",")]
        sweep_part(hs, sys.argv[3])
    elif mode == "sweep_plot":
        sweep_plot(sys.argv[2].split(","), sys.argv[3] if len(sys.argv) > 3 else ".")
