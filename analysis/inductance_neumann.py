#!/usr/bin/env python3
"""Numerical Neumann (partial-inductance) calculation of L_A, L_D, M, k
for the two interleaved planar windings, using the real Gerber polylines.

Run:  python3 -I analysis/inductance_neumann.py
Needs work/theory/geometry.json (from gerber_geometry.py).

Method (magnetostatic, uniform current, air core, mu_r = 1 everywhere):
  L_loop = sum_i sum_j Mp_ij,   Mp_ij = mu0/(4 pi) * Int_i Int_j dl_i . dl_j / r
over straight segments i, j of the closed-loop current path (each winding is
treated from its outer pad to its inner pad; the terminal gap is not closed,
i.e. partial inductance of the open path, which is what a 2-terminal LCR meter
sees apart from the fixture loop).

Pair evaluation:
  (a) i == j, or i, j collinear pieces of the same bar: exact Campbell/Grover
      formula for parallel filaments, lateral distance = self-GMD of the
      w x t rectangle (Maxwell's exact rectangle GMD).  Reproduces Grover's
      straight-bar formula L = (mu0 l / 2 pi)[ln(2l/(w+t)) + 1/2 + ...].
  (b) exactly parallel (or anti-parallel) segments: exact Campbell formula;
      lateral distance replaced by the mutual GMD of the two w x t rectangles
      if they are coplanar and closer than 3 mm (numerical 4-D Gauss GMD).
  (c) non-parallel segments sharing an end point: exact Grover formula for
      two filaments meeting at a point at an angle.
  (d) all other pairs: composite Gauss-Legendre over segment i (n_g nodes per
      piece, piece length <= alpha * min-distance(i, j)) of the closed-form
      line integral over segment j:  Int_j ds/r = 2 atanh(L_j/(r_a + r_b)).
      Symmetrised: (M_ij + M_ji)/2.
Convergence is studied over alpha, n_g and the polyline simplification
(Douglas-Peucker tolerance; raw = every Gerber vertex).
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.special import ellipe, ellipk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gerber_geometry as gg  # noqa: E402

MU0 = 4e-7 * math.pi
IN = 0.0254
WORK = gg.WORK
OUT_JSON = os.path.join(WORK, "neumann_results.json")


# --------------------------------------------------------------------------
# GMD helpers
# --------------------------------------------------------------------------
def gmd_rect_self(a, b):
    """Maxwell's exact GMD of a rectangle a x b from itself (Grover 1946)."""
    if b == 0:
        return a * math.exp(-1.5)
    lng = (0.5 * math.log(a * a + b * b)
           - (b * b / (12 * a * a)) * math.log(1 + a * a / (b * b))
           - (a * a / (12 * b * b)) * math.log(1 + b * b / (a * a))
           + (2 * b / (3 * a)) * math.atan(a / b)
           + (2 * a / (3 * b)) * math.atan(b / a)
           - 25.0 / 12.0)
    return math.exp(lng)


_gl12 = np.polynomial.legendre.leggauss(12)
_gmd_cache = {}


def gmd_rect_pair(d, w, t, n=12):
    """GMD of two equal w x t rectangles side by side (centre distance d, same plane)."""
    key = (round(d, 9), w, t)
    if key in _gmd_cache:
        return _gmd_cache[key]
    x, wx = np.polynomial.legendre.leggauss(n)
    xs = 0.5 * w * x
    ys = 0.5 * t * x
    W = np.einsum("a,b,c,d->abcd", wx, wx, wx, wx) / 16.0
    X1, Y1, X2, Y2 = np.meshgrid(xs, ys, xs, ys, indexing="ij")
    r = np.sqrt((X1 - X2 + d) ** 2 + (Y1 - Y2) ** 2)
    g = math.exp(float(np.sum(W * np.log(r))))
    _gmd_cache[key] = g
    return g


# --------------------------------------------------------------------------
# filament formulas
# --------------------------------------------------------------------------
def _H2(u, d):
    return u * np.arcsinh(u / d) - np.sqrt(u * u + d * d)


def campbell(a1, b1, a2, b2, d):
    """Int_{a1}^{b1} Int_{a2}^{b2} dx dy / sqrt((x-y)^2 + d^2)  (parallel filaments)."""
    return _H2(b1 - a2, d) - _H2(a1 - a2, d) - _H2(b1 - b2, d) + _H2(a1 - b2, d)


def meeting(l, m, cos_e):
    """Int_0^l Int_0^m dx dy / sqrt(x^2 + y^2 - 2 x y cos_e) (filaments from a common point)."""
    R = np.sqrt(np.maximum(l * l + m * m - 2 * l * m * cos_e, 0.0))
    return 2.0 * (l * np.arctanh(m / (l + R)) + m * np.arctanh(l / (m + R)))


def inner_line_integral(P, A, B):
    """Int over segment AB of ds / |P - r(s)| for points P (...,3)."""
    L = np.linalg.norm(B - A, axis=-1)
    ra = np.linalg.norm(P - A, axis=-1)
    rb = np.linalg.norm(P - B, axis=-1)
    x = L / (ra + rb)
    x = np.minimum(x, 1 - 1e-15)
    return 2.0 * np.arctanh(x)


def seg_dist3(A0, A1, B0, B1):
    """approximate min distance between 3-D segments (pairs, vectorised):
    min of endpoint-to-segment distances, then lower-bounded by the planar
    crossing case (returns |dz| if xy-projections cross)."""
    def ps(P, S0, S1):
        d = S1 - S0
        L2 = np.sum(d * d, -1)
        t = np.clip(np.sum((P - S0) * d, -1) / np.where(L2 == 0, 1, L2), 0, 1)
        return np.linalg.norm(P - (S0 + t[..., None] * d), axis=-1)
    dm = np.minimum.reduce([ps(A0, B0, B1), ps(A1, B0, B1), ps(B0, A0, A1), ps(B1, A0, A1)])

    def cr(o, p, q):
        return (p[..., 0] - o[..., 0]) * (q[..., 1] - o[..., 1]) - (p[..., 1] - o[..., 1]) * (q[..., 0] - o[..., 0])
    cross_xy = (cr(A0, A1, B0) * cr(A0, A1, B1) < 0) & (cr(B0, B1, A0) * cr(B0, B1, A1) < 0)
    dz = np.abs(0.5 * (A0[..., 2] + A1[..., 2]) - 0.5 * (B0[..., 2] + B1[..., 2]))
    return np.where(cross_xy, np.minimum(dm, np.maximum(dz, 1e-12)), dm)


# --------------------------------------------------------------------------
# main solver
# --------------------------------------------------------------------------
class Conductor:
    def __init__(self, P0, P1, wid, bar, w, t, gmd_self_override=None):
        self.P0 = np.asarray(P0, float)
        self.P1 = np.asarray(P1, float)
        self.wid = np.asarray(wid, int)        # winding id per segment
        self.bar = np.asarray(bar, int)        # original straight-bar id
        d = self.P1 - self.P0
        self.L = np.linalg.norm(d, axis=1)
        self.T = d / self.L[:, None]
        self.w = np.asarray(w, float)          # width per segment
        self.t = np.asarray(t, float)          # thickness per segment
        self.gself = np.array([gmd_rect_self(a, b) for a, b in zip(self.w, self.t)])
        if gmd_self_override is not None:
            for k, v in gmd_self_override.items():
                self.gself[k] = v


def _closest_param(P0p, P1p, P0q, P1q):
    """parameter s* (metres along p) of the point of p closest to segment q, and the distance.
    Planar-aware: if xy projections cross, s* = crossing parameter and distance = |dz|."""
    d = P1p - P0p
    Lp = np.linalg.norm(d, axis=1)
    tp = d / Lp[:, None]
    e = P1q - P0q
    Lq2 = np.sum(e * e, 1)

    def proj_on_p(X):
        s = np.clip(np.sum((X - P0p) * tp, 1), 0, Lp)
        return s, np.linalg.norm(X - (P0p + s[:, None] * tp), axis=1)

    def dist_to_q(X):
        t = np.clip(np.sum((X - P0q) * e, 1) / np.where(Lq2 == 0, 1, Lq2), 0, 1)
        return np.linalg.norm(X - (P0q + t[:, None] * e), axis=1)
    s1, d1 = proj_on_p(P0q)
    s2, d2 = proj_on_p(P1q)
    d3 = dist_to_q(P0p)
    d4 = dist_to_q(P1p)
    S = np.stack([s1, s2, np.zeros_like(s1), Lp], 1)
    D = np.stack([d1, d2, d3, d4], 1)
    k = np.argmin(D, 1)
    s_star = S[np.arange(len(k)), k]
    dmin = D[np.arange(len(k)), k]
    # xy crossing (different layers or genuine crossing)
    def cr(o, a, b):
        return (a[:, 0] - o[:, 0]) * (b[:, 1] - o[:, 1]) - (a[:, 1] - o[:, 1]) * (b[:, 0] - o[:, 0])
    c1, c2 = cr(P0p, P1p, P0q), cr(P0p, P1p, P1q)
    c3, c4 = cr(P0q, P1q, P0p), cr(P0q, P1q, P1p)
    crossing = (c1 * c2 < 0) & (c3 * c4 < 0)
    if crossing.any():
        lam = c3[crossing] / (c3[crossing] - c4[crossing])   # fraction along p
        s_star[crossing] = lam * Lp[crossing]
        dz = np.abs(0.5 * (P0p[crossing, 2] + P1p[crossing, 2]) - 0.5 * (P0q[crossing, 2] + P1q[crossing, 2]))
        dmin[crossing] = np.maximum(dz, 1e-9)
    return s_star, dmin


def _inner(X, A, B, geff):
    """Int_AB ds / sqrt(|X - r(s)|^2 + geff^2)  (geff = 0: plain Neumann kernel)."""
    e = B - A
    L = np.linalg.norm(e, axis=-1)
    t = e / L[..., None]
    u = np.sum((X - A) * t, -1)
    rho2 = np.maximum(np.sum((X - A) ** 2, -1) - u * u, 0.0) + geff * geff
    rho = np.sqrt(rho2)
    plain = geff == 0
    out = np.arcsinh((L - u) / np.maximum(rho, 1e-300)) + np.arcsinh(u / np.maximum(rho, 1e-300))
    # numerically robust plain-kernel form away from the line
    ra = np.linalg.norm(X - A, axis=-1)
    rb = np.linalg.norm(X - B, axis=-1)
    alt = 2.0 * np.arctanh(np.minimum(L / (ra + rb), 1 - 1e-15))
    return np.where(plain & (rho > 1e-3 * L), alt, out)


def _graded_rule(Lp, s_star, h0, ng, K):
    """breakpoints graded geometrically (ratio 2) away from s_star; K pieces per side
    (padding with zero-length pieces).  returns nodes s (n, 2K*ng) and weights (n, 2K*ng)."""
    gx, gw = np.polynomial.legendre.leggauss(ng)
    k = np.arange(K + 1)
    off = h0[:, None] * (2.0 ** k[None, :] - 1.0)          # (n, K+1)
    right = np.minimum(s_star[:, None] + off, Lp[:, None])
    left = np.maximum(s_star[:, None] - off, 0.0)
    right[:, -1] = Lp
    left[:, -1] = 0.0
    a = np.concatenate([right[:, :-1], left[:, 1:]], 1)
    b = np.concatenate([right[:, 1:], left[:, :-1]], 1)
    hl = 0.5 * (b - a)                                     # may be zero or negative-free
    mid = 0.5 * (a + b)
    s = (mid[:, :, None] + hl[:, :, None] * gx[None, None, :]).reshape(len(Lp), -1)
    w = (np.abs(hl)[:, :, None] * gw[None, None, :]).reshape(len(Lp), -1)
    return s, w


def pair_sum(C, alpha=0.25, ng=8, gmd_pair_max=3e-3, use_pair_gmd=True, r_loc=0.8e-3,
             use_reduced_local=True):
    """Return S[wi, wj] = sum_{i in wi, j in wj} Mp_ij (H) and statistics.

    Kernel per pair:
      exactly parallel      -> Campbell, d_eff = self-GMD (collinear), pair-GMD (other winding or
                               own winding >= r_loc, coplanar, < 3 mm), sqrt(d^2+g^2) (own winding < r_loc)
      same winding, d<r_loc -> reduced kernel 1/sqrt(r^2+g^2) (thin-wire GMD regularisation;
                               identical to the Campbell/GMD bar formula for collinear pieces)
      otherwise             -> plain Neumann kernel 1/r
    Outer integral: Gauss-Legendre (ng nodes/piece); single piece if L_p <= alpha*max(d,g),
    else geometrically graded pieces (ratio 2) from the closest point, first piece alpha*max(d,g).
    """
    n = len(C.L)
    nw = int(C.wid.max()) + 1
    S = np.zeros((nw, nw))
    gx, gw = np.polynomial.legendre.leggauss(ng)
    stats = {"parallel": 0, "perpendicular": 0, "gauss_single": 0, "gauss_graded": 0,
             "reduced_kernel_pairs": 0, "evals": 0}
    K = MU0 / (4 * math.pi)
    I, J = np.triu_indices(n)
    CH = 300000
    for c0 in range(0, len(I), CH):
        ii = I[c0:c0 + CH]
        jj = J[c0:c0 + CH]
        Ti, Tj = C.T[ii], C.T[jj]
        dot = np.sum(Ti * Tj, 1)
        crs = np.linalg.norm(np.cross(Ti, Tj), axis=1)
        mult = np.where(ii == jj, 1.0, 2.0)
        val = np.zeros(len(ii))
        same_w = C.wid[ii] == C.wid[jj]
        perp = np.abs(dot) < 1e-12
        stats["perpendicular"] += int(perp.sum())
        par = (crs < 1e-9) & ~perp
        if par.any():
            a, b = ii[par], jj[par]
            u = C.T[a]
            r0 = C.P0[b] - C.P0[a]
            r1 = C.P1[b] - C.P0[a]
            c_ = np.sum(r0 * u, 1)
            e_ = np.sum(r1 * u, 1)
            a2, b2 = np.minimum(c_, e_), np.maximum(c_, e_)
            d = np.linalg.norm(r0 - c_[:, None] * u, axis=1)
            g = np.sqrt(C.gself[a] * C.gself[b])
            same_plane = np.abs(C.P0[a, 2] - C.P0[b, 2]) < 1e-9
            sw = same_w[par]
            deff = d.copy()
            coll = d < 1e-9
            deff[coll] = g[coll]
            loc = (~coll) & sw & (d < r_loc) & use_reduced_local
            deff[loc] = np.sqrt(d[loc] ** 2 + g[loc] ** 2)
            if use_pair_gmd:
                close = (~coll) & (~loc) & same_plane & (d < gmd_pair_max)
                for k in np.where(close)[0]:
                    wk = 0.5 * (C.w[a[k]] + C.w[b[k]])
                    tk = 0.5 * (C.t[a[k]] + C.t[b[k]])
                    if d[k] > wk * 1.0001:
                        deff[k] = gmd_rect_pair(d[k], wk, tk)
            val[par] = K * np.sign(dot[par]) * campbell(0.0, C.L[a], a2, b2, deff)
            stats["parallel"] += int(par.sum())
        rest = np.where(~(par | perp))[0]
        if len(rest):
            a, b = ii[rest], jj[rest]
            res = np.zeros(len(rest))
            for (p, q) in ((a, b), (b, a)):
                s_star, dmin = _closest_param(C.P0[p], C.P1[p], C.P0[q], C.P1[q])
                red = (C.wid[p] == C.wid[q]) & (dmin < r_loc) & use_reduced_local
                geff = np.where(red, np.sqrt(C.gself[p] * C.gself[q]), 0.0)
                scale = np.maximum(dmin, geff)
                scale = np.maximum(scale, 1e-7)
                h0 = alpha * scale
                single = C.L[p] <= h0
                # single piece
                if single.any():
                    ps, qs = p[single], q[single]
                    X = C.P0[ps][:, None, :] + (0.5 * (gx + 1))[None, :, None] * (C.P1[ps] - C.P0[ps])[:, None, :]
                    F = _inner(X, C.P0[qs][:, None, :], C.P1[qs][:, None, :], geff[single][:, None])
                    res[single] += 0.5 * C.L[ps] * np.sum(F * (0.5 * gw)[None, :], 1)
                    stats["gauss_single"] += int(single.sum())
                    stats["evals"] += int(single.sum()) * ng
                gr = ~single
                if gr.any():
                    idx = np.where(gr)[0]
                    Kp = int(np.ceil(np.log2(C.L[p[idx]].max() / h0[idx].min() + 1))) + 1
                    for k0 in range(0, len(idx), 20000):
                        sl = idx[k0:k0 + 20000]
                        ps, qs = p[sl], q[sl]
                        sN, wN = _graded_rule(C.L[ps], s_star[sl], h0[sl], ng, Kp)
                        X = C.P0[ps][:, None, :] + (sN / C.L[ps][:, None])[:, :, None] * (C.P1[ps] - C.P0[ps])[:, None, :]
                        F = _inner(X, C.P0[qs][:, None, :], C.P1[qs][:, None, :], geff[sl][:, None])
                        res[sl] += 0.5 * np.sum(F * wN, 1)
                        stats["evals"] += sN.size
                    stats["gauss_graded"] += int(gr.sum())
                stats["reduced_kernel_pairs"] += int(red.sum())
            val[rest] = K * dot[rest] * res
        np.add.at(S, (C.wid[ii], C.wid[jj]), mult * val * 0.5)
        np.add.at(S, (C.wid[jj], C.wid[ii]), mult * val * 0.5)
    return S, stats


# --------------------------------------------------------------------------
# validation tests
# --------------------------------------------------------------------------
def polygon_loop(R, N, z=0.0, phase=0.0):
    ang = phase + 2 * np.pi * np.arange(N + 1) / N
    P = np.stack([R * np.cos(ang), R * np.sin(ang), np.full(N + 1, z)], 1)
    return P[:-1], P[1:]


def validation():
    out = {}
    # (1) Grover meeting formula vs brute-force double Gauss
    l, m, eps = 0.010, 0.006, math.radians(50)
    xg, wg = np.polynomial.legendre.leggauss(400)
    # graded: integrate in polar-like split to handle corner singularity -> use many nodes
    x = 0.5 * l * (xg + 1)
    y = 0.5 * m * (xg + 1)
    Xg, Yg = np.meshgrid(x, y, indexing="ij")
    Wg = np.outer(wg, wg) * 0.25 * l * m
    brute = float(np.sum(Wg / np.sqrt(Xg ** 2 + Yg ** 2 - 2 * Xg * Yg * math.cos(eps))))
    out["meeting_formula"] = float(meeting(l, m, math.cos(eps)))
    out["meeting_bruteforce_400x400"] = brute
    # (2) square loop of round wire vs Grover closed form L = (2 mu0 a/pi)[ln(a/r) - 0.774]
    #     (Grover's -0.774 is external inductance; uniform-current GMD r*e^-1/4 adds +0.25)
    a, r = 0.010, 1e-4
    sq = np.array([[0, 0, 0], [a, 0, 0], [a, a, 0], [0, a, 0], [0, 0, 0]], float)
    C = Conductor(sq[:-1], sq[1:], np.zeros(4), np.arange(4), [r] * 4, [r] * 4,
                  gmd_self_override={k: r * math.exp(-0.25) for k in range(4)})
    S, _ = pair_sum(C, use_pair_gmd=False)
    out["square_loop_code_nH"] = S[0, 0] * 1e9
    out["square_loop_grover_uniform_current_nH"] = (2 * MU0 * a / math.pi) * (math.log(a / r) - 0.774 + 0.25) * 1e9
    # (3) circular loop (polygon) vs mu0 R (ln(8R/r) - 2 + 1/4)  [round wire, uniform current]
    R = 0.010
    res = []
    for N in (16, 64, 256):
        P0, P1 = polygon_loop(R, N)
        C = Conductor(P0, P1, np.zeros(N), np.arange(N), [r] * N, [r] * N,
                      gmd_self_override={k: r * math.exp(-0.25) for k in range(N)})
        S, _ = pair_sum(C, use_pair_gmd=False)
        res.append((N, S[0, 0] * 1e9))
    out["circle_loop_polygon_nH"] = res
    out["circle_loop_exact_uniform_current_nH"] = MU0 * R * (math.log(8 * R / r) - 1.75) * 1e9
    # (4) two coplanar concentric circles, mutual (Maxwell) vs polygons
    R1, R2 = 0.010, 0.006
    k2 = 4 * R1 * R2 / (R1 + R2) ** 2
    kk = math.sqrt(k2)
    Mex = MU0 * math.sqrt(R1 * R2) * ((2 / kk - kk) * ellipk(k2) - (2 / kk) * ellipe(k2))
    res = []
    for N in (16, 64, 256):
        a0, a1 = polygon_loop(R1, N)
        b0, b1 = polygon_loop(R2, N, phase=0.37)
        C = Conductor(np.vstack([a0, b0]), np.vstack([a1, b1]), np.r_[np.zeros(N), np.ones(N)],
                      np.arange(2 * N), [r] * (2 * N), [r] * (2 * N))
        S, _ = pair_sum(C, use_pair_gmd=False)
        res.append((N, S[0, 1] * 1e9))
    out["coplanar_circles_M_polygon_nH"] = res
    out["coplanar_circles_M_exact_nH"] = Mex * 1e9
    # (5) GMD sanity: rectangle self GMD thin limit, pair GMD vs thin-strip closed form
    w, t, d = 0.254e-3, 35e-6, 0.508e-3
    F = lambda u: 0.5 * u * u * math.log(abs(u)) - 0.75 * u * u
    lng_thin = (F(d + w) + F(d - w) - 2 * F(d)) / w ** 2
    out["gmd_pair_thin_strip_closed_form_um"] = math.exp(lng_thin) * 1e6
    out["gmd_pair_numeric_t=0_um"] = gmd_rect_pair(d, w, 1e-12) * 1e6
    out["gmd_pair_numeric_t=35um_um"] = gmd_rect_pair(d, w, t) * 1e6
    out["gmd_self_w254_t35_um"] = gmd_rect_self(w, t) * 1e6
    out["gmd_self_0.2235(w+t)_um"] = 0.2235 * (w + t) * 1e6
    return out


# --------------------------------------------------------------------------
# build the card's conductor model
# --------------------------------------------------------------------------
def build_card(G, dp_tol_in=0.001, w=0.254e-3, t=35e-6, h_board=0.8e-3,
               include_bottom=True, include_via=True, via_r=0.15e-3, max_piece=None):
    P0s, P1s, wids, bars, ws, ts = [], [], [], [], [], []
    gover = {}
    bar = 0
    for k, name in enumerate(("AC", "DC")):
        W = G["windings"][name]
        top = np.asarray(W["top_xy_in"], float)
        if dp_tol_in > 0:
            top = gg.dp_simplify(top, dp_tol_in)
        else:
            top = gg.dp_simplify(top, 0.0)  # only removes duplicate vertices
        top3 = np.c_[top * IN, np.zeros(len(top))]
        path = [top3]
        if include_via:
            v = top3[-1]
            path.append(np.array([v, v + [0, 0, -h_board]]))
        if include_bottom:
            bt = np.asarray(W["bottom_xy_in"], float) * IN
            bt3 = np.c_[bt, np.full(len(bt), -h_board if include_via else 0.0)]
            path.append(bt3)
        pts = path[0]
        for p in path[1:]:
            pts = np.vstack([pts, p[1:]])
        for i in range(len(pts) - 1):
            seglen = np.linalg.norm(pts[i + 1] - pts[i])
            if seglen < 1e-12:
                continue
            npc = 1 if not max_piece else int(np.ceil(seglen / max_piece))
            vertical = abs(pts[i + 1][2] - pts[i][2]) > 1e-12
            for q in range(npc):
                P0s.append(pts[i] + (pts[i + 1] - pts[i]) * q / npc)
                P1s.append(pts[i] + (pts[i + 1] - pts[i]) * (q + 1) / npc)
                wids.append(k)
                bars.append(bar)
                if vertical:
                    ws.append(2 * via_r)
                    ts.append(2 * via_r)
                    gover[len(P0s) - 1] = via_r  # thin tube GMD ~ radius
                else:
                    ws.append(w)
                    ts.append(t)
            bar += 1
    return Conductor(P0s, P1s, wids, bars, ws, ts, gmd_self_override=gover)


def summarise(S):
    LA, LD, M = S[0, 0], S[1, 1], S[0, 1]
    return {"L_A_uH": LA * 1e6, "L_D_uH": LD * 1e6, "M_uH": M * 1e6,
            "k": M / math.sqrt(LA * LD),
            "L_aid_uH": (LA + LD + 2 * M) * 1e6, "L_opp_uH": (LA + LD - 2 * M) * 1e6,
            "Lleak_A_uH": LA * (1 - (M * M) / (LA * LD)) * 1e6}


def run():
    G = gg.load_geometry()
    res = {}
    t0 = time.time()
    res["validation"] = validation()
    print("validation", json.dumps(res["validation"], indent=1))
    # ---------------- convergence study (baseline geometry: DP 1.0 mil, drawn w, 35 um)
    conv = []
    C = build_card(G, dp_tol_in=0.001)
    print("segments (DP 1.0 mil):", len(C.L))
    for alpha, ng in ((4, 2), (2, 2), (1, 2), (1, 4), (0.5, 4), (0.25, 4), (0.25, 8), (0.125, 8)):
        t1 = time.time()
        S, st = pair_sum(C, alpha=alpha, ng=ng)
        r = summarise(S)
        r.update({"alpha": alpha, "ng": ng, "n_segments": len(C.L), "seconds": time.time() - t1,
                  "stats": st})
        conv.append(r)
        print("conv", alpha, ng, {k: round(v, 5) for k, v in r.items() if isinstance(v, float)}, st)
    res["convergence_quadrature"] = conv
    ALPHA, NG = 0.25, 8
    # ---------------- geometry-discretisation study (DP tolerance)
    geo = []
    for tol in (0.0, 0.0006, 0.001, 0.002, 0.003):
        Cg = build_card(G, dp_tol_in=tol)
        S, st = pair_sum(Cg, alpha=ALPHA, ng=NG)
        r = summarise(S)
        r.update({"dp_tol_mil": tol * 1000, "n_segments": len(Cg.L)})
        geo.append(r)
        print("geo", tol, len(Cg.L), {k: round(v, 5) for k, v in r.items() if isinstance(v, float)})
    res["convergence_geometry"] = geo
    # ---------------- uniform subdivision of every straight bar (should be invariant)
    sub = []
    for mp in (None, 4e-3, 2e-3, 1e-3):
        Cs = build_card(G, dp_tol_in=0.001, max_piece=mp)
        t1 = time.time()
        S, st = pair_sum(Cs, alpha=0.5, ng=4)
        r = summarise(S)
        r.update({"max_piece_mm": (mp or 0) * 1e3, "n_segments": len(Cs.L), "seconds": time.time() - t1})
        sub.append(r)
        print("sub", mp, len(Cs.L), {k: round(v, 5) for k, v in r.items() if isinstance(v, float)})
    res["convergence_subdivision"] = sub
    # ---------------- model variants (baseline numerics)
    var = {}
    S, _ = pair_sum(C, alpha=ALPHA, ng=NG)
    var["baseline"] = summarise(S)
    S, _ = pair_sum(C, alpha=ALPHA, ng=NG, use_pair_gmd=False)
    var["no_pair_gmd_filament_centrelines"] = summarise(S)
    Cn = build_card(G, dp_tol_in=0.001, include_bottom=False, include_via=False)
    S, _ = pair_sum(Cn, alpha=ALPHA, ng=NG)
    var["top_layer_only"] = summarise(S)
    for h in (1.0e-3, 1.2e-3, 1.6e-3):
        Ch = build_card(G, dp_tol_in=0.001, h_board=h)
        S, _ = pair_sum(Ch, alpha=ALPHA, ng=NG)
        var["board_%.1fmm" % (h * 1e3)] = summarise(S)
    res["variants"] = var
    # ---------------- sensitivity: width x thickness
    sens = []
    for w in (0.254e-3, 0.20e-3):
        for t in (35e-6, 30e-6):
            Cw = build_card(G, dp_tol_in=0.001, w=w, t=t)
            S, _ = pair_sum(Cw, alpha=ALPHA, ng=NG)
            r = summarise(S)
            r.update({"w_mm": w * 1e3, "t_um": t * 1e6})
            sens.append(r)
            print("sens", w, t, {k: round(v, 5) for k, v in r.items() if isinstance(v, float)})
    res["sensitivity_w_t"] = sens
    res["numerics"] = {"alpha": ALPHA, "ng": NG, "dp_tol_mil": 1.0, "h_board_mm": 0.8,
                       "mu0": MU0, "runtime_s": time.time() - t0}
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1, default=float)
    write_values(res)
    return res


def write_values(res):
    b = res["variants"]["baseline"]
    sens = {(round(r["w_mm"], 3), round(r["t_um"])): r for r in res["sensitivity_w_t"]}
    conv = res["convergence_quadrature"]
    geo = res["convergence_geometry"]
    # numerical uncertainty: spread of last two quadrature levels + DP 0.6..2 mil spread
    qerr = abs(conv[-1]["L_A_uH"] - conv[-2]["L_A_uH"])
    gvals = [g["L_A_uH"] for g in geo if 0.5 <= g["dp_tol_mil"] <= 2.01]
    gerr = (max(gvals) - min(gvals)) / 2
    meta = ("[MODEL] inductance_neumann.py: Neumann partial-inductance sum over Gerber polylines "
            "(DP 1.0 mil), exact Campbell/Grover filament formulas + GMD (Maxwell rectangle self-GMD, "
            "numerical pair GMD), Gauss alpha=%g ng=%d; top+via+bottom route, board 0.8 mm (JLCPCB, author-confirmed); "
            "LF, uniform current, air core" % (res["numerics"]["alpha"], res["numerics"]["ng"]))
    e = {
        "theory.L_A_uH": (gg.sig(b["L_A_uH"], 4), meta + "; drawn w 0.254 mm, t 35 um"),
        "theory.L_D_uH": (gg.sig(b["L_D_uH"], 4), meta + "; drawn w 0.254 mm, t 35 um"),
        "theory.M_uH": (gg.sig(b["M_uH"], 4), meta + "; positive with dotted ends = outer ends (DC_IN, AC_RTN)"),
        "theory.k": (gg.sig(b["k"], 4), meta),
        "theory.L_series_aiding_uH": (gg.sig(b["L_aid_uH"], 4), meta + "; L_A+L_D+2M"),
        "theory.L_series_opposing_uH": (gg.sig(b["L_opp_uH"], 4), meta + "; L_A+L_D-2M"),
        "theory.L_leak_A_uH": (gg.sig(b["Lleak_A_uH"], 4), meta + "; L_A(1-k^2), D shorted"),
        "theory.L_A_numerical_uncertainty_uH": (gg.sig(qerr + gerr, 2),
                                               "[CALCULATED] quadrature (last two levels) + half spread over DP 0.6-2 mil"),
        "theory.L_A_top_only_uH": (gg.sig(res["variants"]["top_layer_only"]["L_A_uH"], 4), meta + "; top copper only"),
        "theory.k_top_only": (gg.sig(res["variants"]["top_layer_only"]["k"], 4), meta + "; top copper only"),
        "theory.L_A_w0p20_t30_uH": (gg.sig(sens[(0.2, 30)]["L_A_uH"], 4), meta + "; w 0.20 mm, t 30 um"),
        "theory.L_D_w0p20_t30_uH": (gg.sig(sens[(0.2, 30)]["L_D_uH"], 4), meta + "; w 0.20 mm, t 30 um"),
        "theory.M_w0p20_t30_uH": (gg.sig(sens[(0.2, 30)]["M_uH"], 4), meta + "; w 0.20 mm, t 30 um"),
        "theory.k_w0p20_t30": (gg.sig(sens[(0.2, 30)]["k"], 4), meta + "; w 0.20 mm, t 30 um"),
    }
    gg.update_values(e)


if __name__ == "__main__":
    run()
