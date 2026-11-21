#!/usr/bin/env python3
"""inter-winding capacitance C_AB of the interleaved (A-D-A-D...) coplanar windings.

Run:  python3 -I analysis/capacitance.py
Needs work/theory/geometry.json (adjacent length) and, for the resonance estimates,
work/theory/neumann_results.json.

1) Conformal mapping (primary closed form).  Periodic coplanar strips, alternating
   potentials, zero metal thickness, half-spaces of permittivity eps_up / eps_dn.
   The half-cell (strip centre -> gap centre) maps (zeta = sin(pi z/(2a)), a=(w+s)/2,
   then squaring) onto the quarter-plane coplanar-waveguide problem, giving per gap
   and per unit length
        C'_gap = eps0 (eps_up + eps_dn)/2 * K(k)/K(k'),   k = sin(pi eta / 2),
        eta = w/(w+s)  (metallisation ratio).
   This is the infinite-layer interior-cell modulus used for interdigital electrodes
   by Igreja & Dias, Sens. Actuators A 112 (2004) 291-301 (full text not opened here;
   the expression above is derived independently in A2_theory_report.md and checked
   numerically below).  The fabricated board is 0.8 mm FR-4 (JLCPCB), so the
   FD solver models a 0.8 mm substrate with air below; JLCPCB states FR-4 eps_r 4.5 (2-layer),
   LPI solder mask eps_r 3.8 and 10-15 um thick.
2) 2-D finite-volume Laplace solver of the same periodic half-cell (own code) to
   (a) verify (1), (b) add copper thickness 35 um and a solder-mask layer
   (conformal ~10 um, eps 3.8, JLCPCB; and a 'gap filled' 15 um upper bound).
3) C_AB = C'_gap * L_adj, with L_adj from gerber_geometry.py.
4) Resonance estimates with stated C (see report).
"""
import json
import math
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import ellipk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gerber_geometry as gg  # noqa: E402

EPS0 = 8.8541878128e-12
C0 = 299792458.0


def KK(k):
    """K(k)/K(k') with modulus k (scipy uses parameter m = k^2)."""
    return ellipk(k * k) / ellipk(1 - k * k)


def cgap_conformal(w, pitch, eps_up, eps_dn):
    eta = w / pitch
    k = math.sin(math.pi * eta / 2)
    return EPS0 * (eps_up + eps_dn) / 2 * KK(k)


def cps_isolated(w, s, eps_up, eps_dn):
    """isolated coplanar strip pair (for comparison): C' = eps0 (eps_up+eps_dn)/2 K(k')/K(k), k = s/(s+2w)."""
    k = s / (s + 2 * w)
    return EPS0 * (eps_up + eps_dn) / 2 / KK(k)


SUB_T = 0.8e-3   # fabricated board thickness (JLCPCB 0.8 mm, author)


def fd_halfcell(w, pitch, t_cu, eps_sub, mask_t=0.0, mask_eps=1.0, mask_mode="none",
                h=None, H=1.2e-3, sub_t=SUB_T, H_below=1.6e-3):
    """finite-volume solution of the periodic half-cell; returns C'_gap (F/m)."""
    if h is None:
        h = pitch / 200
    nx = int(round(pitch / h))
    h = pitch / nx
    ns = int(round(max(H_below, sub_t) / h))
    na = int(round(H / h))
    ny = ns + na + 1                       # node rows; y index ns <-> y = 0 (substrate top)
    yc = (np.arange(ny - 1) - ns + 0.5) * h          # cell centre y
    xc = (np.arange(nx) + 0.5) * h                   # cell centre x
    Xc, Yc = np.meshgrid(xc, yc, indexing="ij")      # (nx, ny-1)
    eps = np.where((Yc < 0) & (Yc >= -sub_t), eps_sub, 1.0)
    nt = int(round(t_cu / h))
    tcu = nt * h
    half = w / 2
    inA = (Xc <= half) & (Yc >= 0) & (Yc <= tcu)
    inD = (Xc >= pitch - half) & (Yc >= 0) & (Yc <= tcu)
    cu = inA | inD
    if mask_mode != "none" and mask_t > 0:
        if mask_mode == "conformal":
            # distance to copper (rectangles) or to FR-4 surface, above y=0
            dxA = np.maximum(Xc - half, 0)
            dyA = np.maximum(Yc - tcu, 0)
            dA = np.hypot(dxA, dyA)
            dxD = np.maximum((pitch - half) - Xc, 0)
            dD = np.hypot(dxD, dyA)
            m = (Yc >= 0) & ((Yc <= mask_t) | (dA <= mask_t) | (dD <= mask_t)) & ~cu
        elif mask_mode == "filled":
            m = (Yc >= 0) & (Yc <= tcu + mask_t) & ~cu
        else:
            raise ValueError(mask_mode)
        eps = np.where(m, mask_eps, eps)
    # nodes
    Xn, Yn = np.meshgrid(np.arange(nx + 1) * h, (np.arange(ny) - ns) * h, indexing="ij")
    tol = 1e-12
    nodeA = (Xn <= half + tol) & (Yn >= -tol) & (Yn <= tcu + tol)
    nodeD = (Xn >= pitch - half - tol) & (Yn >= -tol) & (Yn <= tcu + tol)
    N = (nx + 1) * ny
    idx = np.arange(N).reshape(nx + 1, ny)
    # edge coefficients: horizontal edges (i,j)-(i+1,j): cells (i,j-1),(i,j)
    epad = np.zeros((nx, ny + 1))
    epad[:, 1:-1] = eps
    ch = 0.5 * (epad[:, :-1] + epad[:, 1:])          # (nx, ny) horizontal edge coeffs
    epad2 = np.zeros((nx + 2, ny - 1))
    epad2[1:-1, :] = eps
    cv = 0.5 * (epad2[:-1, :] + epad2[1:, :])         # (nx+1, ny-1) vertical edge coeffs
    a = idx[:-1, :].ravel()
    b = idx[1:, :].ravel()
    c = ch.ravel()
    rows = np.concatenate([a, b, a, b])
    cols = np.concatenate([a, b, b, a])
    vals = np.concatenate([c, c, -c, -c])
    a2 = idx[:, :-1].ravel()
    b2 = idx[:, 1:].ravel()
    c2 = cv.ravel()
    rows = np.concatenate([rows, a2, b2, a2, b2])
    cols = np.concatenate([cols, a2, b2, b2, a2])
    vals = np.concatenate([vals, c2, c2, -c2, -c2])
    A = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    V = np.zeros(N)
    fixed = np.zeros(N, bool)
    fixed[idx[nodeA]] = True
    V[idx[nodeA]] = 1.0
    fixed[idx[nodeD]] = True
    V[idx[nodeD]] = 0.0
    free = ~fixed
    Aff = A[free][:, free]
    rhs = -A[free][:, fixed] @ V[fixed]
    V[free] = spla.spsolve(Aff.tocsc(), rhs)
    Vg = V.reshape(nx + 1, ny)
    W2 = np.sum(ch * (Vg[1:, :] - Vg[:-1, :]) ** 2) + np.sum(cv * (Vg[:, 1:] - Vg[:, :-1]) ** 2)
    return EPS0 * W2      # = 2 W / V^2 with V = 1


COARSE_NX = {0.254: 100, 0.2: 127}   # coarse grids on which w/2 is an integer number of cells


def fd_rich(w, pitch, t_cu, eps_sub, mask_t=0.0, **kw):
    """Richardson extrapolation (first order, observed O(h) from the strip-edge singularities):
    C = 2 C(h) - C(2h).  All dimensions are snapped to the coarse grid so that both grids
    describe the identical geometry (snapped values are recorded by the caller)."""
    nxc = COARSE_NX[round(w * 1e3, 3)]
    hc = pitch / nxc
    t_s = round(t_cu / hc) * hc
    m_s = round(mask_t / hc) * hc
    c2 = fd_halfcell(w, pitch, t_s, eps_sub, h=hc, mask_t=m_s, **kw)
    c1 = fd_halfcell(w, pitch, t_s, eps_sub, h=hc / 2, mask_t=m_s, **kw)
    return 2 * c1 - c2


def snapped(w, pitch, x):
    hc = pitch / COARSE_NX[round(w * 1e3, 3)]
    return round(x / hc) * hc


def main():
    G = gg.load_geometry()
    pitch = G["pitch_AD_mm"] * 1e-3
    Ladj = G["adjacent_length_AD_mm"] * 1e-3
    res = {"pitch_m": pitch, "L_adj_m": Ladj}
    # ---------- conformal mapping
    conf = {}
    for w in (0.254e-3, 0.20e-3):
        for er in (4.2, 4.7):
            conf["w%.3f_er%.1f" % (w * 1e3, er)] = cgap_conformal(w, pitch, 1.0, er)
        conf["w%.3f_air" % (w * 1e3)] = cgap_conformal(w, pitch, 1.0, 1.0)
    res["conformal_Cgap_F_per_m"] = conf
    res["cps_isolated_w0.254_er4.4_F_per_m"] = cps_isolated(0.254e-3, pitch - 0.254e-3, 1.0, 4.4)
    # ---------- FD verification (zero thickness, air), grid refinement
    ver = {}
    for w in (0.254e-3, 0.20e-3):
        base_nx = COARSE_NX[round(w * 1e3, 3)]
        for nxc in (base_nx, 2 * base_nx, 4 * base_nx):
            c = fd_halfcell(w, pitch, 0.0, 1.0, h=pitch / nxc)
            ver["w%.3f_air_nx%d" % (w * 1e3, nxc)] = c
        ver["w%.3f_air_conformal" % (w * 1e3)] = cgap_conformal(w, pitch, 1.0, 1.0)
    ver["w0.254_er4.4_richardson"] = fd_rich(0.254e-3, pitch, 0.0, 4.4)
    for w in (0.254e-3, 0.20e-3):
        ver["w%.3f_air_richardson" % (w * 1e3)] = fd_rich(w, pitch, 0.0, 1.0)
    ver["w0.254_er4.4_conformal"] = cgap_conformal(0.254e-3, pitch, 1.0, 4.4)
    res["fd_verification_F_per_m"] = ver
    print("verification", {k: round(v * 1e12, 3) for k, v in ver.items()})
    # ---------- FD cases with copper thickness and mask
    cases = []
    air_cache = {}
    for w in (0.254e-3, 0.20e-3):
        for t in (35e-6, 30e-6):
            for er in (4.2, 4.5, 4.7):
                for (mode, mt, me) in (("none", 0, 1.0), ("conformal", 10e-6, 3.8),
                                       ("filled", 15e-6, 3.8)):
                    c = fd_rich(w, pitch, t, er, mask_t=mt, mask_eps=me, mask_mode=mode)
                    if (w, t) not in air_cache:
                        air_cache[(w, t)] = fd_rich(w, pitch, t, 1.0)
                    cair = air_cache[(w, t)]
                    cases.append({"w_mm": w * 1e3, "t_um": t * 1e6, "eps_FR4": er, "mask": mode,
                                  "t_um_snapped": snapped(w, pitch, t) * 1e6,
                                  "mask_t_um_snapped": snapped(w, pitch, mt) * 1e6,
                                  "mask_t_um": mt * 1e6, "mask_eps": me, "Cgap_pF_per_m": c * 1e12,
                                  "Cgap_air_pF_per_m": cair * 1e12, "eps_eff": c / cair,
                                  "C_AB_pF": c * Ladj * 1e12})
                    print(cases[-1])
    res["fd_cases"] = cases
    # ---------- summary ranges
    bare = [c for c in cases if c["mask"] == "none"]
    masked = [c for c in cases if c["mask"] == "conformal"]
    filled = [c for c in cases if c["mask"] == "filled"]
    drawn = [c for c in cases if abs(c["w_mm"] - 0.254) < 1e-6 and abs(c["t_um"] - 35) < 1e-6]
    zero_t = {er: cgap_conformal(0.254e-3, pitch, 1.0, er) * Ladj * 1e12 for er in (4.2, 4.7)}
    summ = {
        "C_AB_conformal_zero_t_bare_pF": [zero_t[4.2], zero_t[4.7]],
        "C_AB_bare_pF_range": [min(c["C_AB_pF"] for c in bare), max(c["C_AB_pF"] for c in bare)],
        "C_AB_masked_pF_range": [min(c["C_AB_pF"] for c in masked), max(c["C_AB_pF"] for c in masked)],
        "C_AB_filled_upper_pF": max(c["C_AB_pF"] for c in filled),
        "C_AB_overall_pF_range": [min(c["C_AB_pF"] for c in cases), max(c["C_AB_pF"] for c in cases)],
        "C_AB_drawn_bare_pF": [c["C_AB_pF"] for c in drawn if c["mask"] == "none"],
        "C_AB_drawn_masked_pF": [c["C_AB_pF"] for c in drawn if c["mask"] == "conformal"],
        "eps_eff_range": [min(c["eps_eff"] for c in cases), max(c["eps_eff"] for c in cases)],
    }
    res["summary"] = summ
    print(json.dumps(summ, indent=1))
    with open(os.path.join(gg.WORK, "capacitance_results.json"), "w") as f:
        json.dump(res, f, indent=1)
    allC = [c["C_AB_pF"] for c in cases if c["mask"] != "filled"]
    cmid = [c["C_AB_pF"] for c in drawn if c["mask"] == "conformal" and c["mask_t_um"] == 10 and c["eps_FR4"] == 4.5]
    meta = ("[MODEL] capacitance.py: C'_gap (periodic coplanar strips; conformal mapping "
            "k=sin(pi*eta/2), verified + extended with own 2-D finite-volume solver for 30-35 um copper and "
            "solder mask; 0.8 mm FR-4, JLCPCB eps 4.5, mask 10 um eps 3.8) x adjacent length %.3f m from gerber_geometry.py" % Ladj)
    gg.update_values({
        "theory.C_AB_conformal_bare_zero_t_pF_lo": (gg.sig(zero_t[4.2], 3), meta + "; eps_r 4.2, zero thickness, no mask"),
        "theory.C_AB_conformal_bare_zero_t_pF_hi": (gg.sig(zero_t[4.7], 3), meta + "; eps_r 4.7, zero thickness, no mask"),
        "theory.C_AB_pF_lo": (gg.sig(min(allC), 3), meta + "; minimum over w{0.254,0.20} mm, t{30,35} um, eps_r{4.2,4.5,4.7}, bare/conformal mask"),
        "theory.C_AB_pF_hi": (gg.sig(max(allC), 3), meta + "; maximum over same set (excluding 'gap filled' bound)"),
        "theory.C_AB_pF_filled_upper": (gg.sig(summ["C_AB_filled_upper_pF"], 3), meta + "; mask fills gaps to copper top +15 um, eps 3.8 (upper bound)"),
        "theory.C_AB_pF_nominal": (gg.sig(cmid[0], 3), meta + "; nominal: drawn w, 35 um, eps_r 4.5, 10 um mask eps 3.8 (JLCPCB)"),
        "theory.Cgap_conformal_pF_per_m_er4p4": (gg.sig(cgap_conformal(0.254e-3, pitch, 1, 4.4) * 1e12, 4), meta + "; per gap, eps_r 4.4, zero thickness"),
        "theory.C_AB_bare_pF_lo": (gg.sig(summ["C_AB_bare_pF_range"][0], 3), meta + "; no mask, over w, t, eps_r"),
        "theory.C_AB_bare_pF_hi": (gg.sig(summ["C_AB_bare_pF_range"][1], 3), meta + "; no mask, over w, t, eps_r"),
        "theory.C_AB_masked_pF_lo": (gg.sig(summ["C_AB_masked_pF_range"][0], 3), meta + "; conformal mask 10 um/eps 3.8, over w, t, eps_r 4.2-4.7"),
        "theory.C_AB_masked_pF_hi": (gg.sig(summ["C_AB_masked_pF_range"][1], 3), meta + "; conformal mask 10 um/eps 3.8, over w, t, eps_r 4.2-4.7"),
        "theory.eps_eff_lo": (gg.sig(summ["eps_eff_range"][0], 3), meta + "; C'/C'_air"),
        "theory.eps_eff_hi": (gg.sig(summ["eps_eff_range"][1], 3), meta + "; C'/C'_air"),
    })
    return res


if __name__ == "__main__":
    main()
