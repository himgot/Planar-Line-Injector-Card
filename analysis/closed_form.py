#!/usr/bin/env python3
"""closed-form self-inductance estimates (Mohan et al., IEEE JSSC 34(10):1419-1424, 1999)
for each winding of the card, with explicit adaptations of the rectangular,
chamfered spiral to the square / octagonal formulas.

Run:  python3 -I analysis/closed_form.py   (needs work/theory/geometry.json)

Formulas (Mohan 1999; coefficient values cross-checked against the reproductions in
PMC8309821 Tables 1-2, IEEE SSC Magazine Fall 2020 p.89 (K1,K2) and arXiv:2507.16982 (square)):
  modified Wheeler : L = K1 mu0 n^2 d_avg / (1 + K2 rho)
  current sheet    : L = mu0 n^2 d_avg c1/2 [ln(c2/rho) + c3 rho + c4 rho^2]
  monomial         : L[nH] = beta d_out^a1 w^a2 d_avg^a3 n^a4 s^a5   (lengths in um)
  d_avg = (d_out + d_in)/2,  rho = (d_out - d_in)/(d_out + d_in);
  d_out, d_in = outer / inner diameters (outer edge / inner edge).

Rectangular adaptations used (D1 x D2 outer edges, d1 x d2 inner edges of one winding):
  R1 'mean side'      : D = (D1+D2)/2, d = (d1+d2)/2
  R2 'power mean'     : D = PM_p(D1,D2), d = D - 2n(w+s) + 2s  (Papadopoulos & Antonopoulos,
                        arXiv:2507.16982, 2025; p = -0.16 WH, -0.05 current sheet, -1.13 monomial;
                        published for SQUARE coefficients only)
  R3 'equal perimeter': d_avg = (mean turn length)/(perimeter per unit diameter of the shape:
                        4 for square, 8 tan 22.5 deg for octagon), rho taken from R1
Series-aiding estimate: the two windings carry the same current in the same sense,
which is the current distribution of one spiral with 2n turns at the A-D pitch.
"""
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gerber_geometry as gg  # noqa: E402

MU0 = 4e-7 * math.pi
IN2MM = 25.4

WHEELER = {"square": (2.34, 2.75), "octagonal": (2.25, 3.55)}
CSHEET = {"square": (1.27, 2.07, 0.18, 0.13), "octagonal": (1.07, 2.29, 0.00, 0.19)}
MONO = {"square": (1.62e-3, -1.21, -0.147, 2.40, 1.78, -0.030),
        "octagonal": (1.33e-3, -1.21, -0.163, 2.43, 1.75, -0.049)}
PERIM = {"square": 4.0, "octagonal": 8 * math.tan(math.radians(22.5))}
P_PA = {"wheeler": -0.16, "current_sheet": -0.05, "monomial": -1.13}


def L_wheeler(n, dout, din, shape):
    K1, K2 = WHEELER[shape]
    davg, rho = 0.5 * (dout + din), (dout - din) / (dout + din)
    return K1 * MU0 * n * n * davg / (1 + K2 * rho)


def L_csheet(n, dout, din, shape):
    c1, c2, c3, c4 = CSHEET[shape]
    davg, rho = 0.5 * (dout + din), (dout - din) / (dout + din)
    return MU0 * n * n * davg * c1 / 2 * (math.log(c2 / rho) + c3 * rho + c4 * rho * rho)


def L_mono(n, dout, din, w, s, shape):
    b, a1, a2, a3, a4, a5 = MONO[shape]
    davg = 0.5 * (dout + din)
    um = 1e6
    return b * (dout * um) ** a1 * (w * um) ** a2 * (davg * um) ** a3 * n ** a4 * (s * um) ** a5 * 1e-9


def power_mean(a, b, p):
    if abs(p) < 1e-12:
        return math.sqrt(a * b)
    return ((a ** p + b ** p) / 2) ** (1 / p)


def inner_span(scan, letter):
    pos = np.array(scan["positions_in"])
    order = scan["order"]
    g = int(np.argmax(np.diff(pos)))
    below = [i for i in range(g + 1) if order[i] == letter][-1]
    above = [i for i in range(g + 1, len(order)) if order[i] == letter][0]
    return (pos[above] - pos[below]) * IN2MM


def winding_params(G):
    w_mm = G["trace_width_mm"]
    pitch = G["pitch_AD_mm"]
    out = {}
    for name, letter in (("AC", "A"), ("DC", "D")):
        W = G["windings"][name]
        xy = np.array(W["top_xy_in"])
        inb = (xy[:, 0] > 0.30) & (xy[:, 0] < 3.05)
        k = int(np.argmax(inb))
        body = xy[k:]
        n_body = abs(gg.winding_angle(body, G["body_centre_in"])) / (2 * math.pi)
        lbody = gg.polyline_len(gg.dp_simplify(body, 0.001)) * IN2MM
        bb = W["body_bbox_in"]
        D1 = (bb[2] - bb[0]) * IN2MM + w_mm
        D2 = (bb[3] - bb[1]) * IN2MM + w_mm
        d1 = inner_span(G["scans"]["horizontal_through_centre"], letter) - w_mm
        d2 = inner_span(G["scans"]["vertical_through_centre"], letter) - w_mm
        out[name] = {"n_total": W["turns"], "n_body": n_body, "len_body_mm": lbody,
                     "mean_turn_mm": lbody / n_body, "D1_mm": D1, "D2_mm": D2, "d1_mm": d1, "d2_mm": d2,
                     "w_mm": w_mm, "s_mm": 2 * pitch - w_mm}
    # series-aiding composite (both windings, 2n turns at A-D pitch)
    bbx = G["body_centreline_bbox_in"]
    out["SERIES"] = {"n_body": out["AC"]["n_body"] + out["DC"]["n_body"],
                     "len_body_mm": out["AC"]["len_body_mm"] + out["DC"]["len_body_mm"],
                     "D1_mm": (bbx[2] - bbx[0]) * IN2MM + w_mm, "D2_mm": (bbx[3] - bbx[1]) * IN2MM + w_mm,
                     "d1_mm": G["inner_hole_centreline_W_mm"] - w_mm,
                     "d2_mm": G["inner_hole_centreline_H_mm"] - w_mm,
                     "w_mm": w_mm, "s_mm": pitch - w_mm}
    out["SERIES"]["mean_turn_mm"] = out["SERIES"]["len_body_mm"] / out["SERIES"]["n_body"]
    return out


def estimates(p):
    """return list of (formula, shape, adaptation, L_H) for one winding parameter set."""
    n = p["n_body"]
    w, s = p["w_mm"] * 1e-3, p["s_mm"] * 1e-3
    D1, D2, d1, d2 = (p[k] * 1e-3 for k in ("D1_mm", "D2_mm", "d1_mm", "d2_mm"))
    res = []
    Dm, dm = 0.5 * (D1 + D2), 0.5 * (d1 + d2)
    rho1 = (Dm - dm) / (Dm + dm)
    for shape in ("square", "octagonal"):
        # R1
        res.append(("wheeler", shape, "R1_mean_side", L_wheeler(n, Dm, dm, shape)))
        res.append(("current_sheet", shape, "R1_mean_side", L_csheet(n, Dm, dm, shape)))
        res.append(("monomial", shape, "R1_mean_side", L_mono(n, Dm, dm, w, s, shape)))
        # R3 equal mean perimeter
        davg = p["mean_turn_mm"] * 1e-3 / PERIM[shape]
        Do, Di = davg * (1 + rho1), davg * (1 - rho1)
        res.append(("wheeler", shape, "R3_equal_perimeter", L_wheeler(n, Do, Di, shape)))
        res.append(("current_sheet", shape, "R3_equal_perimeter", L_csheet(n, Do, Di, shape)))
        res.append(("monomial", shape, "R3_equal_perimeter", L_mono(n, Do, Di, w, s, shape)))
    # R2 power mean (square only, as published)
    for f, fn in (("wheeler", L_wheeler), ("current_sheet", L_csheet), ("monomial", None)):
        D = power_mean(D1, D2, P_PA[f])
        d = D - 2 * n * (w + s) + 2 * s
        if d <= 0:
            res.append((f, "square", "R2_power_mean", float("nan")))
            continue
        if fn is None:
            res.append((f, "square", "R2_power_mean", L_mono(n, D, d, w, s, "square")))
        else:
            res.append((f, "square", "R2_power_mean", fn(n, D, d, "square")))
    return res


def main():
    G = gg.load_geometry()
    P = winding_params(G)
    rows = []
    allres = {}
    for name in ("AC", "DC", "SERIES"):
        r = estimates(P[name])
        allres[name] = r
        for f, shape, ad, L in r:
            rows.append([name, f, shape, ad, L * 1e6])
    # implied M and k for each (formula, shape, adaptation) combination
    kk = []
    for i, (f, shape, ad, LA) in enumerate(allres["AC"]):
        LD = allres["DC"][i][3]
        LS = allres["SERIES"][i][3]
        M = 0.5 * (LS - LA - LD)
        kk.append([f, shape, ad, LA * 1e6, LD * 1e6, LS * 1e6, M * 1e6, M / math.sqrt(LA * LD)])
    os.makedirs(gg.WORK, exist_ok=True)
    with open(os.path.join(gg.WORK, "closed_form_table.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["formula", "shape_coeffs", "adaptation", "L_A_uH", "L_D_uH", "L_series2n_uH",
                     "M_implied_uH", "k_implied"])
        for r in kk:
            cw.writerow([r[0], r[1], r[2]] + ["%.4g" % v for v in r[3:]])
    with open(os.path.join(gg.WORK, "closed_form_params.json"), "w") as f:
        import json
        json.dump(P, f, indent=1)
    LAs = np.array([r[3] for r in kk if np.isfinite(r[3])])
    sq_r1 = {r[0]: r for r in kk if r[1] == "square" and r[2] == "R1_mean_side"}
    meta = ("[MODEL] closed_form.py: Mohan et al. 1999 (IEEE JSSC 34(10):1419-1424) %s; "
            "n=%.2f body turns, w=%.3f mm, s=%.3f mm; outer %.2f x %.2f mm, inner %.2f x %.2f mm (edges)" %
            ("%s", P["AC"]["n_body"], P["AC"]["w_mm"], P["AC"]["s_mm"], P["AC"]["D1_mm"], P["AC"]["D2_mm"],
             P["AC"]["d1_mm"], P["AC"]["d2_mm"]))
    e = {
        "theory.cf.L_A_wheeler_sq_meanside_uH": (gg.sig(sq_r1["wheeler"][3], 3), meta % "modified Wheeler, square coeffs, mean-side rectangle"),
        "theory.cf.L_A_csheet_sq_meanside_uH": (gg.sig(sq_r1["current_sheet"][3], 3), meta % "current sheet, square coeffs, mean-side rectangle"),
        "theory.cf.L_A_monomial_sq_meanside_uH": (gg.sig(sq_r1["monomial"][3], 3), meta % "monomial, square coeffs, mean-side rectangle"),
        "theory.cf.L_A_min_uH": (gg.sig(LAs.min(), 3), meta % "minimum over 3 formulas x {square, octagonal} x {R1,R3} + R2 (square)"),
        "theory.cf.L_A_max_uH": (gg.sig(LAs.max(), 3), meta % "maximum over 3 formulas x {square, octagonal} x {R1,R3} + R2 (square)"),
        "theory.cf.L_A_median_uH": (gg.sig(float(np.median(LAs)), 3), meta % "median over all combinations"),
        "theory.cf.n_body_turns_A": (gg.sig(P["AC"]["n_body"], 4), "[CALCULATED] closed_form.py: winding angle of AC body (outer lead excluded)/360"),
        "theory.cf.mean_turn_length_A_mm": (gg.sig(P["AC"]["mean_turn_mm"], 4), "[CALCULATED] closed_form.py: body length (DP 1 mil) / body turns"),
    }
    gg.update_values(e)
    print("params", {k: {a: round(b, 3) for a, b in v.items()} for k, v in P.items()})
    for r in kk:
        print("%-14s %-9s %-19s LA=%7.3f LD=%7.3f Lser=%8.3f M=%7.3f k=%.3f" % tuple(r))
    print("L_A range %.2f .. %.2f uH, median %.2f" % (LAs.min(), LAs.max(), np.median(LAs)))


if __name__ == "__main__":
    main()
