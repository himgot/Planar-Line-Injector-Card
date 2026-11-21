#!/usr/bin/env python3
"""Values that use the fabrication data (JLCPCB, 0.8 mm FR-4, black mask, white silkscreen).

Run:  python3 -I analysis/fabrication_values.py
Reads data/values.json (theory.*, meas.*) and data/bench_reconciled.csv; writes data/values_fab.json.
Every quantity here is [CALCULATED] from values already in the paper or from data-sheet figures stated in
the comments.
"""
import csv
import json
import math
import os
import sys

import numpy as np

ROOT = "."
sys.path.insert(0, os.path.join(ROOT, "analysis"))
V = json.load(open(os.path.join(ROOT, "data", "values.json")))
V = V.get("values", V)


def g(k):
    return float(V[k])


out = {}
RL, Rg50 = 111.6, 50.0
# ---------------------------------------------------------------- measured-parameter model (= M2)
LA, LB, Mm, R = g("meas.L_A_open_10k") * 1e-6, g("meas.L_B_open_10k") * 1e-6, g("meas.M_pair_mean_10k") * 1e-6, g("meas.R_dc_ohm")


def H2(f, Rg=Rg50, RA=None, RB=None, LA_=LA, LB_=LB, M_=Mm, RL_=RL):
    RA = R if RA is None else RA
    RB = R if RB is None else RB
    s = 2j * math.pi * f
    za, zb = Rg + RA + s * LA_, RB + RL_ + s * LB_
    return abs(s * M_ * RL_ / (za * zb - (s * M_) ** 2))


# 1) Run X deviation from M2 when normalised at its ~1 MHz point (as plotted), vs fitted-gain convention
rows = list(csv.DictReader(open(os.path.join(ROOT, "data", "bench_reconciled.csv"))))
X = [r for r in rows if r["run"] == "X" and r["rel_dB_re_1MHz"] not in ("", "nan") and float(r["f_kHz"]) <= 1100]
fref = [float(r["f_kHz"]) for r in X if abs(float(r["rel_dB_re_1MHz"])) < 1e-9][0] * 1e3
dev = {}
for r in X:
    f = float(r["f_kHz"]) * 1e3
    m = 20 * math.log10(H2(f) / H2(fref))
    dev[float(r["f_kHz"])] = float(r["rel_dB_re_1MHz"]) - m
inband = {f: d for f, d in dev.items() if f >= 100}
out["fab.X_norm1MHz_max_dB"] = round(max(abs(d) for d in inband.values()), 1)
out["fab.X_norm1MHz_f_kHz"] = round(max(inband, key=lambda f: abs(inband[f])), 1)
d3 = np.array(list(inband.values()))
out["fab.X_fitgain_max_dB"] = round(float(np.max(np.abs(d3 - d3.mean()))), 1)

# 2) skin depth and a sensitivity bound for AC resistance (both windings)
rho, mu0 = 1.724e-8, 4e-7 * math.pi
w, t = 0.254e-3, 35e-6
for f, tag in ((1e6, "1M"), (5.18e6, "5M")):
    d = math.sqrt(rho / (math.pi * f * mu0))
    out[f"fab.skin_um_{tag}"] = round(d * 1e6)
    out[f"fab.w_over_delta_{tag}"] = round(w / d, 1)
    out[f"fab.t_over_delta_{tag}"] = round(t / d, 1)
for fac in (2, 3):
    for f, tag in ((1e6, "1M"), (5e6, "5M")):
        out[f"fab.dH_R{fac}x_{tag}_dB"] = round(20 * math.log10(H2(f, RA=fac * R, RB=fac * R) / H2(f)), 1)

# 3) theory response referred to AC_IN (Rg = 0): denominator roots and -3 dB points re peak
LAt, LBt, Mt, Rt = g("theory.L_A_uH") * 1e-6, g("theory.L_D_uH") * 1e-6, g("theory.M_uH") * 1e-6, g("theory.R_A_w0p254_t35_ohm")
ra, rb = Rt, Rt + RL
a2, a1, a0 = LAt * LBt - Mt ** 2, ra * LBt + rb * LAt, ra * rb
roots = sorted(abs(r) / (2 * math.pi) for r in np.roots([a2, a1, a0]))
out["fab.Rg0_root_lo_kHz"] = round(roots[0] / 1e3, 1)
out["fab.Rg0_root_hi_MHz"] = round(roots[1] / 1e6, 2)
ff = np.logspace(3, 8, 20000)
h = np.array([H2(f, Rg=0.0, RA=Rt, RB=Rt, LA_=LAt, LB_=LBt, M_=Mt) for f in ff])
ipk = int(np.argmax(h)); thr = h[ipk] / math.sqrt(2)
lo = ff[np.where(h[:ipk] < thr)[0][-1] + 1]; hi = ff[ipk + np.where(h[ipk:] < thr)[0][0]]
out["fab.Rg0_3dB_lo_kHz"] = round(lo / 1e3, 1)
out["fab.Rg0_3dB_hi_MHz"] = round(hi / 1e6, 2)
out["fab.Rg0_peak_dB"] = round(20 * math.log10(h[ipk]), 1)

# 4) closed-form (Mohan) rectangular-adaptation dimensions, one winding
import closed_form as cf  # noqa: E402
import gerber_geometry as gg  # noqa: E402
P = cf.winding_params(gg.load_geometry())["AC"]
Dm, dm = 0.5 * (P["D1_mm"] + P["D2_mm"]), 0.5 * (float(P["d1_mm"]) + float(P["d2_mm"]))
out["fab.cf_dout_mm"] = round(Dm, 1)
out["fab.cf_din_mm"] = round(dm, 1)
out["fab.cf_davg_mm"] = round(0.5 * (Dm + dm), 1)
out["fab.cf_rho"] = round((Dm - dm) / (Dm + dm), 3)
out["fab.cf_davg_sq_mm"] = round(P["mean_turn_mm"] / cf.PERIM["square"], 1)
out["fab.cf_davg_oct_mm"] = round(P["mean_turn_mm"] / cf.PERIM["octagonal"], 1)
out["fab.cf_s_mm"] = round(P["s_mm"], 3)

# 5) safety arithmetic for the printed (untested) 30 V / 500 mA notes and a lead-disconnect transient
out["fab.Isc_30V_A"] = round(30.0 / R, 1)
out["fab.Psc_30V_W"] = round(30.0 ** 2 / R)
out["fab.kick_LB_V"] = round(LB * 0.1 / 100e-9)       # 100 mA interrupted in 100 ns
out["fab.kick_M_V"] = round(Mm * 0.1 / 100e-9)

# 6) LCR999A envelope (catalogue: L measuring range from 4 uH; D display range 0.0000-9.999)
out["fab.L_opp_min_uH"] = min(g("meas.L_opp_c7_10k"), g("meas.L_opp_c8_10k"), g("meas.L_opp_c7_1k"), g("meas.L_opp_c8_1k"))
out["fab.L_opp_max_uH"] = max(g("meas.L_opp_c7_10k"), g("meas.L_opp_c8_10k"), g("meas.L_opp_c7_1k"), g("meas.L_opp_c8_1k"))

# 7) DC-loop cross-check of the supply's 0.03 A reading against the displayed DC level
Rtot_003 = 5.0 / 0.03
out["fab.Rextra_ohm"] = round(Rtot_003 - (RL + R))
VL_alt = 5.0 * RL / Rtot_003
out["fab.VL_if_Rextra_V"] = round(VL_alt, 1)
out["fab.ratio_if_Rextra"] = round(g("bench.dc_level_displayed_V") / VL_alt)

# 8) k gap and the eps_eff range of the measured C_AB (air C for drawn and narrowed width)
out["fab.k_gap"] = round(g("theory.k") - g("meas.k_pair_mean_10k"), 3)
out["fab.k_gap_pct"] = round(100 * (g("theory.k") - g("meas.k_pair_mean_10k")) / g("theory.k"), 1)
_cr = json.load(open(os.path.join(ROOT, "work", "theory", "capacitance_results.json")))
_air = {(round(c["w_mm"], 3), round(c["t_um"])): c["Cgap_air_pF_per_m"] * _cr["L_adj_m"] for c in _cr["fd_cases"]}
out["fab.Cair_lo_pF"] = round(min(_air.values()))
out["fab.Cair_hi_pF"] = round(max(_air.values()))
out["fab.eps_eff_meas_lo"] = round(274 / max(_air.values()), 1)
out["fab.eps_eff_meas_hi"] = round(274 / min(_air.values()), 1)
out["fab.Cgap_er45_pFm"] = round(8.8541878128e-12 * (1 + 4.5) / 2 * 1e12, 1)
out["fab.D_aid_10k"] = round((2 * R) / (2 * math.pi * 1e4 * 0.5 * (g("meas.L_aid_c5_10k") + g("meas.L_aid_c6_10k")) * 1e-6), 1)
out["fab.D_aid_1k"] = round((2 * R) / (2 * math.pi * 1e3 * 0.5 * (g("meas.L_aid_c5_1k") + g("meas.L_aid_c6_1k")) * 1e-6))
out["fab.D_single_1k"] = round(R / (2 * math.pi * 1e3 * LA))
out["fab.dRel5_lo_dB"] = round(min(out["fab.dH_R2x_5M_dB"] - out["fab.dH_R2x_1M_dB"], out["fab.dH_R3x_5M_dB"] - out["fab.dH_R3x_1M_dB"]), 1)
out["fab.dRel5_hi_dB"] = round(max(out["fab.dH_R2x_5M_dB"] - out["fab.dH_R2x_1M_dB"], out["fab.dH_R3x_5M_dB"] - out["fab.dH_R3x_1M_dB"]), 1)
out["fab.cf_D1_mm"] = round(P["D1_mm"], 1); out["fab.cf_D2_mm"] = round(P["D2_mm"], 1)
out["fab.cf_d1_mm"] = round(float(P["d1_mm"]), 1); out["fab.cf_d2_mm"] = round(float(P["d2_mm"]), 1)

# 9) Run X slope 104-1030 kHz vs M2 and capacitive feed-through; 274 pF vs predicted bounds
_xr = {float(r["f_kHz"]): float(r["rel_dB_re_1MHz"]) for r in X}
_flo = min(f for f in _xr if f >= 100); _fhi = max(_xr)
out["fab.X_slope_dB"] = round(_xr[_fhi] - _xr[_flo], 1)
out["fab.X_slope_M2_dB"] = round(20 * math.log10(H2(_fhi * 1e3) / H2(_flo * 1e3)), 1)
_cap = lambda f: (lambda wc: wc / math.sqrt(1 + wc * wc))(2 * math.pi * f * 274e-12 * (50 + RL))
out["fab.X_slope_cap_dB"] = round(20 * math.log10(_cap(_fhi * 1e3) / _cap(_flo * 1e3)), 1)
out["fab.C274_over_nom_pct"] = round(100 * (274 / g("theory.C_AB_pF_nominal") - 1))
out["fab.C274_over_filled_pct"] = round(100 * (274 / g("theory.C_AB_pF_filled_upper") - 1))

# 10) session spread vs the prediction: October records read 0.35-0.7 uH above November
_Lnov = 0.5 * (g("meas.L_A_open_10k") + g("meas.L_B_open_10k"))
out["fab.Loct_lo_uH"] = round(_Lnov + g("meas.L_oct_minus_nov_min_uH"), 2)
out["fab.Loct_hi_uH"] = round(_Lnov + g("meas.L_oct_minus_nov_max_uH"), 2)
_Lth = 0.5 * (g("theory.L_A_uH") + g("theory.L_D_uH"))
out["fab.Loct_vs_th_lo_pct"] = round(100 * (out["fab.Loct_lo_uH"] / _Lth - 1))
out["fab.Loct_vs_th_hi_pct"] = round(100 * (out["fab.Loct_hi_uH"] / _Lth - 1))

# 11) Run X slope with larger generator source resistance
for _rg in (200.0, 600.0):
    out[f"fab.X_slope_Rg{int(_rg)}_dB"] = round(20 * math.log10(H2(_fhi * 1e3, Rg=_rg) / H2(_flo * 1e3, Rg=_rg)), 1)

# 12) quasi-TEM cross-check of the bifilar pair: L_opp * C_air ~ l^2/c^2 (odd mode), l = mean winding length
_l = 0.5 * (g("theory.geom.A_len_total_mm") + g("theory.geom.D_len_total_mm")) * 1e-3 if "theory.geom.A_len_total_mm" in V else 3.0
_u = (_l / 2.99792458e8) ** 2
out["fab.tem_len_m"] = round(_l, 2)
out["fab.tem_model_lo"] = round(g("theory.L_series_opposing_uH") * 1e-6 * out["fab.Cair_lo_pF"] * 1e-12 / _u, 1)
out["fab.tem_model_hi"] = round(g("theory.L_series_opposing_uH") * 1e-6 * out["fab.Cair_hi_pF"] * 1e-12 / _u, 1)
_Lopp_m = 0.5 * (g("meas.L_opp_c7_10k") + g("meas.L_opp_c8_10k"))
out["fab.tem_meas"] = round(_Lopp_m * 1e-6 * 274e-12 / _u, 1)

json.dump(out, open(os.path.join(ROOT, "data", "values_fab.json"), "w"), indent=1)
for k, v in out.items():
    print(k, v)
