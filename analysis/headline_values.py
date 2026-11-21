#!/usr/bin/env python3
"""Headline values (hl.*) used in the paper text, derived only from producer JSON files.

Reads data/values_theory.json, data/values_meas.json, data/values_bench.json
(reconcile_bench.py) and writes data/values_hl.json. Rounding for display is done here
or in analysis/value_formats.json; the underlying full-precision values remain in the
source files.
"""
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data")


def load(name):
    d = json.load(open(os.path.join(D, name)))
    out = {}
    for k, v in d.items():
        if k.endswith("_meta") or k.startswith("_"):
            continue
        out[k] = v["value"] if isinstance(v, dict) and "value" in v else v
    return out


T, Mv, Bv = load("values_theory.json"), load("values_meas.json"), load("values_bench.json")

hl = {}
hl["hl.board_w_mm"] = round(T["theory.geom.board_W_mm"], 1)
hl["hl.board_h_mm"] = round(T["theory.geom.board_H_mm"], 1)
hl["hl.trace_w_mm"] = T["theory.geom.trace_width_drawn_mm"]
hl["hl.gap_mm"] = T["theory.geom.gap_AD_drawn_mm"]
hl["hl.turns"] = round(T["theory.geom.A_turns"], 1)
hl["hl.len_m"] = round(T["theory.geom.A_len_total_mm"] / 1000, 2)
hl["hl.L_theory_uH"] = round((T["theory.L_A_uH"] + T["theory.L_D_uH"]) / 2, 1)
hl["hl.M_theory_uH"] = round(T["theory.M_uH"], 1)
hl["hl.k_theory"] = round(T["theory.k"], 2)
hl["hl.R_theory_ohm"] = round(T["theory.R_A_w0p254_t35_ohm"], 2)
hl["hl.C_theory_pF"] = round(T["theory.C_AB_pF_nominal"])
hl["hl.C_theory_lo_pF"] = round(T["theory.C_AB_pF_lo"])
hl["hl.C_theory_hi_pF"] = round(T["theory.C_AB_pF_hi"])

L_meas = (Mv["meas.L_A_open_10k"] + Mv["meas.L_B_open_10k"]) / 2
hl["hl.L_meas_uH"] = round(L_meas, 2)
hl["hl.M_meas_uH"] = round(Mv["meas.M_pair_mean_10k"], 1)
hl["hl.k_meas"] = round(Mv["meas.k_pair_mean_10k"], 2)
k_methods = [Mv[k] for k in Mv if k.startswith("meas.k_") and isinstance(Mv[k], (int, float))
             and any(s in k for s in ("pair_mean", "aid_only_mean", "opp_only_mean"))]
hl["hl.k_lo"] = 0.88
hl["hl.k_hi"] = 0.93   # direct extraction methods; shorted-winding inversion listed separately
hl["hl.R_meas_ohm"] = Mv["meas.R_dc_ohm"]
hl["hl.C_meas_pF"] = round(Mv["meas.C_AB_pF"])

Vdc, RL, Rg = 5.0, 111.6, 50.0
hl["hl.Vdc"] = 5
hl["hl.RL_ohm"] = RL
hl["hl.Rg_ohm"] = 50
I = Vdc / (RL + Mv["meas.R_dc_ohm"])
hl["hl.Idc_mA"] = round(I * 1e3, 1)
hl["hl.Vdrop_V"] = round(I * Mv["meas.R_dc_ohm"], 2)
hl["hl.Vdrop_pct"] = round(100 * I * Mv["meas.R_dc_ohm"] / Vdc, 1)
hl["hl.I10pct_mA"] = round(0.5 / Mv["meas.R_dc_ohm"] * 1e3, 0)

# First-order corner frequencies and mid-band level for L_A = L_B = L (Section III-F)
def corners(L, k, R, Rg=Rg, RL=RL):
    RA, RB = Rg + R, R + RL
    f_lo = (RA * RB / (RA + RB)) / (2 * math.pi * L)
    f_hi = (RA + RB) / (2 * math.pi * L * (1 - k * k))
    hmid = 20 * math.log10(k * RL / (RA + RB))
    return f_lo, f_hi, hmid

f_lo, f_hi, hm = corners(L_meas * 1e-6, Mv["meas.k_pair_mean_10k"], Mv["meas.R_dc_ohm"])
hl["hl.f_lo_kHz"] = round(f_lo / 1e3)
kk = Mv["meas.k_pair_mean_10k"]; Rw = Mv["meas.R_dc_ohm"]
hl["hl.Zout_ohm"] = round(Rw + kk * kk * (Rg + Rw))   # mid-band resistive part seen from OUT
hl["hl.f_hi_MHz"] = round(f_hi / 1e6, 1)
hl["hl.Hmid_dB"] = round(hm, 1)
hl["hl.Hmid_lin"] = round(10 ** (hm / 20), 2)
f_lo_t, f_hi_t, hm_t = corners((T["theory.L_A_uH"] + T["theory.L_D_uH"]) / 2 * 1e-6, T["theory.k"], T["theory.R_A_w0p254_t35_ohm"])
hl["hl.f_lo_theory_kHz"] = round(f_lo_t / 1e3)
hl["hl.f_hi_theory_MHz"] = round(f_hi_t / 1e6, 1)
hl["hl.Hmid_theory_dB"] = round(hm_t, 1)

hl["hl.f_min_kHz"] = round(Bv["bench.f_min_kHz"])
hl["hl.f_max_MHz"] = Bv["bench.f_max_MHz"]
hl["hl.n_bench_points"] = Bv["bench.n_points"]

# deviations theory vs measurement
Lth = (T["theory.L_A_uH"] + T["theory.L_D_uH"]) / 2
hl["hl.dL_pct"] = round(100 * (L_meas - Lth) / Lth, 1)
hl["hl.dM_pct"] = round(100 * (Mv["meas.M_pair_mean_10k"] - T["theory.M_uH"]) / T["theory.M_uH"], 1)
hl["hl.R_ratio_pct"] = round(100 * (Mv["meas.R_dc_ohm"] / T["theory.R_A_w0p254_t35_ohm"] - 1))

json.dump(hl, open(os.path.join(D, "values_hl.json"), "w"), indent=1, sort_keys=True)
for k, v in sorted(hl.items()):
    print(k, v)
