#!/usr/bin/env python3
"""Derived values for the tables and prose of Sections III and IV.

Reads only producer files and writes data/values_tables.json (keys "tab.*").
Inputs:
  data/values_theory.json            data/values_meas.json
  data/values_peec_check.json      (PEEC check)   data/values_hl.json   (headline_values.py)
  work/theory/neumann_results.json           (Neumann model variants)
  work/theory/capacitance_results.json       (capacitance cases: assumed eps / mask)
  work/theory/closed_form_table.csv          (closed-form variants)
  data/photo_index_oct24.csv             (LCR-meter photos P24-53..59, text notes)
Run:  python3 -I analysis/table_values.py
"""
import csv, json, math, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data")
W = os.path.join(ROOT, "work", "theory")


def load_flat(path):
    d = json.load(open(path))
    out = {}
    for k, v in d.items():
        if k.endswith("_meta") or k.startswith("_"):
            continue
        out[k] = v["value"] if isinstance(v, dict) and "value" in v else v
    return out


T = load_flat(os.path.join(D, "values_theory.json"))
Mv = load_flat(os.path.join(D, "values_meas.json"))
HL = json.load(open(os.path.join(D, "values_hl.json")))
PC = json.load(open(os.path.join(D, "values_peec_check.json")))
NR = json.load(open(os.path.join(W, "neumann_results.json")))
CR = json.load(open(os.path.join(W, "capacitance_results.json")))

out, meta = {}, {}


def put(key, value, src):
    out[key] = value
    meta[key] = src


# ---------------------------------------------------------------- legend (design notes)
spec = T["theory.legend.top_spec_line"]
bot = T["theory.legend.bot_lines"]
tol = T["theory.legend.top_tolerance_line"]
m = re.search(r"k\s*\S\s*([0-9.]+)\s*-\s*([0-9.]+)", spec)
put("tab.legend.k_lo", float(m.group(1)), "theory.legend.top_spec_line (Gerber legend, read)")
put("tab.legend.k_hi", float(m.group(2)), "theory.legend.top_spec_line (Gerber legend, read)")
m = re.search(r"DCR\s*([0-9.]+)\s*-\s*([0-9.]+)", spec)
put("tab.legend.DCR_lo_ohm", float(m.group(1)), "theory.legend.top_spec_line")
put("tab.legend.DCR_hi_ohm", float(m.group(2)), "theory.legend.top_spec_line")
put("tab.legend.Cu_um", float(re.search(r"Copper thickness = ([0-9.]+)", bot).group(1)), "theory.legend.bot_lines")
put("tab.legend.len_DC_mm", float(re.search(r"DC coil length = ([0-9.]+)", bot).group(1)), "theory.legend.bot_lines")
put("tab.legend.len_AC_mm", float(re.search(r"AC coil length = ([0-9.]+)", bot).group(1)), "theory.legend.bot_lines")
put("tab.legend.tol_pct", float(re.search(r"([0-9.]+)\s*%", tol).group(1)), "theory.legend.top_tolerance_line")

# ---------------------------------------------------------------- geometry helpers
put("tab.geom.pitch_same_mm", 2 * T["theory.geom.pitch_AD_mm"], "2 x theory.geom.pitch_AD_mm")
put("tab.geom.adj_len_m", T["theory.geom.adjacent_length_AD_mm"] / 1000, "theory.geom.adjacent_length_AD_mm / 1000")

# ---------------------------------------------------------------- DC resistance
rho = T["theory.rho_Cu_ohm_m"]
put("tab.rho_nOhm_m", rho * 1e9, "theory.rho_Cu_ohm_m in nOhm m (IACS, 20 C)")
R_meas = Mv["meas.R_dc_ohm"]
lA = T["theory.geom.A_len_total_mm"] * 1e-3
Aeff = rho * lA / R_meas                       # m^2
w_d, t_d = T["theory.geom.trace_width_drawn_mm"] * 1e-3, 35e-6
put("tab.R.A_drawn_um2", w_d * t_d * 1e12, "drawn 0.254 mm x 35 um")
put("tab.R.Aeff_um2", Aeff * 1e12, "rho*l_A/R_meas (R_meas = meas.R_dc_ohm; includes any lead/contact R)")
put("tab.R.Aeff_ratio_pct", 100 * Aeff / (w_d * t_d), "Aeff / drawn cross-section")
put("tab.R.w_eff_t35_mm", Aeff / 35e-6 * 1e3, "Aeff / 35 um")
put("tab.R.t_eff_w254_um", Aeff / w_d * 1e6, "Aeff / 0.254 mm")
put("tab.R.excess_pct", 100 * (R_meas / T["theory.R_A_w0p254_t35_ohm"] - 1), "meas.R_dc_ohm / theory.R_A_w0p254_t35_ohm - 1")
put("tab.R.meas_vs_w020t35_pct", 100 * (R_meas / T["theory.R_A_w0p20_t35_ohm"] - 1), "meas.R_dc_ohm / theory.R_A_w0p20_t35_ohm - 1")
put("tab.R.meas_vs_w020t35_abs_pct", abs(out["tab.R.meas_vs_w020t35_pct"]), "absolute value")

# ---------------------------------------------------------------- closed forms
with open(os.path.join(W, "closed_form_table.csv")) as fh:
    cf = list(csv.DictReader(fh))
put("tab.cf.n_variants", len(cf), "rows of work/theory/closed_form_table.csv")
put("tab.cf.dev_median_vs_neumann_pct", 100 * (T["theory.cf.L_A_median_uH"] / T["theory.L_A_uH"] - 1),
    "closed-form median / Neumann L_A - 1")

# ---------------------------------------------------------------- Neumann variants
top = NR["variants"]["top_layer_only"]
put("tab.nm.top.L_A_uH", top["L_A_uH"], "neumann_results.json variants.top_layer_only (AC)")
put("tab.nm.top.L_B_uH", top["L_D_uH"], "neumann_results.json variants.top_layer_only (DC)")
put("tab.nm.top.M_uH", top["M_uH"], "neumann_results.json variants.top_layer_only")
put("tab.nm.top.k", top["k"], "neumann_results.json variants.top_layer_only")
for row in NR["sensitivity_w_t"]:
    if abs(row["w_mm"] - 0.2) < 1e-9 and abs(row["t_um"] - 35) < 1e-9:
        put("tab.nm.w020t35.L_A_uH", row["L_A_uH"], "neumann_results.json sensitivity_w_t w=0.20 t=35")
        put("tab.nm.w020t35.k", row["k"], "neumann_results.json sensitivity_w_t w=0.20 t=35")
put("tab.nm.L_num_unc_nH", T["theory.L_A_numerical_uncertainty_uH"] * 1e3, "theory.L_A_numerical_uncertainty_uH in nH")
put("tab.nm.dk_w020t30", T["theory.k_w0p20_t30"] - T["theory.k"], "k(w0.20,t30) - k(drawn)")
put("tab.nm.dk_w020t30_abs", abs(T["theory.k_w0p20_t30"] - T["theory.k"]), "|k(w0.20,t30) - k(drawn)|")
put("tab.nm.dL_w020t30_pct", 100 * (T["theory.L_A_w0p20_t30_uH"] / T["theory.L_A_uH"] - 1), "L_A(w0.20,t30)/L_A - 1")
for row in NR["sensitivity_w_t"]:
    if abs(row["w_mm"] - 0.2) < 1e-9 and abs(row["t_um"] - 30) < 1e-9:
        put("tab.nm.w020t30.L_opp_uH", row["L_opp_uH"], "neumann_results.json sensitivity_w_t w=0.20 t=30")
b16 = NR["variants"]["board_1.6mm"]
put("tab.nm.board_sens_pct", abs(100 * (b16["L_A_uH"] / NR["variants"]["baseline"]["L_A_uH"] - 1)),
    "|L_A(1.6 mm board)/L_A(0.8 mm, fabricated) - 1|")
put("tab.nm.L_minus_M_uH", T["theory.L_A_uH"] - T["theory.M_uH"], "theory.L_A_uH - theory.M_uH")

# ---------------------------------------------------------------- PEEC check mapped to paper naming
# Paper: A = AC winding (AC_IN-AC_RTN), B = DC winding (DC_IN-OUT).  The PEEC check labels its windings A/B by
# file order.  Map by the raw top-layer length, which the geometry script assigns to the AC winding (theory.geom.A_len_top_raw_mm).
wv = PC["windings"]
ac_key = min(wv, key=lambda kk: abs(wv[kk]["path_length_top_mm"] - T["theory.geom.A_len_top_raw_mm"]))
dc_key = [kk for kk in wv if kk != ac_key][0]
assert abs(wv[dc_key]["path_length_top_mm"] - T["theory.geom.D_len_top_raw_mm"]) < 0.1
put("tab.peec.ac_label", ac_key, "PEEC winding label that is the AC winding (length match)")
for tag, blk in (("top", "inductance_top_layer_only"),
                 ("p2p", "inductance_pad_to_pad_incl_bottom_leads_and_vias_1.6mm_board")):
    b = PC[blk]
    put(f"tab.peec.{tag}.L_A_uH", b[f"L{ac_key}_uH"], f"values_peec_check.json {blk} L{ac_key} (AC)")
    put(f"tab.peec.{tag}.L_B_uH", b[f"L{dc_key}_uH"], f"values_peec_check.json {blk} L{dc_key} (DC)")
    put(f"tab.peec.{tag}.M_uH", b["M_uH"], f"values_peec_check.json {blk}")
    put(f"tab.peec.{tag}.k", b["k"], f"values_peec_check.json {blk}")
# agreement Neumann vs PEEC
pairs_top = [(out["tab.peec.top.L_A_uH"], top["L_A_uH"]), (out["tab.peec.top.L_B_uH"], top["L_D_uH"]),
             (out["tab.peec.top.M_uH"], top["M_uH"])]
pairs_p2p = [(out["tab.peec.p2p.L_A_uH"], T["theory.L_A_uH"]), (out["tab.peec.p2p.L_B_uH"], T["theory.L_D_uH"]),
             (out["tab.peec.p2p.M_uH"], T["theory.M_uH"])]
put("tab.peec.maxdev_top_pct", max(abs(100 * (a / b - 1)) for a, b in pairs_top), "max |PEEC/Neumann-1| over L_A, L_B, M, top layer")
put("tab.peec.maxdev_p2p_pct", max(abs(100 * (a / b - 1)) for a, b in pairs_p2p), "max |PEEC/Neumann-1| over L_A, L_B, M, pad to pad")
put("tab.peec.dk_top", abs(out["tab.peec.top.k"] - top["k"]), "|k_V2 - k_A2|, top layer")
put("tab.peec.dk_p2p", abs(out["tab.peec.p2p.k"] - T["theory.k"]), "|k_V2 - k_A2|, pad to pad")
put("tab.peec.R_A_p2p_ohm", wv[ac_key]["R_dc"]["w0.254mm_t35um"]["pad_to_pad_incl_bottom_lead_ohm"], "values_peec_check.json AC winding R pad-to-pad")
ic = PC["interwinding_capacitance"]
put("tab.peec.C_lo_pF", ic["estimate_pF_range"][0], "values_peec_check.json interwinding_capacitance")
put("tab.peec.C_hi_pF", ic["estimate_pF_range"][1], "values_peec_check.json interwinding_capacitance")
put("tab.peec.C_central_pF", ic["central_estimate_pF"], "values_peec_check.json interwinding_capacitance")

# ---------------------------------------------------------------- capacitance assumptions (cases)
cases = CR["fd_cases"]
put("tab.C.epsFR4_lo", min(c["eps_FR4"] for c in cases), "capacitance_results.json fd_cases")
put("tab.C.epsFR4_hi", max(c["eps_FR4"] for c in cases), "capacitance_results.json fd_cases")
mk = [c for c in cases if c["mask"] == "conformal"]
put("tab.C.mask_t_lo_um", min(c["mask_t_um"] for c in mk), "capacitance_results.json fd_cases (conformal mask)")
put("tab.C.mask_t_hi_um", max(c["mask_t_um"] for c in mk), "capacitance_results.json fd_cases (conformal mask)")
put("tab.C.mask_eps_lo", min(c["mask_eps"] for c in mk), "capacitance_results.json fd_cases (conformal mask)")
put("tab.C.mask_eps_hi", max(c["mask_eps"] for c in mk), "capacitance_results.json fd_cases (conformal mask)")
put("tab.C.Cgap_eta05_er44_pFm", T["theory.Cgap_conformal_pF_per_m_er4p4"], "theory.Cgap_conformal_pF_per_m_er4p4")
put("tab.C.Cnom_over3_pF", T["theory.C_AB_pF_nominal"] / 3, "theory.C_AB_pF_nominal / 3")

# ---------------------------------------------------------------- Section III-F: injector response
Rg, RL = float(HL["hl.Rg_ohm"]), float(HL["hl.RL_ohm"])
R = T["theory.R_A_w0p254_t35_ohm"]
LA, LB, Mt = T["theory.L_A_uH"] * 1e-6, T["theory.L_D_uH"] * 1e-6, T["theory.M_uH"] * 1e-6
k = T["theory.k"]
Ra, Rb = Rg + R, R + RL
put("tab.F.Ra_prime_ohm", Ra, "R_g + R_A (hl.Rg_ohm + theory.R_A_w0p254_t35_ohm)")
put("tab.F.Rb_prime_ohm", Rb, "R_B + R_L (theory.R_D... = R_A drawn; hl.RL_ohm)")
put("tab.F.pole_interaction_pct", 100 * Ra * Rb * (1 - k * k) / (Ra + Rb) ** 2, "R_A'R_B'(1-k^2)/(R_A'+R_B')^2")
# exact roots of Z_A Z_B + w^2 M^2 = 0 with unequal L_A, L_B (C = 0)
a2, a1, a0 = LA * LB - Mt * Mt, Ra * LB + Rb * LA, Ra * Rb
disc = math.sqrt(a1 * a1 - 4 * a2 * a0)
s1, s2 = (a1 - disc) / (2 * a2), (a1 + disc) / (2 * a2)
put("tab.F.f_lo_exact_kHz", s1 / (2 * math.pi) / 1e3, "exact real pole, eq. (H), unequal L_A, L_B, C = 0")
put("tab.F.f_hi_exact_MHz", s2 / (2 * math.pi) / 1e6, "exact real pole, eq. (H), unequal L_A, L_B, C = 0")


def H_node(f, C):
    """V(OUT)/V_g for the Fig. 2 circuit, optional capacitor C between AC_IN and OUT."""
    import numpy as np
    w = 2 * math.pi * f
    jw = 1j * w
    ZA0, ZB0 = R + jw * LA, R + jw * LB
    Y = jw * C
    # unknowns x = [V_in, V_out, I_A, I_B]; V_g = 1
    A = np.array([
        [1, 0, -ZA0, -jw * Mt],                # V_in = ZA0 I_A + jwM I_B
        [0, 1, -jw * Mt, -ZB0],                # V_out = ZB0 I_B + jwM I_A
        [1 / Rg + Y, -Y, 1, 0],                # KCL AC_IN
        [-Y, 1 / RL + Y, 0, 1],                # KCL OUT
    ], dtype=complex)
    b = np.array([0, 0, 1 / Rg, 0], dtype=complex)
    return np.linalg.solve(A, b)[1]


def corners3dB(C, fmax=2e7):
    """-3 dB points re the peak, first crossings on each side of the peak (as in
    theory_bench_model.py), sweep 100 Hz .. fmax (default 20 MHz, as the geometry scripts)."""
    import numpy as np
    f = np.logspace(2, math.log10(fmax), 6000)
    g = np.array([abs(H_node(x, C)) for x in f])
    i = int(np.argmax(g))
    tgt = g[i] / math.sqrt(2)
    j = np.where(g[:i] < tgt)[0][-1]
    lo = math.exp(np.interp(tgt, [g[j], g[j + 1]], [math.log(f[j]), math.log(f[j + 1])]))
    j = i + np.where(g[i:] < tgt)[0][0]
    hi = math.exp(np.interp(tgt, [g[j], g[j - 1]], [math.log(f[j]), math.log(f[j - 1])]))
    return lo, hi, 20 * math.log10(g[i]), f[i]


def notch(C, f1=1e6, f2=1e8):
    import numpy as np
    f = np.logspace(math.log10(f1), math.log10(f2), 8000)
    g = np.array([abs(H_node(x, C)) for x in f])
    i = int(np.argmin(g))
    return f[i], 20 * math.log10(g[i])


# sanity: node model equals eq. (H) at C = 0
_f = 1e6
_H = 1j * 2 * math.pi * _f * Mt * RL / ((Ra + 1j * 2 * math.pi * _f * LA) * (Rb + 1j * 2 * math.pi * _f * LB) + (2 * math.pi * _f * Mt) ** 2)
assert abs(H_node(_f, 0.0) - _H) < 1e-9
lo0, hi0, pk0, fpk0 = corners3dB(0.0)
loC, hiC, pkC, fpkC = corners3dB(out["tab.C.Cnom_over3_pF"] * 1e-12)
fn, gn = notch(out["tab.C.Cnom_over3_pF"] * 1e-12)
put("tab.F.notch_Cnom3_MHz", fn / 1e6, "|H| minimum (inductive/capacitive cancellation), C_AB,nom/3 AC_IN-OUT, 1-100 MHz")
put("tab.F.f_lo_3dB_C0_kHz", lo0 / 1e3, "numerical -3 dB re peak, C = 0")
put("tab.F.f_hi_3dB_C0_MHz", hi0 / 1e6, "numerical -3 dB re peak, C = 0")
put("tab.F.peak_C0_dB", pk0, "numerical peak |H|, C = 0")
put("tab.F.f_hi_3dB_Cnom3_MHz", hiC / 1e6, "numerical -3 dB re peak, C_AB,nom/3 between AC_IN and OUT")
put("tab.F.f_lo_3dB_Cnom3_kHz", loC / 1e3, "numerical -3 dB re peak, C_AB,nom/3 between AC_IN and OUT")

# terminal-referred (V_OUT / V_AC_IN): R_g removed from R_A'
Ra_t = R
a1t, a0t = Ra_t * LB + Rb * LA, Ra_t * Rb
disc_t = math.sqrt(a1t * a1t - 4 * a2 * a0t)
put("tab.F.term.f_lo_exact_kHz", (a1t - disc_t) / (2 * a2) / (2 * math.pi) / 1e3, "exact pole, V_OUT/V_AC_IN, C = 0")
put("tab.F.term.f_hi_exact_MHz", (a1t + disc_t) / (2 * a2) / (2 * math.pi) / 1e6, "exact pole, V_OUT/V_AC_IN, C = 0")
put("tab.F.term.Hmid_dB", 20 * math.log10(Mt * RL / a1t), "M R_L/(R_A L_B + R_B' L_A), dB re V(AC_IN)")

# ---------------------------------------------------------------- Section IV: other L records
recs = []
with open(os.path.join(D, "photo_index_oct24.csv")) as fh:
    for r in csv.DictReader(fh):
        n = int(r["id"].split("-")[1]) if r["id"].startswith("P24-") else -1
        if 53 <= n <= 59 and "LCR" in r["subject"]:
            mm = re.search(r"([0-9]+\.[0-9]+)\s*uH", r["notes"])
            if mm:
                recs.append(float(mm.group(1)))
assert len(recs) == 7, recs
put("tab.Lrec.photo_n", len(recs), "photo_index_oct24.csv P24-53..59 (LCR999A display, 1 kHz)")
put("tab.Lrec.photo_min_uH", min(recs), "photo_index_oct24.csv P24-53..59")
put("tab.Lrec.photo_max_uH", max(recs), "photo_index_oct24.csv P24-53..59 (18.9 read from an inverted display, 9/6 ambiguity)")
put("tab.Lrec.photo_n_18p6", sum(1 for x in recs if abs(x - 18.6) < 1e-9), "count of 18.6 uH photos")

# ---------------------------------------------------------------- Section IV: coupling derived
for f in ("1k", "10k"):
    La, Lb = Mv[f"meas.L_A_open_{f}"], Mv[f"meas.L_B_open_{f}"]
    g = math.sqrt(La * Lb)
    Ma5 = (Mv[f"meas.L_aid_c5_{f}"] - La - Lb) / 2
    Ma6 = (Mv[f"meas.L_aid_c6_{f}"] - La - Lb) / 2
    Mo7 = (La + Lb - Mv[f"meas.L_opp_c7_{f}"]) / 2
    Mo8 = (La + Lb - Mv[f"meas.L_opp_c8_{f}"]) / 2
    assert abs(Ma5 / g - Mv[f"meas.k_aid_only_c5_{f}"]) < 1e-4
    assert abs(Mo8 / g - Mv[f"meas.k_opp_only_c8_{f}"]) < 1e-4
    put(f"tab.k.M_aid_c5_{f}", Ma5, "(L_aid,c5 - L_A - L_B)/2")
    put(f"tab.k.M_aid_c6_{f}", Ma6, "(L_aid,c6 - L_A - L_B)/2")
    put(f"tab.k.M_aid_mean_{f}", (Ma5 + Ma6) / 2, "mean")
    put(f"tab.k.M_opp_c7_{f}", Mo7, "(L_A + L_B - L_opp,c7)/2")
    put(f"tab.k.M_opp_c8_{f}", Mo8, "(L_A + L_B - L_opp,c8)/2")
    put(f"tab.k.M_opp_mean_{f}", (Mo7 + Mo8) / 2, "mean")
    put(f"tab.k.L_aid_mean_{f}", (Mv[f"meas.L_aid_c5_{f}"] + Mv[f"meas.L_aid_c6_{f}"]) / 2, "mean of c5, c6")
    put(f"tab.k.L_opp_mean_{f}", (Mv[f"meas.L_opp_c7_{f}"] + Mv[f"meas.L_opp_c8_{f}"]) / 2, "mean of c7, c8")
put("tab.k.Lp_1k_mH", Mv["meas.Lp_if_parallel_1k_uH"] / 1e3, "meas.Lp_if_parallel_1k_uH in mH")
put("tab.k.Lp_min_1k_mH", Mv["meas.Lp_min_possible_1k_uH"] / 1e3, "meas.Lp_min_possible_1k_uH in mH")
sd = [abs(Mv[f"meas.sum_dev_pair{p_}_{f}_pct"]) for p_ in (1, 2) for f in ("1k", "10k")]
put("tab.k.sumdev_min_pct", min(sd), "smallest |sum-identity deviation|, pairs 1-2, 1 and 10 kHz")
put("tab.k.sumdev_max_pct", max(sd), "largest |sum-identity deviation|, pairs 1-2, 1 and 10 kHz")
put("tab.k.pair_mean_absdiff", abs(Mv["meas.k_pair_mean_diff_1k_10k"]), "|k_pair_mean(1k) - k_pair_mean(10k)|")
put("tab.k.w2LC_aid_10k_e5", Mv["meas.w2LC_aid_10k"] * 1e5, "meas.w2LC_aid_10k x 1e5")
put("tab.k.X_opp_min_ohm", min(Mv["meas.X_opp_1k_ohm"], Mv["meas.X_opp_10k_ohm"]), "opposing reactance, 1 kHz")
put("tab.k.L_opp_min_uH", min(Mv[f"meas.L_opp_c{c}_{f}"] for c in (7, 8) for f in ("1k", "10k")), "min opposing reading")
put("tab.k.L_opp_max_uH", max(Mv[f"meas.L_opp_c{c}_{f}"] for c in (7, 8) for f in ("1k", "10k")), "max opposing reading")
put("tab.k.dL_short_meas_lo_uH", min(Mv["meas.dL_short_meas_A_10k_uH"], Mv["meas.dL_short_meas_B_10k_uH"]), "")

# ---------------------------------------------------------------- Table V: theory vs measurement (10 kHz)
Lt = (T["theory.L_A_uH"] + T["theory.L_D_uH"]) / 2
Lm = (Mv["meas.L_A_open_10k"] + Mv["meas.L_B_open_10k"]) / 2
rows = {
    "L": (Lt, Lm),
    "M": (T["theory.M_uH"], Mv["meas.M_pair_mean_10k"]),
    "k": (T["theory.k"], Mv["meas.k_pair_mean_10k"]),
    "R": (T["theory.R_A_w0p254_t35_ohm"], Mv["meas.R_dc_ohm"]),
    "C": (T["theory.C_AB_pF_nominal"], Mv["meas.C_AB_pF"]),
    "Laid": (T["theory.L_series_aiding_uH"], out["tab.k.L_aid_mean_10k"]),
    "Lopp": (T["theory.L_series_opposing_uH"], out["tab.k.L_opp_mean_10k"]),
    "Lleak": (T["theory.L_leak_A_uH"], Lm * (1 - Mv["meas.k_pair_mean_10k"] ** 2)),
}
for name, (t, mm) in rows.items():
    put(f"tab.cmp.{name}_th", t, "theory (values_theory.json)")
    put(f"tab.cmp.{name}_meas", mm, "measured (values_meas.json, 10 kHz where applicable)")
    put(f"tab.cmp.{name}_dev_pct", 100 * (mm / t - 1), "measured/theory - 1")
put("tab.cmp.k_diff_abs", abs(Mv["meas.k_pair_mean_10k"] - T["theory.k"]), "|k_meas - k_theory|")
put("tab.cmp.L_dev_abs_pct", abs(out["tab.cmp.L_dev_pct"]), "absolute value")
put("tab.cmp.Lopp_dev_abs_pct", abs(out["tab.cmp.Lopp_dev_pct"]), "absolute value")
put("tab.cmp.C_dev_abs_pct", abs(out["tab.cmp.C_dev_pct"]), "absolute value")
put("tab.cmp.Lleak_dev_abs_pct", abs(out["tab.cmp.Lleak_dev_pct"]), "absolute value")
put("tab.cmp.M_dev_abs_pct", abs(out["tab.cmp.M_dev_pct"]), "absolute value")
put("tab.cmp.Laid_dev_abs_pct", abs(out["tab.cmp.Laid_dev_pct"]), "absolute value")
put("tab.cmp.k_diff", Mv["meas.k_pair_mean_10k"] - T["theory.k"], "k_meas - k_theory")
put("tab.cmp.L_minus_M_meas_uH", Lm - Mv["meas.M_pair_mean_10k"], "L_meas - M_meas (10 kHz)")
put("tab.cmp.dL_abs_uH", Lm - Lt, "L_meas - L_theory")
put("tab.cmp.C_dev_hi_pct", 100 * (Mv["meas.C_AB_pF"] / T["theory.C_AB_pF_hi"] - 1), "C_meas / C_theory,hi - 1")
put("tab.cmp.C_dev_V2hi_pct", 100 * (Mv["meas.C_AB_pF"] / ic["estimate_pF_range"][1] - 1), "C_meas / C_V2,hi - 1")
put("tab.cmp.Lrec_oct_minus_nov_lo_uH", Mv["meas.L_oct_minus_nov_min_uH"], "meas.L_oct_minus_nov_min_uH")

# ---------------------------------------------------------------- LCR999A published figures (not computed)
# Source: MECO catalogue sheet LCR999A, https://www.meco.in/uploads/products/catalog/LCR999A.pdf
# (opened 17 Nov 2026; see data/instrument_specs.md). The L-range lower limit's unit
# prefix is illegible in the extracted text: 4 uH is an inference (flagged as a gap in the paper).
SPEC = {"tab.spec.lcr.Vrms": 0.6, "tab.spec.lcr.Rsrc_ohm": 100, "tab.spec.lcr.acc_pct": 0.1,
        "tab.spec.lcr.Zmin_ohm": 0.4, "tab.spec.lcr.Lmin_uH_inferred": 4}
for kk, vv in SPEC.items():
    put(kk, vv, "MECO LCR999A catalogue sheet via data/instrument_specs.md (published figure, not computed)")

# ================================================================ additional derived values
# Fix 8: closed-form spread without the power-mean (R2) adaptation (citation not available)
cf2 = sorted(float(r["L_A_uH"]) for r in cf if r["adaptation"] != "R2_power_mean")
n2 = len(cf2)
med2 = cf2[n2 // 2] if n2 % 2 else (cf2[n2 // 2 - 1] + cf2[n2 // 2]) / 2
put("tab.cf2.n_variants", n2, "closed_form_table.csv rows excluding R2_power_mean")
put("tab.cf2.L_A_min_uH", cf2[0], "min, excluding R2_power_mean")
put("tab.cf2.L_A_max_uH", cf2[-1], "max, excluding R2_power_mean")
put("tab.cf2.L_A_median_uH", med2, "median, excluding R2_power_mean")
put("tab.cf2.dev_median_vs_neumann_pct", 100 * (med2 / T["theory.L_A_uH"] - 1), "median / Neumann L_A - 1")
# Fix 8: separately-written-code agreement computed from the values as displayed in Table II (0.01 uH)
disp = [(round(out["tab.peec.top.L_A_uH"], 2), round(top["L_A_uH"], 2)),
        (round(out["tab.peec.top.L_B_uH"], 2), round(top["L_D_uH"], 2)),
        (round(out["tab.peec.top.M_uH"], 2), round(top["M_uH"], 2))]
put("tab.peec.maxdev_top_disp_pct", max(abs(100 * (a / b - 1)) for a, b in disp),
    "max |PEEC/Neumann-1| over L_A, L_B, M (top layer), from the Table II values rounded to 0.01 uH")
# Fix 4: k aliases for tables (3 decimals) vs text (2 decimals)
put("tab.tbl.k_p2p", T["theory.k"], "theory.k (table alias, 3 decimals)")
put("tab.tbl.k_w020t30", T["theory.k_w0p20_t30"], "theory.k_w0p20_t30 (table alias)")
put("tab.txt.k_short_10k", Mv["meas.k_from_short_B_10k"], "meas.k_from_short_B_10k (text alias, 2 decimals)")
# Fix 4: effective resolution at 1 kHz
r1k = [Mv[k_] for k_ in ("meas.L_A_open_1k", "meas.L_A_shortB_1k", "meas.L_B_open_1k", "meas.L_B_shortA_1k",
                          "meas.L_aid_c5_1k", "meas.L_aid_c6_1k", "meas.L_opp_c7_1k", "meas.L_opp_c8_1k")]
put("tab.k.n_1k", len(r1k), "number of 1 kHz readings")
put("tab.k.n_1k_last0", sum(1 for x in r1k if round(x * 100) % 10 == 0), "1 kHz readings whose 0.01 uH digit is 0")
# Fix 4: dissipation factor of the opposing readings at 10 kHz
w10 = 2 * math.pi * 1e4
Rser = 2 * R_meas
Xo = [w10 * Mv[f"meas.L_opp_c{c}_10k"] * 1e-6 for c in (7, 8)]
put("tab.k.X_opp_10k_max_ohm", max(Xo), "omega L_opp, 10 kHz, largest opposing reading")
put("tab.k.D_opp_10k_lo", Rser / max(Xo), "D = 2R/(omega L_opp), 10 kHz, c7")
put("tab.k.D_opp_10k_hi", Rser / min(Xo), "D = 2R/(omega L_opp), 10 kHz, c8")
put("tab.k.Z_opp_10k_ohm", math.hypot(Rser, sum(Xo) / 2), "|Z| of the opposing connection, 10 kHz")
# Fix 5: direct-method range = span of the series-connection estimates listed in Table IV(b)
direct = [Mv[k_] for f in ("1k", "10k") for k_ in (
    f"meas.k_pair1_{f}", f"meas.k_pair2_{f}", f"meas.k_pair_mean_{f}", f"meas.k_aid_only_c5_{f}",
    f"meas.k_aid_only_c6_{f}", f"meas.k_aid_only_mean_{f}", f"meas.k_opp_only_mean_{f}")]
put("tab.k.direct_lo", min(direct), "min of the direct (series-connection) estimates in Table IV(b)")
put("tab.k.direct_hi", max(direct), "max of the direct (series-connection) estimates in Table IV(b)")
# Fix 10: data-driven sensitivities of the pair-method k at 10 kHz
LA10, LB10 = Mv["meas.L_A_open_10k"], Mv["meas.L_B_open_10k"]
g10 = math.sqrt(LA10 * LB10)
k10 = Mv["meas.k_pair_mean_10k"]
put("tab.k.dk_dsingle_per_uH", -k10 * (1 / (2 * LA10) + 1 / (2 * LB10)), "dk/d(delta) for a common offset on L_A and L_B, per uH")
put("tab.k.dk_dLopp_per_uH", -1 / (4 * g10), "dk/dL_opp per uH (pair method)")
put("tab.k.dk_dLaid_per_uH", 1 / (4 * g10), "dk/dL_aid per uH (pair method)")
# Fix 2: split of k_meas - k_theory by substituting measured inputs one at a time (10 kHz)
gth = math.sqrt(T["theory.L_A_uH"] * T["theory.L_D_uH"])
Lai_t, Lop_t = T["theory.L_series_aiding_uH"], T["theory.L_series_opposing_uH"]
Lai_m, Lop_m = out["tab.k.L_aid_mean_10k"], out["tab.k.L_opp_mean_10k"]
k_pt = (Lai_t - Lop_t) / (4 * gth)
put("tab.split.k_single", (Lai_t - Lop_t) / (4 * g10) - k_pt, "k change from measured L_A, L_B only")
put("tab.split.k_opp", (Lai_t - Lop_m) / (4 * gth) - k_pt, "k change from measured L_opp only")
put("tab.split.k_aid", (Lai_m - Lop_t) / (4 * gth) - k_pt, "k change from measured L_aid only")
for kk in ("single", "opp", "aid"):
    put(f"tab.split.k_{kk}_abs", abs(out[f"tab.split.k_{kk}"]), "absolute value")
# Fix 3: common-offset identity (L_aid+L_opp) - 2(L_A+L_B) = -2 delta, per pair; delta-corrected k
for f in ("1k", "10k"):
    La, Lb = Mv[f"meas.L_A_open_{f}"], Mv[f"meas.L_B_open_{f}"]
    for p_, (ca, co) in (("1", (5, 7)), ("2", (6, 8))):
        Laid, Lopp = Mv[f"meas.L_aid_c{ca}_{f}"], Mv[f"meas.L_opp_c{co}_{f}"]
        dlt = -0.5 * ((Laid + Lopp) - 2 * (La + Lb))
        Mp = (Laid - Lopp) / 4
        put(f"tab.delta.pair{p_}_{f}_uH", dlt, "-1/2[(L_aid+L_opp)-2(L_A+L_B)]")
        put(f"tab.delta.k_corr_pair{p_}_{f}", Mp / math.sqrt((La - dlt) * (Lb - dlt)), "M/sqrt((L_A-delta)(L_B-delta))")
# Fix 6: reproducibility floor
put("tab.floor.L_nov_uH", Lm, "mean open-circuit L, 10 kHz, November sheet")
# Fix 10: quarter-wave estimate with the measured C_AB
Cm = Mv["meas.C_AB_pF"] * 1e-12
put("tab.srf.qw_meas_Lopp_Cmeas_MHz", 1 / (4 * math.sqrt(Lop_m * 1e-6 * Cm)) / 1e6, "1/(4 sqrt(L_opp,meas C_AB,meas))")
put("tab.srf.qw_th_Lopp_Cmeas_MHz", 1 / (4 * math.sqrt(Lop_t * 1e-6 * Cm)) / 1e6, "1/(4 sqrt(L_opp,theory C_AB,meas))")

res = {k_: {"value": v_, "src": meta[k_]} for k_, v_ in out.items()}
res["_meta"] = {"producer": "analysis/table_values.py", "date": "2026-11-17"}
json.dump(res, open(os.path.join(D, "values_tables.json"), "w"), indent=1, sort_keys=True)
for k_, v_ in sorted(out.items()):
    print(f"{k_} = {v_}")
