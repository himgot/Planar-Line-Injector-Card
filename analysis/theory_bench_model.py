#!/usr/bin/env python3
"""DC resistance (task B), resonance estimates and low-frequency model of the bench
circuit (task F).

Run:  python3 -I analysis/theory_bench_model.py
Needs work/theory/geometry.json, neumann_results.json, capacitance_results.json.

DCR: R = rho * l / (w * t), rho_Cu = 1.724e-8 Ohm m at 20 C (International Annealed
Copper Standard, IACS: 1/58 Ohm mm^2/m = 0.017241 Ohm mm^2/m); temperature coefficient
0.00393 /K (IACS) used only for the 25 C column.

Bench circuit (assumed topology, see report):
  generator EMF Vg with Rg = 50 Ohm -> AC_IN (inner, undotted) ; AC_RTN (outer, dotted) = ground
  DC_IN (outer, dotted) = AC ground (ideal DC supply) ; OUT (inner, undotted) -> R_L -> ground
  R_L = 111.5 Ohm ; windings R_A, R_D (DCR), L_A, L_D, M (Neumann).
  Optional inter-winding capacitance lumped as C_AB/3 between AC_IN and OUT
  (energy-equivalent of the linear voltage-difference profile with both dotted ends grounded).
  Equations (dot convention, i_a into AC_IN, i_d from DC_IN to OUT):
     V1 = Z_A i_a - jwM i_d ;  V2 = -Z_D i_d + jwM i_a
     (Vg - V1)/Rg = i_a + jwC (V1 - V2) ;  i_d = V2/R_L + jwC (V2 - V1)
With C = 0:  V2/Vg = s M R_L / [s^2 L_A L_D (1-k^2) + s (R1 L_D + R2 L_A) + R1 R2],
             R1 = Rg + R_A, R2 = R_L + R_D.
"""
import csv
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gerber_geometry as gg  # noqa: E402

RHO20 = 1.724e-8
ALPHA_CU = 0.00393
RG = 50.0
RL = 111.5


def dcr(length_m, w, t, T=20.0):
    return RHO20 * (1 + ALPHA_CU * (T - 20)) * length_m / (w * t)


def solve(f, LA, LD, M, RA, RD, C=0.0, Rg=RG, RL_=RL):
    w = 2 * np.pi * f
    ZA = RA + 1j * w * LA
    ZD = RD + 1j * w * LD
    H_emf = np.zeros(len(f), complex)
    H_term = np.zeros(len(f), complex)
    for n, wn in enumerate(w):
        # unknowns x = [V1, V2, ia, id]
        jwM = 1j * wn * M
        jwC = 1j * wn * C
        A = np.array([
            [1, 0, -ZA[n], jwM],
            [0, 1, -jwM, ZD[n]],
            [1 / Rg + jwC, -jwC, 1, 0],
            [-jwC, 1 / RL_ + jwC, 0, -1]], complex)
        b = np.array([0, 0, 1 / Rg, 0], complex)       # Vg = 1
        V1, V2, ia, id_ = np.linalg.solve(A, b)
        H_emf[n] = V2
        H_term[n] = V2 / V1
    return H_emf, H_term


def corners(f, H):
    m = np.abs(H)
    i = int(np.argmax(m))
    lvl = m[i] / math.sqrt(2)
    lo = hi = float("nan")
    below = np.where(m[:i] < lvl)[0]
    if len(below):
        j = below[-1]
        lo = float(np.exp(np.interp(lvl, [m[j], m[j + 1]], [np.log(f[j]), np.log(f[j + 1])])))
    above = np.where(m[i:] < lvl)[0]
    if len(above):
        j = i + above[0]
        hi = float(np.exp(np.interp(lvl, [m[j], m[j - 1]], [np.log(f[j]), np.log(f[j - 1])])))
    return lo, hi, float(m[i]), float(f[i])


def analytic(LA, LD, M, RA, RD, ref="emf"):
    R1 = (RG + RA) if ref == "emf" else RA
    R2 = RL + RD
    a = LA * LD * (1 - M * M / (LA * LD))
    b = R1 * LD + R2 * LA
    c = R1 * R2
    disc = b * b - 4 * a * c
    p1 = (-b + math.sqrt(disc)) / (2 * a)
    p2 = (-b - math.sqrt(disc)) / (2 * a)
    return {"f_pole_lo_Hz": abs(p1) / (2 * math.pi), "f_pole_hi_Hz": abs(p2) / (2 * math.pi),
            "midband_gain": M * RL / b, "midband_gain_dB": 20 * math.log10(M * RL / b)}


def main():
    G = gg.load_geometry()
    NR = json.load(open(os.path.join(gg.WORK, "neumann_results.json")))
    CR = json.load(open(os.path.join(gg.WORK, "capacitance_results.json")))
    out = {}
    # ---------------- B: DC resistance
    lens = {"A": G["windings"]["AC"]["top_len_dp1.0mil_mm"] + G["windings"]["AC"]["bottom_len_mm"],
            "D": G["windings"]["DC"]["top_len_dp1.0mil_mm"] + G["windings"]["DC"]["bottom_len_mm"]}
    R = {}
    rows = []
    for key, Lmm in lens.items():
        for w in (0.254e-3, 0.20e-3):
            for t in (35e-6, 30e-6):
                r20 = dcr(Lmm * 1e-3, w, t)
                r25 = dcr(Lmm * 1e-3, w, t, 25.0)
                R[(key, round(w * 1e3, 3), round(t * 1e6))] = r20
                rows.append([key, Lmm, w * 1e3, t * 1e6, r20, r25])
    # via barrel (plating thickness assumed 20 um) for completeness
    r_via = RHO20 * 1.6e-3 / (math.pi * 0.3e-3 * 20e-6)
    out["dcr_rows"] = rows
    out["r_via_ohm_assumed_20um_plating"] = r_via
    # ---------------- inductances
    base = NR["variants"]["baseline"]
    sens = {(round(r["w_mm"], 3), round(r["t_um"])): r for r in NR["sensitivity_w_t"]}
    sets = {
        "drawn_w0.254_t35": (sens[(0.254, 35)], R[("A", 0.254, 35)], R[("D", 0.254, 35)]),
        "etched_w0.20_t30": (sens[(0.2, 30)], R[("A", 0.2, 30)], R[("D", 0.2, 30)]),
    }
    Cnom = None
    for c in CR["fd_cases"]:
        if abs(c["w_mm"] - 0.254) < 1e-6 and abs(c["t_um"] - 35) < 1e-6 and c["eps_FR4"] == 4.5 \
                and c["mask"] == "conformal" and c["mask_t_um"] == 10:
            Cnom = c["C_AB_pF"] * 1e-12
    # v2.2: the board is mask-covered -> range = conformal-mask minimum .. gap-filled bound
    Clo = min(c["C_AB_pF"] for c in CR["fd_cases"] if c["mask"] == "conformal") * 1e-12
    Chi = max(c["C_AB_pF"] for c in CR["fd_cases"] if c["mask"] == "filled") * 1e-12
    # ---------------- resonance estimates
    LA = base["L_A_uH"] * 1e-6
    Lopp = base["L_opp_uH"] * 1e-6
    ell = lens["A"] * 1e-3
    eps_eff = CR["summary"]["eps_eff_range"]
    srf = {
        "f_LA_with_CAB_lo_C_Hz": 1 / (2 * math.pi * math.sqrt(LA * Chi)),
        "f_LA_with_CAB_hi_C_Hz": 1 / (2 * math.pi * math.sqrt(LA * Clo)),
        "f_LA_with_CAB_over3_lo_Hz": 1 / (2 * math.pi * math.sqrt(LA * Chi / 3)),
        "f_LA_with_CAB_over3_hi_Hz": 1 / (2 * math.pi * math.sqrt(LA * Clo / 3)),
        "f_quarterwave_Lopp_CAB_lo_Hz": 1 / (4 * math.sqrt(Lopp * Chi)),
        "f_quarterwave_Lopp_CAB_hi_Hz": 1 / (4 * math.sqrt(Lopp * Clo)),
        "f_quarterwave_v_over_4l_lo_Hz": 299792458.0 / math.sqrt(eps_eff[1]) / (4 * ell),
        "f_quarterwave_v_over_4l_hi_Hz": 299792458.0 / math.sqrt(eps_eff[0]) / (4 * ell),
        "C_used_pF": [Clo * 1e12, Chi * 1e12], "L_A_uH": LA * 1e6, "L_opp_uH": Lopp * 1e6, "line_length_m": ell,
    }
    # D floating (bifilar tracking): energy-equivalent capacitance across A from
    # (i) half-turn offset between neighbours, C_AB/(4 N^2), and (ii) imperfect tracking (1-k)^2 C_AB/3
    N = G["windings"]["AC"]["turns"]
    kb = base["k"]
    cfl = Cnom / (4 * N * N) + (1 - kb) ** 2 * Cnom / 3
    srf["C_eff_A_with_D_floating_pF"] = cfl * 1e12
    srf["f_LA_C_eff_floating_Hz"] = 1 / (2 * math.pi * math.sqrt(LA * cfl))
    srf["f_halfwave_Lopp_CAB_lo_Hz"] = 2 * srf["f_quarterwave_Lopp_CAB_lo_Hz"]
    srf["f_halfwave_Lopp_CAB_hi_Hz"] = 2 * srf["f_quarterwave_Lopp_CAB_hi_Hz"]
    out["resonance"] = srf
    # ---------------- bench LF model
    f = np.logspace(2, math.log10(2e7), 600)
    bench = {}
    curves = {"f_Hz": f}
    for name, (Ls, RA, RD) in sets.items():
        LAx, LDx, Mx = Ls["L_A_uH"] * 1e-6, Ls["L_D_uH"] * 1e-6, Ls["M_uH"] * 1e-6
        d = {"R_A_ohm": RA, "R_D_ohm": RD, "L_A_uH": Ls["L_A_uH"], "L_D_uH": Ls["L_D_uH"],
             "M_uH": Ls["M_uH"], "k": Ls["k"]}
        d["analytic_emf"] = analytic(LAx, LDx, Mx, RA, RD, "emf")
        d["analytic_terminal"] = analytic(LAx, LDx, Mx, RA, RD, "terminal")
        for cname, C in (("C0", 0.0), ("Cnom_over3", Cnom / 3), ("Chi_over3", Chi / 3)):
            He, Ht = solve(f, LAx, LDx, Mx, RA, RD, C)
            lo, hi, pk, fpk = corners(f, He)
            lo2, hi2, pk2, fpk2 = corners(f, Ht)
            d["numeric_%s" % cname] = {"emf": {"f_lo_3dB_Hz": lo, "f_hi_3dB_Hz": hi, "peak_gain": pk,
                                                "peak_gain_dB": 20 * math.log10(pk), "f_peak_Hz": fpk},
                                        "terminal": {"f_lo_3dB_Hz": lo2, "f_hi_3dB_Hz": hi2, "peak_gain": pk2,
                                                     "peak_gain_dB": 20 * math.log10(pk2), "f_peak_Hz": fpk2},
                                        "C_F": C}
            curves["%s_%s_emf_dB" % (name, cname)] = 20 * np.log10(np.abs(He))
            curves["%s_%s_term_dB" % (name, cname)] = 20 * np.log10(np.abs(Ht))
            curves["%s_%s_emf_phase_deg" % (name, cname)] = np.degrees(np.angle(He))
        bench[name] = d
    out["bench"] = bench
    with open(os.path.join(gg.WORK, "bench_model_curve.csv"), "w", newline="") as fh:
        cw = csv.writer(fh)
        keys = list(curves)
        cw.writerow(keys)
        for i in range(len(f)):
            cw.writerow(["%.6g" % curves[k][i] for k in keys])
    with open(os.path.join(gg.WORK, "bench_model_results.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float)[:8000])
    # ---------------- values
    rho_meta = "rho_Cu = 1.724e-8 Ohm m at 20 C (IACS annealed copper)"
    e = {}
    for key in ("A", "D"):
        for (w, t) in ((0.254, 35), (0.254, 30), (0.2, 35), (0.2, 30)):
            k = "theory.R_%s_w%s_t%d_ohm" % (key, ("0p254" if w == 0.254 else "0p20"), t)
            e[k] = (gg.sig(R[(key, w, t)], 3), "[CALCULATED] theory_bench_model.py: R = rho l/(w t), l = %.1f mm "
                    "(top DP1.0 mil + bottom route, gerber_geometry.py), w = %.3f mm, t = %d um, %s, 20 C" %
                    (lens[key], w, t, rho_meta))
    m = ("[MODEL] theory_bench_model.py: coupled-inductor LF model, Rg 50 Ohm, R_L 111.5 Ohm, ideal DC supply "
         "(AC ground), AC_RTN grounded; L, M from inductance_neumann.py; R = DCR; no capacitance unless stated")
    for name, tag in (("drawn_w0.254_t35", "drawn"), ("etched_w0.20_t30", "w0p20t30")):
        d = bench[name]
        e["theory.bench.%s.f_corner_lo_emf_kHz" % tag] = (gg.sig(d["analytic_emf"]["f_pole_lo_Hz"] / 1e3, 3), m + "; lower pole, V_OUT/V_emf")
        e["theory.bench.%s.f_corner_hi_emf_MHz" % tag] = (gg.sig(d["analytic_emf"]["f_pole_hi_Hz"] / 1e6, 3), m + "; upper (leakage) pole, V_OUT/V_emf")
        e["theory.bench.%s.midband_emf_dB" % tag] = (gg.sig(d["analytic_emf"]["midband_gain_dB"], 3), m + "; M R_L/(R1 L_D + R2 L_A), dB re generator EMF")
        e["theory.bench.%s.f_corner_lo_term_kHz" % tag] = (gg.sig(d["analytic_terminal"]["f_pole_lo_Hz"] / 1e3, 3), m + "; lower pole, V_OUT/V_AC_IN (generator terminal)")
        e["theory.bench.%s.f_corner_hi_term_MHz" % tag] = (gg.sig(d["analytic_terminal"]["f_pole_hi_Hz"] / 1e6, 3), m + "; upper pole, V_OUT/V_AC_IN")
        e["theory.bench.%s.midband_term_dB" % tag] = (gg.sig(d["analytic_terminal"]["midband_gain_dB"], 3), m + "; dB re V_AC_IN")
        e["theory.bench.%s.f_3dB_lo_emf_Cnom_kHz" % tag] = (gg.sig(d["numeric_Cnom_over3"]["emf"]["f_lo_3dB_Hz"] / 1e3, 3), m + "; numeric -3 dB, with C_AB/3 (nominal) AC_IN-OUT")
    rs = srf
    e["theory.srf.f_LA_CAB_MHz_lo"] = (gg.sig(rs["f_LA_with_CAB_lo_C_Hz"] / 1e6, 2), "[MODEL] 1/(2 pi sqrt(L_A C)), C = C_AB_hi (lumped, pessimistic)")
    e["theory.srf.f_LA_CAB_MHz_hi"] = (gg.sig(rs["f_LA_with_CAB_hi_C_Hz"] / 1e6, 2), "[MODEL] 1/(2 pi sqrt(L_A C)), C = C_AB_lo")
    e["theory.srf.f_LA_CABover3_MHz_lo"] = (gg.sig(rs["f_LA_with_CAB_over3_lo_Hz"] / 1e6, 2), "[MODEL] 1/(2 pi sqrt(L_A C_AB/3)), C_AB_hi")
    e["theory.srf.f_LA_CABover3_MHz_hi"] = (gg.sig(rs["f_LA_with_CAB_over3_hi_Hz"] / 1e6, 2), "[MODEL] 1/(2 pi sqrt(L_A C_AB/3)), C_AB_lo")
    e["theory.srf.f_quarterwave_MHz_lo"] = (gg.sig(rs["f_quarterwave_Lopp_CAB_lo_Hz"] / 1e6, 2), "[MODEL] bifilar line lambda/4: 1/(4 sqrt(L_opp C_AB)), C_AB_hi")
    e["theory.srf.f_quarterwave_MHz_hi"] = (gg.sig(rs["f_quarterwave_Lopp_CAB_hi_Hz"] / 1e6, 2), "[MODEL] bifilar line lambda/4: 1/(4 sqrt(L_opp C_AB)), C_AB_lo")
    e["theory.srf.C_eff_A_D_floating_pF"] = (gg.sig(rs["C_eff_A_with_D_floating_pF"], 2), "[MODEL] C_AB_nominal/(4N^2) + (1-k)^2 C_AB_nominal/3 (energy argument, linear profiles)")
    e["theory.srf.f_LA_Ceff_floating_MHz"] = (gg.sig(rs["f_LA_C_eff_floating_Hz"] / 1e6, 2), "[MODEL] 1/(2 pi sqrt(L_A C_eff)); beyond the lambda/4 estimate, so not meaningful as a lumped SRF")
    e["theory.rho_Cu_ohm_m"] = (RHO20, "[DESIGN] IACS annealed copper resistivity at 20 C (1/58 Ohm mm^2/m)")
    gg.update_values(e)


if __name__ == "__main__":
    main()
