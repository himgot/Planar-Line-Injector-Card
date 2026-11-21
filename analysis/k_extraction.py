#!/usr/bin/env python3
"""LCR records, coupling (M, k) extraction, Q / Ls-Lp check, shorted-winding analysis,
lead-inductance (delta) fit, uncertainty budget and DC-drop numbers.

Run from the package root:
    python3 -I analysis/k_extraction.py

Reads:
    data/lcr_sheet_2026-11-17.csv   (the 16 readings of the handwritten LCR sheet, as written)
    data/other_meter_records.csv    (multimeter readings and earlier inductance notes)
Writes:
    data/lcr_raw.csv, data/k_table.csv, data/values_meas.json

Conventions
    Winding A = AC winding (AC_IN - AC_RTN). Winding B = DC winding (DC_IN - OUT).
    All L in uH, R in ohm, f in Hz. Labels: [MEASURED] = as recorded by the author,
    [CALCULATED] = arithmetic in this script, [MODEL] = fitted / assumed-model result.
"""
import csv
import json
import math
import re
import sys

import numpy as np
from scipy.optimize import least_squares

OUT = "."
SHEET = "data/lcr_sheet_2026-11-17.csv"
OTHER = "data/other_meter_records.csv"

SRC_SHEET = "handwritten LCR sheet, photographed 17 Nov 2026"
R_DC = 7.31            # ohm per winding [MEASURED], Keysight U1282A
C_AB_PF = 274.0        # pF interwinding [MEASURED], Keysight U1282A, leads nulled
R_LOAD = 111.6         # ohm, Keysight U1282A
R_LOAD_ALT = 111.5
V_RAIL = 5.0
LSD = 0.01             # uH, last digit written at 10 kHz (the 1 kHz readings were written to 0.1 uH)
BASIC_ACC = 0.001      # 0.1 % MECO LCR999A catalogue "Basic Accuracy" (best case)
L_RANGE_MIN = 4.0      # uH, catalogue "L 4 [u]H ~ 1000H" (mu glyph lost in extraction)
Z_RANGE_MIN = 0.4      # ohm, catalogue "R / Z 0.4 ~ 10M"

V = {}                 # flat values dict -> values_meas.json


def put(key, value, meta, nd=None):
    if isinstance(value, float) and nd is not None:
        value = round(value, nd)
    V[key] = value
    V[key + "_meta"] = meta


def read_text(p):
    with open(p, "r", encoding="utf-8") as fh:
        return fh.read()


# --------------------------------------------------------------------------------------
# 1. Readings of the handwritten LCR sheet (MECO LCR999A, 1 kHz and 10 kHz)
# --------------------------------------------------------------------------------------
with open(SHEET, newline="", encoding="utf-8") as fh:
    sheet = list(csv.DictReader(fh))
if len(sheet) != 8:
    sys.exit(f"expected 8 connections in {SHEET}, got {len(sheet)}")
CONF = []
R = {1000: {}, 10000: {}}      # R[f][config_no] = value (uH)
for row in sheet:
    c = dict(no=int(row["config_no"]), expect=row["description"], link=row["link"], terms=row["measured_terminals"],
             wind=row["measured_winding"], other=row["other_winding"], kind=row["connection"],
             str1k=row["L_1kHz_uH"], str10k=row["L_10kHz_uH"])
    CONF.append(c)
    R[1000][c["no"]] = float(c["str1k"])
    R[10000][c["no"]] = float(c["str10k"])
CONF.sort(key=lambda c: c["no"])
with open(OTHER, newline="", encoding="utf-8") as fh:
    OTHERS = {r["record_id"]: r for r in csv.DictReader(fh)}

# --------------------------------------------------------------------------------------
# 2. lcr_raw.csv (all records in one table)
# --------------------------------------------------------------------------------------
raw_rows = []
for c in CONF:
    for f, s_, res in [(1000, c["str1k"], "0.1 uH (as written)"), (10000, c["str10k"], "0.01 uH (as written)")]:
        raw_rows.append(dict(
            record_id=f"C{c['no']}-{f//1000}k", quantity="L (Ls inferred; mode not recorded)",
            config_no=c["no"], link=c["link"], measured_terminals=c["terms"],
            measured_winding=c["wind"], configuration=c["kind"], other_winding=c["other"],
            f_Hz=f, value=s_, value_min="", value_max="", unit="uH", recorded_resolution=res,
            instrument="MECO LCR999A", date="Nov 2026 (sheet photographed 17 Nov 2026)",
            source=SRC_SHEET, label="[MEASURED]"))
for r in OTHERS.values():
    raw_rows.append(dict(
        record_id=r["record_id"], quantity=r["quantity"], config_no="", link="", measured_terminals=r["measured_terminals"],
        measured_winding=r["measured_winding"], configuration=r["condition"], other_winding="",
        f_Hz=r["f_Hz"] or "not recorded", value=r["value"], value_min=r["value_min"], value_max=r["value_max"],
        unit=r["unit"], recorded_resolution=r["resolution"] or "not recorded", instrument=r["instrument"],
        date=r["date"], source=r["source"], label="[MEASURED]"))
with open(f"{OUT}/data/lcr_raw.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(raw_rows[0].keys()))
    w.writeheader()
    w.writerows(raw_rows)

# raw values into JSON
for f in (1000, 10000):
    tag = "1k" if f == 1000 else "10k"
    names = {1: "L_B_open", 2: "L_B_shortA", 3: "L_A_open", 4: "L_A_shortB",
             5: "L_aid_c5", 6: "L_aid_c6", 7: "L_opp_c7", 8: "L_opp_c8"}
    for n, nm in names.items():
        put(f"meas.{nm}_{tag}", R[f][n],
            f"[MEASURED] MECO LCR999A, {f} Hz, config {n} ({CONF[n-1]['expect']}); {SRC_SHEET}")
put("meas.R_dc_ohm", R_DC, "[MEASURED] per winding; Keysight U1282A; lead nulling not recorded")
put("meas.C_AB_pF", C_AB_PF, "[MEASURED] interwinding, leads nulled; Keysight U1282A")
put("meas.L_xlsx_note_uH", float(OTHERS["BENCH-NOTE"]["value"]), "[MEASURED] " + OTHERS["BENCH-NOTE"]["source"] + "; f, mode, other-winding state not recorded")
put("meas.L_first_test_uH", float(OTHERS["FIRST-TEST"]["value"]), "[MEASURED] first bench test; author's project notes; conditions not recorded")
put("meas.L_notes_oct_min_uH", float(OTHERS["OCT-NOTES"]["value_min"]), "[MEASURED] author's project notes, Oct 2026")
put("meas.L_notes_oct_max_uH", float(OTHERS["OCT-NOTES"]["value_max"]), "[MEASURED] author's project notes, Oct 2026")
for nm, rid in [("A_open", "NOTE-1KHZ-A-open"), ("A_shortB", "NOTE-1KHZ-A-short"), ("B_open", "NOTE-1KHZ-B-open"), ("B_shortA", "NOTE-1KHZ-B-short")]:
    put(f"meas.L_{nm}_1k_earlier", float(OTHERS[rid]["value"]), "[MEASURED] author's project notes, early Nov 2026 (0.1 uH resolution)")

# --------------------------------------------------------------------------------------
# 3. Coupling extraction
# --------------------------------------------------------------------------------------
ktab = []


def addk(f, method, variant, LA, LB, M, k, note):
    ktab.append(dict(f_Hz=f, method=method, variant=variant,
                     L_A_uH=None if LA is None else round(LA, 4),
                     L_B_uH=None if LB is None else round(LB, 4),
                     M_uH=None if M is None else round(M, 4),
                     k=None if k is None else round(k, 4), note=note))


EXT = {}
for f in (1000, 10000):
    tag = "1k" if f == 1000 else "10k"
    r = R[f]
    LA, LB = r[3], r[1]                 # open-other-winding readings
    G = math.sqrt(LA * LB)
    aid = {5: r[5], 6: r[6]}
    opp = {7: r[7], 8: r[8]}
    pairs = {"pair1(c5,c7)": (5, 7), "pair2(c6,c8)": (6, 8),
             "cross(c5,c8)": (5, 8), "cross(c6,c7)": (6, 7)}
    kp = {}
    for pn, (a, o) in pairs.items():
        M = (r[a] - r[o]) / 4.0
        k = M / G
        kp[pn] = (M, k)
        dev = (r[a] + r[o]) / (2 * (LA + LB)) - 1.0
        addk(f, "pair M=(Laid-Lopp)/4", pn, LA, LB, M, k,
             f"sum identity: Laid+Lopp={r[a]+r[o]:.2f} vs 2(LA+LB)={2*(LA+LB):.2f} uH ({100*dev:+.2f} %)")
    M1, k1 = kp["pair1(c5,c7)"]
    M2, k2 = kp["pair2(c6,c8)"]
    Mm, km = (M1 + M2) / 2, (k1 + k2) / 2
    addk(f, "pair M=(Laid-Lopp)/4", "mean of pair1,pair2", LA, LB, Mm, km, "headline method (standard series aiding/opposing)")
    # aiding-only and opposing-only
    ka = {a: (r[a] - LA - LB) / (2 * G) for a in aid}
    ko = {o: (LA + LB - r[o]) / (2 * G) for o in opp}
    for a, kk in ka.items():
        addk(f, "aiding-only k=(Laid-LA-LB)/(2sqrt(LA LB))", f"c{a}", LA, LB, kk * G, kk,
             "assumes opposing readings are the faulty ones")
    for o, kk in ko.items():
        addk(f, "opposing-only k=(LA+LB-Lopp)/(2sqrt(LA LB))", f"c{o}", LA, LB, kk * G, kk,
             "assumes aiding readings are the faulty ones")
    kam = sum(ka.values()) / 2
    kom = sum(ko.values()) / 2
    addk(f, "aiding-only", "mean", LA, LB, kam * G, kam, "")
    addk(f, "opposing-only", "mean", LA, LB, kom * G, kom, "")
    s1 = r[5] + r[7]
    s2 = r[6] + r[8]
    ideal = 2 * (LA + LB)
    put(f"meas.M_pair1_{tag}", M1, f"[CALCULATED] (L_c5-L_c7)/4 at {f} Hz ({SRC_SHEET})", 4)
    put(f"meas.M_pair2_{tag}", M2, f"[CALCULATED] (L_c6-L_c8)/4 at {f} Hz", 4)
    put(f"meas.M_pair_mean_{tag}", Mm, f"[CALCULATED] mean of pair1/pair2 M at {f} Hz", 4)
    put(f"meas.k_pair1_{tag}", k1, f"[CALCULATED] M_pair1/sqrt(L_A_open*L_B_open) at {f} Hz", 4)
    put(f"meas.k_pair2_{tag}", k2, f"[CALCULATED] M_pair2/sqrt(L_A_open*L_B_open) at {f} Hz", 4)
    put(f"meas.k_pair_mean_{tag}", km, f"[CALCULATED] mean of pair1/pair2 k at {f} Hz; inputs {SRC_SHEET}", 4)
    put(f"meas.k_cross_c5c8_{tag}", kp["cross(c5,c8)"][1], "[CALCULATED] cross-pairing (c5 aiding, c8 opposing)", 4)
    put(f"meas.k_cross_c6c7_{tag}", kp["cross(c6,c7)"][1], "[CALCULATED] cross-pairing (c6 aiding, c7 opposing)", 4)
    put(f"meas.k_aid_only_c5_{tag}", ka[5], "[CALCULATED] (L_c5-L_A-L_B)/(2 sqrt(L_A L_B))", 4)
    put(f"meas.k_aid_only_c6_{tag}", ka[6], "[CALCULATED] (L_c6-L_A-L_B)/(2 sqrt(L_A L_B))", 4)
    put(f"meas.k_aid_only_mean_{tag}", kam, "[CALCULATED] mean aiding-only k", 4)
    put(f"meas.M_aid_only_c6_{tag}", ka[6] * G, "[CALCULATED] (L_c6-L_A-L_B)/2", 4)
    put(f"meas.k_opp_only_c7_{tag}", ko[7], "[CALCULATED] (L_A+L_B-L_c7)/(2 sqrt(L_A L_B))", 4)
    put(f"meas.k_opp_only_c8_{tag}", ko[8], "[CALCULATED] (L_A+L_B-L_c8)/(2 sqrt(L_A L_B))", 4)
    put(f"meas.k_opp_only_mean_{tag}", kom, "[CALCULATED] mean opposing-only k", 4)
    put(f"meas.sum_ideal_{tag}_uH", ideal, "[CALCULATED] 2(L_A_open+L_B_open)", 3)
    put(f"meas.sum_pair1_{tag}_uH", s1, "[CALCULATED] L_c5+L_c7", 3)
    put(f"meas.sum_pair2_{tag}_uH", s2, "[CALCULATED] L_c6+L_c8", 3)
    put(f"meas.sum_dev_pair1_{tag}_pct", 100 * (s1 / ideal - 1), "[CALCULATED] (Laid+Lopp)/(2(LA+LB))-1, pair1", 2)
    put(f"meas.sum_dev_pair2_{tag}_pct", 100 * (s2 / ideal - 1), "[CALCULATED] same, pair2", 2)
    put(f"meas.aid_spread_{tag}_uH", abs(r[5] - r[6]), "[CALCULATED] |L_c5-L_c6|: two nominally identical aiding circuits", 3)
    put(f"meas.opp_spread_{tag}_uH", abs(r[7] - r[8]), "[CALCULATED] |L_c7-L_c8|: two nominally identical opposing circuits", 3)
    EXT[f] = dict(LA=LA, LB=LB, G=G, Mm=Mm, km=km, kam=kam, kom=kom, ka=ka, ko=ko, k1=k1, k2=k2)

# frequency agreement and polarity inference
put("meas.k_pair_mean_diff_1k_10k", EXT[10000]["km"] - EXT[1000]["km"], "[CALCULATED] k_pair_mean_10k - k_pair_mean_1k", 4)
put("meas.polarity_inference", "like-polarity (dotted) ends: DC_IN and AC_RTN (equivalently OUT and AC_IN)",
    "[CALCULATED] from which links gave the large (aiding) readings: c5 DC_IN<->AC_IN and c6 OUT<->AC_RTN aiding; c7/c8 opposing; all four consistent")

# --------------------------------------------------------------------------------------
# 4. Q and Ls / Lp
# --------------------------------------------------------------------------------------
for f in (1000, 10000):
    tag = "1k" if f == 1000 else "10k"
    w = 2 * math.pi * f
    L = (EXT[f]["LA"] + EXT[f]["LB"]) / 2 * 1e-6
    Q = w * L / R_DC
    Lp = L * (1 + 1 / Q**2)
    Lp_min = 2 * R_DC / w           # minimum of Ls + Rs^2/(w^2 Ls) over Ls (at w Ls = Rs)
    Qa = w * (sum(R[f][a] for a in (5, 6)) / 2e6) / (2 * R_DC)
    Qo = w * (sum(R[f][o] for o in (7, 8)) / 2e6) / (2 * R_DC)
    put(f"meas.Q_single_{tag}", Q, f"[CALCULATED] w*L/R, L=mean open reading, R={R_DC} ohm (dc value)", 5)
    put(f"meas.D_single_{tag}", 1 / Q, "[CALCULATED] D=1/Q", 2)
    put(f"meas.phase_single_{tag}_deg", math.degrees(math.atan(Q)), "[CALCULATED] atan(Q)", 3)
    put(f"meas.Q_aid_{tag}", Qa, "[CALCULATED] w*mean(L_aid)/(2R)", 5)
    put(f"meas.Q_opp_{tag}", Qo, "[CALCULATED] w*mean(L_opp)/(2R)", 5)
    put(f"meas.Lp_if_parallel_{tag}_uH", Lp * 1e6, "[CALCULATED] Lp=Ls(1+1/Q^2) that a parallel-mode reading would show", 1)
    put(f"meas.Lp_min_possible_{tag}_uH", Lp_min * 1e6, "[CALCULATED] min over Ls of Lp=Ls+Rs^2/(w^2 Ls) = 2Rs/w; any Lp reading must exceed this", 1)
    put(f"meas.XL_single_{tag}_ohm", w * L, "[CALCULATED] w*L (single winding reactance)", 4)
    put(f"meas.X_opp_{tag}_ohm", w * sum(R[f][o] for o in (7, 8)) / 2e6, "[CALCULATED] w*mean(L_opp)", 4)

# --------------------------------------------------------------------------------------
# 5. Shorted-secondary analysis
# --------------------------------------------------------------------------------------


def dL_short(f, k, L1, L2, R2):
    w = 2 * math.pi * f
    M2 = k * k * L1 * L2
    return w * w * M2 * L2 / (R2 * R2 + w * w * L2 * L2)          # units follow L (uH ok: ratio form)


def dL_short_SI(f, k, L1, L2, R2):
    w = 2 * math.pi * f
    L1s, L2s = L1 * 1e-6, L2 * 1e-6
    return 1e6 * (w * w * k * k * L1s * L2s * L2s / (R2 * R2 + w * w * L2s * L2s))


def k_from_dL(f, dL, L1, L2, R2):
    w = 2 * math.pi * f
    L1s, L2s, dLs = L1 * 1e-6, L2 * 1e-6, dL * 1e-6
    return math.sqrt(dLs * (R2 * R2 + w * w * L2s * L2s) / (w * w * L1s * L2s * L2s))


for f in (1000, 10000):
    tag = "1k" if f == 1000 else "10k"
    r = R[f]
    kk = EXT[f]["km"]
    # DC (B) measured with AC (A) shorted: L1=L_B, L2=L_A ; AC measured with DC shorted: L1=L_A, L2=L_B
    for nm, (op, sh, L1, L2) in {"B": (1, 2, r[1], r[3]), "A": (3, 4, r[3], r[1])}.items():
        meas = r[op] - r[sh]
        pred = dL_short_SI(f, kk, L1, L2, R_DC)
        pred_k1 = dL_short_SI(f, 1.0, L1, L2, R_DC)
        put(f"meas.dL_short_meas_{nm}_{tag}_uH", meas, f"[CALCULATED] L_open-L_short for winding {nm} at {f} Hz (readings {SRC_SHEET})", 3)
        put(f"meas.dL_short_pred_{nm}_{tag}_uH", pred, f"[MODEL] w^2 M^2 L2/(R2^2+w^2 L2^2), k=k_pair_mean_{tag}, R2={R_DC} ohm", 4)
        put(f"meas.dL_short_pred_k1_{nm}_{tag}_uH", pred_k1, "[MODEL] same at k=1 (upper bound)", 4)
        tb = L1 * (1 - kk * kk)
        put(f"meas.L_short_textbook_{nm}_{tag}_uH", tb, "[MODEL] lossless textbook L1(1-k^2) at k_pair_mean (does NOT apply: Q2<<1)", 3)
        k_tb = math.sqrt(max(0.0, 1 - r[sh] / r[op]))
        put(f"meas.k_textbook_misuse_{nm}_{tag}", k_tb, "[CALCULATED] sqrt(1-Lshort/Lopen): invalid estimator at Q2<<1, shown only to warn", 4)
        if f == 10000:
            k0 = k_from_dL(f, meas, L1, L2, R_DC)
            put(f"meas.k_from_short_{nm}_10k", k0, "[MODEL] k implied by measured dL alone, R2=7.31 ohm, L from open readings", 4)
            for d, lab in [(-0.01, "m001"), (0.01, "p001"), (-0.02, "m002"), (0.02, "p002")]:
                put(f"meas.k_from_short_{nm}_10k_dL{lab}", k_from_dL(f, meas + d, L1, L2, R_DC),
                    f"[MODEL] same with dL {d:+.2f} uH (resolution sensitivity)", 4)
            for Rv, lab in [(7.01, "R7p01"), (7.61, "R7p61")]:
                put(f"meas.k_from_short_{nm}_10k_{lab}", k_from_dL(f, meas, L1, L2, Rv),
                    f"[MODEL] same with R2={Rv} ohm (R2 sensitivity, +/-0.3 ohm)", 4)
            # M implied by shorted test
            w = 2 * math.pi * f
            Mshort = math.sqrt(meas * 1e-6 * (R_DC**2 + (w * L2 * 1e-6) ** 2) / (w * w * L2 * 1e-6)) * 1e6
            put(f"meas.M_from_short_{nm}_10k_uH", Mshort, "[MODEL] M implied by dL alone (delta-immune)", 3)
            # R2 needed for consistency with pair M
            Mp = EXT[f]["Mm"] * 1e-6
            R2need = math.sqrt(max(0.0, w * w * Mp * Mp * L2 * 1e-6 / (meas * 1e-6) - (w * L2 * 1e-6) ** 2))
            put(f"meas.R2_for_pairM_consistency_{nm}_10k_ohm", R2need, "[MODEL] R2 that would make shorted-test M equal pair-mean M", 3)
            dLneed = dL_short_SI(f, EXT[f]["km"], L1, L2, R_DC)
            put(f"meas.dL_gap_short_vs_pair_{nm}_10k_uH", meas - dLneed, "[CALCULATED] measured dL minus dL predicted at k_pair_mean", 3)
    # attenuation factor of reflected reactance
    w = 2 * math.pi * f
    Q2 = w * r[3] * 1e-6 / R_DC
    put(f"meas.reflect_factor_{tag}", Q2**2 / (1 + Q2**2), "[MODEL] Q2^2/(1+Q2^2): fraction of the ideal k^2 L drop that survives with R2", 6)
    Rrefl = (w * EXT[f]["Mm"] * 1e-6) ** 2 * R_DC / (R_DC**2 + (w * r[3] * 1e-6) ** 2)
    put(f"meas.R_reflected_{tag}_ohm", Rrefl, "[MODEL] w^2 M^2 R2/(R2^2+w^2 L2^2) added to primary Rs when secondary shorted (not recorded on the sheet)", 4)
put("meas.f_Q1_kHz", R_DC / (2 * math.pi * 17.95e-6) / 1e3, "[CALCULATED] f where w L = R (Q=1), L=17.95 uH, R=7.31 ohm", 1)
# earlier 1 kHz notes: 17.9 -> 17.8 is 0.1 uH, vs k=1 maximum
put("meas.dL_short_earlier_A_1k_uH", 0.1, "[CALCULATED] 17.9-17.8 (earlier project notes); exceeds k=1 physical maximum -> reading scatter, not coupling")

# --------------------------------------------------------------------------------------
# 6. Lead / fixture inductance (delta) least-squares fit on 6 two-terminal readings
# --------------------------------------------------------------------------------------
Adesign = np.array([
    # LA  LB   M   delta
    [0, 1, 0, 1],    # c1  L_B open
    [1, 0, 0, 1],    # c3  L_A open
    [1, 1, 2, 1],    # c5  aiding
    [1, 1, 2, 1],    # c6  aiding
    [1, 1, -2, 1],   # c7  opposing
    [1, 1, -2, 1],   # c8  opposing
], dtype=float)
DELTA = {}
for f in (1000, 10000):
    tag = "1k" if f == 1000 else "10k"
    r = R[f]
    y = np.array([r[1], r[3], r[5], r[6], r[7], r[8]])
    sol, res, rank, sv = np.linalg.lstsq(Adesign, y, rcond=None)
    LAc, LBc, Mc, dc = sol
    resid = y - Adesign @ sol
    kc = Mc / math.sqrt(LAc * LBc)
    DELTA[f] = dict(LA=LAc, LB=LBc, M=Mc, d=dc, k=kc)
    put(f"meas.delta_fit_{tag}_uH", dc, "[MODEL] common additive series inductance in every 2-terminal reading; LSQ on c1,c3,c5-c8", 4)
    put(f"meas.delta_fit_LA_{tag}_uH", LAc, "[MODEL] delta-corrected L_A", 4)
    put(f"meas.delta_fit_LB_{tag}_uH", LBc, "[MODEL] delta-corrected L_B", 4)
    put(f"meas.delta_fit_M_{tag}_uH", Mc, "[MODEL] delta-fit M (identical to pair-mean M by construction)", 4)
    put(f"meas.k_delta_corrected_{tag}", kc, "[MODEL] M/sqrt(LA LB) with delta removed", 4)
    put(f"meas.delta_fit_rank_{tag}", int(rank), "[CALCULATED] design-matrix rank (4 = identifiable)")
    put(f"meas.delta_fit_dof_{tag}", int(len(y) - rank), "[CALCULATED] residual degrees of freedom")
    put(f"meas.delta_fit_resid_rms_{tag}_uH", float(np.sqrt(np.mean(resid**2))), "[CALCULATED] RMS residual; residuals equal +/- half the within-pair spreads", 4)
    put(f"meas.delta_fit_resid_{tag}_uH", [round(float(x), 4) for x in resid], "[CALCULATED] residuals order c1,c3,c5,c6,c7,c8")
    addk(f, "delta (lead-inductance) LSQ fit", "LA,LB,M,delta on c1,c3,c5-c8", LAc, LBc, Mc, kc, f"delta={dc:.4f} uH")
    if f == 10000:
        # cross-check: shorted-secondary k with delta-corrected L (delta cancels in dL)
        for nm, (op, sh, L1, L2) in {"B": (1, 2, LBc, LAc), "A": (3, 4, LAc, LBc)}.items():
            dL = r[op] - r[sh]
            put(f"meas.k_from_short_deltacorr_{nm}_10k", k_from_dL(f, dL, L1, L2, R_DC),
                "[MODEL] shorted-test k using delta-corrected L1,L2 (>1 is unphysical)", 4)
            put(f"meas.k_from_short_deltacorr_{nm}_10k_dLm002", k_from_dL(f, dL - 0.02, L1, L2, R_DC),
                "[MODEL] same with dL-0.02 uH", 4)
        # largest delta compatible with k_short <= 1 (nominal dL and dL-0.02)
        for dd, lab in [(0.0, "nom"), (-0.02, "dLm002")]:
            lo, hi = 0.0, 3.0
            for _ in range(60):
                mid = (lo + hi) / 2
                kk = k_from_dL(f, r[1] - r[2] + dd, r[1] - mid, r[3] - mid, R_DC)
                lo, hi = (mid, hi) if kk <= 1 else (lo, mid)
            put(f"meas.delta_max_for_kshort_le1_{lab}_10k_uH", lo, f"[MODEL] largest common delta keeping shorted-test k<=1 (winding B, dL{dd:+.2f})", 3)

# Joint nonlinear fits at 10 kHz on all 8 readings (sensitivity, not headline)
f = 10000
w = 2 * math.pi * f
r = R[f]
yobs = np.array([r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]])


def model(p, kind):
    LA, LB, M = p[0], p[1], p[2]
    d = p[3] if kind in ("delta", "both") else 0.0
    b = p[-1] if kind in ("beta", "both") else 0.0
    dA = w * w * (M * 1e-6) ** 2 * (LA * 1e-6) / (R_DC**2 + (w * LA * 1e-6) ** 2) * 1e6   # drop of B when A shorted
    dB = w * w * (M * 1e-6) ** 2 * (LB * 1e-6) / (R_DC**2 + (w * LB * 1e-6) ** 2) * 1e6
    return np.array([LB + d, LB - dA + d, LA + d, LA - dB + d,
                     LA + LB + 2 * M + d, LA + LB + 2 * M + d,
                     LA + LB - 2 * M + d + b, LA + LB - 2 * M + d + b])


for kind, p0 in [("none", [18, 18, 16]), ("delta", [18, 18, 16, 0.3]), ("beta", [18, 18, 16, -0.5])]:
    sol = least_squares(lambda p: model(p, kind) - yobs, p0)
    p = sol.x
    rr = model(p, kind) - yobs
    kk = p[2] / math.sqrt(p[0] * p[1])
    extra = "" if kind == "none" else f"; {'delta' if kind=='delta' else 'beta(opp bias)'}={p[3]:.4f} uH"
    addk(f, "joint unweighted LSQ, 8 readings incl. shorted (reflected-impedance model)", kind,
         p[0], p[1], p[2], kk, f"RMS resid={np.sqrt(np.mean(rr**2)):.4f} uH{extra}")
    put(f"meas.joint_{kind}_k_10k", kk, f"[MODEL] joint fit '{kind}' on c1-c8 at 10 kHz, R2={R_DC}", 4)
    put(f"meas.joint_{kind}_rms_10k_uH", float(np.sqrt(np.mean(rr**2))), "[MODEL] RMS residual of that fit", 4)
    if kind != "none":
        put(f"meas.joint_{kind}_param_10k_uH", float(p[3]), f"[MODEL] fitted {'delta' if kind=='delta' else 'opposing-only bias beta'}", 4)

with open(f"{OUT}/data/k_table.csv", "w", newline="", encoding="utf-8") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(ktab[0].keys()))
    wr.writeheader()
    wr.writerows(ktab)

# --------------------------------------------------------------------------------------
# 7. Uncertainty budget
# --------------------------------------------------------------------------------------
# Per-reading best-case accuracy: 0.1 % (LCR999A basic accuracy) applied to |Z|, i.e. the
# |Z|*sqrt(1+D^2)-type convention used by other handheld meters (DE-5000 manual); LCR999A
# detailed accuracy table = GAP, so this is a LOWER bound on the real accuracy limit.
budget = {}
for f in (1000, 10000):
    tag = "1k" if f == 1000 else "10k"
    w = 2 * math.pi * f
    for cls, Lval, Rser in [("single", (R[f][1] + R[f][3]) / 2, R_DC),
                            ("aid", (R[f][5] + R[f][6]) / 2, 2 * R_DC),
                            ("opp", (R[f][7] + R[f][8]) / 2, 2 * R_DC)]:
        X = w * Lval * 1e-6
        Z = math.hypot(Rser, X)
        uacc = BASIC_ACC * Z / w * 1e6
        budget[(tag, cls)] = uacc
        put(f"meas.u_acc_{cls}_{tag}_uH", uacc, f"[CALCULATED] 0.1 % of |Z|/w (|Z|={Z:.3f} ohm): best-case basic-accuracy bound; LCR999A table GAP", 4)
        put(f"meas.u_acc_{cls}_{tag}_pct", 100 * uacc / Lval, "[CALCULATED] same as % of reading", 2)
        put(f"meas.Z_{cls}_{tag}_ohm", Z, "[CALCULATED] |Z| = sqrt(Rs^2+X^2) with Rs from 7.31 ohm/winding", 3)
    put(f"meas.u_res_{tag}_uH", LSD / 2 / math.sqrt(3), "[CALCULATED] 0.01 uH LSD, rectangular half-width 0.005 uH", 4)
    # phase-offset sensitivity: dL per 0.1 deg
    for cls, Rser in [("single", R_DC), ("series", 2 * R_DC)]:
        put(f"meas.dL_per_0p1deg_{cls}_{tag}_uH", Rser * math.radians(0.1) / w * 1e6,
            "[CALCULATED] Rs*dtheta/w for a 0.1 deg phase error (sensitivity only; LCR999A phase accuracy GAP)", 3)
    # propagate accuracy to M and k (pair-mean), independent errors
    us, ua, uo = budget[(tag, "single")], budget[(tag, "aid")], budget[(tag, "opp")]
    ures = LSD / 2 / math.sqrt(3)
    E = EXT[1000 if tag == "1k" else 10000]
    # two pairs averaged: aid mean and opp mean have u/sqrt(2) if independent; keep conservative (no averaging gain)
    uM = math.sqrt(ua**2 + uo**2 + 2 * ures**2) / 4
    uLrel = math.sqrt(us**2 + ures**2) / ((E["LA"] + E["LB"]) / 2)
    uk = E["km"] * math.sqrt((uM / E["Mm"]) ** 2 + 0.5 * uLrel**2)
    put(f"meas.u_M_acc_{tag}_uH", uM, "[CALCULATED] sqrt(u_aid^2+u_opp^2+2u_res^2)/4, best-case accuracy only", 3)
    put(f"meas.u_k_acc_{tag}", uk, "[CALCULATED] pair-mean k standard-uncertainty-like bound from best-case accuracy + resolution only (excludes systematic method spread)", 4)
    put(f"meas.u_M_pairspread_{tag}_uH", abs(EXT[1000 if tag=='1k' else 10000]['k1'] - EXT[1000 if tag=='1k' else 10000]['k2']) / 2 * E["G"], "[CALCULATED] half the pair1-pair2 M difference", 3)
    put(f"meas.opp_below_Lmin_{tag}", bool(max(R[1000 if tag=='1k' else 10000][o] for o in (7, 8)) < L_RANGE_MIN),
        f"[CALCULATED] opposing readings < {L_RANGE_MIN} uH catalogue L lower limit (mu glyph inferred)")
    put(f"meas.X_opp_below_Zmin_{tag}", bool(w * max(R[1000 if tag=='1k' else 10000][o] for o in (7, 8)) * 1e-6 < Z_RANGE_MIN),
        f"[CALCULATED] opposing reactance below the {Z_RANGE_MIN} ohm catalogue R/Z lower limit (interpretation: limit applied to reactive part)")
    put(f"meas.X_single_below_Zmin_{tag}", bool(w * (E['LA'] + E['LB']) / 2 * 1e-6 < Z_RANGE_MIN),
        "[CALCULATED] single-winding reactance below 0.4 ohm (same interpretation)")

# headline k and interval: method envelope over both frequencies
cands = {
    "aid_only_min": min(EXT[1000]["ka"][5], EXT[1000]["ka"][6], EXT[10000]["ka"][5], EXT[10000]["ka"][6]),
    "pair_mean_10k": EXT[10000]["km"],
    "opp_only_max": max(EXT[1000]["ko"][7], EXT[1000]["ko"][8], EXT[10000]["ko"][7], EXT[10000]["ko"][8]),
    "delta_max": max(DELTA[1000]["k"], DELTA[10000]["k"]),
    "short_B": V["meas.k_from_short_B_10k"], "short_A": V["meas.k_from_short_A_10k"],
}
k_lo = cands["aid_only_min"]
k_hi = max(cands["delta_max"], cands["short_B"], cands["short_A"])
put("meas.k_headline", round(EXT[10000]["km"], 2), "[CALCULATED] pair (series aiding/opposing) method, 10 kHz, rounded to 2 dp; 1 kHz agrees within 0.002", 2)
put("meas.k_interval_lo", round(k_lo, 2), "[CALCULATED] lower edge of method envelope (aiding-only, lowest of 1/10 kHz), rounded to 2 dp; NOT a confidence interval", 2)
put("meas.k_interval_hi", round(k_hi, 2), "[CALCULATED] upper edge of method envelope (max of delta-corrected and shorted-secondary estimates), rounded to 2 dp; NOT a confidence interval", 2)
put("meas.k_interval_lo_raw", k_lo, "[CALCULATED] unrounded lower envelope", 4)
put("meas.k_interval_hi_raw", k_hi, "[CALCULATED] unrounded upper envelope", 4)
put("meas.M_headline_uH", round(EXT[10000]["Mm"], 1), "[CALCULATED] pair-mean M at 10 kHz, 1 dp", 1)
Mlo = min(min(EXT[f]['ka'][a] * EXT[f]['G'] for a in (5, 6)) for f in (1000, 10000))
Mhi = max(V["meas.M_pair1_10k"], V["meas.M_pair1_1k"])
put("meas.M_interval_lo_uH", Mlo, "[CALCULATED] lowest M in envelope (aiding-only); shorted-test M excluded (see report)", 3)
put("meas.M_interval_hi_uH", Mhi, "[CALCULATED] highest pair-method M", 3)
allL = [R[f][n] for f in (1000, 10000) for n in (1, 3)] + [17.9, 18.0, 18.4, 18.6, 18.3]
put("meas.L_records_min_uH", min(allL), "[CALCULATED] min over all open-other / unspecified single-winding L records (incl. 17.9/18.0 earlier, 18.3-18.6, 18.4)", 2)
put("meas.L_records_max_uH", max(allL), "[CALCULATED] max over the same records", 2)
put("meas.L_headline_uH", round((EXT[10000]["LA"] + EXT[10000]["LB"]) / 2, 2), "[CALCULATED] mean open-other reading at 10 kHz", 2)
put("meas.L_records_spread_pct", 100 * (max(allL) - min(allL)) / min(allL), "[CALCULATED] (max-min)/min over L records", 1)

# inter-frequency consistency: a frequency-independent phase offset dtheta adds R*dtheta/w to Ls,
# i.e. 10x more at 1 kHz than at 10 kHz; solve for the dtheta implied by mean open readings
w1, w10 = 2 * math.pi * 1e3, 2 * math.pi * 1e4
L1m = (R[1000][1] + R[1000][3]) / 2 * 1e-6
L10m = (R[10000][1] + R[10000][3]) / 2 * 1e-6
dth = (L1m - L10m) / (R_DC * (1 / w1 - 1 / w10))
put("meas.phase_offset_implied_deg", math.degrees(dth), "[MODEL] constant phase offset that would explain L_open(1k)-L_open(10k) (if that were the only error)", 4)
put("meas.L_open_1k_minus_10k_uH", (L1m - L10m) * 1e6, "[CALCULATED] mean open reading 1 kHz minus 10 kHz", 3)
put("meas.L_oct_minus_nov_min_uH", 18.3 - (R[10000][1] + R[10000][3]) / 2, "[CALCULATED] 18.3 (lowest Oct record) minus 10 kHz Nov open mean", 3)
put("meas.L_oct_minus_nov_max_uH", 18.6 - (R[1000][1] + R[1000][3]) / 2, "[CALCULATED] 18.6 minus 1 kHz Nov open mean", 3)

# capacitance effect check
for f in (1000, 10000):
    tag = "1k" if f == 1000 else "10k"
    w = 2 * math.pi * f
    put(f"meas.w2LC_aid_{tag}", w * w * (R[f][5] * 1e-6) * C_AB_PF * 1e-12, "[CALCULATED] w^2 L_aid C_AB: scale of capacitive correction (negligible if <<1e-3)", 8)
put("meas.XC_AB_10k_kohm", 1 / (2 * math.pi * 1e4 * C_AB_PF * 1e-12) / 1e3, "[CALCULATED] 1/(w C_AB) at 10 kHz", 1)

# --------------------------------------------------------------------------------------
# 8. DC-drop numbers (CALCULATED, not ratings)
# --------------------------------------------------------------------------------------
I42 = V_RAIL / (R_LOAD + R_DC)
put("meas.dc_I_calc_mA", 1e3 * I42, f"[CALCULATED] 5 V/({R_LOAD}+{R_DC}) ohm; not a rating", 2)
put("meas.dc_I_calc_alt_mA", 1e3 * V_RAIL / (R_LOAD_ALT + R_DC), f"[CALCULATED] 5 V/({R_LOAD_ALT}+{R_DC}) ohm", 2)
for I, lab in [(0.030, "30mA"), (I42, "Icalc")]:
    put(f"meas.dc_drop_{lab}_V", I * R_DC, "[CALCULATED] I*7.31 ohm copper drop; not a rating", 4)
    put(f"meas.dc_drop_{lab}_pct", 100 * I * R_DC / V_RAIL, "[CALCULATED] drop as % of 5 V", 2)
    put(f"meas.dc_Pcu_{lab}_mW", 1e3 * I * I * R_DC, "[CALCULATED] I^2*7.31 ohm copper dissipation; not a thermal rating", 2)
I10 = 0.1 * V_RAIL / R_DC
put("meas.dc_I_10pct_drop_mA", 1e3 * I10, "[CALCULATED] 0.5 V/7.31 ohm: 10 % drop criterion at 5 V; NOT a current rating", 1)
put("meas.dc_Pcu_10pct_mW", 1e3 * I10 * I10 * R_DC, "[CALCULATED] copper dissipation at that current; not thermally verified", 1)
put("meas.dc_printed500mA_drop_V", 0.5 * R_DC, "[CALCULATED] printed 500 mA claim x 7.31 ohm (shows the claim is not a usable rating)", 3)
put("meas.dc_printed500mA_P_W", 0.25 * R_DC, "[CALCULATED] (0.5 A)^2 x 7.31 ohm", 3)
put("meas.dc_display_note", "supply display 0.03 A vs calculated ~42 mA: unreconciled (display resolution/accuracy, actual rail voltage, lead resistance not recorded)",
    "[OPEN] bench spreadsheet: 5 V rail, load 111.6 ohm, supply display 0.03 A")

with open(f"{OUT}/data/values_meas.json", "w", encoding="utf-8") as fh:
    json.dump(V, fh, indent=1, ensure_ascii=False)

# console summary
for k in ["meas.k_pair_mean_1k", "meas.k_pair_mean_10k", "meas.k_aid_only_mean_10k", "meas.k_opp_only_mean_10k",
          "meas.delta_fit_10k_uH", "meas.k_delta_corrected_10k", "meas.delta_fit_1k_uH", "meas.k_delta_corrected_1k",
          "meas.k_from_short_B_10k", "meas.k_from_short_deltacorr_B_10k", "meas.delta_max_for_kshort_le1_nom_10k_uH",
          "meas.joint_none_k_10k", "meas.joint_delta_k_10k", "meas.joint_beta_k_10k",
          "meas.k_headline", "meas.k_interval_lo", "meas.k_interval_hi", "meas.u_k_acc_10k"]:
    print(f"{k:50s} {V[k]}")
print(f"wrote {len(raw_rows)} raw rows, {len(ktab)} k rows, {sum(1 for k in V if not k.endswith('_meta'))} values")
