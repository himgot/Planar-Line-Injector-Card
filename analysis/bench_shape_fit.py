#!/usr/bin/env python3
"""M2 with R_g FIXED at 50 ohm and one fitted gain (least squares in dB)
over 0.1-1 MHz; compare with the free-R_g fit; capacitive feed-through slope test;
display-scale check against the load DC voltage; quantisation estimate."""
import csv, json, math, os
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M = json.load(open(os.path.join(ROOT, "data", "values_meas.json")))
g = lambda k: M[k]["value"] if isinstance(M[k], dict) else M[k]
LA, LB, Mm, R = g("meas.L_A_open_10k")*1e-6, g("meas.L_B_open_10k")*1e-6, g("meas.M_pair_mean_10k")*1e-6, g("meas.R_dc_ohm")
RL = 111.6
def H(f, Rg=50.0, RLx=RL):
    s = 2j*math.pi*f; ZA = Rg+R+s*LA; ZB = R+RLx+s*LB
    return abs(s*Mm*RLx/(ZA*ZB - (s*Mm)**2))
rows = list(csv.DictReader(open(os.path.join(ROOT, "data", "bench_reconciled.csv"))))
out = {}
def fit(run, fmin, fmax, Rg=50.0, include45=True):
    pts = [(float(r["f_kHz"])*1e3, float(r["Vpp_displayed"])) for r in rows
           if r["run"] == run and fmin <= float(r["f_kHz"])*1e3 <= fmax*1.04 and (include45 or float(r["f_kHz"]) > 50)]
    f = np.array([p[0] for p in pts]); y = 20*np.log10([p[1] for p in pts])
    m = 20*np.log10([H(x, Rg) for x in f]); gain = np.mean(y-m); res = y-m-gain
    return dict(n=len(f), gain_dB=gain, rms=float(np.sqrt(np.mean(res**2))), maxabs=float(np.max(np.abs(res))), f=f, res=res)
for run, fmin in [("R", 0.09e6), ("X", 0.09e6)]:
    r = fit(run, fmin, 1.0e6)
    out[f"fit.{run}_fixed50_n"] = r["n"]; out[f"fit.{run}_fixed50_rms_dB"] = round(r["rms"], 2); out[f"fit.{run}_fixed50_max_dB"] = round(r["maxabs"], 2)
    # implied generator EMF (Vpp) with display /10: Vdisp/10 = |H| * Vg * 10^(0)  -> Vg = 10^(gain/20)/10
    out[f"fit.{run}_EMF_div10_Vpp"] = round(10**(r["gain_dB"]/20)/10, 2)
    # free Rg (2 parameters)
    best = min(((fit(run, fmin, 1.0e6, Rg)["rms"], Rg) for Rg in np.arange(20, 400, 1.0)))
    out[f"fit.{run}_freeRg_ohm"] = round(best[1]); out[f"fit.{run}_freeRg_rms_dB"] = round(best[0], 2)
    print(run, "fixed 50:", r["n"], "pts rms %.2f max %.2f dB; residuals" % (r["rms"], r["maxabs"]), [f"{a/1e3:.0f}k:{b:+.2f}" for a, b in zip(r["f"], r["res"])], "| free Rg", best)
# Run R if R_L were ~160 ohm total-consistent (0.03 A): R_L_eff = 5/0.03 - 7.31
RLalt = 5/0.03 - R
f_R = [float(r["f_kHz"])*1e3 for r in rows if r["run"] == "R" and float(r["f_kHz"]) <= 1050]
y_R = 20*np.log10([float(r["Vpp_displayed"]) for r in rows if r["run"] == "R" and float(r["f_kHz"]) <= 1050])
m_alt = 20*np.log10([H(x, 50.0, RLalt) for x in f_R]); ga = np.mean(y_R-m_alt)
out["fit.R_fixed50_RLalt_ohm"] = round(RLalt); out["fit.R_fixed50_RLalt_rms_dB"] = round(float(np.sqrt(np.mean((y_R-m_alt-ga)**2))), 2)
# slope test 0.5 -> 1 MHz
def cap_only(f, C=274e-12, Rs=50.0+RL):
    w = 2*math.pi*f; x = w*C*Rs; return x/math.sqrt(1+x*x)
out["fit.slope_cap_0p5to1_dB"] = round(20*math.log10(cap_only(1e6)/cap_only(0.5e6)), 1)
out["fit.slope_cap_corner_MHz"] = round(1/(2*math.pi*274e-12*(50+RL))/1e6, 1)
out["fit.slope_M2_0p5to1_dB"] = round(20*math.log10(H(1.0032e6)/H(0.5008e6)), 1)
v = {r["f_kHz"]: float(r["Vpp_displayed"]) for r in rows if r["run"] == "R"}
out["fit.slope_runR_0p5to1_dB"] = round(20*math.log10(v["1003.2"]/v["500.8"]), 1)
# display scale
Vload = 5.0*RL/(RL+R); out["fit.Vload_dc_V"] = round(Vload, 2); out["fit.display_ratio"] = round(48.3/Vload, 1)
out["fit.Ialt_total_ohm"] = round(5/0.03)
# quantisation: 8-bit over 8 div -> 25 levels/div; 100 kHz trace 0.317 div pk-pk
lv = 256/10.24  # DSO1000B: ~25 levels/div (8 bits over ~10 div incl. off-screen); conservative 25
out["fit.quant_levels_100k"] = round(0.317*25, 1)
out["fit.quant_err_100k_dB"] = round(20*math.log10((0.317*25+1)/(0.317*25)), 1)
# M2 -3 dB band (relative to its mid-band peak), Rg = 50 ohm, R_L = 111.6 ohm
ff = np.logspace(4, 8, 20000); hh = np.array([H(x) for x in ff]); pk = hh.max(); inb = ff[hh >= pk/math.sqrt(2)]
out["fit.M2_peak_dB"] = round(20*math.log10(pk), 1); out["fit.M2_band_lo_kHz"] = round(inb[0]/1e3); out["fit.M2_band_hi_MHz"] = round(inb[-1]/1e6, 1)
# worked example: DUT input = R_L in parallel with an input capacitor C_in
def Hz(f, Cin, Rg=50.0, RLx=RL):
    s_ = 2j*math.pi*f; ZL = 1/(1/RLx + s_*Cin) if Cin else RLx
    ZA = Rg+R+s_*LA; ZB = R+ZL+s_*LB
    return abs(s_*Mm*ZL/(ZA*ZB - (s_*Mm)**2))
for Cin, tag in [(1e-6, "1u"), (10e-6, "10u")]:
    for f0, ftag in [(100e3, "100k"), (1e6, "1M")]:
        out[f"fit.Cin{tag}_{ftag}_dB"] = round(20*math.log10(Hz(f0, Cin)), 0)
# full-decade slope test 0.1 -> 1 MHz (Run R)
out["fit.slope_runR_decade_dB"] = round(20*math.log10(v["1003.2"]/v["99.8"]), 1)
out["fit.slope_M2_decade_dB"] = round(20*math.log10(H(1.0032e6)/H(99.8e3)), 1)
out["fit.slope_cap_decade_dB"] = round(20*math.log10(cap_only(1.0032e6)/cap_only(99.8e3)), 1)
# M1 (geometry-only: theory L, M, R) vs Run R, R_g fixed 50 ohm, one gain
T = json.load(open(os.path.join(ROOT, "data", "values_theory.json")))
tg = lambda k: T[k]["value"] if isinstance(T[k], dict) else T[k]
LA1, LB1, M1, R1 = tg("theory.L_A_uH")*1e-6, tg("theory.L_D_uH")*1e-6, tg("theory.M_uH")*1e-6, tg("theory.R_A_w0p254_t35_ohm")
def H1(f, Rg=50.0):
    s_ = 2j*math.pi*f; ZA = Rg+R1+s_*LA1; ZB = R1+RL+s_*LB1
    return abs(s_*M1*RL/(ZA*ZB - (s_*M1)**2))
m1 = 20*np.log10([H1(x) for x in f_R]); g1 = np.mean(y_R-m1)
out["fit.R_M1_fixed50_rms_dB"] = round(float(np.sqrt(np.mean((y_R-m1-g1)**2))), 2)
json.dump(out, open(os.path.join(ROOT, "data", "values_fit.json"), "w"), indent=1)
for k, v in out.items(): print(k, v)
