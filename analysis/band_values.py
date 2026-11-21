#!/usr/bin/env python3
"""Band, output-impedance and corner computations.
All inputs from data/*.json and data/bench_reconciled.csv; outputs data/values_band.json."""
import csv, json, math, os
import numpy as np
from scipy.optimize import minimize_scalar
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = lambda n: json.load(open(os.path.join(ROOT, "data", n)))
M, T = J("values_meas.json"), J("values_theory.json")
g = lambda d, k: d[k]["value"] if isinstance(d[k], dict) else d[k]
LA, LB, Mm, R = g(M,"meas.L_A_open_10k")*1e-6, g(M,"meas.L_B_open_10k")*1e-6, g(M,"meas.M_pair_mean_10k")*1e-6, g(M,"meas.R_dc_ohm")
RL, Rg0 = 111.6, 50.0
out = {}
def H(f, Rg=Rg0, Cab=0.0, k=None, scaleLM=1.0, Rw=R, ZLx=None):
    """Exact two-loop response with optional lumped C_AB/3 between AC_IN and OUT (like-polarity inner ends)."""
    w = 2*math.pi*f; s = 1j*w
    la, lb = LA*scaleLM, LB*scaleLM
    mm = (k*math.sqrt(la*lb)) if k is not None else Mm*scaleLM
    ZL = RL if ZLx is None else ZLx(s)
    Y = s*Cab/3.0
    # unknowns Va (AC_IN), Vo (OUT), I1 (into AC_IN), I2 (into OUT terminal of winding B)
    A = np.array([[1/Rg+Y, -Y, 1, 0], [-Y, Y+1/ZL, 0, 1], [1, 0, -(Rw+s*la), -s*mm], [0, 1, -s*mm, -(Rw+s*lb)]], dtype=complex)
    b = np.array([1/Rg, 0, 0, 0], dtype=complex)
    return abs(np.linalg.solve(A, b)[1])
ff = np.logspace(4, 8, 6000)
def band(**kw):
    """Lower/upper -3 dB points relative to the mid-band peak (searched 0.1-5 MHz);
    upper edge = FIRST crossing below peak/sqrt2 above the peak (the capacitive path
    makes the response rise again far above the band, which is not part of the band)."""
    h = np.array([H(x, **kw) for x in ff]); sel = (ff >= 1e5) & (ff <= 5e6)
    ipk = np.where(sel)[0][np.argmax(h[sel])]; pk = h[ipk]; thr = pk/math.sqrt(2)
    lo = ff[np.where(h[:ipk] < thr)[0][-1] + 1]
    above = np.where(h[ipk:] < thr)[0]; hi = ff[ipk + above[0]] if len(above) else float("nan")
    return lo, hi
# predicted C_AB taken from capacitance.py (0.8 mm FR-4, JLCPCB eps 4.5, mask 10 um eps 3.8)
# board is mask-covered: baseline range = conformal-mask minimum .. gap-filled bound 
_Cset = [("C0", 0.0), ("Clo", g(T, "theory.C_AB_masked_pF_lo")*1e-12), ("Cnom", g(T, "theory.C_AB_pF_nominal")*1e-12),
         ("Chi", g(T, "theory.C_AB_pF_filled_upper")*1e-12), ("C274", 274e-12)]
for tag, C in _Cset:
    lo, hi = band(Cab=C)
    out[f"band.band_hi_{tag}_MHz"] = round(hi/1e6, 1); out[f"band.band_lo_{tag}_kHz"] = round(lo/1e3)
# upper edge depends mainly on k -- delta-corrected k 0.93 / 0.95 with nominal C_AB
for _k, _t in ((0.93, "k093"), (0.95, "k095")):
    _lo, _hi = band(Cab=g(T, "theory.C_AB_pF_nominal")*1e-12, k=_k)
    out[f"band.band_hi_{_t}_Cnom_MHz"] = round(_hi/1e6, 1)
# output impedance seen from OUT (DC winding) with generator connected
def Zout(f, Rg=Rg0):
    s = 2j*math.pi*f; return (R + s*LB) + (s*Mm)**2*(-1)/(Rg + R + s*LA)
for f, t in [(1e5,"100k"), (1e6,"1M"), (5e6,"5M")]:
    z = Zout(f); out[f"band.Zout_{t}_mag"] = round(abs(z)); out[f"band.Zout_{t}_re"] = round(z.real)
# max DUT input capacitance (parallel with RL) for a given level at 1 MHz
for lvl in (-20, -30, -40):
    res = minimize_scalar(lambda lc: (20*math.log10(H(1e6, ZLx=lambda s, C=10**lc: 1/(1/RL + s*C))) - lvl)**2, bounds=(-11, -5), method="bounded")
    out[f"band.Cmax_{abs(lvl)}dB_1MHz_nF"] = round(10**res.x*1e9)
# low-impedance driver (SG-01): corner and Zout with Rg = 5 ohm
lo5, hi5 = band(Rg=5.0); out["band.lowZ_Rg5_lo_kHz"] = round(lo5/1e3); out["band.lowZ_Rg5_Zout_1M"] = round(abs(Zout(1e6, Rg=5.0)))
# k uncertainty: phase-error model dL = |Z| dtheta / w at 10 kHz (assumed 1 and 5 mrad)
w = 2*math.pi*1e4
LAm, LBm = g(M,"meas.L_A_open_10k")*1e-6, g(M,"meas.L_B_open_10k")*1e-6
La1, Lo1 = 68.19e-6, 2.79e-6
Zs, Za, Zo = abs(R+1j*w*LAm), abs(2*R+1j*w*La1), abs(2*R+1j*w*Lo1)
for dth, t in [(1e-3, "1m"), (5e-3, "5m")]:
    dLs, dLa, dLo = Zs*dth/w, Za*dth/w, Zo*dth/w
    k0 = (La1-Lo1)/4/math.sqrt(LAm*LBm)
    # uncorrelated propagation
    dk_dLa = 1/(4*math.sqrt(LAm*LBm)); dk_dLo = -dk_dLa; dk_dLs = -k0/(2*LAm)
    sk = math.sqrt((dk_dLa*dLa)**2 + (dk_dLo*dLo)**2 + 2*(dk_dLs*dLs)**2)
    out[f"band.dL_single_{t}_uH"] = round(dLs*1e6, 2); out[f"band.dL_opp_{t}_uH"] = round(dLo*1e6, 2)
    out[f"band.sigma_k_{t}"] = round(sk, 3)
# Run R shape-test variants (R_g fixed 50 ohm, one fitted gain)
rows = list(csv.DictReader(open(os.path.join(ROOT, "data", "bench_reconciled.csv"))))
pts = [(float(r["f_kHz"])*1e3, float(r["Vpp_displayed"]), r["photos"].split()[0]) for r in rows if r["run"] == "R" and float(r["f_kHz"]) <= 1050]
f_R = np.array([p[0] for p in pts]); y_R = 20*np.log10([p[1] for p in pts]); ids = [p[2] for p in pts]
def rms(model, mask=None):
    m = np.array([20*math.log10(model(x)) for x in f_R]); y = y_R
    if mask is not None: m, y = m[mask], y[mask]
    d = y - m; d = d - d.mean(); return float(np.sqrt(np.mean(d**2)))
out["band.fit_M2"] = round(rms(lambda x: H(x)), 2)
out["band.fit_k030"] = round(rms(lambda x: H(x, k=0.30)), 2)
out["band.fit_k099"] = round(rms(lambda x: H(x, k=0.99)), 2)
out["band.fit_LM08"] = round(rms(lambda x: H(x, scaleLM=0.8)), 2)
out["band.fit_cap_only"] = round(rms(lambda x: (lambda wc: wc/math.sqrt(1+wc*wc))(2*math.pi*x*274e-12*(50+RL))), 2)
best = min(((rms(lambda x, fc=fc: (x/fc)/math.sqrt(1+(x/fc)**2)), fc) for fc in np.arange(200e3, 600e3, 1e3)))
out["band.fit_generic_hp"] = round(best[0], 2); out["band.fit_generic_hp_corner_kHz"] = round(best[1]/1e3)
mask_noT = np.array([i not in ("R09", "R10") for i in ids]); mask_no10 = np.array([i != "R10" for i in ids])
out["band.fit_M2_noTiebreak"] = round(rms(lambda x: H(x), mask_noT), 2); out["band.fit_M2_noR10"] = round(rms(lambda x: H(x), mask_no10), 2)
# winding resistance doubled (proxy for R_ac at MHz)
out["band.R2x_dB_1MHz"] = round(20*math.log10(H(1e6, Rw=2*R)/H(1e6)), 1)
# HF-path bound at 1 MHz: second-order skirt giving +3 dB at 3 MHz (in phase)
a3 = 10**(3/20) - 1; a1 = a3*(1/3)**2; out["band.hf_bound_1MHz_dB"] = round(20*math.log10(1+a1), 1)
# effective permittivity check of the measured C_AB
Lgap = 5.633; e0 = 8.854e-12
Cair0 = e0*Lgap*1.0   # zero-thickness conformal map (air)
_cr = json.load(open(os.path.join(ROOT, "work", "theory", "capacitance_results.json")))
_air = [c for c in _cr["fd_cases"] if abs(c["w_mm"] - 0.254) < 1e-6 and abs(c["t_um"] - 35) < 1e-6][0]
Cair_fv = _air["Cgap_air_pF_per_m"]*1e-12*_cr["L_adj_m"]   # own finite-volume solver, 35 um copper, air
out["band.Cair_pF_zero_t"] = round(e0*_cr["L_adj_m"]*1e12); out["band.Cair_pF_t35"] = round(Cair_fv*1e12)
out["band.eps_eff_meas"] = round(274e-12/Cair_fv, 1)
# Keysight U1282A data sheet (5992-0847EN), specified after Math Null: C 10 nF range 1 % + 5 counts (1 pF resolution);
# R 60 ohm range 0.15 % + 20 counts (1 mOhm resolution); R 600 ohm range 0.05 % + 10 counts (10 mOhm resolution)
out["band.U1282A_dC_pF"] = round(0.01*274 + 5*1.0, 1); out["band.U1282A_dR_mohm"] = round((0.0015*7.31 + 20*0.001)*1e3)
out["band.U1282A_dRL_ohm"] = round(0.0005*111.6 + 10*0.01, 2)
# corner range allowed by Run R (free first-order corner, free gain) and corner vs generator source resistance
import csv as _csv
_rows = [r for r in _csv.DictReader(open(os.path.join(ROOT, "data", "bench_reconciled.csv"))) if r["run"] == "R" and float(r["f_kHz"]) <= 1100]
_f = np.array([float(r["f_kHz"]) for r in _rows])*1e3; _y = np.array([float(r["rel_dB_re_1MHz"]) for r in _rows])
def _rms(fc):
    m = 20*np.log10((_f/fc)/np.sqrt(1 + (_f/fc)**2)); g = np.mean(_y - m); return math.sqrt(np.mean((_y - m - g)**2))
_fcs = np.arange(150e3, 1000e3, 1e3); _r = np.array([_rms(x) for x in _fcs]); _ok = _fcs[_r <= 0.5]
out["band.corner_fit05_lo_kHz"] = round(_ok.min()/1e3); out["band.corner_fit05_hi_kHz"] = round(_ok.max()/1e3)
_L = g(M, "meas.L_A_open_10k")*1e-6
for _Rg in (50, 200, 600):
    _a, _b = _Rg + R, RL + R; _fc = (_a*_b/(_a + _b))/(2*math.pi*_L)
    out[f"band.corner_Rg{_Rg}_kHz"] = round(_fc/1e3); out[f"band.corner_Rg{_Rg}_rms_dB"] = round(_rms(_fc), 2)
json.dump(out, open(os.path.join(ROOT, "data", "values_band.json"), "w"), indent=1)
for k, v in out.items(): print(k, v)
