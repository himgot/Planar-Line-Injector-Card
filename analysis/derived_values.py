#!/usr/bin/env python3
"""Derived display values from data/values_sim.json (unit conversions, differences)."""
import json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S = json.load(open(os.path.join(ROOT, "data", "values_sim.json")))
g = lambda k: S[k]["value"] if isinstance(S[k], dict) else S[k]
out = {
 "sim.M0_peak_f_MHz": round(g("sim.M0_peak_f_Hz") / 1e6, 2),
 "sim.M2_f3dB_kHz": round(g("sim.M2_f_minus3dB_low_Hz") / 1e3),
 "sim.M3_vs_M2_absdiff_1MHz_dB": round(abs(g("sim.M3_dB_1MHz") - g("sim.M2_dB_1MHz")), 2),
 "sim.M3_first_pole_MHz": round(g("sim.M3_first_pole_f_Hz") / 1e6, 2),
}
json.dump(out, open(os.path.join(ROOT, "data", "values_simderived.json"), "w"), indent=1)
print(out)
