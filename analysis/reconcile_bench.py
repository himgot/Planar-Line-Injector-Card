#!/usr/bin/env python3
"""Reconcile independent bench-photo readings into one table (G5).

Readers: A  = reading A (22 Oct raw R* and workbook X* images; data/bench_readings_A_oct22.csv)
         A2 = reading A2 (independent read of the same workbook images, P24-* copies; data/bench_readings_A_oct24.csv)
         B  = reading B, made without access to the other readings ( data/bench_readings_B.csv)
         T  = additional pixel-extent reading for R09/R10 (green-pixel extent minus line width; see logs)
Consensus = median of available centre-line pk-pk readings (divisions); spread = max-min.
Displayed V_pp = div x V/div (as shown on screen). Interpreted V_pp = displayed / 10
(trace DC level reads ~48 V with a 5.0 V supply -> screen shows ~10x tip voltage; inference).
Relative dB per run, normalised to that run's ~1 MHz point.
"""
import csv, json, math, os, statistics as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data")

A = {r["id"]: r for r in csv.DictReader(open(os.path.join(D, "bench_readings_A_oct22.csv")))}
A2 = {r["id"]: r for r in csv.DictReader(open(os.path.join(D, "bench_readings_A_oct24.csv")))}
B = {r["id"]: r for r in csv.DictReader(open(os.path.join(D, "bench_readings_B.csv")))}
T = {"R09": 0.921, "R10": 1.117}   # additional pixel-extent reading (centre-line pk-pk, divisions)

def f(x):
    try:
        return float(x)
    except Exception:
        return None

# point definitions: run, label freq (Hz), V/div, list of (reader, value)
pts = []
def add(run, fHz, vdiv, readings, photos, note=""):
    vals = [v for _, v in readings if v is not None]
    pts.append(dict(run=run, f_Hz=fHz, V_per_div=vdiv, readings=readings, photos=photos,
                    n_readers=len(vals), pkpk_div=st.median(vals),
                    spread_div=(max(vals) - min(vals)) if len(vals) > 1 else None, note=note))

for rid, fHz in [("R05", 99.8e3), ("R06", 199.6e3), ("R07", 300.5e3), ("R08", 399.0e3), ("R09", 500.8e3),
                 ("R10", 1003.2e3), ("R12", 2980e3), ("R14", 3993e3), ("R16", 5180e3)]:
    rd = [("A", f(A[rid]["pkpk_div_centre"])), ("B", f(B[rid]["pp_center_div"]))]
    if rid in T:
        rd.append(("T", T[rid]))
    add("R", fHz, 10.0, rd, [rid, "R04 (duplicate screen of R05)" if rid == "R05" else ""],
        "R04 shows the same 100 kHz screen" if rid == "R05" else "")

add("X", 45.4e3, 10.0, [("A", f(A["X10"]["pkpk_div_centre"])), ("A2", f(A2["P24-03"]["pkpk_div_best"])), ("B", f(B["X10"]["pp_center_div"]))], ["X10", "P24-03"], "about one trace width tall; low confidence")
add("X", 104.3e3, 10.0, [("A", f(A["X09"]["pkpk_div_centre"])), ("A2", f(A2["P24-05"]["pkpk_div_best"]))], ["X09", "P24-05"])
add("X", 500e3, 10.0, [("A", f(A["X08"]["pkpk_div_centre"])), ("A2", f(A2["P24-14"]["pkpk_div_best"]))], ["X08", "P24-14"])
add("X", 1030e3, 20.0, [("A", f(A["X05"]["pkpk_div_centre"])), ("A2", f(A2["P24-27"]["pkpk_div_best"])), ("A2b", f(A2["P24-26"]["pkpk_div_best"]))], ["X05", "P24-27", "P24-26"])
add("X", 2052e3, 20.0, [("A2", f(A2["P24-23"]["pkpk_div_best"]))], ["P24-23"], "single reader; not in workbook (sheet reuses 1017 kHz image)")
add("X", 3040e3, 20.0, [("A", f(A["X06"]["pkpk_div_centre"])), ("A2", f(A2["P24-20"]["pkpk_div_best"]))], ["X06", "P24-20"])
add("X", 4027e3, 20.0, [("A", f(A["X07"]["pkpk_div_centre"])), ("A2", f(A2["P24-17"]["pkpk_div_best"]))], ["X07", "P24-17"], "P24-17 was filed at 3006 kHz by reading A2; screen period gives ~4 MHz")
add("X", 5010e3, 20.0, [("A", f(A["X04"]["pkpk_div_centre"])), ("A2", f(A2["P24-32"]["pkpk_div_best"])), ("B", f(B["X04"]["pp_center_div"]))], ["X04", "P24-32"])

ref = {"R": next(p for p in pts if p["run"] == "R" and abs(p["f_Hz"] - 1.003e6) < 1e4),
       "X": next(p for p in pts if p["run"] == "X" and abs(p["f_Hz"] - 1.03e6) < 2e4)}
rows = []
for p in pts:
    disp = p["pkpk_div"] * p["V_per_div"]
    r = ref[p["run"]]
    rel = 20 * math.log10(disp / (r["pkpk_div"] * r["V_per_div"]))
    rel_spread = None
    if p["spread_div"] is not None:
        rel_spread = 20 * math.log10((p["pkpk_div"] + p["spread_div"] / 2) / max(p["pkpk_div"] - p["spread_div"] / 2, 1e-9)) / 2
    rows.append(dict(run=p["run"], f_kHz=round(p["f_Hz"] / 1e3, 1), V_per_div=p["V_per_div"],
                     readings=";".join(f"{n}={v}" for n, v in p["readings"]), n_readers=p["n_readers"],
                     pkpk_div=round(p["pkpk_div"], 3), spread_div=None if p["spread_div"] is None else round(p["spread_div"], 3),
                     Vpp_displayed=round(disp, 2), Vpp_interpreted_div10=round(disp / 10, 3),
                     rel_dB_re_1MHz=round(rel, 2), rel_dB_halfspread=None if rel_spread is None else round(rel_spread, 2),
                     photos=" ".join(x for x in p["photos"] if x), note=p["note"]))

with open(os.path.join(D, "bench_reconciled.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

vals = {}
for r in rows:
    k = f"bench.{r['run']}_{str(r['f_kHz']).replace('.', 'p')}k"
    vals[k + "_Vpp_disp"] = r["Vpp_displayed"]
    vals[k + "_Vpp_int"] = r["Vpp_interpreted_div10"]
    vals[k + "_reldB"] = r["rel_dB_re_1MHz"]
vals["bench.f_min_kHz"] = min(r["f_kHz"] for r in rows)
vals["bench.f_max_MHz"] = round(max(r["f_kHz"] for r in rows) / 1e3, 2)
vals["bench.n_points"] = len(rows)
vals["bench.dc_level_displayed_V"] = 48.3
vals["bench.dc_level_note"] = "trace DC level ~48 V displayed (Reader B; ground marker check on X04 48.7 V) vs 5.0 V supply"
json.dump(vals, open(os.path.join(D, "values_bench.json"), "w"), indent=1)
for r in rows:
    print(r["run"], r["f_kHz"], r["readings"], "->", r["pkpk_div"], r["Vpp_displayed"], r["Vpp_interpreted_div10"], r["rel_dB_re_1MHz"], r["rel_dB_halfspread"])
