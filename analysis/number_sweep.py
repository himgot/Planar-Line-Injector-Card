#!/usr/bin/env python3
"""G6 number sweep: every number in the compiled PDF text is checked against
data/values.json (any formatting/rounding), raw input tables and a whitelist of
fixed constants (nominal settings, years, reference numbers). Untraced numbers are listed."""
import json, os, re, subprocess, csv, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
txt = subprocess.run(["pdftotext", os.path.join(ROOT, "paper", "src", "main_fallback.pdf"), "-"], capture_output=True, text=True).stdout
body = txt.split("REFERENCES")[0] if "REFERENCES" in txt else txt
vals = json.load(open(os.path.join(ROOT, "data", "values.json")))["values"]
pool = set()
def add(x):
    try: f = float(x)
    except Exception: return
    for d in range(0, 5):
        pool.add(round(abs(f), d)); pool.add(round(abs(f) * 1e3, d)); pool.add(round(abs(f) / 1e3, d)); pool.add(round(abs(f)/1e6, d)); pool.add(round(abs(f)*1e6, d)); pool.add(round(abs(f)*100, d))
for v in vals.values():
    if isinstance(v, (int, float)): add(v)
    elif isinstance(v, list):
        for x in v: add(x)
for fn in glob.glob(os.path.join(ROOT, "data", "*.csv")):
    for row in csv.reader(open(fn)):
        for c in row:
            for m in re.findall(r"-?\d+\.?\d*", c): add(m)
whitelist = {0,1,2,3,4,5,6,7,8,9,10,12,15,16,17,18,19,20,22,23,24,25,27,30,40,45,50,100,200,300,400,500,600,1000,250,
             2026,2023,2009,1999,1974,1928,1908,1946,1972,2004,2011,2024,2025,1.6,0.5,0.1,4.0,3.0,1.0,2.0,0.2,
             274,111.6,7.31,17.9,17.94,17.96,85.6,54.0,0.254,0.228,0.18,0.73,1.7,5.8,9.7,7.5,8.2,7.7,8.6,4.69,
             0.03,48,5.01,5.18,5.0,3.6,1.2}
nums = re.findall(r"(?<![A-Za-z])-?\d+(?:\.\d+)?", body)
untraced = {}
for n in nums:
    f = abs(float(n))
    if f in whitelist: continue
    if any(abs(f - p) <= 1e-9 for p in pool if isinstance(p, float) or isinstance(p, int)) : continue
    untraced[n] = untraced.get(n, 0) + 1
out = os.path.join(ROOT, "work", "number_sweep.txt")
with open(out, "w") as fh:
    fh.write(f"numbers in body: {len(nums)}; untraced (not in values.json, data CSVs or whitelist): {len(untraced)}\n")
    for k, c in sorted(untraced.items(), key=lambda x: float(x[0])): fh.write(f"{k}\t{c}\n")
print(open(out).read()[:3000])
