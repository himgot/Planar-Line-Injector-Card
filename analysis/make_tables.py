#!/usr/bin/env python3
"""Generate LaTeX tables from data files (no hand-typed numbers)."""
import csv, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(os.path.join(ROOT, "data", "bench_reconciled.csv"))))
out = [r"\begin{table}[!t]", r"\centering",
       r"\caption{Bench injection results (22 Oct 2026), $V_\mathrm{DC}=\SI{5}{V}$, $R_L\approx\SI{111.6}{\ohm}$}",
       r"\label{tab:bench}", r"\footnotesize", r"\setlength{\tabcolsep}{3pt}",
       r"\begin{tabular}{@{}c r c c r r r l@{}}", r"\toprule",
       r"Run & $f$ (kHz) & V/div & $n$ & $V_\mathrm{pp}$ disp. & $V_\mathrm{pp}$ int.\textsuperscript{a} & dB re 1\,MHz\textsuperscript{b} & Photo \\",
       r"\midrule"]
prev = None
for r in rows:
    if prev and r["run"] != prev:
        out.append(r"\midrule")
    prev = r["run"]
    hs = r["rel_dB_halfspread"]
    rel = f"{float(r['rel_dB_re_1MHz']):+.1f}".replace("-", "$-$").replace("+", "$+$")
    photo = r["photos"].split()[0]
    rd = [x.split("=")[0] for x in r["readings"].split(";") if x]
    indep = len({x.rstrip("b") for x in rd if x != "T"})
    tb = "+T" if "T" in rd else ""
    flag = r"\gap{*}" if (r["note"] and ("low confidence" in r["note"] or "single reader" in r["note"])) else ""
    out.append(f"{r['run']} & {float(r['f_kHz']):.1f} & {float(r['V_per_div']):.0f} & {indep}{tb} & {float(r['Vpp_displayed']):.2f} & \\gap{{{round(float(r['Vpp_displayed']), 2) / 10:.3f}}} & {rel} & {photo}{flag} \\\\")
out += [r"\midrule",
        r"X & 501 & 20 & -- & \multicolumn{3}{c}{no injected sine visible; excluded} & X02 \\",
        r"X & 1007 & -- & -- & \multicolumn{3}{c}{no injected sine visible; excluded} & X03 \\",
        r"\bottomrule", r"\end{tabular}", r"\par\vspace{2pt}\parbox{\columnwidth}{\scriptsize",
        r"$n$: independent automated readings (no manual reading); +T: non-blind tie-break reading; consensus = median. Run R: raw 22 Oct photos; Run X: workbook images (the ``24 Oct'' phone files are copies; the 2052\,kHz image exists only among the phone files). \gap{*}~low confidence or single reader. ",
        r"\textsuperscript{a}Displayed value $\div 10$: the trace DC level reads $\approx\SI{48}{V}$ on screen with a \SI{5.0}{V} supply; \gap{probe and scope probe-factor settings were not recorded}. ",
        r"\textsuperscript{b}Each run normalised to its own $\approx$1\,MHz point; computed from unrounded readings. Uncertainty: display-pixel quantisation (about 25 px/div) $\approx\pm$1\,dB at 100\,kHz (trace $\approx$0.3\,div), less at higher amplitude; the Run R reference (1003\,kHz) carries +0.8/$-$0.5\,dB reader spread.}",
        r"\end{table}"]
open(os.path.join(ROOT, "paper", "src", "tables", "tab_bench.tex"), "w").write("\n".join(out) + "\n")
print("tab_bench.tex written,", len(rows), "rows")
