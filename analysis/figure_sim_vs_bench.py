#!/usr/bin/env python3
"""Figures (IEEE single column, 3.5 in wide, >= 8 pt, greyscale-readable).
Inputs: work/models/curves.npz (written by circuit_models.py), data/bench_reconciled.csv, data/values_sim.json,
        data/ltspice_initial_model_curve_digitized.csv.
Outputs: paper/src/figures/fig_sim_vs_bench.pdf, paper/src/figures/fig_m0_audit.pdf (+ .png previews in work/models)."""
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
sys.path.insert(0, 'analysis')
import bench_data  # noqa: E402

ROOT = ''
OUT = ROOT + 'paper/src/figures/'
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'legend.fontsize': 8, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix',
                     'pdf.fonttype': 42, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'grid.linewidth': 0.4, 'grid.color': '0.85', 'lines.linewidth': 1.0})
z = np.load(ROOT + 'work/models/curves.npz')
V = json.load(open(ROOT + 'data/values_sim.json'))
F = z['f']


def ndb(h, f0=1e6):
    return 20 * np.log10(np.abs(h)) - 20 * np.log10(np.abs(np.interp(np.log10(f0), np.log10(F), np.abs(h))))


GRID = pd.read_csv(ROOT + 'simulation/results_M4_grid.csv')


def _lab(run):
    """legend text for the best in-range M4 case of a run (from simulation/results_M4_grid.csv)"""
    b = GRID.sort_values(f'rms_all_{run}').iloc[0]
    parts = [f'{b.Rg_ohm:g} \u03a9']
    if b.probe != 'none':
        parts.append({'x10': '\u00d710 probe', 'x1_65pF': '\u00d71 probe (65 pF)', 'x1_115pF': '\u00d71 probe'}[b.probe])
    if b.Lw_uH > 0:
        parts.append(f'$L_w$ {b.Lw_uH:g} \u00b5H')
    if b.L_RL_uH > 0:
        parts.append(f'$L_{{RL}}$ {b.L_RL_uH:g} \u00b5H')
    if b.base == 'M3':
        parts.append('dist. C')
    return ', '.join(parts)


# ------------------------------------------------------------------ fig_sim_vs_bench
fig, ax = plt.subplots(figsize=(3.5, 3.45))
_d24 = pd.read_csv(ROOT + 'data/bench_readings_A_oct24.csv')
_ring = list(_d24[_d24.gen_freq_setting_Hz == 100.1].freq_from_period_Hz.dropna())
_ring.append(json.load(open(ROOT + 'data/bench_readings_A_oct22.json'))['X01']['freq_from_screen_hz'])
ax.axvspan(min(_ring), max(_ring), color='0.88', lw=0, zorder=0)   # ringing bursts read from X01 and P24-43..48
ax.text(8.6e6, -27.0, 'bench\nringing', ha='left', va='bottom', fontsize=8, color='0.3')
ax.axvline(5e6, color='0.5', lw=0.6, ls=(0, (1, 2)), zorder=1)
ax.text(4.6e6, -27.0, 'generator\n5 MHz limit', ha='right', va='bottom', fontsize=8, color='0.3')
lines = [('M1_Rgeom', 'M1 geometry ($L,M$ theory, $R$=%.2f Ω)' % json.load(open(ROOT + 'data/values_theory.json'))['theory.R_A_w0p254_t35_ohm'],
          dict(color='0.55', ls=(0, (1, 1.2)), lw=1.2)),
         ('M2', 'M2 measured $L,M,R$ (magnetic)', dict(color='k', ls='-', lw=1.0)),
         ('M3', 'M3 = M2 + 274 pF distributed', dict(color='k', ls=(0, (4, 2)), lw=1.0)),
         ('M4_R', 'M4 hyp. R: ' + _lab('R'), dict(color='0.45', ls='-', lw=1.8)),
         ('M4_X', 'M4 hyp. X: ' + _lab('X'), dict(color='0.45', ls=(0, (5, 1.5, 1, 1.5)), lw=1.4))]
for key, lab, st in lines:
    ax.semilogx(F, ndb(z[key]), label=lab, zorder=3, **st)
runs = bench_data.load()
for run, mk, lab in (('R', dict(marker='o', mfc='white', mec='k', ms=4.5, mew=0.9), 'Bench run R (22 Oct photos)'),
                     ('X', dict(marker='^', mfc='k', mec='k', ms=4.5, mew=0.6), 'Bench run X (workbook images)')):
    pts = runs[run]['points']
    f = np.array([p['f'] for p in pts]); d = np.array([p['dB_norm'] for p in pts])
    e = np.nan_to_num(np.array([p['dB_half'] for p in pts]))
    ax.errorbar(f, d, yerr=e, ls='none', ecolor='k', elinewidth=0.6, capsize=1.5, zorder=5, label=lab, **mk)
ax.set_xlim(1e4, 3e7); ax.set_ylim(-28, 16)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Normalised $|V_{\\mathrm{OUT}}|$ (dB)')
ax.grid(True, which='major'); ax.grid(True, which='minor', lw=0.25, color='0.93')
ax.set_yticks(np.arange(-25, 20, 5))
h_, l_ = ax.get_legend_handles_labels()
ax.legend(h_, l_, loc='upper center', bbox_to_anchor=(0.42, -0.2), ncol=1, frameon=False, handlelength=2.6,
          borderaxespad=0, labelspacing=0.25)
fig.subplots_adjust(left=0.14, right=0.97, top=0.98, bottom=0.43)
fig.savefig(OUT + 'fig_sim_vs_bench.pdf')
fig.savefig(ROOT + 'work/models/fig_sim_vs_bench_preview.png', dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ fig_m0_audit
dig = pd.read_csv(ROOT + 'data/ltspice_initial_model_curve_digitized.csv')
fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.5, 4.6), gridspec_kw=dict(height_ratios=[1, 1.05], hspace=0.42))
a1.semilogx(dig.freq_Hz, dig.mag_dB_digitized, color='0.7', lw=3.2, label='17 Nov screenshot (digitised)', solid_capstyle='butt')
a1.semilogx(z['f_lt'], 20 * np.log10(np.abs(z['h_lt'])), color='k', lw=0.8, marker='.', ms=3,
            label='M0 reproduced (.ac dec 10 points)')
a1.semilogx(z['f_full'], 20 * np.log10(np.abs(z['h0_full'])), color='k', lw=0.8, ls=(0, (3, 2)), label='M0 reproduced (fine grid)')
a1.set_xlim(1, 1e12); a1.set_ylim(-100, 20)
a1.set_xticks([1, 1e3, 1e6, 1e9, 1e12])
a1.set_ylabel('dB re 1 V source')
a1.set_xlabel('Frequency (Hz)')
a1.grid(True)
a1.legend(loc='lower center', frameon=False, labelspacing=0.2, handlelength=2.4)
a1.text(0.02, 0.95, f"(a) RMS {V['sim.M0_rms_vs_digitized_dB']:.2f} dB", transform=a1.transAxes, va='top', fontsize=8)
for key, lab, st in (('M0', 'M0 as drawn (no load, $R_s$ = 1 Ω)', dict(color='k', ls=(0, (3, 2)), lw=0.9)),
                     ('audit_S2', 'M0 + $R_L$ 111.6 Ω + $R_g$ 50 Ω', dict(color='0.5', ls=(0, (5, 1.5, 1, 1.5)), lw=1.2)),
                     ('M2', 'M2 (measured, bench topology)', dict(color='k', ls='-', lw=1.1)),
                     ('M3', 'M3 (M2 + distributed C)', dict(color='0.55', ls='-', lw=1.8))):
    a2.semilogx(F, 20 * np.log10(np.abs(z[key])), label=lab, **st)
a2.set_xlim(1e4, 3e7); a2.set_ylim(-30, 20)
a2.set_xlabel('Frequency (Hz)'); a2.set_ylabel('dB re generator EMF')
a2.grid(True)
a2.legend(loc='upper left', frameon=False, labelspacing=0.2, handlelength=2.4, fontsize=8)
a2.text(0.98, 0.05, '(b)', transform=a2.transAxes, ha='right', fontsize=8)
fig.subplots_adjust(left=0.16, right=0.97, top=0.98, bottom=0.09)
fig.savefig(OUT + 'fig_m0_audit.pdf')
fig.savefig(ROOT + 'work/models/fig_m0_audit_preview.png', dpi=200)
plt.close(fig)
print('figures written')
