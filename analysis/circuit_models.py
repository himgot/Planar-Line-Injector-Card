#!/usr/bin/env python3
"""Circuit models: M0 audit, M1 geometry, M2 measured, M3 distributed, M4 hypotheses, bench comparison.
Run:  python3 -I analysis/circuit_models.py
Requires data/model_inputs/peec_sections_N24.npz and _N48.npz (analysis/a5_peec_sections.py 24 48).
Writes: simulation/*.cir, simulation/results_*.csv, simulation/audit_of_17Nov_model.md,
        data/values_sim.json, work/models/curves.npz (for the figure script)."""
import sys
import json
import os
import itertools
import numpy as np
import pandas as pd
sys.path.insert(0, 'analysis')
from mna_sim import (parse_spice, bench_lumped, bench_sectioned, scaled_section_matrix,  # noqa: E402
                     uniform_section_matrix, transfer, ring_poles, Circuit)
import bench_data  # noqa: E402

ROOT = ''
SIM = ROOT + 'simulation/'
W5 = ROOT + 'work/models/'
os.makedirs(W5 or '.', exist_ok=True)
V = {}          # values_sim.json


def put(key, val, meta):
    if isinstance(val, (np.floating, np.integer)):
        val = val.item()
    if isinstance(val, float):
        val = None if not np.isfinite(val) else float(f'{val:.6g}')
    V['sim.' + key] = val
    V['sim.' + key + '_meta'] = meta


def dB(h):
    return 20 * np.log10(np.abs(h))


# ------------------------------------------------------------------ inputs
meas = json.load(open(ROOT + 'data/values_meas.json'))
pc = json.load(open(ROOT + 'data/values_peec_check.json'))
theory = json.load(open(ROOT + 'data/values_theory.json'))
LA_m, LB_m = meas['meas.L_A_open_10k'] * 1e-6, meas['meas.L_B_open_10k'] * 1e-6
M_m, k_m, R_m = meas['meas.M_pair_mean_10k'] * 1e-6, meas['meas.k_pair_mean_10k'], meas['meas.R_dc_ohm']
C_m = meas['meas.C_AB_pF'] * 1e-12
# PEEC check: polyline 'B' = AC winding (starts at AC_RTN pad), 'A' = DC winding (starts at DC_IN pad)
ind = pc['inductance_top_layer_only']
LA_v2, LB_v2, M_v2 = ind['LB_uH'] * 1e-6, ind['LA_uH'] * 1e-6, ind['M_uH'] * 1e-6
RA_v2 = pc['windings']['B']['R_dc']['w0.254mm_t35um']['top_spiral_ohm']
RB_v2 = pc['windings']['A']['R_dc']['w0.254mm_t35um']['top_spiral_ohm']
# Neumann theory (preferred, values_theory.json; A = AC winding, D = DC winding; top + via + bottom route)
LA_g, LB_g, M_g = theory['theory.L_A_uH'] * 1e-6, theory['theory.L_D_uH'] * 1e-6, theory['theory.M_uH'] * 1e-6
RA_g, RB_g = theory['theory.R_A_w0p254_t35_ohm'], theory['theory.R_D_w0p254_t35_ohm']
RL = 111.6
SRC_M = ('L_A=meas.L_A_open_10k, L_B=meas.L_B_open_10k, M=meas.M_pair_mean_10k, R=meas.R_dc_ohm '
         '(data/values_meas.json)')
SRC_G = ('data/values_theory.json theory.L_A_uH, theory.L_D_uH, theory.M_uH (Neumann, top+via+bottom route) and '
         'theory.R_A/R_D_w0p254_t35_ohm')
SRC_V2 = ('data/values_peec_check.json inductance_top_layer_only (PEEC check; its B = AC winding, A = DC winding) and '
          'R_dc w0.254mm_t35um top_spiral')
assert theory['theory.geom.A_inner_terminal'] == 'AC_IN' and theory['theory.geom.D_inner_terminal'] == 'OUT'

F = np.logspace(4, np.log10(30e6), int(np.ceil(np.log10(3e3) * 250)) + 1)    # 10 kHz-30 MHz, 250 pts/dec
put('sweep_points_per_decade', 250, '[MODEL] AC sweep 10 kHz-30 MHz, log spaced, 871 points')
runs = bench_data.load()


def bench_err(h_fun, run):
    """normalised model at bench points; returns dict with residuals and RMS (dB)"""
    pts = runs[run]['points']
    f = np.array([p['f'] for p in pts]); b = np.array([p['dB_norm'] for p in pts])
    h = h_fun(np.r_[f, runs[run]['f_ref']])
    m = dB(h[:-1]) - dB(h[-1:])
    r = m - b
    isref = np.abs(f - runs[run]['f_ref']) < 1
    lo = (f >= 0.099e6) & (f < 0.99e6) & ~isref
    hi = (f > 1.05e6) & (f <= 5.3e6) & ~isref
    al = (f >= 0.099e6) & (f <= 5.3e6) & ~isref
    off = np.mean(r[(f >= 0.099e6) & (f <= 5.3e6)])
    return dict(model=m, resid=r, f=f,
                lo=float(np.sqrt(np.mean(r[lo] ** 2))), hi=float(np.sqrt(np.mean(r[hi] ** 2))),
                all=float(np.sqrt(np.mean(r[al] ** 2))),
                freegain=float(np.sqrt(np.mean((r[(f >= 0.099e6) & (f <= 5.3e6)] - off) ** 2))))


# ================================================================== M0
m0_txt = open(SIM + 'm0_initial_model.cir').read()
c0, dots0, prm0 = parse_spice(m0_txt, title='M0 as reconstructed')
c0.meas = 'OUTNODE'
dig = pd.read_csv(ROOT + 'data/ltspice_initial_model_curve_digitized.csv')
f_lt = 10 ** np.arange(0, 12.0001, 0.1)               # .ac dec 10 1 1000G sample points
h_lt = transfer(c0, f_lt)
interp_lt = np.interp(np.log10(dig.freq_Hz), np.log10(f_lt), dB(h_lt))   # LTspice draws straight lines
rms0 = np.sqrt(np.mean((interp_lt - dig.mag_dB_digitized) ** 2))
c0v, _, _ = parse_spice(m0_txt, params={'RAC_SER': 1e-9, 'RDC_SER': 1e-9}); c0v.meas = 'OUTNODE'
rms0v = np.sqrt(np.mean((np.interp(np.log10(dig.freq_Hz), np.log10(f_lt), dB(transfer(c0v, f_lt)))
                         - dig.mag_dB_digitized) ** 2))
f_full = np.logspace(0, 12, 12 * 100 + 1)
h0_full = transfer(c0, f_full)
h0 = transfer(c0, F)
i0 = np.argmax(dB(h0_full))
put('M0_rms_vs_digitized_dB', rms0, '[MODEL] simulation/m0_initial_model.cir solved with analysis/mna_sim.py at the '
    '.ac dec 10 points, linear in log f between points (as LTspice draws), vs data/ltspice_initial_model_curve_digitized.csv (1389 pts)')
put('M0_rms_vs_digitized_visible_only_dB', rms0v, '[MODEL] same, with the inductor series resistances RAC_SER=RDC_SER=0 (values visible on the schematic only)')
put('M0_peak_f_Hz', f_full[i0], '[MODEL] M0 fine grid (100 pts/dec) peak of |V(OUTNODE)|')
put('M0_peak_dB', dB(h0_full[i0]), '[MODEL] dB re 1 V AC source, M0 fine grid')
put('M0_dB_1MHz', dB(transfer(c0, [1e6]))[0], '[MODEL] M0 |V(OUTNODE)| at 1 MHz, dB re 1 V source')

# ================================================================== M1 / M2 lumped
models = {}           # name -> function f -> complex H (for bench comparison)
curves = {'f': F}


def lumped(LA, LB, M, RA, RB, **kw):
    return bench_lumped(LA, LB, M, RA, RB, **kw)


cM1a = lumped(LA_g, LB_g, M_g, RA_g, RB_g, Rg=50, RL=RL, title='M1 geometry R from geometry')
cM1b = lumped(LA_g, LB_g, M_g, R_m, R_m, Rg=50, RL=RL, title='M1 geometry R=7.31')
cM1v = lumped(LA_v2, LB_v2, M_v2, RA_v2, RB_v2, Rg=50, RL=RL, title='M1 cross-check with the PEEC check values')
cM2 = lumped(LA_m, LB_m, M_m, R_m, R_m, Rg=50, RL=RL, title='M2 measured')
cM2_600 = lumped(LA_m, LB_m, M_m, R_m, R_m, Rg=600, RL=RL, title='M2 measured Rg=600')
for nm, c in (('M1_Rgeom', cM1a), ('M1_R7p31', cM1b), ('M1_peec', cM1v), ('M2', cM2), ('M2_Rg600', cM2_600)):
    models[nm] = (lambda cc: (lambda f: transfer(cc, f)))(c)

# ================================================================== M3 sectioned
z24 = np.load(ROOT + 'data/model_inputs/peec_sections_N24.npz'); z48 = np.load(ROOT + 'data/model_inputs/peec_sections_N48.npz')
Ls24, sc = scaled_section_matrix(z24['Lsec'], LA_m, LB_m, M_m)
Ls48, _ = scaled_section_matrix(z48['Lsec'], LA_m, LB_m, M_m)
cM3 = bench_sectioned(Ls24, z24['Cn'], C_m, R_m, R_m, Rg=50, RL=RL, title='M3 N=24')
cM3_600 = bench_sectioned(Ls24, z24['Cn'], C_m, R_m, R_m, Rg=600, RL=RL, title='M3 N=24 Rg=600')
cM3_48 = bench_sectioned(Ls48, z48['Cn'], C_m, R_m, R_m, Rg=50, RL=RL, title='M3 N=48')
cM3_u = bench_sectioned(uniform_section_matrix(24, LA_m, LB_m, M_m, 0.5), z24['Cn'], C_m, R_m, R_m, Rg=50, RL=RL,
                        title='M3 uniform-coupling fallback N=24')
cM3_0 = bench_sectioned(Ls24, z24['Cn'], 0.0, R_m, R_m, Rg=50, RL=RL, title='M3 with C=0')
for nm, c in (('M3', cM3), ('M3_Rg600', cM3_600)):
    models[nm] = (lambda cc: (lambda f: transfer(cc, f)))(c)
for nm, c in (('M0', c0), ('M1_Rgeom', cM1a), ('M1_R7p31', cM1b), ('M1_peec', cM1v), ('M2', cM2), ('M2_Rg600', cM2_600), ('M3', cM3),
              ('M3_Rg600', cM3_600), ('M3_N48', cM3_48), ('M3_uniform', cM3_u), ('M3_C0', cM3_0)):
    curves[nm] = transfer(c, F)

for key, val in sc.items():
    put(f'M3_scale_{key}', val, '[MODEL] scaling of the PEEC section blocks (data/model_inputs/peec_sections_N24.npz) so the block sums '
        'equal measured L_A, L_B, M; rho = gamma/sqrt(alpha*beta) <= 1 keeps the matrix positive definite')
put('M3_N', 24, '[MODEL] sections per winding (~1 turn each); N=48 used as convergence check')
put('M3_min_eig_nH', np.linalg.eigvalsh(Ls24).min() * 1e9, '[MODEL] smallest eigenvalue of the scaled 48x48 inductance matrix')
peec_sum = json.load(open(ROOT + 'data/model_inputs/peec_sections_summary.json'))['24']
put('M3_PEEC_sum_L_AC_uH', peec_sum['L_AC_uH'], '[MODEL] block sum of the unscaled section matrix; equals the PEEC check LB_uH (AC winding)')
put('M3_PEEC_sum_M_uH', peec_sum['M_uH'], '[MODEL] block sum of the unscaled section matrix; equals the PEEC check M_uH')
put('M3_facing_length_mm', peec_sum['facing_length_mm'], '[CALCULATED] AC-DC facing length from 0.1 mm resampling, nearest neighbour per side '
    '(PEEC check coupled_gap_length 5633 mm, geometry adjacent_length_AD 5694 mm)')
put('M3_C_fraction_direct_ACIN_OUT', peec_sum['frac_C_ACIN_OUT_direct'], '[MODEL] fraction of the 274 pF that lands directly between the AC_IN and OUT nodes (N=24 hat-function lumping)')
mask30 = F <= 30e6
put('M3_N24_vs_N48_maxdiff_dB', np.max(np.abs(dB(curves['M3']) - dB(curves['M3_N48']))[mask30]),
    '[MODEL] max |dB difference| N=24 vs N=48 over 10 kHz-30 MHz (convergence check)')
m5 = F <= 5.3e6
put('M3_N24_vs_N48_maxdiff_to5MHz_dB', np.max(np.abs(dB(curves['M3']) - dB(curves['M3_N48']))[m5]), '[MODEL] same, 10 kHz-5.3 MHz')
put('M3_PEEC_vs_uniform_maxdiff_to5MHz_dB', np.max(np.abs(dB(curves['M3']) - dB(curves['M3_uniform']))[m5]),
    '[MODEL] PEEC-sectioned vs uniform-coupling fallback (cw=0.5), 10 kHz-5.3 MHz')
put('M3_PEEC_vs_uniform_maxdiff_dB', np.max(np.abs(dB(curves['M3']) - dB(curves['M3_uniform']))[mask30]),
    '[MODEL] PEEC-sectioned vs uniform-coupling fallback (cw=0.5), 10 kHz-30 MHz')
put('M3_C0_vs_M2_maxdiff_dB', np.max(np.abs(dB(curves['M3_C0']) - dB(curves['M2']))),
    '[MODEL] check: sectioned model with C=0 reproduces lumped M2 (should be ~0)')

for nm in ('M1_Rgeom', 'M1_R7p31', 'M1_peec', 'M2', 'M2_Rg600', 'M3', 'M3_Rg600'):
    h = transfer({'M1_Rgeom': cM1a, 'M1_R7p31': cM1b, 'M1_peec': cM1v, 'M2': cM2, 'M2_Rg600': cM2_600, 'M3': cM3,
                  'M3_Rg600': cM3_600}[nm],
                 [1e5, 1e6, 5e6])
    for fv, lab in zip(h, ('100kHz', '1MHz', '5MHz')):
        put(f'{nm}_dB_{lab}', dB(fv), f'[MODEL] |V(OUT)/V_g| in dB re generator EMF, bench topology, R_L=111.6 ohm, '
            f'{"R_g=600" if "600" in nm else "R_g=50"} ohm, no probe loading; '
            + ((SRC_V2 if nm == 'M1_peec' else SRC_G) if nm.startswith('M1') else SRC_M + ('; C_AB=274 pF distributed (N=24)' if nm.startswith('M3') else '')))
# low-frequency corner (-3 dB below the maximum of the magnetic-only curve, low side)
for nm in ('M2', 'M2_Rg600', 'M1_Rgeom'):
    d = dB(curves[nm]); im = np.argmax(d)
    fl = F[:im][np.argmin(np.abs(d[:im] - (d[im] - 3)))]
    put(f'{nm}_f_minus3dB_low_Hz', fl, '[MODEL] lower -3 dB point relative to the curve maximum (10 kHz-30 MHz)')
    put(f'{nm}_max_dB', d[im], '[MODEL] maximum of |H| over 10 kHz-30 MHz')
    put(f'{nm}_f_max_Hz', F[im], '[MODEL] frequency of the maximum')
for nm, c in (('M3', cM3), ('M3_Rg600', cM3_600)):
    rp = ring_poles(c, 1e6, 60e6)
    put(f'{nm}_first_pole_f_Hz', rp[0][0], '[MODEL] lowest complex natural frequency (damped, Hz) of the network, generalized eigenvalues of (G, -C)')
    put(f'{nm}_first_pole_Q', rp[0][1], '[MODEL] Q = 1/(2 zeta) of that pole')

# ================================================================== M4 grid
grid = []
probes = [('none', None, 0.0), ('x10', 'x10', 0.0), ('x1_65pF', 'x1', 50e-12), ('x1_115pF', 'x1', 100e-12)]
for base, Rg, (pname, pk, cc), Lw, Lrl in itertools.product(('M2', 'M3'), (50.0, 600.0), probes,
                                                            (0.0, 0.25e-6, 0.5e-6, 1e-6),
                                                            (0.0, 0.5e-6, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6)):
    if pk is None and Lw > 0:
        continue
    if base == 'M2':
        c = bench_lumped(LA_m, LB_m, M_m, R_m, R_m, Rg=Rg, RL=RL, probe=pk, Lw=Lw, L_RL=Lrl, c_cable=cc)
    else:
        c = bench_sectioned(Ls24, z24['Cn'], C_m, R_m, R_m, Rg=Rg, RL=RL, probe=pk, Lw=Lw, L_RL=Lrl, c_cable=cc)
    hf = (lambda cc_: (lambda f: transfer(cc_, f)))(c)
    row = dict(base=base, Rg_ohm=Rg, probe=pname, Lw_uH=Lw * 1e6, L_RL_uH=Lrl * 1e6)
    for run in runs:
        e = bench_err(hf, run)
        row.update({f'rms_lo_{run}': e['lo'], f'rms_hi_{run}': e['hi'], f'rms_all_{run}': e['all'],
                    f'rms_freegain_{run}': e['freegain']})
    rp = [p for p in ring_poles(c, 2e6, 30e6) if p[1] >= 0.5]
    near = [p for p in rp if 7.0e6 <= p[0] <= 8.7e6]
    row['pole_f_MHz'] = round(rp[0][0] / 1e6, 3) if rp else np.nan
    row['pole_Q'] = round(rp[0][1], 2) if rp else np.nan
    row['ring_7to8p7MHz_Q2'] = bool(any(p[1] >= 2 for p in near))
    row['H_1MHz'] = abs(hf([1e6])[0])
    grid.append(row)
G = pd.DataFrame(grid)
G.to_csv(SIM + 'results_M4_grid.csv', index=False, float_format='%.4g')
best = {}
for run in runs:
    b = G.sort_values(f'rms_all_{run}').iloc[0]
    best[run] = b
    tag = f'run{run}'
    desc = f"base={b.base}, R_g={b.Rg_ohm:g} ohm, probe={b.probe}, L_w={b.Lw_uH:g} uH, L_RL={b.L_RL_uH:g} uH"
    put(f'M4_best_{tag}_params', desc, f'[MODEL] hypothesis: grid point (simulation/results_M4_grid.csv, {len(G)} cases) with the '
        f'lowest RMS over 0.1-5.3 MHz for Run {run}; NOT a conclusion')
    for k_ in ('lo', 'hi', 'all'):
        put(f'M4_best_{tag}_rms_{k_}', b[f'rms_{k_}_{run}'], f'[MODEL] RMS dB error of that hypothesis vs Run {run} ({k_} band)')
    put(f'M4_best_{tag}_pole_f_Hz', b.pole_f_MHz * 1e6, '[MODEL] lowest complex natural frequency of that hypothesis (Q>=0.5)')
    put(f'M4_best_{tag}_pole_Q', b.pole_Q, '[MODEL] its Q')
    ring = G[G.ring_7to8p7MHz_Q2].sort_values(f'rms_all_{run}')
    put(f'M4_n_cases_ring_7to8p7MHz', int(G.ring_7to8p7MHz_Q2.sum()), '[MODEL] grid cases with a pole 7.0-8.7 MHz and Q>=2')
    if len(ring):
        r0 = ring.iloc[0]
        put(f'M4_best_ringing_{tag}_params', f"base={r0.base}, R_g={r0.Rg_ohm:g}, probe={r0.probe}, L_w={r0.Lw_uH:g} uH, "
            f"L_RL={r0.L_RL_uH:g} uH; pole {r0.pole_f_MHz:.2f} MHz Q={r0.pole_Q:.1f}; rms lo/hi "
            f"{r0[f'rms_lo_{run}']:.2f}/{r0[f'rms_hi_{run}']:.2f} dB",
            '[MODEL] best-fitting grid case that ALSO has a 7.0-8.7 MHz pole with Q>=2 (hypothesis only)')


def build_case(b):
    pk, cc = {'none': (None, 0), 'x10': ('x10', 0), 'x1_65pF': ('x1', 50e-12), 'x1_115pF': ('x1', 100e-12)}[b['probe']]
    if b['base'] == 'M2':
        return bench_lumped(LA_m, LB_m, M_m, R_m, R_m, Rg=b['Rg_ohm'], RL=RL, probe=pk, Lw=b['Lw_uH'] * 1e-6,
                            L_RL=b['L_RL_uH'] * 1e-6, c_cable=cc, title=f"M4 hypothesis {b['base']}")
    return bench_sectioned(Ls24, z24['Cn'], C_m, R_m, R_m, Rg=b['Rg_ohm'], RL=RL, probe=pk, Lw=b['Lw_uH'] * 1e-6,
                           L_RL=b['L_RL_uH'] * 1e-6, c_cable=cc, title=f"M4 hypothesis {b['base']}")


cM4 = {run: build_case(best[run]) for run in runs}
for run in runs:
    models[f'M4_{run}'] = (lambda cc: (lambda f: transfer(cc, f)))(cM4[run])
    curves[f'M4_{run}'] = transfer(cM4[run], F)

# continuous R_g fit (magnetic-only M2, no probe), low band
for run in runs:
    rg = np.logspace(1, 3.3, 200)
    e = [bench_err(lambda f, r=r: transfer(lumped(LA_m, LB_m, M_m, R_m, R_m, Rg=r, RL=RL), f), run)['lo'] for r in rg]
    put(f'Rg_fit_run{run}_lowband_ohm', rg[int(np.argmin(e))], f'[MODEL] source resistance minimising the 0.1-1 MHz RMS of M2 '
        f'(magnetic only, no probe) vs Run {run}; diagnostic only (true SM5072 switch state UNKNOWN)')
    put(f'Rg_fit_run{run}_lowband_rms_dB', min(e), '[MODEL] RMS at that R_g')

# ================================================================== extended-range check (OUTSIDE the brief's ranges)
ext = []
for base, Rg, (pname, pk, cc), Lw, Lrl in itertools.product(('M2', 'M3'), (50.0, 600.0), probes,
                                                            (0.0, 1e-6, 1.5e-6, 2e-6, 3e-6), (5e-6, 7e-6, 10e-6, 15e-6)):
    if pk is None and Lw > 0:
        continue
    b_ = dict(base=base, Rg_ohm=Rg, probe=pname, Lw_uH=Lw * 1e6, L_RL_uH=Lrl * 1e6)
    c = build_case(b_)
    hf = (lambda cc_: (lambda f: transfer(cc_, f)))(c)
    for run in runs:
        e = bench_err(hf, run)
        b_.update({f'rms_lo_{run}': e['lo'], f'rms_hi_{run}': e['hi'], f'rms_all_{run}': e['all']})
    rp = [p for p in ring_poles(c, 2e6, 30e6) if p[1] >= 0.5]
    b_['pole_f_MHz'] = round(rp[0][0] / 1e6, 3) if rp else np.nan
    b_['pole_Q'] = round(rp[0][1], 2) if rp else np.nan
    ext.append(b_)
E = pd.DataFrame(ext)
E.to_csv(SIM + 'results_M4_grid_extended.csv', index=False, float_format='%.4g')
best_ext = {}
for run in runs:
    b = E.sort_values(f'rms_all_{run}').iloc[0]
    best_ext[run] = b
    put(f'M4ext_best_run{run}_params', f"base={b.base}, R_g={b.Rg_ohm:g} ohm, probe={b.probe}, L_w={b.Lw_uH:g} uH, L_RL={b.L_RL_uH:g} uH; "
        f"pole {b.pole_f_MHz:.2f} MHz Q={b.pole_Q:.1f}; rms lo/hi {b[f'rms_lo_{run}']:.2f}/{b[f'rms_hi_{run}']:.2f} dB",
        '[MODEL] hypothesis with L_w up to 3 uH and L_RL up to 15 uH, i.e. OUTSIDE the ranges set in the task brief '
        '(simulation/results_M4_grid_extended.csv); shows whether the in-range optimum sits on a range boundary')
    put(f'M4ext_best_run{run}_rms_all', b[f'rms_all_{run}'], '[MODEL] RMS 0.1-5.3 MHz of that extended-range case')
cM4ext = {run: build_case(best_ext[run]) for run in runs}
for run in runs:
    models[f'M4ext_{run}'] = (lambda cc: (lambda f: transfer(cc, f)))(cM4ext[run])
    curves[f'M4ext_{run}'] = transfer(cM4ext[run], F)

# ================================================================== bench comparison table
rows = []
for nm, hf in models.items():
    for run in runs:
        e = bench_err(hf, run)
        tag = f'run{run}'
        put(f'rms_err_{nm}_{tag}_lowband', e['lo'], f'[MODEL] RMS of (model - bench) in dB, both normalised to the Run {run} '
            f'reference point ({runs[run]["f_ref"] / 1e3:g} kHz), bench points 0.1-1 MHz excluding the reference; '
            'bench = data/bench_reconciled.csv rel_dB_re_1MHz')
        put(f'rms_err_{nm}_{tag}_highband', e['hi'], f'[MODEL] same, bench points 1.05-5.3 MHz')
        put(f'rms_err_{nm}_{tag}_freegain', e['freegain'], '[MODEL] RMS after removing the mean offset over 0.1-5.3 MHz '
            '(insensitive to the reference-point reading)')
        for p, mdl, r in zip(runs[run]['points'], e['model'], e['resid']):
            rows.append(dict(model=nm, run=run, point=p['id'], f_Hz=p['f'], bench_dB=p['dB_norm'], bench_halfspread_dB=p['dB_half'],
                             model_dB_norm=mdl, resid_dB=r))
pd.DataFrame(rows).to_csv(SIM + 'results_bench_compare.csv', index=False, float_format='%.4g')

# ================================================================== absolute (EMF) consistency check
erows = []
for nm in ('M2', 'M2_Rg600', 'M3', 'M3_Rg600', 'M4_R', 'M4_X', 'M4ext_R', 'M4ext_X'):
    for run in runs:
        pts = [p for p in runs[run]['points'] if 0.099e6 <= p['f'] <= 5.3e6]
        H = np.abs(models[nm](np.array([p['f'] for p in pts])))
        emf_as = np.array([p['Vpp'] for p in pts]) / H
        emf_10 = np.array([p['Vpp_div10'] for p in pts]) / H
        for p, h_, a_, t_ in zip(pts, H, emf_as, emf_10):
            erows.append(dict(model=nm, run=run, point=p['id'], f_Hz=p['f'], H_abs=h_, Vpp_displayed=p['Vpp'],
                              EMF_Vpp_as_displayed=a_, EMF_Vpp_div10=t_))
        put(f'EMF_{nm}_run{run}_div10_median_Vpp', np.median(emf_10), f'[MODEL]/[CALCULATED] generator EMF implied by {nm} if the '
            f'screen shows 10x the tip voltage (Vpp_interpreted_div10 / |H|), median over Run {run} points 0.1-5.3 MHz')
        put(f'EMF_{nm}_run{run}_div10_range_Vpp', [float(f'{emf_10.min():.3g}'), float(f'{emf_10.max():.3g}')], '[CALCULATED] min, max')
        put(f'EMF_{nm}_run{run}_asis_median_Vpp', np.median(emf_as), f'[CALCULATED] same, displayed V_pp taken as the tip voltage')
pd.DataFrame(erows).to_csv(SIM + 'results_emf_check.csv', index=False, float_format='%.4g')
vdc = 5.0 * RL / (RL + R_m)
put('DC_out_expected_V', vdc, '[CALCULATED] 5 V x R_L/(R_L + R_DC) with R_L=111.6, R_DC=7.31 ohm (ideal 5 V rail)')
put('DC_out_expected_displayed_if_x10_V', 10 * vdc, '[CALCULATED] x10 of the above; compare with the displayed DC level '
    '46-56 V inferred from the screen readings (data/bench_readings_A_oct22.json _dc_level_inference) and the 48.8 V trigger level')

# ================================================================== k sweep at 1 MHz
ks = np.round(np.arange(0.30, 0.9701, 0.01), 2)
krows = []
for kk in ks:
    Mk = kk * np.sqrt(LA_m * LB_m)
    h50 = abs(transfer(lumped(LA_m, LB_m, Mk, R_m, R_m, Rg=50, RL=RL), [1e6])[0])
    h600 = abs(transfer(lumped(LA_m, LB_m, Mk, R_m, R_m, Rg=600, RL=RL), [1e6])[0])
    krows.append(dict(k=kk, H_1MHz_Rg50=h50, dB_1MHz_Rg50=20 * np.log10(h50), H_1MHz_Rg600=h600, dB_1MHz_Rg600=20 * np.log10(h600)))
K = pd.DataFrame(krows)
K.to_csv(SIM + 'results_k_sweep.csv', index=False, float_format='%.5g')
for kk in (0.3, 0.5, 0.7, 0.9, 0.97):
    r = K[np.isclose(K.k, kk)].iloc[0]
    put(f'k_sweep_dB_1MHz_Rg50_k{str(kk).replace(".", "p")}', r.dB_1MHz_Rg50, '[MODEL] M2 (L_A, L_B, R fixed; M = k sqrt(L_A L_B)), '
        'R_g=50, R_L=111.6, |V(OUT)/V_g| at 1 MHz')
for kk, lab in ((k_m, 'kmeas'), (M_g / np.sqrt(LA_g * LB_g), 'ktheory'), (M_v2 / np.sqrt(LA_v2 * LB_v2), 'kpeec')):
    Mk = kk * np.sqrt(LA_m * LB_m)
    put(f'k_sweep_dB_1MHz_Rg50_{lab}', dB(transfer(lumped(LA_m, LB_m, Mk, R_m, R_m, Rg=50, RL=RL), [1e6]))[0],
        f'[MODEL] as above at k={kk:.4f}')
s = np.polyfit(K.k[(K.k >= 0.85)], K.dB_1MHz_Rg50[(K.k >= 0.85)], 1)[0]
put('k_sweep_slope_dB_per_0p01k_near_0p9', s * 0.01, '[MODEL] local slope of |H(1 MHz)| in dB per 0.01 change of k (k 0.85-0.97)')

# ================================================================== M0 audit waterfall
steps = []


def m0_variant(Rs=1.0, RL_=None, dots='as_drawn', caps='as_drawn', Lval=(18e-6, 18e-6), k=0.9, R=(7.31, 7.31)):
    c = Circuit('M0 variant')
    c.add_V('V1', 'src', '0', 1.0)
    c.add_R('Rs', 'src', 'AC_IN', Rs)
    sgn = 1 if dots == 'as_drawn' else -1
    # relative phasing as in the schematic: +k for (0->AC_IN, DCDOT->0), equivalent to AC dot at AC_IN with DC dot at ground
    c.add_Lmatrix(['LAC', 'LDC'], [('0', 'AC_IN'), ('DCDOT', '0')],
                  [[Lval[0], sgn * k * np.sqrt(Lval[0] * Lval[1])], [sgn * k * np.sqrt(Lval[0] * Lval[1]), Lval[1]]], rser=list(R))
    if caps == 'as_drawn':
        c.add_C('C1', 'AC_IN', 'DCDOT', 274e-12); c.add_C('C2', '0', 'DCDOT', 274e-12); c.add_C('C3', 'AC_IN', '0', 274e-12)
    if RL_:
        c.add_R('RL', 'DCDOT', '0', RL_)
    c.meas = 'DCDOT'
    return c


variants = [
    ('S0', 'M0 initial model: Rser=1, winding R 7.31/7.31 (in the inductor properties), L 18/18 uH, k 0.9 (phasing opposite to the measured polarity), C1-C3 274 pF, no load',
     dict()),
    ('S1', '+ load R_L = 111.6 ohm at the probed node', dict(RL_=RL)),
    ('S2', '+ generator source resistance 50 ohm (instead of Rser = 1)', dict(RL_=RL, Rs=50)),
    ('S3', '+ like-polarity dots as measured (AC_IN with OUT)', dict(RL_=RL, Rs=50, dots='fixed')),
    ('S4', '- remove the three 274 pF lumped capacitors (magnetic only)', dict(RL_=RL, Rs=50, dots='fixed', caps='none')),
    ('S5', '+ measured L_A, L_B, k, R (= M2)', dict(RL_=RL, Rs=50, dots='fixed', caps='none', Lval=(LA_m, LB_m), k=M_m / np.sqrt(LA_m * LB_m),
                                                       R=(R_m, R_m))),
]
for sid, desc, kw in variants:
    c = m0_variant(**kw)
    h = transfer(c, F)
    d = dB(h); im = np.argmax(d)
    eR = bench_err(lambda f, cc=c: transfer(cc, f), 'R')
    eX = bench_err(lambda f, cc=c: transfer(cc, f), 'X')
    steps.append(dict(step=sid, change=desc, dB_100kHz=dB(transfer(c, [1e5]))[0], dB_1MHz=dB(transfer(c, [1e6]))[0],
                      dB_5MHz=dB(transfer(c, [5e6]))[0], max_dB_10k_30M=d[im], f_max_MHz=F[im] / 1e6,
                      rms_all_runR=eR['all'], rms_all_runX=eX['all']))
    curves['audit_' + sid] = h
h_m3 = curves['M3']
eR, eX = bench_err(models['M3'], 'R'), bench_err(models['M3'], 'X')
steps.append(dict(step='S6', change='+ distributed 274 pF ladder, N=24 PEEC sections (= M3)', dB_100kHz=dB(transfer(cM3, [1e5]))[0],
                  dB_1MHz=dB(transfer(cM3, [1e6]))[0], dB_5MHz=dB(transfer(cM3, [5e6]))[0], max_dB_10k_30M=dB(h_m3).max(),
                  f_max_MHz=F[np.argmax(dB(h_m3))] / 1e6, rms_all_runR=eR['all'], rms_all_runX=eX['all']))
S = pd.DataFrame(steps)
S.to_csv(SIM + 'results_m0_audit_steps.csv', index=False, float_format='%.4g')
assert abs(S.iloc[0].dB_1MHz - V['sim.M0_dB_1MHz']) < 0.05, 'audit S0 must equal parsed M0'
assert np.max(np.abs(dB(curves['audit_S5']) - dB(curves['M2']))) < 1e-6, 'audit S5 must equal M2'
for _, r in S.iterrows():
    put(f'audit_{r.step}_dB_1MHz', r.dB_1MHz, f'[MODEL] audit step {r.step}: {r.change}; |H(1 MHz)| dB re 1 V source')
    put(f'audit_{r.step}_max_dB', r.max_dB_10k_30M, f'[MODEL] audit step {r.step}: maximum over 10 kHz-30 MHz')
    put(f'audit_{r.step}_f_max_MHz', r.f_max_MHz, f'[MODEL] audit step {r.step}: frequency of that maximum')

# ================================================================== curves CSV
out = pd.DataFrame({'f_Hz': F})
for nm in ('M0', 'M1_Rgeom', 'M1_R7p31', 'M1_peec', 'M2', 'M2_Rg600', 'M3', 'M3_Rg600', 'M3_N48', 'M3_uniform', 'M4_R', 'M4_X',
           'M4ext_R', 'M4ext_X'):
    out[nm + '_dB'] = dB(curves[nm])
    out[nm + '_deg'] = np.degrees(np.angle(curves[nm]))
out.to_csv(SIM + 'results_models.csv', index=False, float_format='%.5g')
pd.DataFrame({'f_Hz': f_full, 'M0_dB': dB(h0_full), 'M0_deg': np.degrees(np.angle(h0_full))}).to_csv(
    SIM + 'results_M0.csv', index=False, float_format='%.5g')
np.savez(W5 + 'curves.npz', f_lt=f_lt, h_lt=h_lt, f_full=f_full, h0_full=h0_full, **{k: v for k, v in curves.items()})

# ================================================================== netlists
HDR_BENCH = ['Bench topology: V_g (AC 1 = generator EMF) -> R_g -> AC_IN; AC_RTN = 0; DC_IN = 0 (5 V rail is an AC ground);',
             'OUT -> R_L = 111.6 ohm -> 0. Measured quantity: V(<prefix>out) / V_g.',
             'Inductor dot = FIRST node; positive K = like-polarity first nodes (measured: AC_IN with OUT, AC_RTN with DC_IN).',
             'Winding resistance is an explicit resistor (R*_ser). LTspice adds its default 1 mohm inductor Rser; negligible here.',
             'Generated by analysis/circuit_models.py with analysis/mna_sim.py. Label: [MODEL].']


def merged(circs, tags):
    m = Circuit(' + '.join(tags))
    for c, t in zip(circs, tags):
        rn = lambda n, t=t: n if n == '0' else f'{t}_{n}'  # noqa: E731
        for nm, a, b, v in c.V: m.add_V(f'{nm}_{t}', rn(a), rn(b), v)
        for nm, a, b, v in c.R: m.add_R(f'{nm}_{t}', rn(a), rn(b), v)
        for nm, a, b, v in c.C: m.add_C(f'{nm}_{t}', rn(a), rn(b), v)
        off = len(m.L)
        for nm, a, b, v in c.L: m.L.append((f'{nm}_{t}', rn(a), rn(b), v)); m.n(rn(a)); m.n(rn(b))
        for (i, j), k in c.Kc.items(): m.Kc[(i + off, j + off)] = k
    return m


def write_net(path, circ, title_lines, probes, check):
    meas_l = [f'.meas AC {p.upper()}_1MHz FIND V({p}) AT 1Meg' for p in probes]
    txt = circ.to_spice(header_lines=title_lines, tail_lines=['.ac dec 250 10k 30Meg'] + meas_l + ['.backanno', '.end'])
    open(path, 'w').write(txt)
    c_back, _, _ = parse_spice(txt)
    worst = 0.0
    for p, ref in zip(probes, check):
        worst = max(worst, float(np.max(np.abs(dB(c_back.v(c_back.solve(F), p)) - dB(ref)))))
    return worst


w1 = write_net(SIM + 'm1_geometry.cir', merged([cM1a, cM1b], ['m1a', 'm1b']),
               [f'M1 geometry-only prediction (magnetic only). m1a: R_AC = {RA_g:g}, R_DC = {RB_g:g} ohm (values_theory.json, drawn 0.254 mm x 35 um);',
                f'm1b: R = {R_m:g} ohm (measured). L_AC {LA_g*1e6:g}, L_DC {LB_g*1e6:g}, M {M_g*1e6:g} uH (Neumann, top+via+bottom route), R_g = 50 ohm.'] + HDR_BENCH,
               ['m1a_out', 'm1b_out'], [curves['M1_Rgeom'], curves['M1_R7p31']])
w2 = write_net(SIM + 'corrected_injector.cir', merged([cM2, cM3], ['m2', 'm3']),
               ['Corrected injector model (supersedes the 17 Nov schematic). Two independent circuits in one file:',
                f'm2 = M2 lumped, measured L_A {LA_m*1e6:g}, L_B {LB_m*1e6:g}, M {M_m*1e6:g} uH (k {M_m/np.sqrt(LA_m*LB_m):.4f}), '
                f'R {R_m:g} ohm each, R_g 50 ohm, R_L {RL:g} ohm (data/values_meas.json).',
                f'm3 = M3: M2 + {C_m*1e12:g} pF inter-winding C distributed over N=24 sections per winding. Section inductance matrix =',
                f"PEEC partial inductances (data/model_inputs/peec_sections_N24.npz) scaled to the measured totals (AC x{sc['alpha']:.4f}, "
                f"DC x{sc['beta']:.4f}, cross x{sc['gamma']:.4f}).",
                'Section i runs outer->inner. AC nodes: 0 (AC_RTN), m3_a1..m3_a23, m3_ac_in (AC_IN); DC nodes: 0 (DC_IN), m3_b1..m3_b23, m3_out (OUT).',
                'Change Rg_m2 / Rg_m3 to 600 for the 600-ohm generator setting.'] + HDR_BENCH,
               ['m2_out', 'm3_out'], [curves['M2'], curves['M3']])
w3 = write_net(SIM + 'm4_bench_hypotheses.cir', merged([cM4['R'], cM4['X']], ['hr', 'hx']),
               ['M4 HYPOTHESES (not conclusions): parasitic sets that best reproduce the normalised bench curves.',
                f"hr (Run R): {V['sim.M4_best_runR_params']}", f"hx (Run X): {V['sim.M4_best_runX_params']}",
                'Probe x1 = 1 Mohm || (15 pF + cable C); x10 = 10 Mohm || 12 pF; Lw = probe loop inductance between OUT and the probe;',
                f"LRL = series inductance of the load resistor. Measured nodes: V(hr_{cM4['R'].meas}) and V(hx_{cM4['X'].meas})."] + HDR_BENCH,
               [f"hr_{cM4['R'].meas}", f"hx_{cM4['X'].meas}"], [curves['M4_R'], curves['M4_X']])
put('netlist_roundtrip_maxdiff_dB', max(w1, w2, w3), '[MODEL] max |dB| difference between the written .cir files re-parsed and solved '
    'and the in-memory circuits (10 kHz-30 MHz); checks the LTspice netlists encode the same circuits')
put('netlist_n_K_cards_corrected_injector', len(cM3.Kc) + len(cM2.Kc), '[MODEL] number of K statements in corrected_injector.cir')
put('netlist_n_C_M3', len(cM3.C), '[MODEL] number of inter-winding capacitors in the M3 ladder (fractions < 1e-5 of 274 pF dropped)')

# figure normalisation check: models are normalised at 1.000 MHz in the figure, bench at 1003.2 / 1030 kHz
dn = max(abs(dB(transfer_f([runs[r]['f_ref']])[0]) - dB(transfer_f([1e6])[0]))
         for transfer_f in [models[m] for m in ('M1_Rgeom', 'M2', 'M3', 'M4_R', 'M4_X')] for r in runs)
put('fig_norm_offset_max_dB', dn, '[MODEL] largest change of any plotted model if normalised at the bench reference '
    'frequency (1003.2 or 1030 kHz) instead of 1.000 MHz')

_d24 = pd.read_csv(ROOT + 'data/bench_readings_A_oct24.csv')
_ring = list(_d24[_d24.gen_freq_setting_Hz == 100.1].freq_from_period_Hz.dropna())
_ring.append(json.load(open(ROOT + 'data/bench_readings_A_oct22.json'))['X01']['freq_from_screen_hz'])
put('bench_ringing_f_range_Hz', [float(f'{min(_ring):.4g}'), float(f'{max(_ring):.4g}')], '[MEASURED] ringing-burst frequencies '
    'from screen periods: X01 (bench_readings_A_oct22.json) and P24-43/45/48 (bench_readings_A_oct24.csv); used as the '
    'target band 7.0-8.7 MHz (with margin) for the M4 ringing test')

# ================================================================== JSON
json.dump(V, open(ROOT + 'data/values_sim.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in V.items() if not k.endswith('_meta')}, indent=0)[:6000])
print(S.to_string())
