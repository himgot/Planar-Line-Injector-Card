"""bench data for the simulation comparison.
Source: data/bench_reconciled.csv, written by
analysis/reconcile_bench.py (consensus median of three independent photo readers).
  Run R = 22 Oct raw photos (10 V/div), normalised to its 1003.2 kHz point.
  Run X = workbook images (10 V/div below 1 MHz, 20 V/div above), normalised to its 1030 kHz point.
Columns used: f_kHz (screen/generator frequency as reconciled), rel_dB_re_1MHz, rel_dB_halfspread,
Vpp_displayed, Vpp_interpreted_div10. The two runs are never merged."""
import numpy as np
import pandas as pd

SRC = 'data/bench_reconciled.csv'


def load():
    d = pd.read_csv(SRC)
    runs = {}
    for run, g in d.groupby('run'):
        g = g.sort_values('f_kHz')
        pts = []
        for _, r in g.iterrows():
            pts.append(dict(id=str(r.photos).split()[0], photos=r.photos, f=float(r.f_kHz) * 1e3,
                            dB_norm=float(r.rel_dB_re_1MHz),
                            dB_half=float(r.rel_dB_halfspread) if np.isfinite(r.rel_dB_halfspread) else np.nan,
                            Vpp=float(r.Vpp_displayed), Vpp_div10=float(r.Vpp_interpreted_div10),
                            n_readers=int(r.n_readers), note='' if pd.isna(r.note) else str(r.note)))
        ref = [p for p in pts if abs(p['dB_norm']) < 1e-9]
        assert len(ref) == 1, run
        runs[run] = dict(points=pts, f_ref=ref[0]['f'], ref=ref[0]['id'], source=SRC)
    return runs


if __name__ == '__main__':
    for k, v in load().items():
        print(k, v['ref'], v['f_ref'])
        for p in v['points']:
            print('   %-6s %9.0f  %6.2f V  %+7.2f dB (+-%.2f)' % (p['id'], p['f'], p['Vpp'], p['dB_norm'], p['dB_half']))
