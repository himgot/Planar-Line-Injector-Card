#!/usr/bin/env python3
"""mna_sim.py -- small-signal AC modified-nodal-analysis (MNA) solver.

No SPICE engine is used; this file provides:
  * Circuit: R, C, L (with an arbitrary inductance matrix, i.e. coupled inductors), independent
    voltage sources (AC magnitude), descriptor form (G + s*Cm) x = b, AC sweep, natural poles.
  * parse_spice(): reader for the simple SPICE/LTspice netlists used in this project
    (R, C, L, K, V cards, .param with {expressions}, '+' continuation, '*' and ';' comments).
  * Circuit.to_spice(): LTspice-compatible netlist writer (explicit resistors, K cards).
  * Model builders for the bench topology (lumped M1/M2, sectioned M3, M4 parasitics).

Formulation (unknown vector x = [node voltages | inductor branch currents | V-source currents]):
  resistor  : conductance stamp in G
  capacitor : capacitance stamp in Cm (multiplied by s = j*2*pi*f)
  inductor  : KCL +-1 in G; KVL row  v_a - v_b - s * sum_j Lmat[i,j] * i_j = 0 (Lmat in Cm)
  V source  : KCL +-1 in G; constraint row v_a - v_b = V_ac
Natural frequencies are the finite generalised eigenvalues s of G x = -s Cm x.
Labels: everything computed here is [MODEL].
"""
import re
import numpy as np

_SUF = [('meg', 1e6), ('mil', 25.4e-6), ('t', 1e12), ('g', 1e9), ('k', 1e3), ('m', 1e-3),
        ('u', 1e-6), ('µ', 1e-6), ('n', 1e-9), ('p', 1e-12), ('f', 1e-15)]


def spice_val(tok, params=None):
    """Convert a SPICE value token ('18u', '274p', '1Meg', '1e-6', '{RAC_HID}') to float."""
    tok = tok.strip()
    if tok.startswith('{') and tok.endswith('}'):
        expr = tok[1:-1]
        return float(eval(expr, {'__builtins__': {}}, dict(params or {})))  # own netlists only
    if params and tok in params:
        return float(params[tok])
    m = re.match(r'^([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)([a-zA-Zµ]*)$', tok)
    if not m:
        raise ValueError(f'cannot parse value {tok!r}')
    x = float(m.group(1))
    s = m.group(2).lower()
    for suf, mul in _SUF:
        if s.startswith(suf):
            return x * mul
    return x


def fmt(x):
    """compact number for netlists"""
    return f'{x:.6g}'


class Circuit:
    def __init__(self, title=''):
        self.title = title
        self.nodes = {}            # name -> index (ground '0' not stored)
        self.R, self.C, self.V = [], [], []
        self.L = []                # (name, a, b, Lself)
        self.Kc = {}               # (i, j) i<j -> coupling coefficient k_ij
        self.notes = []

    # ---------------------------------------------------------------- building
    def n(self, name):
        name = str(name)
        if name.lower() in ('0', 'gnd'):
            return -1
        if name not in self.nodes:
            self.nodes[name] = len(self.nodes)
        return self.nodes[name]

    def add_R(self, name, a, b, val):
        if val == 0:
            raise ValueError('zero resistor; short the nodes instead')
        self.n(a); self.n(b); self.R.append((name, str(a), str(b), float(val)))

    def add_C(self, name, a, b, val):
        self.n(a); self.n(b); self.C.append((name, str(a), str(b), float(val)))

    def add_L(self, name, a, b, val, rser=0.0):
        """inductor from a to b; the dot (like-polarity terminal) is the FIRST node a.
        rser > 0 inserts an explicit series resistor on the a side (internal node name+'_r')."""
        if rser and rser > 0:
            mid = name + '_r'
            self.add_R('R' + name[1:] + '_ser' if name[0] in 'Ll' else 'R' + name + '_ser', a, mid, rser)
            a = mid
        self.n(a); self.n(b)
        self.L.append((name, str(a), str(b), float(val)))
        return len(self.L) - 1

    def add_K(self, l1, l2, k):
        names = [x[0] for x in self.L]
        i, j = names.index(l1), names.index(l2)
        if i == j:
            raise ValueError('K on same inductor')
        self.Kc[(min(i, j), max(i, j))] = float(k)

    def add_Lmatrix(self, names, ends, Lmat, rser=None):
        """add len(names) inductors with full inductance matrix Lmat (henry).
        ends[i] = (a_i, b_i) with a_i the dotted (reference) end."""
        Lmat = np.asarray(Lmat, float)
        idx = []
        for i, nm in enumerate(names):
            idx.append(self.add_L(nm, ends[i][0], ends[i][1], Lmat[i, i], 0.0 if rser is None else rser[i]))
        d = np.sqrt(np.diag(Lmat))
        for p in range(len(idx)):
            for q in range(p + 1, len(idx)):
                if Lmat[p, q] != 0.0:
                    self.Kc[(idx[p], idx[q])] = Lmat[p, q] / (d[p] * d[q])

    def add_V(self, name, a, b, ac=1.0):
        self.n(a); self.n(b); self.V.append((name, str(a), str(b), complex(ac)))

    # ---------------------------------------------------------------- matrices
    def Lmatrix(self):
        nl = len(self.L)
        Ls = np.array([x[3] for x in self.L])
        Lm = np.diag(Ls)
        for (i, j), k in self.Kc.items():
            Lm[i, j] = Lm[j, i] = k * np.sqrt(Ls[i] * Ls[j])
        return Lm if nl else np.zeros((0, 0))

    def build(self):
        nn, nl, nv = len(self.nodes), len(self.L), len(self.V)
        N = nn + nl + nv
        G = np.zeros((N, N)); Cm = np.zeros((N, N)); b = np.zeros(N, complex)
        ix = lambda nm: self.n(nm)  # noqa: E731

        def stamp(M, a, c, y):
            if a >= 0: M[a, a] += y
            if c >= 0: M[c, c] += y
            if a >= 0 and c >= 0:
                M[a, c] -= y; M[c, a] -= y
        for _, a, c, r in self.R:
            stamp(G, ix(a), ix(c), 1.0 / r)
        for _, a, c, cv in self.C:
            stamp(Cm, ix(a), ix(c), cv)
        Lm = self.Lmatrix()
        for k, (_, a, c, _) in enumerate(self.L):
            r = nn + k
            ia, ic = ix(a), ix(c)
            if ia >= 0: G[ia, r] += 1; G[r, ia] += 1
            if ic >= 0: G[ic, r] -= 1; G[r, ic] -= 1
        if nl:
            Cm[nn:nn + nl, nn:nn + nl] -= Lm
        for k, (_, a, c, v) in enumerate(self.V):
            r = nn + nl + k
            ia, ic = ix(a), ix(c)
            if ia >= 0: G[ia, r] += 1; G[r, ia] += 1
            if ic >= 0: G[ic, r] -= 1; G[r, ic] -= 1
            b[r] = v
        return G, Cm, b

    def solve(self, freqs, chunk=64):
        """AC sweep; returns complex solution array (nf, N)"""
        G, Cm, b = self.build()
        f = np.atleast_1d(np.asarray(freqs, float))
        out = np.empty((len(f), len(b)), complex)
        for c0 in range(0, len(f), chunk):
            s = 2j * np.pi * f[c0:c0 + chunk]
            A = G[None, :, :] + s[:, None, None] * Cm[None, :, :]
            out[c0:c0 + chunk] = np.linalg.solve(A, np.broadcast_to(b, (len(s), len(b)))[..., None])[..., 0]
        return out

    def v(self, X, node):
        i = self.n(node)
        return np.zeros(X.shape[0], complex) if i < 0 else X[:, i]

    def poles(self, fmin=1e3, fmax=1e9):
        """finite natural frequencies (complex s, rad/s) with fmin <= |s|/2pi <= fmax"""
        import scipy.linalg as sla
        G, Cm, _ = self.build()
        w = sla.eigvals(G, -Cm)
        w = w[np.isfinite(w)]
        a = np.abs(w) / (2 * np.pi)
        return w[(a >= fmin) & (a <= fmax)]

    # ---------------------------------------------------------------- export
    def to_spice(self, header_lines=(), tail_lines=()):
        L = []
        L.append('* ' + self.title)
        for h in header_lines:
            L.append('* ' + h if not h.startswith(('*', '.')) else h)
        for nm, a, c, v in self.V:
            L.append(f'{nm} {a} {c} AC {fmt(abs(v))}')
        for nm, a, c, r in self.R:
            L.append(f'{nm} {a} {c} {fmt(r)}')
        for nm, a, c, cv in self.C:
            L.append(f'{nm} {a} {c} {fmt(cv)}')
        for nm, a, c, lv in self.L:
            L.append(f'{nm} {a} {c} {fmt(lv)}')
        nk = 0
        for (i, j), k in sorted(self.Kc.items()):
            nk += 1
            L.append(f'K{nk} {self.L[i][0]} {self.L[j][0]} {k:.9f}')
        L += list(tail_lines)
        return '\n'.join(L) + '\n'


def parse_spice(text, params=None, title=''):
    """Parse a simple SPICE / LTspice netlist into a Circuit. Returns (circuit, dot_commands, params)."""
    prm = {}
    lines = []
    for raw in text.splitlines():
        s = raw.split(';')[0].rstrip()
        if not s.strip() or s.lstrip().startswith('*'):
            continue
        if s.lstrip().startswith('+') and lines:
            lines[-1] += ' ' + s.lstrip()[1:]
        else:
            lines.append(s.strip())
    dots = []
    for s in lines:
        if s.lower().startswith('.param'):
            for m in re.finditer(r'(\w+)\s*=\s*(\{[^}]*\}|[^\s]+)', s[6:]):
                prm[m.group(1)] = spice_val(m.group(2), prm)
    if params:
        prm.update(params)
    ckt = Circuit(title)
    kcards = []
    for s in lines:
        tok = s.split()
        c0 = tok[0][0].upper()
        if c0 == '.':
            dots.append(s)
            continue
        if c0 == 'R':
            ckt.add_R(tok[0], tok[1], tok[2], spice_val(tok[3], prm))
        elif c0 == 'C':
            ckt.add_C(tok[0], tok[1], tok[2], spice_val(tok[3], prm))
        elif c0 == 'L':
            rser = 0.0
            for t in tok[4:]:
                if t.lower().startswith('rser='):
                    rser = spice_val(t.split('=', 1)[1], prm)
            ckt.add_L(tok[0], tok[1], tok[2], spice_val(tok[3], prm), rser)
        elif c0 == 'K':
            kcards.append(tok)
        elif c0 == 'V':
            ac = 0.0
            rest = re.sub(r'\w+\([^)]*\)', '', ' '.join(tok[3:])).split()   # drop SINE(...) etc.
            rser = 0.0
            for i, t in enumerate(rest):
                if t.upper() == 'AC':
                    ac = spice_val(rest[i + 1], prm)
                if t.lower().startswith('rser='):
                    rser = spice_val(t.split('=', 1)[1], prm)
            if rser:
                mid = tok[0] + '_int'
                ckt.add_V(tok[0], mid, tok[2], ac)
                ckt.add_R('R' + tok[0] + '_ser', mid, tok[1], rser)
            else:
                ckt.add_V(tok[0], tok[1], tok[2], ac)
        else:
            raise ValueError(f'unsupported card: {s}')
    for tok in kcards:
        k = spice_val(tok[-1], prm)
        names = tok[1:-1]
        for p in range(len(names)):
            for q in range(p + 1, len(names)):
                ckt.add_K(names[p], names[q], k)
    return ckt, dots, prm


# ====================================================================== bench-topology builders
def probe_params(kind, c_cable=100e-12):
    """[MODEL] probe loading at the measurement node.
    x10 : 10 Mohm || 12 pF (10-15 pF range given in the task brief)
    x1  : 1 Mohm || (15 pF scope input + cable capacitance)."""
    if kind in (None, 'none'):
        return None
    if kind == 'x10':
        return dict(R=10e6, C=12e-12)
    if kind == 'x1':
        return dict(R=1e6, C=15e-12 + c_cable)
    raise ValueError(kind)


def bench_lumped(LA, LB, M, RA, RB, Rg=50.0, RL=111.6, probe=None, Lw=0.0, L_RL=0.0, C_direct=0.0,
                 c_cable=100e-12, title='bench lumped'):
    """Bench topology with a lumped coupled-inductor pair.
    V_g (AC 1) --Rg-- ac_in ; AC winding ac_in->ac_rtn(=0) ; DC winding out->dc_in(=0, 5 V rail = AC ground);
    load RL (+ optional series L_RL) out->0 ; probe at node 'meas' (through optional loop inductance Lw).
    Dots: ac_in and out (like-polarity, from the measured polarity). Output quantity: V(meas)/V_g."""
    c = Circuit(title)
    c.add_V('Vg', 'src', '0', 1.0)
    c.add_R('Rg', 'src', 'ac_in', Rg)
    c.add_Lmatrix(['LAC', 'LDC'], [('ac_in', '0'), ('out', '0')], [[LA, M], [M, LB]], rser=[RA, RB])
    _load_probe(c, RL, L_RL, probe, Lw, c_cable)
    if C_direct:
        c.add_C('Cdir', 'ac_in', 'out', C_direct)
    return c


def _load_probe(c, RL, L_RL, probe, Lw, c_cable):
    if L_RL > 0:
        c.add_R('RL', 'out', 'rl_l', RL)
        c.add_L('LRL', 'rl_l', '0', L_RL)
    else:
        c.add_R('RL', 'out', '0', RL)
    pp = probe_params(probe, c_cable)
    if pp is None:
        c.meas = 'out'
        return
    node = 'out'
    if Lw > 0:
        c.add_L('Lw', 'out', 'meas', Lw)
        node = 'meas'
    c.add_R('Rprobe', node, '0', pp['R'])
    c.add_C('Cprobe', node, '0', pp['C'])
    c.meas = node


def scaled_section_matrix(Lsec, LA, LB, M):
    """scale PEEC section blocks so that block sums equal (LA, LB, M):
    AC block x alpha, DC block x beta, cross blocks x gamma. PSD is preserved when
    gamma <= sqrt(alpha*beta) (proof: convex combination with the block-diagonal part)."""
    N = Lsec.shape[0] // 2
    A, D, X = Lsec[:N, :N], Lsec[N:, N:], Lsec[:N, N:]
    al, be, ga = LA / A.sum(), LB / D.sum(), M / X.sum()
    out = np.block([[al * A, ga * X], [ga * X.T, be * D]])
    return out, dict(alpha=al, beta=be, gamma=ga, rho=ga / np.sqrt(al * be))


def uniform_section_matrix(N, LA, LB, M, cw=0.5):
    """fallback: uniform-coupling sections that preserve the totals LA, LB, M.
    within-winding blocks a[(1-cw) I + cw 11^T], cross block (M/N^2) 11^T."""
    J = np.ones((N, N)); I = np.eye(N)
    base = (1 - cw) * I + cw * J
    sA = LA / base.sum(); sB = LB / base.sum()
    return np.block([[sA * base, M / N**2 * J], [M / N**2 * J, sB * base]])


def bench_sectioned(Lsec_scaled, Cn, Ctot, RA, RB, Rg=50.0, RL=111.6, probe=None, Lw=0.0, L_RL=0.0,
                    c_cable=100e-12, cmin_frac=1e-5, title='bench sectioned'):
    """M3: N sections per winding, full section inductance matrix, distributed inter-winding C.
    AC nodes a0..aN : a0 = AC_RTN (ground), aN = 'ac_in'.  DC nodes b0..bN : b0 = DC_IN (AC ground), bN = 'out'.
    Section i of each winding runs outer->inner (a_i -> a_{i+1}); the outer end of every section is its
    dotted end, so M>0 (PEEC sign: both windings outer->inner) = like-polarity outer ends (AC_RTN, DC_IN).
    Cn[i, j] = fraction of Ctot between AC node i and DC node j (from geometry adjacency)."""
    N = Lsec_scaled.shape[0] // 2
    an = ['0'] + [f'a{i}' for i in range(1, N)] + ['ac_in']
    bn = ['0'] + [f'b{i}' for i in range(1, N)] + ['out']
    c = Circuit(title)
    c.add_V('Vg', 'src', '0', 1.0)
    c.add_R('Rg', 'src', 'ac_in', Rg)
    names = [f'LA{i}' for i in range(N)] + [f'LB{i}' for i in range(N)]
    ends = [(an[i], an[i + 1]) for i in range(N)] + [(bn[i], bn[i + 1]) for i in range(N)]
    c.add_Lmatrix(names, ends, Lsec_scaled, rser=[RA / N] * N + [RB / N] * N)
    nc = 0
    dropped = 0.0
    for i in range(N + 1):
        for j in range(N + 1):
            fr = Cn[i, j]
            if fr <= 0:
                continue
            if an[i] == '0' and bn[j] == '0':
                continue                      # both ends at AC ground
            if fr < cmin_frac:
                dropped += fr
                continue
            nc += 1
            c.add_C(f'Cw{nc}', an[i], bn[j], fr * Ctot)
    c.notes.append(f'{nc} inter-winding capacitors; dropped fraction {dropped:.2e}')
    _load_probe(c, RL, L_RL, probe, Lw, c_cable)
    return c


def transfer(c, freqs, node=None):
    X = c.solve(freqs)
    return c.v(X, node or c.meas)


def ring_poles(c, fmin=1e6, fmax=50e6):
    """complex natural frequencies with positive imaginary part: list of (f_d [Hz], Q)"""
    p = c.poles(fmin, fmax)
    p = p[p.imag > 0]
    out = []
    for s in p:
        wn = abs(s)
        zeta = -s.real / wn
        Q = 1 / (2 * zeta) if zeta > 0 else np.inf
        out.append((s.imag / (2 * np.pi), Q, wn / (2 * np.pi)))
    return sorted(out)


if __name__ == '__main__':
    # self-test: ideal transformer limit and RC low-pass
    c = Circuit('rc')
    c.add_V('V1', 'in', '0', 1); c.add_R('R1', 'in', 'o', 1e3); c.add_C('C1', 'o', '0', 1e-9)
    f = np.array([1e3, 159154.94, 1e7])
    h = c.v(c.solve(f), 'o')
    assert abs(abs(h[1]) - 2 ** -0.5) < 1e-6, h
    p = c.poles(1, 1e9)
    assert abs(abs(p[0]) - 1e6) < 1e-3, p
    c2 = bench_lumped(18e-6, 18e-6, 0.9999 * 18e-6, 1e-6, 1e-6, Rg=1e-3, RL=1e6)
    h2 = transfer(c2, [1e6])
    assert abs(h2[0] - 0.9999) < 1e-3, h2
    txt = c2.to_spice(tail_lines=['.ac dec 10 1k 10Meg', '.end'])
    c3, _, _ = parse_spice(txt)
    c3.meas = 'out'
    assert abs(transfer(c3, [1e6])[0] - h2[0]) < 1e-9
    print('mna_sim self-test OK')
