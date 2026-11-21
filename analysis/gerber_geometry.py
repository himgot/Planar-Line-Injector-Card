#!/usr/bin/env python3
"""Geometry extraction from the line-injector-card Gerber set.

Run with:  python3 -I analysis/gerber_geometry.py
(the Gerber zip is opened read-only with zipfile; nothing in it is executed.)

Own minimal RS-274X parser written for this file set:
  * %FSTAX23Y23*%  -> trailing zeros omitted, 2 integer + 3 decimal digits;
    a coordinate string is padded on the RIGHT to 5 digits, value = int/1000 in.
  * %MOIN*%        -> inches.
  * modal X/Y, modal D01/D02/D03, aperture select G54Dnn, G36/G37 regions,
    %LPD*% / %LPC*% polarity layers (legend files).  Only linear interpolation
    occurs in these files (no G02/G03 -- asserted).

Outputs (all written by this script):
  work/theory/geometry.json        (machine-readable geometry)
  data/geometry_table.csv      (human table)
  work/theory/*.png                (legend / copper renders, OCR crops)
  keys "theory.geom.*" merged into data/values_theory.json
"""
import csv
import json
import math
import os
import re
import subprocess
import zipfile
from collections import defaultdict

import numpy as np

ZIP = "hardware/Planar_Line_Injector_gerbers.zip"
WORK = "work/theory"
DATA = "data"
VALUES = os.path.join(DATA, "values_theory.json")
IN2MM = 25.4
SRC = "gerber.zip (New folder/)"

os.makedirs(WORK, exist_ok=True)


# --------------------------------------------------------------------------
# helpers shared by the other scripts
# --------------------------------------------------------------------------
def update_values(entries):
    """Merge {key: (value, meta_str)} into values_theory.json (flat keys + _meta)."""
    if os.path.exists(VALUES):
        with open(VALUES) as f:
            d = json.load(f)
    else:
        d = {"_meta": {}}
    meta = d.pop("_meta", {})
    for k, (v, m) in entries.items():
        if isinstance(v, (np.floating, np.integer)):
            v = v.item()
        d[k] = v
        meta[k] = m
    out = {k: d[k] for k in sorted(d)}
    out["_meta"] = {k: meta[k] for k in sorted(meta)}
    with open(VALUES, "w") as f:
        json.dump(out, f, indent=1)


def load_geometry():
    with open(os.path.join(WORK, "geometry.json")) as f:
        return json.load(f)


def sig(x, n=4):
    """round to n significant figures (for JSON readability)."""
    x = float(x)
    if x == 0 or not np.isfinite(x):
        return x
    return float(round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1)))


# --------------------------------------------------------------------------
# Gerber parser
# --------------------------------------------------------------------------
def _coord(s):
    """FSTAX23: trailing zeros omitted -> pad right to 5 digits; inches."""
    sign = -1 if s.startswith("-") else 1
    s = s.lstrip("+-")
    if len(s) < 5:
        s = (s + "00000")[:5]
    return sign * int(s) / 1000.0


def parse_gerber(text):
    """Return dict with apertures and a list of polarity layers.

    each layer: {'name','polarity','draws':[(ap,x1,y1,x2,y2)],
                 'flashes':[(ap,x,y)],'regions':[[(x,y),...],...]}
    """
    assert "%FSTAX23Y23*%" in text and "%MOIN*%" in text, "unexpected format"
    assert "G02" not in text and "G03" not in text and "G75" not in text
    apertures = {}
    layers = []
    cur = {"name": None, "polarity": "D", "draws": [], "flashes": [], "regions": []}
    layers.append(cur)
    x = y = 0.0
    dmode = None
    ap = None
    in_region = False
    contour = []
    tokens = re.findall(r"%[^%]*%|[^*%\r\n]+\*", text)
    for tok in tokens:
        tok = tok.strip()
        if tok.startswith("%"):
            body = tok.strip("%").rstrip("*")
            m = re.match(r"ADD(\d+)C,([0-9.]+)", body)
            if m:
                apertures[int(m.group(1))] = ("C", float(m.group(2)))
                continue
            assert not body.startswith("AD"), "non-circular aperture: " + body
            if body.startswith("LN"):
                if cur["draws"] or cur["flashes"] or cur["regions"] or cur["name"]:
                    cur = {"name": body[2:], "polarity": "D", "draws": [],
                           "flashes": [], "regions": []}
                    layers.append(cur)
                else:
                    cur["name"] = body[2:]
                continue
            if body.startswith("LP"):
                cur["polarity"] = body[2]
                continue
            continue  # FS, MO, SF, IP: format asserted above
        tok = tok.rstrip("*")
        if tok == "M02":
            break
        if tok == "G36":
            in_region, contour = True, []
            continue
        if tok == "G37":
            if len(contour) > 2:
                cur["regions"].append(contour)
            in_region, contour = False, []
            continue
        m = re.match(r"G54D(\d+)$", tok)
        if m:
            ap = int(m.group(1))
            continue
        m = re.match(r"^D(\d+)$", tok)
        if m and int(m.group(1)) >= 10:
            ap = int(m.group(1))
            continue
        m = re.match(r"^(?:G01)?(?:X([+-]?\d+))?(?:Y([+-]?\d+))?(?:D0?([123]))?$", tok)
        if not m:
            raise ValueError("unparsed token %r" % tok)
        nx = _coord(m.group(1)) if m.group(1) is not None else x
        ny = _coord(m.group(2)) if m.group(2) is not None else y
        if m.group(3):
            dmode = int(m.group(3))
        if dmode == 1:
            if in_region:
                if not contour:
                    contour.append((x, y))
                contour.append((nx, ny))
            else:
                cur["draws"].append((ap, x, y, nx, ny))
        elif dmode == 2:
            if in_region:
                if len(contour) > 2:
                    cur["regions"].append(contour)
                contour = [(nx, ny)]
        elif dmode == 3:
            cur["flashes"].append((ap, nx, ny))
        x, y = nx, ny
    return {"apertures": apertures, "layers": layers}


# Aperture codes used by the analysis, identified by shape and size (inch). The released files number their
# apertures differently, so they are renumbered by size before parsing. Rectangular apertures belong only to the
# unpopulated "Noise injection" footprint and are set aside (see drop_footprint).
APERTURE_BY_SIZE = {0.062992: 11, 0.010000: 12, 0.023622: 14, 0.039370: 15, 0.001000: 16,
                    0.011811: 9500, 0.043307: 9501}
LAYER_NAMES = {"copper_signal_top": "card_copper_signal_top.gbr", "copper_signal_bot": "card_copper_signal_bot.gbr",
               "legend_top": "card_legend_top.gbr", "legend_bot": "card_legend_bot.gbr",
               "soldermask_top": "card_soldermask_top.gbr", "soldermask_bot": "card_soldermask_bot.gbr",
               "profile": "card_profile.gbr"}


FOOTPRINT_PADS = {}   # layer name -> [(x, y)] of the set-aside rectangular footprint pads


def _renumber(text, legend=False):
    """Renumber circular apertures by size. Rectangular apertures (footprint pads) are set aside: their flashes
    become moves (D02), so the modal current point is unchanged, and their positions are returned."""
    if legend:
        return text, []
    mapping, dropped = {}, set()
    for num, shape, size in re.findall(r"%ADD(\d+)([CR]),([0-9.X]+)\*%", text):
        if shape == "R":
            dropped.add(num)
        else:
            mapping[num] = str(APERTURE_BY_SIZE[round(float(size), 6)])
    out, cur, pads = [], None, []
    x = y = 0.0
    for tok in re.findall(r"%[^%]*%|[^*%\r\n]+\*|\r?\n", text):
        m = re.match(r"%ADD(\d+)([CR]),", tok)
        if m:
            if m.group(1) in dropped:
                continue
            tok = tok.replace("ADD%s%s," % (m.group(1), m.group(2)), "ADD%s%s," % (mapping[m.group(1)], m.group(2)), 1)
            out.append(tok)
            continue
        if tok.startswith("%"):
            out.append(tok)
            continue
        m = re.match(r"^(G54)?D(\d+)\*$", tok)
        if m and int(m.group(2)) >= 10:
            cur = m.group(2)
            if cur not in dropped:
                out.append("%sD%s*" % (m.group(1) or "", mapping[cur]))
            continue
        mx = re.search(r"X([+-]?\d+)", tok)
        my = re.search(r"Y([+-]?\d+)", tok)
        if mx:
            x = _coord(mx.group(1))
        if my:
            y = _coord(my.group(1))
        if cur in dropped and re.search(r"D0?3\*$", tok):
            pads.append((x, y))
            tok = re.sub(r"D0?3\*$", "D02*", tok)
        out.append(tok)
    return "".join(out), pads


def read_zip():
    z = zipfile.ZipFile(ZIP)
    out = {}
    for n in z.namelist():
        if n.endswith("/"):
            continue
        base = n.split("/")[-1]
        text = z.read(n).decode("ascii")
        if base.endswith(".txt"):          # NC drill (Gerber-style)
            out["card.txt"] = _renumber(text)[0]
            continue
        for key, std in LAYER_NAMES.items():
            if base.endswith(key + ".gbr"):
                out[std], pads = _renumber(text, legend=key.startswith("legend"))
                if pads:
                    FOOTPRINT_PADS[std] = pads
    return out


def drop_footprint(layer, pad_xy, tol=1e-6):
    """Remove the copper draws of the unpopulated footprint from the top layer.

    Starting at the footprint pads, draws are followed through shared end points until a round edge pad
    (pad_xy) is reached. Those draws connect the footprint to AC_IN / AC_RTN and are not part of a winding.
    Returns the number of draws removed."""
    fp = layer.get("footprint_pads", [])
    if not fp:
        return 0
    key = lambda x, y: (round(x, 4), round(y, 4))
    pads = {key(x, y) for x, y in pad_xy}
    adj = {}
    for i, d in enumerate(layer["draws"]):
        for a in (key(d[1], d[2]), key(d[3], d[4])):
            adj.setdefault(a, []).append(i)
    seen, stack, drop = set(), [key(x, y) for x, y in fp], set()
    while stack:
        p = stack.pop()
        if p in seen:
            continue
        seen.add(p)
        if p in pads:
            continue
        for i in adj.get(p, []):
            if i in drop:
                continue
            drop.add(i)
            d = layer["draws"][i]
            for q in (key(d[1], d[2]), key(d[3], d[4])):
                stack.append(q)
    layer["draws"] = [d for i, d in enumerate(layer["draws"]) if i not in drop]
    return len(drop)


# --------------------------------------------------------------------------
# rasteriser (PIL)
# --------------------------------------------------------------------------
BOARD_BBOX = (0.0, 0.0, 3.370, 2.125)


def render(parsed, path, dpi=1200, mirror=False, bbox=BOARD_BBOX):
    from PIL import Image, ImageDraw
    x0, y0, x1, y1 = bbox
    W = int(round((x1 - x0) * dpi)) + 1
    H = int(round((y1 - y0) * dpi)) + 1
    img = Image.new("L", (W, H), 255)
    dr = ImageDraw.Draw(img)

    def P(x, y):
        px = (x - x0) * dpi
        if mirror:
            px = W - 1 - px
        return (px, (y1 - y) * dpi)

    aps = parsed["apertures"]
    for layer in parsed["layers"]:
        col = 0 if layer["polarity"] == "D" else 255
        for poly in layer["regions"]:
            dr.polygon([P(*p) for p in poly], fill=col)
        for ap, xa, ya, xb, yb in layer["draws"]:
            w = max(1, int(round(aps[ap][1] * dpi)))
            dr.line([P(xa, ya), P(xb, yb)], fill=col, width=w)
            r = aps[ap][1] * dpi / 2
            for (xx, yy) in ((xa, ya), (xb, yb)):
                cx, cy = P(xx, yy)
                dr.ellipse([cx - r, cy - r, cx + r, cy + r], fill=col)
        for ap, xa, ya in layer["flashes"]:
            r = aps[ap][1] * dpi / 2
            cx, cy = P(xa, ya)
            dr.ellipse([cx - r, cy - r, cx + r, cy + r], fill=col)
    if path:
        img.save(path)
    return img


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------
class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def f(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def u(self, a, b):
        a, b = self.f(a), self.f(b)
        if a != b:
            self.p[a] = b


def seg_point_dist(p, a, b):
    p, a, b = (np.asarray(v, float) for v in (p, a, b))
    ab = b - a
    L2 = float(ab @ ab)
    if L2 == 0:
        return float(np.linalg.norm(p - a))
    t = max(0.0, min(1.0, float((p - a) @ ab) / L2))
    return float(np.linalg.norm(p - (a + t * ab)))


def seg_seg_dist_vec(A0, A1, B0, B1):
    """vectorised min distance between 2-D segments A0A1 (n,2) and B0B1 (m,2) ->(n,m).
    Exact for non-intersecting segments; returns 0 if they properly cross."""
    def pseg(P, S0, S1):
        d = S1 - S0
        L2 = np.sum(d * d, -1)
        L2 = np.where(L2 == 0, 1e-30, L2)
        t = np.clip(np.sum((P - S0) * d, -1) / L2, 0, 1)
        Q = S0 + t[..., None] * d
        return np.sqrt(np.sum((P - Q) ** 2, -1))
    a0, a1 = A0[:, None, :], A1[:, None, :]
    b0, b1 = B0[None, :, :], B1[None, :, :]
    d = np.minimum.reduce([pseg(a0, b0, b1), pseg(a1, b0, b1),
                           pseg(b0, a0, a1), pseg(b1, a0, a1)])

    def cross(o, p, q):
        return (p[..., 0] - o[..., 0]) * (q[..., 1] - o[..., 1]) - \
               (p[..., 1] - o[..., 1]) * (q[..., 0] - o[..., 0])
    c1 = cross(a0, a1, b0) * cross(a0, a1, b1)
    c2 = cross(b0, b1, a0) * cross(b0, b1, a1)
    d = np.where((c1 < 0) & (c2 < 0), 0.0, d)
    return d


def chain_segments(draws, tol=0.0002):
    """Group draw segments into chains by endpoint coincidence (graph walk).

    tol: endpoint snap tolerance in inches (Gerber grid is 0.001 in; must be
    < 0.001 or consecutive arc vertices would be merged).
    Returns node_xy, edges, degree dict, list of node-index paths.
    """
    pts = []
    for (_, xa, ya, xb, yb) in draws:
        pts.append((xa, ya))
        pts.append((xb, yb))
    P = np.array(pts)
    n = len(P)
    dsu = DSU(n)
    order = np.lexsort((P[:, 1], P[:, 0]))
    for ii in range(n):
        i = order[ii]
        for jj in range(ii + 1, n):
            j = order[jj]
            if P[j, 0] - P[i, 0] > tol:
                break
            if abs(P[j, 1] - P[i, 1]) <= tol:
                dsu.u(i, j)
    node_of = [dsu.f(i) for i in range(n)]
    nodes = sorted(set(node_of))
    nid = {r: k for k, r in enumerate(nodes)}
    node_xy = np.zeros((len(nodes), 2))
    cnt = np.zeros(len(nodes))
    for i in range(n):
        node_xy[nid[node_of[i]]] += P[i]
        cnt[nid[node_of[i]]] += 1
    node_xy /= cnt[:, None]
    edges = []
    adj = defaultdict(list)
    for s in range(len(draws)):
        a = nid[node_of[2 * s]]
        b = nid[node_of[2 * s + 1]]
        if a == b:
            continue
        e = len(edges)
        edges.append((a, b))
        adj[a].append(e)
        adj[b].append(e)
    deg = {k: len(v) for k, v in adj.items()}
    used = [False] * len(edges)
    chains = []
    starts = [k for k in adj if deg[k] == 1] + list(adj.keys())
    for s in starts:
        for e0 in adj[s]:
            if used[e0]:
                continue
            path = [s]
            cur, e = s, e0
            while True:
                used[e] = True
                a, b = edges[e]
                nxt = b if a == cur else a
                path.append(nxt)
                cur = nxt
                cand = [x for x in adj[cur] if not used[x]]
                if deg[cur] != 2 or not cand:
                    break
                e = cand[0]
            chains.append(path)
    return node_xy, edges, deg, chains


def polyline_len(xy):
    d = np.diff(np.asarray(xy, float), axis=0)
    return float(np.sum(np.hypot(d[:, 0], d[:, 1])))


def winding_angle(xy, centre):
    v = np.asarray(xy, float) - np.asarray(centre, float)
    ang = np.unwrap(np.arctan2(v[:, 1], v[:, 0]))
    return float(ang[-1] - ang[0])


def dp_simplify(xy, tol):
    """Douglas-Peucker (iterative). tol in same units as xy. Keeps endpoints."""
    xy = np.asarray(xy, float)
    keep = np.ones(len(xy), bool)
    keep[1:] = np.any(np.abs(np.diff(xy, axis=0)) > 1e-12, axis=1)
    xy = xy[keep]
    if tol <= 0 or len(xy) < 3:
        return xy
    mark = np.zeros(len(xy), bool)
    mark[0] = mark[-1] = True
    stack = [(0, len(xy) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, b = xy[i], xy[j]
        ab = b - a
        L = np.hypot(*ab)
        seg = xy[i + 1:j]
        if L == 0:
            d = np.hypot(*(seg - a).T)
        else:
            d = np.abs(ab[0] * (seg[:, 1] - a[1]) - ab[1] * (seg[:, 0] - a[0])) / L
        k = int(np.argmax(d))
        if d[k] > tol:
            m = i + 1 + k
            mark[m] = True
            stack.append((i, m))
            stack.append((m, j))
    return xy[mark]


def resample(xy, step):
    """points every `step` along polyline (with unit tangents and length fraction)."""
    xy = np.asarray(xy, float)
    d = np.diff(xy, axis=0)
    L = np.hypot(d[:, 0], d[:, 1])
    s = np.concatenate([[0], np.cumsum(L)])
    S = np.arange(step / 2, s[-1], step)
    idx = np.clip(np.searchsorted(s, S) - 1, 0, len(L) - 1)
    t = (S - s[idx]) / np.where(L[idx] == 0, 1, L[idx])
    P = xy[idx] + t[:, None] * d[idx]
    T = d[idx] / np.where(L[idx] == 0, 1, L[idx])[:, None]
    return P, T, S / s[-1]


# --------------------------------------------------------------------------
# legend label reading (render + connected components + tesseract)
# --------------------------------------------------------------------------
def find_label_blobs(img, dpi, mirror, zone):
    """return list of (board-inch bbox, pixel bbox) of word blobs whose centre is in x-zone."""
    import cv2
    a = (np.asarray(img) < 128).astype(np.uint8)
    kx = int(0.06 * dpi)
    k = cv2.getStructuringElement(cv2.MORPH_RECT, (kx, max(1, int(0.02 * dpi))))
    dil = cv2.dilate(a, k)
    n, lab, stats, _ = cv2.connectedComponentsWithStats(dil, 8)
    out = []
    for i in range(1, n):
        x, y, w, h, area = stats[i]
        if h > 0.12 * dpi or w > 0.6 * dpi:
            continue
        bx0 = x / dpi
        bx1 = (x + w) / dpi
        if mirror:
            bx0, bx1 = 3.370 - bx1, 3.370 - bx0
        by1 = 2.125 - y / dpi
        by0 = 2.125 - (y + h) / dpi
        cx = 0.5 * (bx0 + bx1)
        if zone[0] <= cx <= zone[1]:
            out.append(((bx0, by0, bx1, by1), (int(x), int(y), int(w), int(h))))
    return out


def ocr_crop(img, box_px, path, psm=7):
    x, y, w, h = box_px
    pad = 20
    crop = img.crop((max(0, x - pad), max(0, y - pad), x + w + pad, y + h + pad))
    crop.save(path)
    try:
        r = subprocess.run(["tesseract", path, "stdout", "--psm", str(psm)],
                           capture_output=True, text=True, timeout=60)
        return r.stdout.strip()
    except Exception as e:  # pragma: no cover
        return "OCR-failed:%s" % e


# --------------------------------------------------------------------------
def main():
    files = read_zip()
    parsed = {k: parse_gerber(v) for k, v in files.items()}
    prof = parsed["card_profile.gbr"]
    top = parsed["card_copper_signal_top.gbr"]
    # unpopulated "Noise injection" footprint: not part of either winding
    top["layers"][0]["footprint_pads"] = FOOTPRINT_PADS.get("card_copper_signal_top.gbr", [])
    round_pads = [(x, y) for ap, x, y in top["layers"][0]["flashes"] if ap == 11]
    n_fp = drop_footprint(top["layers"][0], round_pads)
    bot = parsed["card_copper_signal_bot.gbr"]
    drl = parsed["card.txt"]
    mtop = parsed["card_soldermask_top.gbr"]
    G = {"units": "inch unless suffixed _mm", "source": os.path.basename(ZIP),
         "footprint_pads_set_aside": len(FOOTPRINT_PADS.get("card_copper_signal_top.gbr", [])),
         "footprint_draws_set_aside": n_fp}

    # ---------------- outline
    pd = prof["layers"][0]["draws"]
    opts = np.array([(d[1], d[2]) for d in pd] + [(d[3], d[4]) for d in pd])
    G["outline_bbox_in"] = [float(opts[:, 0].min()), float(opts[:, 1].min()),
                            float(opts[:, 0].max()), float(opts[:, 1].max())]
    G["outline_W_mm"] = float(np.ptp(opts[:, 0]) * IN2MM)
    G["outline_H_mm"] = float(np.ptp(opts[:, 1]) * IN2MM)
    bottom_edge = [d for d in pd if d[2] == 0 and d[4] == 0 and abs(d[3] - d[1]) > 1.0]
    xs = sorted([bottom_edge[0][1], bottom_edge[0][3]])
    G["outline_corner_radius_mm"] = float(min(xs[0], 3.370 - xs[1]) * IN2MM)
    G["outline_poly_in"] = [[d[1], d[2]] for d in pd] + [[pd[-1][3], pd[-1][4]]]

    # ---------------- apertures
    G["apertures"] = {}
    for fname in sorted(parsed):
        for k, (shape, dia) in parsed[fname]["apertures"].items():
            G["apertures"]["%s:D%d" % (fname, k)] = {"shape": shape, "dia_in": dia,
                                                      "dia_mm": dia * IN2MM}
    trace_w_in = top["apertures"][12][1]
    G["trace_width_mm"] = trace_w_in * IN2MM

    # ---------------- pads, vias, drills
    pads = [(x, y) for ap, x, y in top["layers"][0]["flashes"] if ap == 11]
    vias = [(x, y) for ap, x, y in top["layers"][0]["flashes"] if ap == 14]
    odd = [(ap, x, y) for ap, x, y in top["layers"][0]["flashes"] if ap not in (11, 14)]
    G["pad_dia_mm"] = top["apertures"][11][1] * IN2MM
    G["via_pad_dia_mm"] = top["apertures"][14][1] * IN2MM
    drills = [(drl["apertures"][ap][1], x, y) for ap, x, y in drl["layers"][0]["flashes"]]
    G["drills"] = [{"dia_mm": d * IN2MM, "x_in": x, "y_in": y} for d, x, y in drills]
    G["mask_openings_top"] = [(x, y) for ap, x, y in mtop["layers"][0]["flashes"]]
    G["mask_opening_dia_mm"] = mtop["apertures"][11][1] * IN2MM
    G["vias_in"] = vias
    G["pads_in"] = pads

    # ---------------- top copper chains
    tdraws = top["layers"][0]["draws"]
    assert {d[0] for d in tdraws} == {12}
    node_xy, edges, deg, chains = chain_segments(tdraws)
    G["top_n_draws"] = len(tdraws)
    G["top_n_zero_length_draws"] = int(sum(1 for d in tdraws if d[1] == d[3] and d[2] == d[4]))
    G["top_n_nodes"] = len(node_xy)
    G["top_deg_hist"] = {str(k): sum(1 for v in deg.values() if v == k)
                         for k in sorted(set(deg.values()))}
    assert len(chains) == 2 and G["top_deg_hist"] == {"1": 4, "2": len(node_xy) - 4}
    G["top_other_flashes"] = []
    for ap, x, y in odd:
        dmin = min(seg_point_dist((x, y), (d[1], d[2]), (d[3], d[4])) for d in tdraws)
        G["top_other_flashes"].append({"aperture": "D%d" % ap, "dia_mm": top["apertures"][ap][1] * IN2MM,
                                       "x_in": x, "y_in": y,
                                       "dist_to_nearest_trace_centreline_mm": dmin * IN2MM})

    # ---------------- bottom copper
    bd = bot["layers"][0]["draws"]
    b12 = [d for d in bd if d[0] == 12]
    b15 = [d for d in bd if d[0] == 15]
    bnode, bedges, bdeg, bchains = chain_segments(b12)
    bottom_routes = [[[float(v) for v in p] for p in bnode[c]] for c in bchains]
    G["bottom_routes_in"] = bottom_routes
    G["bottom_wide_traces"] = [{"from_in": [d[1], d[2]], "to_in": [d[3], d[4]],
                                "width_mm": bot["apertures"][15][1] * IN2MM} for d in b15]

    # ---------------- legend: render, OCR
    leg = {}
    for name, mirror in (("card_legend_top.gbr", False), ("card_legend_bot.gbr", True)):
        tag = "legend_top" if not mirror else "legend_bot_mirrored"
        img = render(parsed[name], os.path.join(WORK, tag + "_1200dpi.png"), dpi=1200, mirror=mirror)
        img.reduce(2).save(os.path.join(WORK, tag + "_600dpi.png"))
        blobs = []
        for zone in ((0.0, 0.40), (2.95, 3.37)):
            for bb, px in find_label_blobs(img, 1200, mirror, zone):
                txt = ocr_crop(img, px, os.path.join(WORK, "ocr_%s_%.3f_%.3f.png" % (tag, bb[0], bb[1])))
                blobs.append({"bbox_in": [round(v, 4) for v in bb], "ocr": txt})
        try:
            r = subprocess.run(["tesseract", os.path.join(WORK, tag + "_1200dpi.png"), "stdout",
                                "--psm", "4"], capture_output=True, text=True, timeout=300)
            full = r.stdout
        except Exception as e:  # pragma: no cover
            full = "OCR failed: %s" % e
        with open(os.path.join(WORK, tag + "_ocr.txt"), "w") as f:
            f.write(full)
        leg[tag] = {"label_blobs": blobs, "ocr_full_text": full}
    G["legend"] = leg

    # ---------------- terminal mapping (label -> nearest pad on the same board edge)
    expected = ["DC_IN", "DC_RTN", "OUT_RTN", "OUT", "AC_IN", "AC_RTN"]

    def norm(s):
        return re.sub(r"[^A-Z_]", "", s.upper().replace(" ", "_"))
    label_pos = {}
    for b in leg["legend_top"]["label_blobs"]:
        t = norm(b["ocr"])
        for e in expected:
            if t == e or t.replace("_", "") == e.replace("_", ""):
                x0, y0, x1, y1 = b["bbox_in"]
                label_pos[e] = (0.5 * (x0 + x1), 0.5 * (y0 + y1))
    G["label_positions_top_in"] = label_pos
    missing = [e for e in expected if e not in label_pos]
    G["label_ocr_missing"] = missing
    left = sorted([p for p in pads if p[0] < 1.0], key=lambda p: -p[1])
    right = sorted([p for p in pads if p[0] > 2.0], key=lambda p: -p[1])
    terminals = {}
    assert not missing, "OCR did not find labels %s" % missing
    for side_pads, names in ((left, ["DC_IN", "DC_RTN", "OUT_RTN", "OUT"]),
                             (right, ["AC_IN", "AC_RTN"])):
        labs = sorted(names, key=lambda e: -label_pos[e][1])
        assert labs == names, "label order differs from expectation"
        for p, e in zip(side_pads, labs):
            dists = [math.hypot(p2[0] - label_pos[e][0], p2[1] - label_pos[e][1]) for p2 in side_pads]
            assert side_pads[int(np.argmin(dists))] == p, "nearest-pad mismatch for " + e
            terminals[e] = list(p)
    G["terminals_in"] = terminals
    pad_name = {tuple(v): k for k, v in terminals.items()}

    # ---------------- winding identification via connectivity
    vset = [tuple(v) for v in vias]
    via_to_pad = {}
    for r in bottom_routes:
        a, b = tuple(r[0]), tuple(r[-1])
        if a in vset:
            via_to_pad[a] = b
        elif b in vset:
            via_to_pad[b] = a
    windings = {}
    for c in chains:
        xy = node_xy[c]
        e0, e1 = tuple(np.round(xy[0], 4)), tuple(np.round(xy[-1], 4))
        if e0 in vset:
            xy = xy[::-1]  # orient pad(outer) -> via(inner)
            e0, e1 = e1, e0
        outer_pad = pad_name[e0]
        inner_pad = pad_name[tuple(via_to_pad[e1])]
        name = "AC" if outer_pad.startswith("AC") else "DC"
        windings[name] = {"outer_terminal": outer_pad, "inner_terminal": inner_pad,
                          "inner_via_in": list(e1), "top_xy_in": xy.tolist()}
        broute = [r for r in bottom_routes if tuple(r[0]) == e1 or tuple(r[-1]) == e1][0]
        if tuple(broute[0]) != e1:
            broute = broute[::-1]
        windings[name]["bottom_xy_in"] = broute
    G["d15_connects"] = [[pad_name.get(tuple(d["from_in"])), pad_name.get(tuple(d["to_in"]))]
                         for d in G["bottom_wide_traces"]]
    for f in G["top_other_flashes"]:
        best = None
        for n, w in windings.items():
            xy = np.array(w["top_xy_in"])
            dm = min(seg_point_dist((f["x_in"], f["y_in"]), xy[i], xy[i + 1]) for i in range(len(xy) - 1))
            if best is None or dm < best[1]:
                best = (n, dm)
        f["on_winding"] = best[0]

    # ---------------- per-winding metrics
    allxy = np.vstack([np.array(w["top_xy_in"]) for w in windings.values()])
    inb_all = (allxy[:, 0] > 0.30) & (allxy[:, 0] < 3.05)
    body = allxy[inb_all]
    bx0, by0 = body.min(0)
    bx1, by1 = body.max(0)
    centre = (0.5 * (bx0 + bx1), 0.5 * (by0 + by1))
    G["body_centreline_bbox_in"] = [float(bx0), float(by0), float(bx1), float(by1)]
    G["body_centre_in"] = [float(centre[0]), float(centre[1])]
    G["body_outer_copper_W_mm"] = float((bx1 - bx0 + trace_w_in) * IN2MM)
    G["body_outer_copper_H_mm"] = float((by1 - by0 + trace_w_in) * IN2MM)

    for name, w in windings.items():
        xy = np.array(w["top_xy_in"])
        bxy = np.array(w["bottom_xy_in"])
        w["n_vertices_raw"] = len(xy)
        w["top_len_raw_mm"] = polyline_len(xy) * IN2MM
        for tol_mil in (0.3, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0):
            w["top_len_dp%.1fmil_mm" % tol_mil] = polyline_len(dp_simplify(xy, tol_mil / 1000)) * IN2MM
        w["n_vertices_dp1.0mil"] = len(dp_simplify(xy, 0.0010))
        w["bottom_len_mm"] = polyline_len(bxy) * IN2MM
        ang = winding_angle(xy, centre)
        w["winding_angle_deg_outer_to_inner"] = math.degrees(ang)
        w["turns"] = abs(ang) / (2 * math.pi)
        w["sense_outer_to_inner"] = "clockwise (viewed from top)" if ang < 0 else "counter-clockwise"
        inb = (xy[:, 0] > 0.30) & (xy[:, 0] < 3.05)
        w["body_bbox_in"] = [float(xy[inb, 0].min()), float(xy[inb, 1].min()),
                             float(xy[inb, 0].max()), float(xy[inb, 1].max())]
        k = int(np.argmax(inb))
        w["outer_lead_len_mm"] = polyline_len(xy[:k + 1]) * IN2MM

    # ---------------- cross-section scans (pitch, gap, interleave order, inner hole)
    segs = []
    for name, w in windings.items():
        xy = np.array(w["top_xy_in"])
        for i in range(len(xy) - 1):
            segs.append((name, xy[i], xy[i + 1]))

    def cross(axis, c0):
        out = []
        u, v = (0, 1) if axis == "v" else (1, 0)
        for n, a, b in segs:
            if (a[u] - c0) * (b[u] - c0) < 0:
                t = (c0 - a[u]) / (b[u] - a[u])
                out.append((a[v] + t * (b[v] - a[v]), n))
        return sorted(out)
    scans = {}
    for nm, axis, c0 in (("vertical_through_centre", "v", centre[0] + 0.0003),
                         ("horizontal_through_centre", "h", centre[1] + 0.0003)):
        cr = cross(axis, c0)
        pos = np.array([c[0] for c in cr])
        order = "".join("A" if c[1] == "AC" else "D" for c in cr)
        dd = np.diff(pos) * IN2MM
        gap_i = int(np.argmax(dd))
        scans[nm] = {"n_crossings": len(cr), "order": order,
                     "pitch_mm_median": float(np.median(dd)),
                     "pitch_mm_min": float(dd[dd < 1.0].min()),
                     "pitch_mm_max_excl_hole": float(dd[dd < 1.0].max()),
                     "inner_hole_centreline_mm": float(dd[gap_i]),
                     "positions_in": pos.tolist()}
    G["scans"] = scans
    pitch = scans["vertical_through_centre"]["pitch_mm_median"]
    G["pitch_AD_mm"] = pitch
    G["pitch_same_winding_mm"] = 2 * pitch
    G["gap_AD_drawn_mm"] = pitch - G["trace_width_mm"]
    G["inner_hole_centreline_W_mm"] = scans["horizontal_through_centre"]["inner_hole_centreline_mm"]
    G["inner_hole_centreline_H_mm"] = scans["vertical_through_centre"]["inner_hole_centreline_mm"]

    # ---------------- clearances (centreline distance -> edge gap = d - w)
    def seg_arrays(xy):
        xy = dp_simplify(np.asarray(xy), 0.0003)
        return xy[:-1], xy[1:]
    A0, A1 = seg_arrays(windings["AC"]["top_xy_in"])
    D0, D1 = seg_arrays(windings["DC"]["top_xy_in"])
    dAD = seg_seg_dist_vec(A0, A1, D0, D1)
    G["min_centreline_dist_AC_DC_mm"] = float(dAD.min() * IN2MM)
    G["min_edge_gap_AC_DC_mm"] = float((dAD.min() - trace_w_in) * IN2MM)
    ia, idd = np.unravel_index(np.argmin(dAD), dAD.shape)
    G["min_gap_AC_DC_location_in"] = [float(A0[ia][0]), float(A0[ia][1])]

    def self_min(S0, S1):
        # exclude pairs closer than 3 mm along the conductor (same corner / arc)
        d = seg_seg_dist_vec(S0, S1, S0, S1)
        Ls = np.hypot(*(S1 - S0).T)
        s0 = np.concatenate([[0], np.cumsum(Ls)[:-1]])
        s1 = s0 + Ls
        along = np.maximum(np.maximum(s0[None, :] - s1[:, None], s0[:, None] - s1[None, :]), 0)
        d = np.where(along < 3.0 / IN2MM, np.inf, d)
        return float(d.min())
    G["min_centreline_dist_AC_self_nonadjacent_mm"] = self_min(A0, A1) * IN2MM
    G["min_centreline_dist_DC_self_nonadjacent_mm"] = self_min(D0, D1) * IN2MM
    padclr = {}
    for pn, p in terminals.items():
        for wn, (S0, S1) in (("AC", (A0, A1)), ("DC", (D0, D1))):
            d = min(seg_point_dist(p, S0[i], S1[i]) for i in range(len(S0)))
            padclr["%s_to_%s_top_edge_gap_mm" % (pn, wn)] = (d - trace_w_in / 2 - top["apertures"][11][1] / 2) * IN2MM
    G["pad_clearances"] = padclr

    # ---------------- A-D adjacency (for inter-winding capacitance)
    step_mm = 0.1
    P, T, frac = resample(windings["AC"]["top_xy_in"], step_mm / IN2MM)
    Nrm = np.stack([-T[:, 1], T[:, 0]], 1)
    Dxy = dp_simplify(np.array(windings["DC"]["top_xy_in"]), 0.0003)
    Ds0, Ds1 = Dxy[:-1], Dxy[1:]
    pitch_in = pitch / IN2MM
    adj_count = np.zeros(len(P))
    for s in (+1, -1):
        Q = P + s * Nrm * pitch_in
        hit = np.zeros(len(P), bool)
        for c0 in range(0, len(Q), 1000):
            d = seg_seg_dist_vec(Q[c0:c0 + 1000], Q[c0:c0 + 1000], Ds0, Ds1).min(1)
            hit[c0:c0 + 1000] = d < 0.25 * pitch_in
        adj_count += hit
    G["adjacent_length_AD_mm"] = float(adj_count.sum() * step_mm)
    G["adjacency_hist"] = {str(k): int((adj_count == k).sum()) for k in (0, 1, 2)}
    G["adjacency_step_mm"] = step_mm

    G["windings"] = windings
    with open(os.path.join(WORK, "geometry.json"), "w") as f:
        json.dump(G, f, indent=1, default=float)
    return G


def write_outputs(G):
    W = G["windings"]
    tbl = []

    def row(param, val, unit, label, method):
        tbl.append([param, val, unit, label, method])
    row("board outline W x H", "%.3f x %.3f" % (G["outline_W_mm"], G["outline_H_mm"]), "mm",
        "[DESIGN]", "card_profile.gbr extents")
    row("outline corner radius", "%.2f" % G["outline_corner_radius_mm"], "mm", "[DESIGN]",
        "card_profile.gbr straight-edge start points")
    row("drawn trace width (D12)", "%.3f" % G["trace_width_mm"], "mm", "[DESIGN]",
        "card_copper_signal_top.gbr %ADD12C,0.010000*%")
    row("centre pitch A-D (adjacent conductors)", "%.3f" % G["pitch_AD_mm"], "mm", "[CALCULATED]",
        "median spacing of centreline crossings, vertical scan through body centre")
    row("pitch within one winding", "%.3f" % G["pitch_same_winding_mm"], "mm", "[CALCULATED]", "2 x A-D pitch")
    row("drawn edge gap A-D", "%.3f" % G["gap_AD_drawn_mm"], "mm", "[CALCULATED]", "pitch - drawn width")
    row("min centreline distance A-D (whole board)", "%.3f" % G["min_centreline_dist_AC_DC_mm"], "mm",
        "[CALCULATED]", "segment-segment distance, DP 0.3 mil polylines")
    row("min drawn edge gap A-D", "%.3f" % G["min_edge_gap_AC_DC_mm"], "mm", "[CALCULATED]",
        "min centreline distance - drawn width")
    row("spiral body outer copper W x H", "%.2f x %.2f" % (G["body_outer_copper_W_mm"], G["body_outer_copper_H_mm"]),
        "mm", "[CALCULATED]", "bbox of centrelines + drawn width")
    row("inner hole (centreline) W x H", "%.2f x %.2f" % (G["inner_hole_centreline_W_mm"], G["inner_hole_centreline_H_mm"]),
        "mm", "[CALCULATED]", "largest crossing gap in horizontal/vertical scans")
    for n in ("AC", "DC"):
        w = W[n]
        row("%s winding terminals (outer -> inner)" % n, "%s -> %s" % (w["outer_terminal"], w["inner_terminal"]),
            "-", "[CALCULATED]", "chain connectivity + via + bottom route; legend label nearest pad")
        row("%s top-layer length (raw Gerber vertices)" % n, "%.1f" % w["top_len_raw_mm"], "mm", "[CALCULATED]",
            "sum of segment lengths, 1-mil-quantised arc vertices")
        row("%s top-layer length (DP 1.0 mil)" % n, "%.1f" % w["top_len_dp1.0mil_mm"], "mm", "[CALCULATED]",
            "Douglas-Peucker 1.0 mil (> 0.71 mil max vertex-rounding jitter) removes 1-mil grid zig-zag on arcs")
        row("%s bottom-layer route length" % n, "%.1f" % w["bottom_len_mm"], "mm", "[CALCULATED]",
            "card_copper_signal_bot.gbr D12 route via->pad")
        row("%s total length (DP 1.0 mil top + bottom)" % n, "%.1f" % (w["top_len_dp1.0mil_mm"] + w["bottom_len_mm"]),
            "mm", "[CALCULATED]", "excludes 1.6 mm via barrel")
        row("%s turns (winding angle outer->inner / 360)" % n, "%.2f" % w["turns"], "turns", "[CALCULATED]",
            "unwrapped polar angle about body centre, pad to via")
        row("%s sense outer->inner" % n, w["sense_outer_to_inner"], "-", "[CALCULATED]", "sign of winding angle")
        bb = w["body_bbox_in"]
        row("%s body centreline extents" % n, "%.2f x %.2f" % ((bb[2] - bb[0]) * IN2MM, (bb[3] - bb[1]) * IN2MM),
            "mm", "[CALCULATED]", "bbox of chain inside x in [0.30,3.05] in")
    row("A-D adjacent length (sum over both sides)", "%.0f" % G["adjacent_length_AD_mm"], "mm", "[CALCULATED]",
        "AC centreline sampled every 0.1 mm; DC conductor found at +-1 pitch")
    row("pad dia / drill", "%.2f / %.2f" % (G["pad_dia_mm"], G["drills"][-1]["dia_mm"]), "mm", "[DESIGN]",
        "D11 / nc_drill D9501")
    row("via pad dia / drill", "%.2f / %.2f" % (G["via_pad_dia_mm"], G["drills"][0]["dia_mm"]), "mm", "[DESIGN]",
        "D14 / nc_drill D9500")
    with open(os.path.join(DATA, "geometry_table.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["parameter", "value", "unit", "label", "method_source"])
        cw.writerows(tbl)

    m = lambda s: "[CALCULATED] gerber_geometry.py from %s: %s" % (SRC, s)
    d = lambda s: "[DESIGN] %s: %s" % (SRC, s)
    e = {
        "theory.geom.board_W_mm": (sig(G["outline_W_mm"], 5), d("card_profile.gbr extent")),
        "theory.geom.board_H_mm": (sig(G["outline_H_mm"], 5), d("card_profile.gbr extent")),
        "theory.geom.trace_width_drawn_mm": (sig(G["trace_width_mm"]), d("top copper aperture D12 = 0.010 in")),
        "theory.geom.pitch_AD_mm": (sig(G["pitch_AD_mm"]), m("median centreline spacing (vertical scan)")),
        "theory.geom.gap_AD_drawn_mm": (sig(G["gap_AD_drawn_mm"]), m("pitch minus drawn width")),
        "theory.geom.min_gap_AD_drawn_mm": (sig(G["min_edge_gap_AC_DC_mm"]), m("global min centreline distance minus width")),
        "theory.geom.body_W_mm": (sig(G["body_outer_copper_W_mm"]), m("outer copper bbox of spiral body")),
        "theory.geom.body_H_mm": (sig(G["body_outer_copper_H_mm"]), m("outer copper bbox of spiral body")),
        "theory.geom.inner_hole_W_mm": (sig(G["inner_hole_centreline_W_mm"]), m("centreline gap, horizontal scan")),
        "theory.geom.inner_hole_H_mm": (sig(G["inner_hole_centreline_H_mm"]), m("centreline gap, vertical scan")),
        "theory.geom.adjacent_length_AD_mm": (sig(G["adjacent_length_AD_mm"]), m("sampled adjacency at +-1 pitch")),
    }
    for n, key in (("AC", "A"), ("DC", "D")):
        w = W[n]
        e["theory.geom.%s_turns" % key] = (sig(w["turns"]), m("winding angle/360 pad->via, %s winding" % n))
        e["theory.geom.%s_len_top_mm" % key] = (sig(w["top_len_dp1.0mil_mm"], 5), m("top length, DP 1.0 mil"))
        e["theory.geom.%s_len_top_raw_mm" % key] = (sig(w["top_len_raw_mm"], 5), m("top length, raw 1-mil vertices"))
        e["theory.geom.%s_len_bottom_mm" % key] = (sig(w["bottom_len_mm"], 4), m("bottom route length"))
        e["theory.geom.%s_len_total_mm" % key] = (sig(w["top_len_dp1.0mil_mm"] + w["bottom_len_mm"], 5),
                                                  m("top (DP 1.0 mil) + bottom"))
        e["theory.geom.%s_outer_terminal" % key] = (w["outer_terminal"], m("connectivity + legend label"))
        e["theory.geom.%s_inner_terminal" % key] = (w["inner_terminal"], m("connectivity via + bottom route"))
    # legend transcription: visual reading of work/theory/legend_top_600dpi.png and
    # legend_bot_mirrored_600dpi.png (rendered by this script); corroborated by the tesseract
    # output in work/theory/legend_top_ocr.txt / legend_bot_mirrored_ocr.txt (OCR renders mu as 'u',
    # Omega as 'O', +- as '+'; the visual reading is authoritative).
    rd = "[DESIGN] text READ visually from rendered %s (gerber.zip), corroborated by tesseract OCR"
    e["theory.legend.top_spec_line"] = ("Two interleaved spirals. __ \u00b5H \u00b7 k \u2248 0.8-0.9 \u00b7 DCR 5-7 \u03a9",
                                         rd % "card_legend_top.gbr")
    e["theory.legend.top_tolerance_line"] = ("Values nominal, etch tolerance \u00b115%", rd % "card_legend_top.gbr")
    e["theory.legend.top_limit_line"] = ("MAX 30 V \u00b7 500 mA \u00b7 USE A CURRENT-LIMITED SUPPLY", rd % "card_legend_top.gbr")
    e["theory.legend.top_terminal_labels"] = ("left edge top->bottom: DC_IN, DC_RTN, OUT_RTN, OUT; right edge: AC_IN, AC_RTN",
                                               rd % "card_legend_top.gbr" + "; positions from connected-component OCR blobs")
    e["theory.legend.bot_lines"] = ("Trace width = 0.254 mm; Trace clearance = 0.254 mm; Copper thickness = 35 \u00b5m (1 oz); "
                                    "DC coil length = 3002 mm; AC coil length = 2999 mm",
                                    rd % "card_legend_bot.gbr (mirrored for reading)")
    e["theory.legend.matches_artwork_text"] = (False, "[CALCULATED] comparison: artwork (FRONT pdf, pdftotext) reads "
                                               "'~20 uH \u00b7 k~0.9 \u00b7 DCR 9 ohm' and 'MAX 30 V \u00b7 50 mA \u00b7 VALUES NOMINAL +/-15%' "
                                               "with GND_S/GND_D labels; Gerber legend has blank uH field, k 0.8-0.9, DCR 5-7 ohm, "
                                               "500 mA, and DC_RTN/OUT_RTN instead of GND_S/GND_D")
    if G["top_other_flashes"]:
        e["theory.geom.stray_flash_D13"] = ("0.254 mm round flash at (%.3f, %.3f) in lies on a DC-winding vertex: redundant, no electrical effect"
                                            % (G["top_other_flashes"][0]["x_in"], G["top_other_flashes"][0]["y_in"]),
                                            m("distance to DC centreline %.3f mm" % G["top_other_flashes"][0]["dist_to_nearest_trace_centreline_mm"]))
    e["theory.geom.rtn_pads_link"] = ("DC_RTN and OUT_RTN joined by a 1.0 mm bottom-layer trace; no connection to either winding",
                                      m("bottom copper D15 draw + connectivity"))
    e["theory.geom.top_soldermask_openings"] = ("only the 6 terminal pads (1.6 mm); spiral is mask-covered in the design",
                                                d("card_soldermask_top.gbr / _bot.gbr flashes"))
    e["theory.geom.dotted_pair"] = ("DC_IN & AC_RTN (outer ends); equivalently OUT & AC_IN (inner ends)",
                                    m("both windings run clockwise outer->inner (viewed from top)"))
    update_values(e)


if __name__ == "__main__":
    G = main()
    write_outputs(G)
    show = {k: v for k, v in G.items() if k not in ("outline_poly_in", "windings", "legend", "scans", "apertures")}
    print(json.dumps(show, indent=1, default=float)[:6000])
    for n, w in G["windings"].items():
        print(n, {k: v for k, v in w.items() if not k.endswith("xy_in")})
    for k, s in G["scans"].items():
        print(k, {a: b for a, b in s.items() if a != "positions_in"})
    for t, l in G["legend"].items():
        print(t, [(b["bbox_in"], b["ocr"]) for b in l["label_blobs"]])
