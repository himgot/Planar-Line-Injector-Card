#!/usr/bin/env python3
"""Fig. windings: the two interleaved windings from the Gerber polylines.

Run:  python3 -I analysis/fig_windings.py
Output: paper/src/figures/fig_windings.pdf (+ .png preview in work/theory)

(a) whole card: board outline, AC winding (blue, solid) and DC winding (orange, solid;
    dashed in the zoom), bottom-layer routes (dotted), vias, pads with terminal names,
    dot markers on the dotted (outer) ends, 10 mm scale bar, zoom window.
(b) zoom of the interleave at true trace width with pitch / gap annotation.
Colours: validated pair (#2a62b0, #f08a4b) that also separates in greyscale
(relative luminance ~0.12 vs ~0.37); identity is also carried by labels and line style.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, Rectangle, ConnectionPatch  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gerber_geometry as gg  # noqa: E402

OUT = "paper/src/figures/fig_windings.pdf"
MM = 25.4
C_AC = "#2a62b0"
C_DC = "#f08a4b"
INK = "#0b0b0b"
INK2 = "#52514e"

plt.rcParams.update({
    "font.size": 8, "font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.5, "legend.fontsize": 8,
})


def main():
    G = gg.load_geometry()
    W = G["windings"]
    T = {k: np.array(v) * MM for k, v in G["terminals_in"].items()}
    outline = np.array(G["outline_poly_in"]) * MM
    fig = plt.figure(figsize=(3.5, 3.75))
    # ---------------- panel (a)
    ax = fig.add_axes([0.0, 0.43, 1.0, 0.57])
    xl, xr = -21.0, 85.6 + 17.5
    ax.set_xlim(xl, xr)
    ax.set_ylim(-6.5, 54.6)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.plot(outline[:, 0], outline[:, 1], color=INK, lw=0.6)
    lw_tr = 0.5  # ~ true 0.254 mm at this scale
    for name, col in (("AC", C_AC), ("DC", C_DC)):
        xy = gg.dp_simplify(np.array(W[name]["top_xy_in"]), 0.0005) * MM
        ax.plot(xy[:, 0], xy[:, 1], color=col, lw=lw_tr, solid_capstyle="round", solid_joinstyle="round", zorder=3)
        b = np.array(W[name]["bottom_xy_in"]) * MM
        ax.plot(b[:, 0], b[:, 1], color=col, lw=1.0, ls=(0, (1.2, 1.2)), zorder=4)
        v = np.array(W[name]["inner_via_in"]) * MM
        ax.add_patch(Circle(v, 0.3 * 2, fc="white", ec=col, lw=0.8, zorder=5))
    # bottom wide trace DC_RTN-OUT_RTN
    for d in G["bottom_wide_traces"]:
        p0 = np.array(d["from_in"]) * MM
        p1 = np.array(d["to_in"]) * MM
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=INK2, lw=1.4, ls=(0, (1.2, 1.2)), zorder=4)
    # pads
    for nm, p in T.items():
        ax.add_patch(Circle(p, 0.8, fc="#b8b8b0", ec=INK, lw=0.5, zorder=6))
        ax.add_patch(Circle(p, 0.55, fc="white", ec="none", zorder=7))
    # labels: left column right-aligned outside board, right column left-aligned
    lab_dot = {"DC_IN", "AC_RTN"}
    for nm, p in T.items():
        if p[0] < 40:
            ax.text(-1.4, p[1], nm, ha="right", va="center", color=INK, fontsize=8)
        else:
            ax.text(85.6 + 1.4, p[1], nm, ha="left", va="center", color=INK, fontsize=8)
        if nm in lab_dot:
            # dot marker beside the pad, on the board side
            dx = 2.3 if p[0] < 40 else -2.3
            ax.add_patch(Circle((p[0] + dx, p[1] + 2.0), 0.75, fc=INK, ec="none", zorder=8))
    # scale bar
    ax.plot([0, 10], [-3.2, -3.2], color=INK, lw=1.2, solid_capstyle="butt")
    ax.plot([0, 0], [-3.9, -2.5], color=INK, lw=0.6)
    ax.plot([10, 10], [-3.9, -2.5], color=INK, lw=0.6)
    ax.text(11.2, -3.2, "10 mm", ha="left", va="center", fontsize=8, color=INK)
    # zoom window
    zx0, zx1, zy0, zy1 = 9.3, 15.3, 25.0, 28.0
    ax.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, fill=False, ec=INK, lw=0.6, zorder=9))
    ax.text(zx1 + 0.8, zy1 + 0.3, "(b)", fontsize=8, ha="left", va="bottom", color=INK, zorder=9,
            bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))
    ax.text(xl + 0.5, 54.0, "(a)", fontsize=8, ha="left", va="top", color=INK)
    # ---------------- panel (b): zoom at true width
    bx = fig.add_axes([0.12, 0.175, 0.80, 0.235])
    bx.set_xlim(zx0, zx1)
    bx.set_ylim(zy0, zy1)
    bx.set_aspect("equal")
    for s in bx.spines.values():
        s.set_color(INK)
        s.set_linewidth(0.6)
    bx.tick_params(labelsize=8, length=2, width=0.5, colors=INK2)
    bx.set_xlabel("x (mm)", fontsize=8, color=INK2, labelpad=1)
    bx.set_ylabel("y (mm)", fontsize=8, color=INK2, labelpad=1)
    fig.canvas.draw()
    # points per mm in the zoom axes
    bb = bx.get_window_extent()
    ppmm = bb.width / (zx1 - zx0) * 72 / fig.dpi
    w_pt = G["trace_width_mm"] * ppmm
    for name, col, ls in (("AC", C_AC, "-"), ("DC", C_DC, (0, (3.0, 1.2)))):
        xy = np.array(W[name]["top_xy_in"]) * MM
        bx.plot(xy[:, 0], xy[:, 1], color=col, lw=w_pt, ls=ls, solid_capstyle="butt", dash_capstyle="butt",
                zorder=3)
    # annotate which line is which: label the first few traces at top
    pos = np.array(G["scans"]["horizontal_through_centre"]["positions_in"]) * MM
    order = G["scans"]["horizontal_through_centre"]["order"]
    for xx, o in zip(pos, order):
        if zx0 + 0.2 < xx < zx1 - 0.2:
            bx.text(xx, zy1 - 0.15, "A" if o == "A" else "B", ha="center", va="top", fontsize=8,
                    color=INK, zorder=6,
                    bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
    # pitch / gap annotation between first two inner traces in window
    sel = [(xx, o) for xx, o in zip(pos, order) if zx0 + 0.5 < xx < zx1 - 0.5]
    x_a, x_b = sel[2][0], sel[3][0]
    ya = zy0 + 0.75
    bx.annotate("", xy=(x_a, ya), xytext=(x_b, ya),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color=INK, shrinkA=0, shrinkB=0), zorder=7)
    bx.text(0.5 * (x_a + x_b), ya - 0.12, "pitch %.3f mm" % G["pitch_AD_mm"], ha="center", va="top", fontsize=8,
            color=INK, zorder=7, bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))
    xg0 = sel[5][0] + G["trace_width_mm"] / 2
    xg1 = sel[6][0] - G["trace_width_mm"] / 2
    yb = zy0 + 1.75
    bx.annotate("", xy=(xg0, yb), xytext=(xg1, yb),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color=INK, shrinkA=0, shrinkB=0), zorder=7)
    bx.text(xg1 + 0.12, yb, "gap %.3f mm" % G["gap_AD_drawn_mm"], ha="left", va="center", fontsize=8, color=INK,
            zorder=7, bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))
    # connectors from window to zoom
    for (xa_, ya_), (xb_, yb_) in (((zx0, zy0), (zx0, zy1)), ((zx1, zy0), (zx1, zy1))):
        con = ConnectionPatch(xyA=(xa_, ya_), coordsA=ax.transData, xyB=(xb_, yb_), coordsB=bx.transData,
                              color=INK2, lw=0.4, ls=":")
        fig.add_artist(con)
    # ---------------- legend (figure-level, below panel b)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=C_AC, lw=1.6, label="AC winding (A)"),
         Line2D([], [], color=C_DC, lw=1.6, ls=(0, (3.0, 1.2)), label="DC winding (B)"),
         Line2D([], [], color=INK2, lw=1.2, ls=(0, (1.2, 1.2)), label="bottom-layer route"),
         Line2D([], [], color=INK, marker="o", ls="none", ms=4, label="dot (outer ends)")]
    fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, fontsize=8, handlelength=2.0,
               columnspacing=1.0, handletextpad=0.5, borderaxespad=0.1, bbox_to_anchor=(0.5, -0.005))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT)
    fig.savefig(os.path.join(gg.WORK, "fig_windings_preview.png"), dpi=300)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
