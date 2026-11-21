"""Figure of four oscilloscope screens from bench Run R.

Crops the screen area from the original photographs (no other image
processing) and arranges them in a 2 x 2 grid.

Usage:
    python analysis/figure_scope_screens.py [photo_dir] [output_pdf]

Defaults: photos/scope_screens and paper/src/figures/fig_scope_screens.pdf
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

PHOTO_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("photos", "scope_screens")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join("paper", "src", "figures", "fig_scope_screens.pdf")

# reading ID, photo file, panel title, crop box (left, top, right, bottom) in pixels
PANELS = [
    ("R06", "E-scope-oct22-unlisted-200kHz-R06.jpg", "200 kHz", (250, 110, 1124, 760)),
    ("R10", "E-scope-oct22-1003kHz.jpg", "1.00 MHz", (293, 250, 1059, 809)),
    ("R14", "E-scope-oct22-4001kHz-s13.jpg", "3.99 MHz", (293, 232, 1047, 790)),
    ("R16", "E-scope-oct22-5182kHz.jpg", "5.18 MHz", (331, 256, 1051, 811)),
]


def main():
    plt.rcParams.update({"font.family": "serif", "font.size": 7})
    fig, axes = plt.subplots(2, 2, figsize=(3.5, 2.75))
    for ax, (rid, fname, title, box) in zip(axes.ravel(), PANELS):
        im = Image.open(os.path.join(PHOTO_DIR, fname)).convert("RGB").crop(box)
        im = im.resize((800, int(800 * im.height / im.width)), Image.LANCZOS)
        ax.imshow(im)
        ax.set_title(f"{title} (Run R) [{rid}]", pad=2)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.01, wspace=0.04, hspace=0.16)
    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    fig.savefig(OUT, dpi=300, metadata={"CreationDate": None})
    print("wrote", OUT)


if __name__ == "__main__":
    main()
