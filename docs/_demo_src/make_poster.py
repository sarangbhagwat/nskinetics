# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Compose the README quickstart-demo poster from real repo assets.

Writes ``docs/source/_static/images/examples/quickstart_demo_poster.png`` -- a
branded, clickable still (logo, title, step sequence, a play badge, and the real
step-2 kinetics figure in a device card) that the README links to the live
interactive demo. The PNG is ``*.png``-gitignored like the other example
figures, so commit it with ``git add -f``.

Run with the IBO_2026 interpreter::

    python docs/_demo_src/make_poster.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

ROOT = pathlib.Path(__file__).resolve().parents[2]
IMG_DIR = ROOT / "docs" / "source" / "_static" / "images"
OUT = IMG_DIR / "examples" / "quickstart_demo_poster.png"

# brand tokens
BG      = "#0B1519"   # deep teal-black
INK     = "#EAF1F3"
MUTED   = "#93A6AD"
ACCENT  = "#5FB6D3"   # bright teal (reads on dark)
ACCENTD = "#459DB9"   # brand teal
SCREEN  = "#FFFFFF"

W, H = 1200.0, 630.0


def make():
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=200)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

    # faint accent hairline down the seam
    ax.plot([636, 636], [70, 560], color=ACCENTD, lw=1.2, alpha=0.28,
            solid_capstyle="round")

    # right: the real step-2 kinetics figure in a white "screen" card
    ax.add_patch(FancyBboxPatch((656, 58), 500, 498,
                 boxstyle="round,pad=0,rounding_size=18", linewidth=0,
                 facecolor="#000000", alpha=0.28, zorder=1))          # shadow
    ax.add_patch(FancyBboxPatch((650, 66), 500, 498,
                 boxstyle="round,pad=0,rounding_size=18", linewidth=0,
                 facecolor=SCREEN, zorder=2))                         # card
    kin = plt.imread(str(IMG_DIR / "examples" / "tutorial_01_quickstart_kinetics.png"))
    ax_img = fig.add_axes([668 / W, 86 / H, 464 / W, 458 / H])
    ax_img.axis("off"); ax_img.set_facecolor(SCREEN)
    ax_img.imshow(kin, aspect="equal", zorder=3)

    # left: logo (light-on-dark variant)
    logo = plt.imread(str(IMG_DIR / "logo" / "logo_nskinetics_dark.png"))
    disp_h = 58.0
    disp_w = disp_h * logo.shape[1] / logo.shape[0]
    ax_logo = fig.add_axes([70 / W, (H - 70 - disp_h) / H, disp_w / W, disp_h / H])
    ax_logo.axis("off"); ax_logo.imshow(logo, zorder=3)

    # left: title (serif, echoing the demo's Spectral)
    ax.text(70, 452, "Quickstart,", fontsize=41, fontweight="bold",
            family="serif", color=INK, va="baseline")
    ax.text(70, 398, "end to end.", fontsize=41, fontweight="bold",
            family="serif", color=INK, va="baseline")

    # subtitle
    ax.text(70, 340, "One factory call → a full fed-batch", fontsize=16.5,
            color=MUTED, va="baseline")
    ax.text(70, 314, "fermentation + its process economics.", fontsize=16.5,
            color=MUTED, va="baseline")

    # step sequence
    x = 70
    steps = ["Build", "Simulate", "Re-tune", "Inspect"]
    for i, s in enumerate(steps):
        ax.text(x, 250, s, fontsize=13, color=ACCENT, va="baseline",
                family="monospace")
        x += 11 * len(s) + 10
        if i < len(steps) - 1:
            ax.text(x, 250, "·", fontsize=15, color=MUTED, va="baseline",
                    family="monospace")
            x += 22

    # play badge pill
    ax.add_patch(FancyBboxPatch((70, 150), 415, 50,
                 boxstyle="round,pad=0,rounding_size=25", linewidth=0,
                 facecolor=ACCENTD, zorder=3))
    ax.add_patch(Polygon([[100, 163], [100, 187], [121, 175]], closed=True,
                 facecolor="white", zorder=4))
    ax.text(138, 175, "Watch the demo  ·  ~1 min", fontsize=15.5,
            fontweight="bold", color="white", va="center", zorder=4)

    fig.savefig(str(OUT), facecolor=BG)
    plt.close(fig)
    return OUT


if __name__ == "__main__":
    out = make()
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} KB)")
