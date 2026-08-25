# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Render the animated "strain-to-TEA" loop GIFs for the docs landing page.

The scene is a captioned three-stage pipeline drawn left to right:

* **Engineered strain** — a teal ovoid microbe with granules and an amber
  plasmid loop, above a control card of three tick-marked sliders (the
  middle, amber one is the parameter the animation turns).
* **Fermentation kinetics** — a glass stirred tank: gradient broth with a
  waving surface, a sparger streaming rising bubbles, a two-tier Rushton
  impeller spinning in side view, and an amber product curve that boosts
  when the parameter changes.
* **Process & TEA** — a plant building with vapor wisps drifting off its
  stack, a cutaway flowsheet panel that lights amber, a trayed
  distillation column, and a large semicircular $ gauge above.

Forward arrows link the stages and a dashed feedback arc curves back
underneath. Amber comet pulses travel those paths, each stage glowing as a
pulse reaches it.

One loop is 10 s at 20 fps (200 frames), rendered 2000 x 720 px. Ambient
motion (impeller, bubbles, surface wave, vapor) runs in every frame with an
integer number of cycles per loop; on top sit three acts — a pulse down the
pipeline, then a slider move that boosts the curve and swings the gauge,
then a feedback pulse right-to-left while everything eases home, masking the
reset. Frame 0 and the last frame differ no more than any adjacent pair, so
the loop wraps seamlessly.

Run with no arguments to (re)build both theme variants:

    python docs/_demo_src/make_loop_gif.py
    ->  docs/source/_static/images/demo/loop_light.gif
        docs/source/_static/images/demo/loop_dark.gif
        docs/source/_static/images/demo/loop_light_still.png
        docs/source/_static/images/demo/loop_dark_still.png

The ``*_still.png`` files are the frame-0 stills the landing page serves
instead of the GIF under ``prefers-reduced-motion: reduce``.

Use ``--stills DIR`` to render one PNG per preset state for visual review.
Every frame is quantized against one shared 256-color palette with dithering
off; pass ``--dither`` to enable Floyd-Steinberg if the gradients band.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import (Arc, Circle, Ellipse, FancyArrowPatch,
                                FancyBboxPatch, Polygon, Rectangle, Wedge)
from PIL import Image

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / 'source' / '_static' / 'images' / 'demo'

DPI = 200    # figsize (10, 3.6) -> 2000 x 720 px (displayed at 200 px height;
             # rendered large for hi-DPI screens and reuse at bigger sizes)
DUR = 10.0   # seconds per loop
FPS = 20

# 'bg' must match the docs page background (pydata-sphinx-theme): GIFs are
# opaque, so an off-theme bg shows as a visible rectangle on the page.
# Palette rule: navy + amber dominate; amber/orange is reserved for
# attention elements (pulses, animated knob, plasmid, product curve, lit
# flowsheet, $, glows); teal/red/steel are quieter supporting hues.
THEMES = {
    'light': dict(bg='#ffffff', stroke='#1f2a63', navy_mid='#3d4d8f',
                  peri='#8a93b8', faint='#8a93b8',
                  card='#eef0f7', card_edge='#c3c9de', shadow='#1f2a63',
                  steel='#7d8799', steel_dark='#5a6478',
                  accent='#f5a623', accent2='#ef7b45',
                  mb0='#58b3a4', mb1='#cdeae4', teal_edge='#1f6f64',
                  broth0='#8ed0c3', broth1='#d3efe9', bubble='#ffffff',
                  red='#d64545', text='#1f2a63', glass_hi='#ffffff',
                  glass0='#e7eaf4', glass1='#f8f9fc',
                  bldg0='#1f2a63', bldg1='#3d4d8f'),
    'dark':  dict(bg='#14181e', stroke='#d7dce8', navy_mid='#38466e',
                  peri='#5a6478', faint='#5a6478',
                  card='#202836', card_edge='#3a4459', shadow='#000000',
                  steel='#8a94a8', steel_dark='#67718a',
                  accent='#f5a623', accent2='#ef7b45',
                  mb0='#1f5b52', mb1='#47a08f', teal_edge='#63c9b8',
                  broth0='#17453e', broth1='#2a7d70', bubble='#8fd8ca',
                  red='#e05c5c', text='#d7dce8', glass_hi='#c3cde8',
                  glass0='#1b222d', glass1='#273043',
                  bldg0='#232d47', bldg1='#3a4a75'),
}

# %% Scene geometry

# Stage centers: microbe/card ~x 1.35, reactor x 5.0, plant ~x 8.3.
# Captions sit at y 0.15; the feedback arc dips to ~y 0.45 mid-canvas,
# clearing the reactor caption below it and the control card above it.
ARROW1 = ((2.32, 2.00), (4.24, 2.00))
ARROW2 = ((5.76, 2.00), (7.00, 2.00))

# Feedback arc: quadratic bezier (right to left) below the pipeline.
# FB_RAD reproduces the same control point via arc3 (offset/chord = -1.0/7.0).
FB_P0, FB_P1, FB_P2 = (8.45, 0.95), (4.95, -0.05), (1.45, 0.80)
FB_RAD = -0.143

CAPTION_Y = 0.15

MICROBE_C = (1.35, 2.55)
CARD = (0.50, 0.88, 1.70, 1.06)          # x, y, w, h
SLIDER_X0, SLIDER_X1 = 0.98, 2.02
SLIDER_YS = (1.68, 1.41, 1.14)
SLIDER_LABELS = (r'$\mu$', r'$K$', r'$Y$')

VESSEL = (4.42, 1.15, 1.16, 1.70)   # x, y, w, h; rounded corners 0.20
LIQ_Y = 2.28                        # nominal broth surface
SPIN_REVS = 8                       # impeller revolutions per loop (integer)
WAVE_CYCLES = 4                     # surface-wave cycles per loop (integer)
# (x, rise cycles per loop, phase offset) per bubble; integer cycles keep
# the loop seamless.
BUBBLES = ((4.62, 2, 0.05), (4.78, 3, 0.55), (4.93, 2, 0.35),
           (5.07, 3, 0.80), (5.22, 2, 0.62), (5.38, 3, 0.15))

PLANT_BLDG = (7.18, 1.05, 1.37, 1.25)    # x, y, w, h
STACK = (7.60, 2.30, 0.20, 0.60)
PANEL = (7.30, 1.28, 1.06, 0.78)         # flowsheet cutaway card
FLOW_BOXES = ((7.38, 1.72, 0.24, 0.20), (7.74, 1.84, 0.24, 0.20),
              (7.74, 1.42, 0.24, 0.20), (8.06, 1.60, 0.24, 0.20))
FLOW_LINES = (((7.62, 1.84), (7.74, 1.92)), ((7.62, 1.80), (7.74, 1.54)),
              ((7.98, 1.92), (8.10, 1.76)), ((7.98, 1.52), (8.10, 1.64)))
COLUMN = (8.76, 1.05, 0.44, 1.85)
VAPOR_CYCLES = 2                          # integer cycles per loop
GAUGE_C = (8.28, 3.03)
GAUGE_R = 0.50                            # ~1.5x the old 0.34


def default_state():
    """Scene state for one frame; every animated quantity in one dict.

    't' is the absolute loop time in seconds — ambient motion (impeller,
    bubbles, wave, vapor) derives from it so it runs in every frame.
    """
    return dict(t=0.0,
                slider_frac=0.35,   # animated (middle) slider knob, 0..1
                needle_frac=0.50,   # gauge needle, 0 (far left) .. 1 (far right)
                curve_boost=0.0,    # 0 = base product curve, 1 = boosted
                pulse1_s=None,      # forward comet path params, None = hidden
                pulse2_s=None,
                fb_s=None,          # feedback comet path param, None = hidden
                glow_microbe=0.0, glow_reactor=0.0, glow_plant=0.0)


# %% Drawing helpers

def gradient_fill(ax, clip_patch, extent, c0, c1, direction='v',
                  alpha=1.0, zorder=2):
    """Linear-gradient fill clipped to an already-added patch.

    ``clip_patch`` must already be on ``ax`` (an invisible clip shape or the
    visible outline itself). ``extent`` = (x0, y0, x1, y1) in data coords;
    ``c0`` maps to the bottom (or left, for direction='h'), ``c1`` to the
    top (or right). aspect='auto' is safe: the axes fills the figure and the
    data aspect already equals the figure aspect (10 : 3.6).
    """
    grad = np.linspace(0.0, 1.0, 256)
    grad = grad[:, None] if direction == 'v' else grad[None, :]
    cmap = LinearSegmentedColormap.from_list('_grad', [c0, c1])
    x0, y0, x1, y1 = extent
    im = ax.imshow(grad, extent=(x0, x1, y0, y1), origin='lower', cmap=cmap,
                   aspect='auto', alpha=alpha, zorder=zorder,
                   interpolation='bilinear')
    im.set_clip_path(clip_patch)
    return im


def soft_shadow(ax, xy, rx, ry, th, alpha=0.18, zorder=1):
    """Soft elliptical drop shadow (stacked low-alpha ellipses)."""
    for k, a in ((1.30, 0.35), (1.12, 0.65), (1.0, 1.0)):
        ax.add_patch(Ellipse(xy, 2*rx*k, 2*ry*k, fc=th['shadow'], ec='none',
                             alpha=alpha*a*0.4, zorder=zorder))


def draw_caption(ax, th, x, text):
    ax.text(x, CAPTION_Y, text, ha='center', va='center', fontsize=11.5,
            color=th['text'], alpha=0.85, zorder=3)


# %% Drawing

def _halo(ax, th, xy, rx, ry, strength):
    """Soft amber glow behind a stage while a pulse passes it."""
    if strength <= 0:
        return
    for k, a in ((1.35, 0.10), (1.0, 0.16)):
        ax.add_patch(Ellipse(xy, 2*rx*k, 2*ry*k, fc=th['accent'], ec='none',
                             alpha=a*min(strength, 1.0), zorder=1))


def draw_microbe(ax, th, glow):
    """Engineered ovoid microbe: teal gradient body, granules, amber plasmid."""
    cx, cy = MICROBE_C
    _halo(ax, th, MICROBE_C, 0.66, 0.48, glow)
    soft_shadow(ax, (cx, cy - 0.38), 0.42, 0.07, th, alpha=0.14)
    body_clip = Ellipse(MICROBE_C, 1.04, 0.66, angle=-14, fc='none',
                        ec='none', zorder=2)
    ax.add_patch(body_clip)
    gradient_fill(ax, body_clip, (cx - 0.6, cy - 0.4, cx + 0.6, cy + 0.4),
                  th['mb0'], th['mb1'], zorder=2.5)
    ax.add_patch(Ellipse(MICROBE_C, 1.04, 0.66, angle=-14, fc='none',
                         ec=th['teal_edge'], lw=3, zorder=3.5))
    # membrane highlight along the upper-left rim
    ax.add_patch(Arc(MICROBE_C, 0.84, 0.46, angle=-14, theta1=95, theta2=175,
                     ec=th['glass_hi'], lw=2.5, alpha=0.55, zorder=3))
    # granules (quiet navy)
    for gx, gy, r in ((1.10, 2.62, 0.05), (1.24, 2.42, 0.045),
                      (1.16, 2.52, 0.028)):
        ax.add_patch(Circle((gx, gy), r, fc=th['stroke'], ec='none',
                            alpha=0.85, zorder=3))
    # engineering cue: amber plasmid loop (nested loops = supercoil)
    ax.add_patch(Circle((1.58, 2.52), 0.10, fc='none', ec=th['accent'],
                        lw=2.5, zorder=3))
    ax.add_patch(Circle((1.58, 2.52), 0.045, fc='none', ec=th['accent'],
                        lw=1.5, alpha=0.7, zorder=3))


def draw_control_card(ax, th, frac):
    """Rounded control card with three tick-marked, labeled sliders.

    The middle (amber) knob is the animated one; its track fills amber up
    to the knob so the change reads at a glance.
    """
    x, y, w, h = CARD
    soft_shadow(ax, (x + w/2, y + 0.02), w*0.52, 0.09, th)
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle='round,pad=0,rounding_size=0.12',
                                fc=th['card'], ec=th['card_edge'], lw=1.5,
                                zorder=2))
    knob_fracs = (0.30, frac, 0.72)
    for yy, lab, kf, animated in zip(SLIDER_YS, SLIDER_LABELS, knob_fracs,
                                     (False, True, False)):
        ax.text(0.76, yy, lab, ha='center', va='center', fontsize=10,
                color=th['text'], zorder=3)
        ax.plot([SLIDER_X0, SLIDER_X1], [yy, yy], color=th['peri'], lw=3.5,
                solid_capstyle='round', alpha=0.75, zorder=3)
        for tk in np.linspace(SLIDER_X0, SLIDER_X1, 5):
            ax.plot([tk, tk], [yy - 0.035, yy + 0.035], color=th['peri'],
                    lw=1.2, alpha=0.6, zorder=3)
        kx = SLIDER_X0 + kf*(SLIDER_X1 - SLIDER_X0)
        if animated:
            ax.plot([SLIDER_X0, kx], [yy, yy], color=th['accent'], lw=3.5,
                    solid_capstyle='round', alpha=0.9, zorder=3.5)
        color = th['accent'] if animated else th['stroke']
        ax.add_patch(Circle((kx, yy), 0.085, fc=color, ec=th['card'],
                            lw=1.8, zorder=4))


def draw_arrow(ax, th, p0, p1):
    """Tapered filled forward arrow (steel, quiet under the amber comets)."""
    ax.add_patch(FancyArrowPatch(
        p0, p1,
        arrowstyle='simple,tail_width=0.18,head_width=0.85,head_length=1.2',
        mutation_scale=16, fc=th['steel'], ec='none', alpha=0.9, zorder=2))


def draw_feedback_arrow(ax, th):
    """Dashed 'iterate' bezier right-to-left with a solid head at the end."""
    ss = np.linspace(0.0, 0.96, 60)
    xy = np.array([fb_xy(s) for s in ss])
    ax.plot(xy[:, 0], xy[:, 1], color=th['steel'], lw=2.2,
            ls=(0, (5, 4)), alpha=0.9, zorder=2)
    ax.add_patch(FancyArrowPatch(fb_xy(0.93), fb_xy(1.0), arrowstyle='-|>',
                                 mutation_scale=20, lw=0, fc=th['steel'],
                                 ec=th['steel'], zorder=2))


def _vessel_patch(fc='none', ec='none', lw=0, zorder=2):
    x, y, w, h = VESSEL
    return FancyBboxPatch((x, y), w, h,
                          boxstyle='round,pad=0,rounding_size=0.20',
                          fc=fc, ec=ec, lw=lw, zorder=zorder)


def _liquid_poly(t):
    """Broth silhouette: chamfered bottom, wavy top (WAVE_CYCLES per loop)."""
    x, y, w, _ = VESSEL
    ph = 2*np.pi*WAVE_CYCLES*t/DUR
    xs = np.linspace(x + 0.04, x + w - 0.04, 40)
    top = [(xi, LIQ_Y + 0.02*np.sin(2*np.pi*xi*2.5 - ph)) for xi in xs]
    c = 0.13    # bottom-corner chamfer, approximating the vessel rounding
    bottom = [(x + w - 0.04, y + 0.04 + c), (x + w - 0.04 - c, y + 0.04),
              (x + 0.04 + c, y + 0.04), (x + 0.04, y + 0.04 + c)]
    return Polygon(top + bottom, closed=True, fc='none', ec='none')


def draw_reactor(ax, th, t, glow, curve_boost):
    x, y, w, h = VESSEL
    cx, top = x + w/2, y + h
    _halo(ax, th, (cx, y + h/2), 0.82, 1.02, glow)
    soft_shadow(ax, (cx, y - 0.03), w*0.55, 0.08, th)
    # glass wall: subtle horizontal gradient clipped to the vessel silhouette
    wall = _vessel_patch(zorder=2)
    ax.add_patch(wall)
    gradient_fill(ax, wall, (x, y, x + w, y + h), th['glass0'], th['glass1'],
                  direction='h', zorder=2)
    # broth: teal gradient clipped to the wavy liquid polygon
    liq = _liquid_poly(t)
    ax.add_patch(liq)
    gradient_fill(ax, liq, (x, y, x + w, LIQ_Y + 0.04),
                  th['broth0'], th['broth1'], zorder=2.4)
    xs = np.linspace(x + 0.04, x + w - 0.04, 40)
    ph = 2*np.pi*WAVE_CYCLES*t/DUR
    ax.plot(xs, LIQ_Y + 0.02*np.sin(2*np.pi*xs*2.5 - ph),
            color=th['teal_edge'], lw=2, alpha=0.7, zorder=2.9)
    # sparger + bubbles
    ax.plot([cx - 0.30, cx + 0.30], [y + 0.10, y + 0.10],
            color=th['steel_dark'], lw=2, solid_capstyle='round', zorder=2.6)
    for bx, cycles, phb in BUBBLES:
        u = (cycles*t/DUR + phb) % 1.0
        by = y + 0.14 + u*(LIQ_Y - 0.24 - (y + 0.14))
        wob = 0.028*np.sin(2*np.pi*(3*u + phb))
        a = 0.55*(1.0 - smooth(0.85, 1.0, u))
        ax.add_patch(Circle((bx + wob, by), 0.018 + 0.030*u, fc='none',
                            ec=th['bubble'], lw=1.4, alpha=a, zorder=2.8))
    # impeller: motor block, shaft, two spinning Rushton tiers. Side-view
    # rotation: two blade pairs 90 deg out of phase; each pair's apparent
    # half-span is L|cos(theta)| (front pair dark, back pair lighter).
    ax.add_patch(FancyBboxPatch((cx - 0.14, top - 0.02), 0.28, 0.24,
                                boxstyle='round,pad=0,rounding_size=0.05',
                                fc=th['navy_mid'], ec=th['stroke'], lw=1.5,
                                zorder=4))
    ax.plot([cx, cx], [top, 1.50], color=th['stroke'], lw=2.2, zorder=3)
    theta = 2*np.pi*SPIN_REVS*t/DUR
    for ty in (1.95, 1.55):
        for phase, front in ((theta, True), (theta + np.pi/2, False)):
            half = 0.26*abs(np.cos(phase))
            if half < 0.02:
                continue
            color = th['stroke'] if front else th['navy_mid']
            a = 0.95 if front else 0.6
            ax.plot([cx - half, cx + half], [ty, ty], color=color, lw=5,
                    solid_capstyle='round', alpha=a, zorder=3.1)
            for sx in (cx - half, cx + half):
                ax.plot([sx, sx], [ty - 0.055, ty + 0.055], color=color,
                        lw=2.2, alpha=a, zorder=3.1)
    # product curve (attention amber) over the broth, with a soft under-glow
    cxs = np.linspace(x + 0.18, x + w - 0.18, 40)
    amp = 0.30 + 0.28*curve_boost
    cys = 1.52 + amp/(1 + np.exp(-(cxs - 4.95)*9))
    ax.plot(cxs, cys, color=th['accent'], lw=7, alpha=0.25,
            solid_capstyle='round', zorder=3.4)
    ax.plot(cxs, cys, color=th['accent'], lw=2.8, solid_capstyle='round',
            zorder=3.5)
    # vessel outline + glass sheen on top of everything inside
    ax.add_patch(_vessel_patch(ec=th['stroke'], lw=3, zorder=4))
    ax.plot([x + 0.16, x + 0.16], [y + 0.28, y + h - 0.38],
            color=th['glass_hi'], lw=5, alpha=0.30, solid_capstyle='round',
            zorder=4.2)


def draw_plant(ax, th, t, glow, flow_lit):
    _halo(ax, th, (8.30, 1.95), 1.12, 1.05, glow)
    x, y, w, h = PLANT_BLDG
    soft_shadow(ax, (x + w/2 + 0.5, y - 0.02), 1.15, 0.08, th)
    # building + stack: filled navy gradient with a stroke outline
    for (bx, by, bw, bh) in (PLANT_BLDG, STACK):
        p = FancyBboxPatch((bx, by), bw, bh,
                           boxstyle='round,pad=0,rounding_size=0.05',
                           fc='none', ec='none', zorder=2)
        ax.add_patch(p)
        gradient_fill(ax, p, (bx, by, bx + bw, by + bh),
                      th['bldg0'], th['bldg1'], zorder=2)
        ax.add_patch(FancyBboxPatch((bx, by), bw, bh,
                                    boxstyle='round,pad=0,rounding_size=0.05',
                                    fc='none', ec=th['stroke'], lw=2.5,
                                    zorder=3))
    # vapor wisps drifting up-right from the stack (ambient, periodic)
    sx, sy = STACK[0] + STACK[2]/2, STACK[1] + STACK[3]
    for i in range(3):
        u = (VAPOR_CYCLES*t/DUR + i/3.0) % 1.0
        ax.add_patch(Circle((sx + 0.10*u + 0.03*np.sin(2*np.pi*u),
                             sy + 0.05 + 0.34*u),
                            0.05 + 0.08*u, fc=th['peri'], ec='none',
                            alpha=0.35*(1.0 - u), zorder=2.5))
    # flowsheet cutaway panel; base pass in steel, lit pass in amber
    px, py, pw, phh = PANEL
    ax.add_patch(FancyBboxPatch((px, py), pw, phh,
                                boxstyle='round,pad=0,rounding_size=0.08',
                                fc=th['card'], ec=th['card_edge'], lw=1.2,
                                alpha=0.95, zorder=3.5))
    for color, alpha, lw in ((th['steel'], 1.0, 1.8),
                             (th['accent'], min(flow_lit, 1.0), 2.4)):
        if alpha <= 0:
            continue
        for bx, by, bw, bh in FLOW_BOXES:
            ax.add_patch(Rectangle((bx, by), bw, bh, fc='none', ec=color,
                                   lw=lw, alpha=alpha, zorder=4))
        for (x0, y0), (x1, y1) in FLOW_LINES:
            ax.plot([x0, x1], [y0, y1], color=color, lw=lw, alpha=alpha,
                    zorder=4)
    # distillation column: glass gradient, tray lines, connecting pipe
    colx, coly, colw, colh = COLUMN
    col = FancyBboxPatch((colx, coly), colw, colh,
                         boxstyle='round,pad=0,rounding_size=0.20',
                         fc='none', ec='none', zorder=2)
    ax.add_patch(col)
    gradient_fill(ax, col, (colx, coly, colx + colw, coly + colh),
                  th['glass0'], th['glass1'], direction='h', zorder=2)
    ax.add_patch(FancyBboxPatch((colx, coly), colw, colh,
                                boxstyle='round,pad=0,rounding_size=0.20',
                                fc='none', ec=th['stroke'], lw=3, zorder=3))
    for ty in np.linspace(1.45, 2.55, 4):
        ax.plot([colx + 0.06, colx + colw - 0.06], [ty, ty],
                color=th['steel'], lw=1.6, alpha=0.8, zorder=3)
    ax.plot([x + w, colx], [2.02, 2.02], color=th['steel_dark'], lw=2.5,
            zorder=2.6)
    ax.plot([(x + w + colx)/2]*2, [1.95, 2.09], color=th['steel_dark'],
            lw=2, zorder=2.6)   # flange tick
    # ground line
    ax.plot([7.05, 9.50], [1.05, 1.05], color=th['stroke'], lw=3,
            solid_capstyle='round', zorder=3)


def draw_gauge(ax, th, needle_frac):
    """Semicircular $-gauge: filled face, ticks, warm-red needle, amber $."""
    cx, cy = GAUGE_C
    ax.add_patch(Wedge((cx, cy), GAUGE_R, 0, 180, fc=th['card'], ec='none',
                       zorder=3))
    ax.add_patch(Arc((cx, cy), 2*GAUGE_R, 2*GAUGE_R, theta1=0, theta2=180,
                     ec=th['stroke'], lw=3.5, zorder=4))
    ax.plot([cx - GAUGE_R, cx + GAUGE_R], [cy, cy], color=th['stroke'],
            lw=3.5, solid_capstyle='round', zorder=4)
    for ang in (150, 120, 90, 60, 30):
        a = np.deg2rad(ang)
        ax.plot([cx + 0.38*np.cos(a), cx + 0.46*np.cos(a)],
                [cy + 0.38*np.sin(a), cy + 0.46*np.sin(a)],
                color=th['steel'], lw=2, zorder=4)
    a = np.deg2rad(160 - 140*needle_frac)
    ax.plot([cx, cx + 0.36*np.cos(a)], [cy, cy + 0.36*np.sin(a)],
            color=th['red'], lw=3, solid_capstyle='round', zorder=4.5)
    ax.add_patch(Circle((cx, cy), 0.06, fc=th['stroke'], ec='none',
                        zorder=5))
    ax.text(cx, cy - 0.17, '$', ha='center', va='center', fontsize=15,
            fontweight='bold', color=th['accent'], zorder=4)


def pulse_xy(s):
    """Forward-comet position at path parameter s in [0, 1].

    Returns None while the comet is "inside" a stage (the stage glows
    instead): s 0.34-0.50 crosses the reactor; s > 0.80 is inside the plant.
    """
    if s < 0.34:
        u = s/0.34
        (x0, y0), (x1, y1) = ARROW1
    elif s < 0.50:
        return None
    elif s < 0.80:
        u = (s - 0.50)/0.30
        (x0, y0), (x1, y1) = ARROW2
    else:
        return None
    return (x0 + u*(x1 - x0), y0 + u*(y1 - y0))


def fb_xy(s):
    """Quadratic-bezier position of the feedback pulse (right to left)."""
    (x0, y0), (x1, y1), (x2, y2) = FB_P0, FB_P1, FB_P2
    a, b = (1 - s)**2, 2*s*(1 - s)
    return (a*x0 + b*x1 + s*s*x2, a*y0 + b*y1 + s*s*y2)


def draw_comet(ax, th, path_fn, s):
    """Amber comet: glowing head + fading trail of diminishing circles."""
    for i in range(9):
        si = s - 0.022*i
        if si < 0:
            break
        xy = path_fn(si)
        if xy is None:
            continue
        if i == 0:
            ax.add_patch(Circle(xy, 0.16, fc=th['accent'], ec='none',
                                alpha=0.30, zorder=6))
            ax.add_patch(Circle(xy, 0.085, fc=th['accent'], ec='none',
                                zorder=6))
        else:
            k = 1.0 - i/9.0
            ax.add_patch(Circle(xy, 0.075*k, fc=th['accent'], ec='none',
                                alpha=0.45*k, zorder=6))


def draw_scene(ax, th, state):
    t = state['t']
    draw_caption(ax, th, 1.35, 'Engineered strain')
    draw_caption(ax, th, 5.0, 'Fermentation kinetics')
    draw_caption(ax, th, 8.35, 'Process & TEA')
    draw_feedback_arrow(ax, th)
    draw_arrow(ax, th, *ARROW1)
    draw_arrow(ax, th, *ARROW2)
    draw_control_card(ax, th, state['slider_frac'])
    draw_microbe(ax, th, state['glow_microbe'])
    draw_reactor(ax, th, t, state['glow_reactor'], state['curve_boost'])
    draw_plant(ax, th, t, state['glow_plant'], min(state['glow_plant'], 1.0))
    draw_gauge(ax, th, state['needle_frac'])
    for s in (state['pulse1_s'], state['pulse2_s']):
        if s is not None:
            draw_comet(ax, th, pulse_xy, s)
    if state['fb_s'] is not None:
        draw_comet(ax, th, fb_xy, state['fb_s'])


def render_frame(th, state):
    """Render one frame to an opaque RGB PIL image (2000 x 720)."""
    fig = plt.figure(figsize=(10, 3.6), dpi=DPI)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.6)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor(th['bg'])
    draw_scene(ax, th, state)
    # imshow (gradient fills) can autoscale the data limits; re-assert.
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.6)
    fig.canvas.draw()
    img = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
    plt.close(fig)
    return img


# %% Animation timeline
#
# One 10 s loop. Ambient motion (impeller spin, bubbles, surface wave, vapor
# wisps) derives from state['t'] and runs in every frame, each with an
# integer number of cycles per loop so the wrap is seamless. On top, three
# acts:
#   Act 1 (~0-3.6 s)   forward comet microbe -> reactor -> plant; each stage
#                      glows as it passes; the needle settles mid-scale.
#   Act 2 (~3.9-7.8 s) the amber slider moves; a second comet propagates;
#                      the product curve boosts; the needle swings to a
#                      lower cost with a small overshoot-and-settle.
#   Act 3 (~8.0-9.8 s) the feedback comet runs the dashed arc right-to-left
#                      while slider, curve, and needle ease home, masking
#                      the reset.


def smooth(a, b, t):
    """Smoothstep from 0 (t <= a) to 1 (t >= b)."""
    if t <= a:
        return 0.0
    if t >= b:
        return 1.0
    x = (t - a)/(b - a)
    return x*x*(3.0 - 2.0*x)


def bump(a, b, t):
    """0 -> 1 -> 0 over [a, b] (glow envelope)."""
    m = 0.5*(a + b)
    return smooth(a, m, t)*(1.0 - smooth(m, b, t))


def timeline(t):
    """Scene state at time t (seconds) within the DUR-second loop."""
    s = default_state()
    s['t'] = t
    if 0.3 <= t < 3.5:
        s['pulse1_s'] = (t - 0.3)/3.2
    if 4.6 <= t < 7.4:
        s['pulse2_s'] = (t - 4.6)/2.8
    if 8.0 <= t < 9.6:
        s['fb_s'] = (t - 8.0)/1.6
    # Glow envelopes are timed to the comets' stage arrivals (path params
    # 0.34/0.50 = reactor, 0.80 = plant).
    s['glow_microbe'] = bump(0.1, 1.2, t) + bump(3.9, 4.9, t)
    s['glow_reactor'] = bump(1.2, 2.2, t) + bump(5.5, 6.4, t)
    s['glow_plant'] = bump(2.7, 3.8, t) + bump(6.6, 7.7, t)
    s['slider_frac'] = (0.35 + 0.40*smooth(3.9, 4.6, t)
                        - 0.40*smooth(8.6, 9.6, t))
    s['curve_boost'] = smooth(5.6, 6.3, t) - smooth(8.6, 9.6, t)
    s['needle_frac'] = (0.50 + 0.20*smooth(3.0, 3.6, t)
                        - 0.46*smooth(7.0, 7.6, t)
                        + 0.06*bump(7.5, 8.4, t)      # overshoot-settle
                        + 0.26*smooth(8.8, 9.8, t))   # ease home
    return s


# Theme colors forced into the palette exactly, in priority order. Both are
# large, flat, low-saturation fills, which is precisely where a nearest-match
# miss is most visible: the eye reads a whole uniform panel as the wrong hue.
PINNED_KEYS = ('bg', 'card')


def _pin_exact(pal, rgb, taken):
    """Overwrite the nearest not-yet-pinned palette slot with ``rgb`` exactly.

    Returns the slot index. Editing the *nearest* slot keeps the disturbance
    minimal — that slot moves by at most its own quantization error — and
    guarantees the pinned color quantizes to itself afterwards (distance 0).
    """
    entries = np.asarray(pal, int).reshape(-1, 3)
    d = ((entries - np.asarray(rgb, int))**2).sum(axis=1).astype(float)
    d[list(taken)] = np.inf
    idx = int(np.argmin(d))
    pal[3*idx:3*idx + 3] = list(rgb)
    return idx


def _shared_palette(frames, th):
    """Build the one palette every frame is quantized against.

    Derived from a probe of *many* frames, not frame 0: at t = 0 no comet,
    glow or boosted curve is on screen, so a frame-0 palette has no entry
    for those colors and nearest-match snaps them to wrong hues.
    MAXCOVERAGE (not MEDIANCUT) because the canvas is mostly background —
    median cut spends its budget subdividing that near-uniform mass and
    starves the sparse accent colors that carry the animation.

    At 2000 x 720 x 200 frames a full-resolution probe is ~0.9 GB of RGB,
    so the probe takes every 2nd frame at half resolution — still covering
    every act (comets, glows, boosted curve, feedback) and every gradient.

    The theme colors in ``PINNED_KEYS`` are then written into the palette
    exactly, because on a big flat fill a nearest-match miss is not a slight
    tint but a wrong color. The GIF is opaque and matted on the docs page, so
    an off-'bg' shows as a visibly bounded off-white (or off-dark) rectangle;
    and 'card' (#eef0f7) otherwise lands on a neighbour whose green exceeds
    its red and blue, flipping the whole control-card face from blue-tinted
    lavender to mint. Each pin takes its own slot so a later pin cannot
    overwrite an earlier one.
    """
    sub = [f.resize((f.width//2, f.height//2), Image.Resampling.BILINEAR)
           for f in frames[::2]]
    w, h = sub[0].size
    probe = Image.new('RGB', (w, h*len(sub)))
    for k, f in enumerate(sub):
        probe.paste(f, (0, k*h))
    base = probe.quantize(colors=256, method=Image.MAXCOVERAGE)
    pal = list(base.getpalette())
    taken = set()
    for key in PINNED_KEYS:
        rgb = tuple(int(th[key][i:i + 2], 16) for i in (1, 3, 5))
        taken.add(_pin_exact(pal, rgb, taken))
    out = Image.new('P', (1, 1))
    out.putpalette(pal)
    out.load()   # refresh the core palette so quantize() sees the edit
    return out


def build_gif(theme_name, out_path, dither=Image.Dither.NONE):
    """Render all frames for one theme and write an infinite-loop GIF."""
    th = THEMES[theme_name]
    n = int(round(DUR*FPS))
    frames = [render_frame(th, timeline(i/FPS)) for i in range(n)]
    # Quantize every frame against one shared palette. Dither stays off by
    # default (flat regions stay clean); pass --dither only if the review
    # stills show visible banding in the gradients.
    base = _shared_palette(frames, th)
    pal_frames = [f.quantize(palette=base, dither=dither) for f in frames]
    pal_frames[0].save(out_path, save_all=True, append_images=pal_frames[1:],
                       duration=int(1000/FPS), loop=0, optimize=True)
    print(f'{out_path} ({out_path.stat().st_size/1024:.0f} KiB, '
          f'{n} rendered frames)')


# %% CLI

# Hand-picked states for reviewing the scene art without the animation.
_PRESETS = {
    'base': {},
    'spin-quarter': dict(t=DUR/32),   # impeller quarter-turn vs 'base'
    'pulse': dict(t=1.0, pulse1_s=0.22, glow_microbe=0.5),
    'reactor-glow': dict(t=1.6, glow_reactor=1.0),
    'shifted': dict(t=6.5, slider_frac=0.75, curve_boost=1.0,
                    needle_frac=0.24, glow_plant=1.0),
    'feedback': dict(t=8.8, fb_s=0.5),
}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Build the landing-page loop GIFs (or review stills).')
    ap.add_argument('--stills', metavar='DIR',
                    help='render one PNG per preset state into DIR and exit')
    ap.add_argument('--dither', action='store_true',
                    help='Floyd-Steinberg dithering (only if gradients band)')
    args = ap.parse_args(argv)
    if args.stills:
        out = Path(args.stills)
        out.mkdir(parents=True, exist_ok=True)
        for theme_name, th in THEMES.items():
            for name, over in _PRESETS.items():
                render_frame(th, default_state() | over).save(
                    out / f'{name}_{theme_name}.png')
        print(f'wrote {len(THEMES)*len(_PRESETS)} stills to {out}')
        return
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dither = (Image.Dither.FLOYDSTEINBERG if args.dither
              else Image.Dither.NONE)
    for theme_name in THEMES:
        build_gif(theme_name, OUT_DIR / f'loop_{theme_name}.gif',
                  dither=dither)
    # Frame-0 stills: what the landing page shows in place of the GIF under
    # prefers-reduced-motion (see the <picture> blocks in docs/source/index.rst).
    for theme_name, th in THEMES.items():
        still = OUT_DIR / f'loop_{theme_name}_still.png'
        render_frame(th, timeline(0.0)).save(still)
        print(f'{still} ({still.stat().st_size/1024:.0f} KiB, frame 0)')


if __name__ == '__main__':
    main()
