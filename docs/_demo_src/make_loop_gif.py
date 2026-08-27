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

* **Metabolic engineering** — a teal ovoid microbe with granules and an amber
  plasmid loop, above a control card of three tick-marked sliders (the
  middle, amber one is the parameter the animation turns).
* **Kinetics & reactor design** — a glass stirred tank: gradient broth with a
  waving surface, a sparger streaming rising bubbles, a two-tier Rushton
  impeller spinning in side view, and an amber product curve that boosts
  when the parameter changes.
* **Facility-scale economics** — a parapet-capped plant building venting
  billowing vapor plumes off its stack and roof vent, a cutaway flowsheet
  panel (pump, exchanger, a mini of the reactor stage with its own teal
  broth and spinning impeller, column, reflux drum, product tank, and a
  recycle back to the feed) that lights amber, and a skirted
  tray column with a domed head, a breathing sump and vapor rising between
  its trays, piped to an overhead condenser drum. Slugs of material ride
  every pipe run and every flowsheet stream. A large semicircular $ gauge
  sits above.

Forward arrows link the stages and a dashed feedback arc curves back
underneath. Amber comet pulses travel those paths, each stage glowing as a
pulse reaches it.

One loop is 10 s at 20 fps (200 frames), rendered 2000 x 720 px. Ambient
motion (impeller, bubbles, surface wave, vapor plumes, pipe and flowsheet
slugs, column vapor, sump/drum levels) runs in every frame with an integer
number of cycles per loop; on top, one iteration plays out: three amber
comets circulate the loop, each departing a stage at the instant that
stage's one widget begins to move — the microbe's slider sends a comet to
the reactor, the reactor's product curve sends one to the plant, and the
plant's $ gauge sends the right-to-left feedback comet back to the microbe.
Each stage lights as its incoming comet arrives, and all three widgets then
ease home together under the feedback comet, masking the reset. Frame 0 and
the last frame differ no more than any adjacent pair, so the loop wraps
seamlessly.

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
# Every value is a color except 'sheen_a', the alpha of the vessel's glass
# highlight: on the dark theme the same 0.30 that reads as a subtle sheen on
# white glass reads as an opaque grey rod, so dark dials it down.
# 'vapor' is 'peri' on light, where a periwinkle plume already separates
# cleanly from white; on dark, 'peri' sits close enough to 'bg' that the
# plume dissolves at the 556 px delivery size, so dark takes a lighter
# periwinkle. It only ever paints small, soft, translucent puffs — no large
# flat region depends on it quantizing exactly.
THEMES = {
    'light': dict(bg='#ffffff', stroke='#1f2a63', navy_mid='#3d4d8f',
                  peri='#8a93b8', vapor='#8a93b8',
                  card='#eef0f7', card_edge='#c3c9de', shadow='#1f2a63',
                  steel='#7d8799', steel_dark='#5a6478',
                  accent='#f5a623',
                  mb0='#58b3a4', mb1='#cdeae4', teal_edge='#1f6f64',
                  broth0='#8ed0c3', broth1='#d3efe9', bubble='#479a8b',
                  red='#d64545', text='#1f2a63', glass_hi='#ffffff',
                  sheen_a=0.30,
                  glass0='#e7eaf4', glass1='#f8f9fc',
                  bldg0='#1f2a63', bldg1='#3d4d8f'),
    'dark':  dict(bg='#14181e', stroke='#d7dce8', navy_mid='#38466e',
                  peri='#5a6478', vapor='#98a3bb',
                  card='#202836', card_edge='#3a4459', shadow='#000000',
                  steel='#8a94a8', steel_dark='#67718a',
                  accent='#f5a623',
                  mb0='#1f5b52', mb1='#47a08f', teal_edge='#63c9b8',
                  broth0='#17453e', broth1='#2a7d70', bubble='#8fd8ca',
                  red='#e05c5c', text='#d7dce8', glass_hi='#c3cde8',
                  sheen_a=0.16,
                  glass0='#1b222d', glass1='#273043',
                  bldg0='#232d47', bldg1='#3a4a75'),
}

# %% Scene geometry

# Stage centers: microbe/card ~x 1.97, reactor x 5.0, plant ~x 8.3.
# The three stages are evenly spaced: the strain stage spans x 1.12-2.82,
# the vessel 4.42-5.58 and the plant block 7.18-, so both inter-stage gaps
# are 1.60 and the two forward arrows are the same length (1.24), each
# inset 0.18 from the stage edges it links.
# Captions sit at y 0.15; the feedback arc dips to ~y 0.45 mid-canvas,
# clearing the reactor caption below it and the control card above it.
ARROW1 = ((3.00, 2.00), (4.24, 2.00))
ARROW2 = ((5.76, 2.00), (7.00, 2.00))

# Feedback arc: quadratic bezier (right to left) below the pipeline.
FB_P0, FB_P1, FB_P2 = (8.45, 0.95), (5.26, -0.05), (2.07, 0.80)

CAPTION_Y = 0.15

MICROBE_C = (1.97, 2.55)
CARD = (1.12, 0.88, 1.70, 1.06)          # x, y, w, h
SLIDER_X0, SLIDER_X1 = 1.60, 2.64
SLIDER_LABEL_X = 1.34
SLIDER_YS = (1.68, 1.41, 1.14)
SLIDER_LABELS = (r'$\mu$', r'$k_{\mathrm{cat}}$', r'$K_{\mathrm{M}}$')

VESSEL = (4.42, 1.15, 1.16, 1.70)   # x, y, w, h; rounded corners 0.20
LIQ_Y = 2.28                        # nominal broth surface
SPIN_REVS = 8                       # impeller revolutions per loop (integer)
WAVE_CYCLES = 4                     # surface-wave cycles per loop (integer)
# (x, rise cycles per loop, phase offset) per bubble; integer cycles keep
# the loop seamless.
BUBBLES = ((4.62, 2, 0.05), (4.78, 3, 0.55), (4.93, 2, 0.35),
           (5.07, 3, 0.80), (5.22, 2, 0.62), (5.38, 3, 0.15))

PLANT_BLDG = (7.18, 1.05, 1.37, 1.25)    # x, y, w, h
PARAPET = (7.12, 2.26, 1.49, 0.11)       # roof slab capping the building
STACK = (7.60, 2.30, 0.20, 0.60)
STACK_BAND = (7.53, 2.71, 0.34, 0.09)    # collar below the stack rim
VENT = (7.28, 2.30, 0.13, 0.30)          # second, smaller roof vent
VENT_BAND = (7.235, 2.495, 0.22, 0.07)
SEAM_YS = (2.16, 1.17)                   # cladding seams on the exposed face
PANEL = (7.26, 1.24, 1.20, 0.84)         # flowsheet cutaway card
PANEL_ALPHA = 0.95                       # panel face, over the building
                                         # gradient (_shared_palette pins the
                                         # resulting composite -- keep in sync)
# The cutaway flowsheet: six *differentiated* unit silhouettes wired by seven
# streams, the last of which recycles the column bottoms into the feed — the
# same feedback idea the scene's own dashed arc carries, one scale down.
# Differentiation matters more than count at the 556 px delivery size, where
# the whole panel is ~67 px wide: a triangle, a circle, a box, a domed
# rectangle, a capsule and a second circle stay tellable apart where six
# boxes would not. Every unit is *stroked, never filled*, so the amber
# 'lit' pass can re-draw the entire sheet in one pass over the same geometry.
# The one exception is the stirred reactor, drawn as a miniature of the
# reactor stage: its teal broth and spinning impeller sit *under* the
# strokes, so the lit pass still re-traces its outline like every other unit.
FS_PUMP, FS_R_PUMP = (7.40, 1.50), 0.055        # feed pump (triangle)
FS_HX, FS_R_HX = (7.60, 1.50), 0.070            # exchanger (circle + duty bar)
FS_REACTOR = (7.68, 1.72, 0.21, 0.20)           # mini stirred reactor
FS_RX_LIQ = 0.62          # broth surface, as a fraction of the mini's height
FS_RX_BLADE = 0.045       # mini impeller half-span at full extension
FS_RX_BUB_CYCLES = 3      # mini bubble rises this many times per loop
FS_COLUMN, FS_DOME = (7.98, 1.45, 0.14, 0.40), 0.055     # domed column
FS_DRUM = (8.22, 1.78, 0.17, 0.095)             # reflux drum (capsule)
FS_TANK, FS_R_TANK = (8.305, 1.50), 0.085       # product tank (circle)
FS_STREAMS = (
    ((7.31, 1.50), (7.345, 1.50)),                                  # feed in
    ((7.455, 1.50), (7.53, 1.50)),                                  # pump->HX
    ((7.67, 1.50), (7.785, 1.50), (7.785, 1.72)),                   # HX->rx
    ((7.89, 1.80), (7.935, 1.80), (7.935, 1.66), (7.98, 1.66)),     # rx->col
    ((8.05, 1.905), (8.05, 1.955), (8.305, 1.955), (8.305, 1.875)), # overhead
    ((8.305, 1.78), (8.305, 1.585)),                                # -> tank
    ((8.05, 1.45), (8.05, 1.36), (7.49, 1.36), (7.49, 1.50)),       # recycle
)
FS_LEVEL_CYCLES = 3       # drum/tank level breathing, integer cycles per loop
# Distillation column: COLUMN is the straight *shell*; elliptical heads are
# added above and below it, so the silhouette tops out at 1.30 + 1.36 + 0.17
# = 2.83 — clear of the gauge, whose face is the semicircle above y 3.03.
COLUMN = (8.76, 1.30, 0.44, 1.36)        # shell only (x, y, w, h)
COL_HEAD_T, COL_HEAD_B = 0.17, 0.10      # elliptical head rise, top / bottom
COL_TRAY_YS = (1.58, 2.48)               # first / last tray line
COL_N_TRAYS = 5
COL_SKIRT = (8.84, 1.12, 0.28, 0.12)     # support skirt under the bottom head
COL_PLINTH = (8.78, 1.05, 0.40, 0.07)    # foundation tick on the ground line
PIPE_Y, PIPE2_Y = 2.02, 1.50             # building -> column runs
OVERHEAD_Y = 2.925                       # overhead vapor run above the column
DRUM = (9.04, 2.84, 0.28, 0.17)          # condenser / reflux drum
RETURN_X, RETURN_Y = 9.38, 1.44          # return leg down to the column base
# Plant ambient motion. Every constant below is an integer number of cycles
# per DUR-second loop (or, for the slug trains, a *speed* that
# ``_flow_marks`` rounds to an integer cycle count), so the wrap stays
# invisible. Speeds are in data units / s: at 200 dpi one data unit is 200 px
# rendered and 55.6 px at the 556 px delivery width, so ~0.3-0.45 data/s is
# 17-25 px/s on the page — plainly moving without strobing at 20 fps.
VAPOR_CYCLES, VAPOR_PUFFS = 3, 5          # stack/vent plumes
COL_SUMP_Y = 1.46                         # nominal column sump level
COL_LEVEL_CYCLES = 4                      # sump/drum level breathing
COL_VAPOR_CYCLES = 3                      # vapor ticks rising between trays
FLOW_SLUG_V = 0.30                        # flowsheet stream slug speed
PIPE_SLUG_V = 0.45                        # plant pipe-run slug speed
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


def _path_pts(pts, step=0.004):
    """Densely resample a polyline; returns (points, cumulative arc length)."""
    xy = [np.asarray(pts[0], float)]
    for p0, p1 in zip(pts[:-1], pts[1:]):
        a, b = np.asarray(p0, float), np.asarray(p1, float)
        k = max(int(np.hypot(*(b - a))/step), 2)
        xy.extend(a + (b - a)*i/k for i in range(1, k + 1))
    xy = np.asarray(xy)
    d = np.r_[0.0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
    return xy, d


def _flow_marks(ax, th, pts, t, speed, length=0.07, lw=2.0, color=None,
                alpha=0.95, zorder=2.68, phase=0.0, max_cycles=12):
    """Slugs of material sliding along a pipe/stream run.

    The plant's counterpart to the reactor's bubbles: a train of short bright
    capsules that ride ``pts`` (a polyline, given in flow direction) and read
    as flow rather than as plumbing. Each capsule is a sub-range of the
    resampled path, so it bends around corners instead of jumping them.

    Rather than fixing a cycle count per run — which would make slugs crawl
    down the 1.7-long return leg and blur along the 0.21 feed runs — the
    caller fixes a *speed* and this solves for the nearest integer number of
    cycles per loop (``max_cycles`` caps the rate on stub-length runs). An
    integer count is what keeps the loop seamless; capsules also enter and
    leave past the path ends, so at u = 0 and u = 1 nothing is drawn and the
    wrap has nothing to show. Mark count follows path length, keeping the
    spacing between slugs roughly constant across runs.
    """
    xy, d = _path_pts(pts)
    L = float(d[-1])
    if L <= 1e-9:
        return
    span = min(length, 0.5*L)
    n = max(1, int(round(L/0.30)))
    cycles = int(min(max(round(speed*DUR/(L + span)), 1), max_cycles))
    for i in range(n):
        u = (cycles*t/DUR + phase + i/n) % 1.0
        s0 = u*(L + span) - span
        a, b = max(s0, 0.0), min(s0 + span, L)
        m = (d >= a) & (d <= b)
        if int(m.sum()) < 2:
            continue
        ax.plot(xy[m, 0], xy[m, 1], color=color or th['stroke'], lw=lw,
                solid_capstyle='round', solid_joinstyle='round', alpha=alpha,
                zorder=zorder)


def draw_caption(ax, th, x, text):
    ax.text(x, CAPTION_Y, text, ha='center', va='center', fontsize=11.5,
            color=th['text'], alpha=0.85, zorder=3)


# %% Drawing

HALO_PEAK = 0.25     # composite opacity of a stage glow at its core; also
                     # the top of the ramp _shared_palette pins (see there)


def _halo(ax, th, xy, rx, ry, strength):
    """Soft amber glow behind a stage while a pulse passes it.

    A stack of nested rings, not two: with only a couple of them the
    outermost steps straight from nothing to its own alpha, and at the 200 px
    delivery size that edge reads as a hard-rimmed brown ellipse on the dark
    theme.

    Rather than hand-tuning per-ring alphas, aim at the *composite* opacity:
    ``cum`` is a smoothstep from HALO_PEAK at the core to 0 at ``k_max``, and
    each ring's alpha is solved so the stack reproduces it exactly
    (1 - a_i = (1 - cum_i)/(1 - cum_{i+1})). Smoothstep has zero slope at
    both ends, which is what removes the visible rim; HALO_PEAK keeps
    roughly the old pair's center weight.

    The widest ring (k_max) sets the footprint, so it also sets the
    clearance: the plant halo's rx = 1.12 about x = 8.30 keeps the outer ring
    inside the right canvas edge (8.30 + 1.12*1.45 = 9.92 < 10), and the
    reactor and microbe halos stay well clear of every edge at the same k.
    """
    if strength <= 0:
        return
    n, k_max = 20, 1.45
    u = np.linspace(0.0, 1.0, n)
    ks = 1.0 + u*(k_max - 1.0)
    cum = HALO_PEAK*(1.0 - u*u*(3.0 - 2.0*u))
    s = min(strength, 1.0)
    for i in range(n - 1, -1, -1):          # outermost ring first
        outside = cum[i + 1] if i + 1 < n else 0.0
        a = 1.0 - (1.0 - cum[i])/(1.0 - outside)
        if a <= 0.0:
            continue
        ax.add_patch(Ellipse(xy, 2*rx*ks[i], 2*ry*ks[i], fc=th['accent'],
                             ec='none', alpha=a*s, zorder=1))


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
    # granules (quiet navy), placed as offsets from the body center so the
    # whole microbe translates with MICROBE_C
    for dx, dy, r in ((-0.25, 0.07, 0.05), (-0.11, -0.13, 0.045),
                      (-0.19, -0.03, 0.028)):
        ax.add_patch(Circle((cx + dx, cy + dy), r, fc=th['stroke'], ec='none',
                            alpha=0.85, zorder=3))
    # engineering cue: amber plasmid loop (nested loops = supercoil)
    plasmid = (cx + 0.23, cy - 0.03)
    ax.add_patch(Circle(plasmid, 0.10, fc='none', ec=th['accent'],
                        lw=2.5, zorder=3))
    ax.add_patch(Circle(plasmid, 0.045, fc='none', ec=th['accent'],
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
        ax.text(SLIDER_LABEL_X, yy, lab, ha='center', va='center', fontsize=10,
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
    # product curve (attention amber) over the broth, with a soft under-glow.
    # The 1.44 baseline sits below the lower impeller tier (y = 1.55) so the
    # flat left end does not slice through the blades; boosted, the curve
    # tops out at 2.02, still under the broth surface.
    cxs = np.linspace(x + 0.18, x + w - 0.18, 40)
    amp = 0.30 + 0.28*curve_boost
    cys = 1.44 + amp/(1 + np.exp(-(cxs - 4.95)*9))
    ax.plot(cxs, cys, color=th['accent'], lw=7, alpha=0.25,
            solid_capstyle='round', zorder=3.4)
    ax.plot(cxs, cys, color=th['accent'], lw=2.8, solid_capstyle='round',
            zorder=3.5)
    # vessel outline + glass sheen on top of everything inside. The sheen
    # stops just under the broth surface: run past it into the headspace and
    # it stops reading as a reflection on the wetted wall and starts reading
    # as a rod hanging in mid-air. Its alpha is per-theme ('sheen_a').
    ax.add_patch(_vessel_patch(ec=th['stroke'], lw=3, zorder=4))
    ax.plot([x + 0.16, x + 0.16], [y + 0.28, LIQ_Y - 0.05],
            color=th['glass_hi'], lw=5, alpha=th['sheen_a'],
            solid_capstyle='round', zorder=4.2)


def _column_poly(fc='none', ec='none', lw=0, zorder=2):
    """Column silhouette: straight shell capped by two elliptical heads.

    Traced as one closed polygon (rather than a rounded box plus separate
    caps) so the outline is a single continuous stroke — no seam line where a
    head meets the shell, which at the 200 px delivery size would read as a
    stray tick rather than a weld.
    """
    x, y, w, h = COLUMN
    cx, rx, yt = x + w/2, w/2, y + h
    top = [(cx + rx*np.cos(a), yt + COL_HEAD_T*np.sin(a))
           for a in np.linspace(0.0, np.pi, 28)]
    bot = [(cx + rx*np.cos(a), y + COL_HEAD_B*np.sin(a))
           for a in np.linspace(np.pi, 2*np.pi, 28)]
    return Polygon(top + bot, closed=True, fc=fc, ec=ec, lw=lw, zorder=zorder)


def _fs_column_poly(ec='none', lw=0, alpha=1.0, zorder=4):
    """Flowsheet column: straight shell closed by one domed head, as a poly.

    Same single-stroke trick as ``_column_poly`` — inside the panel a seam
    tick where a cap meets a shell would be a stray pixel, not a weld.
    """
    x, y, w, h = FS_COLUMN
    cx, rx, yt = x + w/2, w/2, y + h
    dome = [(cx + rx*np.cos(a), yt + FS_DOME*np.sin(a))
            for a in np.linspace(0.0, np.pi, 16)]
    return Polygon(dome + [(x, y), (x + w, y)], closed=True, fc='none',
                   ec=ec, lw=lw, alpha=alpha, zorder=zorder)


def _fs_reactor_patch(**kw):
    x, y, w, h = FS_REACTOR
    return FancyBboxPatch((x, y), w, h,
                          boxstyle='round,pad=0,rounding_size=0.045', **kw)


def _draw_fs_reactor(ax, th, t):
    """The flowsheet's stirred reactor, as a miniature of the reactor stage.

    Same vocabulary as ``draw_reactor`` two scales down: teal broth under a
    waving surface, a rising bubble, and a Rushton tier spinning in side view
    on the same shaft speed (``SPIN_REVS``), so the two reactors visibly turn
    together. The color link is the point -- the only teal inside the plant
    is the unit that *is* the reactor stage, one scale down.

    One tier, not two: at 42 x 40 px a second would fuse into the first. The
    blade half-span and stroke weights are set to the big vessel's
    *proportions* rather than scaled down from its absolute values, which is
    what keeps it reading as the same machine rather than a fat smudge.

    Everything here sits at zorder 3.88-3.95 -- above the panel card (3.5),
    below the flowsheet's stroke passes (4) -- and is drawn once per frame
    rather than once per pass, so the fill is not double-composited and the
    amber ``flow_lit`` pass still re-traces the outline like every other unit.
    """
    x, y, w, h = FS_REACTOR
    cx, liq, by = x + w/2, y + FS_RX_LIQ*h, y + 0.055
    clip = _fs_reactor_patch(fc='none', ec='none', zorder=3.88)
    ax.add_patch(clip)
    gradient_fill(ax, clip, (x, y, x + w, liq), th['broth0'], th['broth1'],
                  zorder=3.9)
    # broth surface, riding the same WAVE_CYCLES as the big vessel. The
    # amplitude is ~1 px rendered: a live edge at full size, and below the
    # noise floor (not a flicker) once downscaled.
    xs = np.linspace(x + 0.012, x + w - 0.012, 24)
    ax.plot(xs, liq + 0.005*np.sin(2*np.pi*(xs - x)*14
                                   - 2*np.pi*WAVE_CYCLES*t/DUR),
            color=th['teal_edge'], lw=1.4, alpha=0.9, zorder=3.92)
    u = (FS_RX_BUB_CYCLES*t/DUR) % 1.0
    ax.add_patch(Circle((cx + 0.062, y + 0.022 + u*(liq - y - 0.05)), 0.010,
                        fc='none', ec=th['bubble'], lw=1.2,
                        alpha=0.9*(1.0 - smooth(0.82, 1.0, u)), zorder=3.93))
    # shaft + one Rushton tier: two blade pairs 90 deg out of phase, front
    # pair dark, back pair lighter -- draw_reactor's side-view trick verbatim
    ax.plot([cx, cx], [y + h, by], color=th['stroke'], lw=1.4, zorder=3.95)
    theta = 2*np.pi*SPIN_REVS*t/DUR
    for phase, front in ((theta, True), (theta + np.pi/2, False)):
        half = FS_RX_BLADE*abs(np.cos(phase))
        if half < 0.006:
            continue
        ax.plot([cx - half, cx + half], [by, by],
                color=th['stroke'] if front else th['navy_mid'], lw=2.2,
                solid_capstyle='round', alpha=0.95 if front else 0.6,
                zorder=3.95)


def _flowsheet_pass(ax, th, levels, color, lw, alpha, zorder=4):
    """Stroke the whole cutaway flowsheet once, in one color.

    Called twice — a steel base pass, then the amber ``flow_lit`` pass — so
    every unit, internal and stream lights together when a pulse reaches the
    plant. Nothing here is filled, which is what lets the second pass simply
    re-trace the first.
    """
    drum_lv, tank_lv = levels
    kw = dict(fc='none', ec=color, lw=lw, alpha=alpha, zorder=zorder)
    line = dict(color=color, lw=lw, alpha=alpha, zorder=zorder,
                solid_capstyle='round', solid_joinstyle='round')
    # feed pump: triangle, flat face to the feed, apex to the exchanger
    pxc, pyc = FS_PUMP
    ax.add_patch(Polygon([(pxc - FS_R_PUMP, pyc - FS_R_PUMP),
                          (pxc - FS_R_PUMP, pyc + FS_R_PUMP),
                          (pxc + FS_R_PUMP, pyc)], closed=True, **kw))
    # exchanger: circle crossed by its duty bar
    hx, hy = FS_HX
    ax.add_patch(Circle(FS_HX, FS_R_HX, **kw))
    dd = FS_R_HX*0.70
    ax.plot([hx - dd, hx + dd], [hy - dd, hy + dd], **line)
    # stirred reactor: only its rounded vessel outline -- the broth and the
    # spinning impeller inside it are drawn once, by _draw_fs_reactor
    ax.add_patch(_fs_reactor_patch(**kw))
    # column: domed shell with three tray decks
    ax.add_patch(_fs_column_poly(ec=color, lw=lw, alpha=alpha, zorder=zorder))
    cx0, cy0, cw, ch = FS_COLUMN
    for ty in np.linspace(cy0 + 0.22*ch, cy0 + 0.82*ch, 3):
        ax.plot([cx0 + 0.20*cw, cx0 + 0.80*cw], [ty, ty],
                color=color, lw=lw*0.75, alpha=alpha, zorder=zorder)
    # reflux drum and product tank, each holding a breathing liquid level
    dx0, dy0, dw, dh = FS_DRUM
    ax.add_patch(FancyBboxPatch((dx0, dy0), dw, dh,
                                boxstyle='round,pad=0,rounding_size=0.045',
                                **kw))
    ax.plot([dx0 + 0.014, dx0 + dw - 0.014], [dy0 + dh*drum_lv]*2,
            color=color, lw=lw*0.75, alpha=alpha*0.85, zorder=zorder)
    tcx, tcy = FS_TANK
    ax.add_patch(Circle(FS_TANK, FS_R_TANK, **kw))
    ly = tcy + FS_R_TANK*(2.0*tank_lv - 1.0)
    half = float(np.sqrt(max(FS_R_TANK**2 - (ly - tcy)**2, 0.0)))
    ax.plot([tcx - half, tcx + half], [ly, ly], color=color, lw=lw*0.75,
            alpha=alpha*0.85, zorder=zorder)
    for pts in FS_STREAMS:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], **line)


def draw_plant(ax, th, t, glow, flow_lit):
    _halo(ax, th, (8.30, 1.95), 1.12, 1.05, glow)
    x, y, w, h = PLANT_BLDG
    soft_shadow(ax, (x + w/2 + 0.5, y - 0.02), 1.15, 0.08, th)
    # building, main stack and the smaller roof vent: navy gradient + stroke
    for (bx, by, bw, bh) in (PLANT_BLDG, STACK, VENT):
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
    # cladding seams across the exposed building face (quiet, low-alpha steel:
    # any stronger and they compete with the flowsheet panel in front)
    for sy_ in SEAM_YS:
        ax.plot([x + 0.09, x + w - 0.09], [sy_, sy_], color=th['steel'],
                lw=1.3, alpha=0.30, zorder=3.05)
    # roof parapet and the two stack collars, in the darkest building navy so
    # they read against the lighter top of the building gradient
    for (bx, by, bw, bh) in (PARAPET, STACK_BAND, VENT_BAND):
        ax.add_patch(Rectangle((bx, by), bw, bh, fc=th['bldg0'],
                               ec=th['stroke'], lw=1.8, zorder=3.1))
    # vapor plumes off the stack and the vent: five puffs each, billowing as
    # they rise and drifting *left* — away from the $ gauge (center 8.28,
    # radius 0.50), which is what buys the room to make them big enough to
    # read. The widest stack puff clears the gauge face by ~0.09 and tops out
    # at y 3.51, inside the 3.6 canvas. Alpha follows sin(pi*u), so a puff
    # fades in at the rim and out at the top: the plume never pops, and the
    # loop wrap has nothing to show since alpha is 0 at both u = 0 and u = 1.
    for (bx, by, bw, bh), rise, drift, r0, r1, aa in (
            (STACK, 0.38, -0.24, 0.060, 0.130, 0.62),
            (VENT, 0.26, -0.08, 0.038, 0.085, 0.48)):
        sx, sy = bx + bw/2, by + bh
        for i in range(VAPOR_PUFFS):
            u = (VAPOR_CYCLES*t/DUR + i/VAPOR_PUFFS) % 1.0
            wob = 0.035*np.sin(2*np.pi*(u + i/VAPOR_PUFFS))
            a = aa*np.sin(np.pi*u)**0.7
            # drift ~ u**0.6: the plume leans away early, where it would
            # otherwise brush the gauge's baseline corner, and lands in the
            # same place at u = 1.
            c = (sx + drift*u**0.6 + wob, sy + 0.04 + rise*u)
            # three nested discs per puff: same footprint as one flat disc,
            # but the rim fades instead of ending on a hard circle, which is
            # what separates vapor from a bubble at this size.
            for k, ka in ((1.0, 0.34), (0.74, 0.34), (0.48, 0.40)):
                ax.add_patch(Circle(c, (r0 + r1*u)*k, fc=th['vapor'],
                                    ec='none', alpha=a*ka, zorder=2.5))
    # flowsheet cutaway panel; base pass in steel, lit pass in amber
    px, py, pw, phh = PANEL
    ax.add_patch(FancyBboxPatch((px, py), pw, phh,
                                boxstyle='round,pad=0,rounding_size=0.08',
                                fc=th['card'], ec=th['card_edge'], lw=1.2,
                                alpha=PANEL_ALPHA, zorder=3.5))
    ph_lv = 2*np.pi*FS_LEVEL_CYCLES*t/DUR
    levels = (0.45 + 0.16*np.sin(ph_lv), 0.50 + 0.18*np.sin(ph_lv + 2.1))
    _draw_fs_reactor(ax, th, t)
    for color, alpha, lw in ((th['steel'], 1.0, 1.8),
                             (th['accent'], min(flow_lit, 1.0), 2.4)):
        if alpha <= 0:
            continue
        _flowsheet_pass(ax, th, levels, color, lw, alpha)
    # slugs riding the streams, drawn over both passes so the sheet keeps
    # moving while it is lit. Navy/pale 'stroke' (never amber — amber is the
    # attention channel) is the highest-contrast bead available against a
    # steel line on the panel card, in both themes.
    for k, pts in enumerate(FS_STREAMS):
        _flow_marks(ax, th, pts, t, FLOW_SLUG_V, length=0.055, lw=1.7,
                    zorder=4.35, phase=0.11*k)
    colx, coly, colw, colh = COLUMN
    cx = colx + colw/2
    # two feed runs from the building into the column shell, one gated by a
    # bowtie valve, the other flanged. Drawn under the column outline (2.6)
    # so they tuck behind the shell rather than crossing it.
    for py_ in (PIPE_Y, PIPE2_Y):
        ax.plot([x + w, colx], [py_, py_], color=th['steel_dark'], lw=2.5,
                zorder=2.6)
    mx = (x + w + colx)/2
    ax.plot([mx, mx], [PIPE2_Y - 0.07, PIPE2_Y + 0.07], color=th['steel_dark'],
            lw=2, zorder=2.6)   # flange tick
    for sgn in (-1, 1):
        ax.add_patch(Polygon([(mx + sgn*0.055, PIPE_Y - 0.048),
                              (mx + sgn*0.055, PIPE_Y + 0.048),
                              (mx, PIPE_Y)], closed=True,
                             fc=th['steel_dark'], ec='none', zorder=2.7))
    # overhead vapor line -> condenser drum -> return leg into the column base
    riser = ((cx, coly + colh + COL_HEAD_T), (cx, OVERHEAD_Y),
             (DRUM[0], OVERHEAD_Y))
    ret = ((DRUM[0] + DRUM[2], OVERHEAD_Y), (RETURN_X, OVERHEAD_Y),
           (RETURN_X, RETURN_Y), (colx + colw, RETURN_Y))
    for run in (riser, ret):
        ax.plot([p[0] for p in run], [p[1] for p in run],
                color=th['steel_dark'], lw=2.2, solid_capstyle='round',
                solid_joinstyle='round', zorder=2.6)
    for fy in (2.30, 1.82):
        ax.plot([RETURN_X - 0.055, RETURN_X + 0.055], [fy, fy],
                color=th['steel_dark'], lw=2, zorder=2.6)   # flange ticks
    # slug trains on all four runs: feed in from the building, vapor up the
    # riser, condensate back down the long return leg. One shared speed, so
    # the whole plant appears to move material at one rate.
    for k, run in enumerate((((x + w, PIPE_Y), (colx, PIPE_Y)),
                             ((x + w, PIPE2_Y), (colx, PIPE2_Y)),
                             riser, ret)):
        _flow_marks(ax, th, run, t, PIPE_SLUG_V, length=0.10, lw=2.0,
                    alpha=0.9, zorder=2.68, phase=0.17*k)
    drum = FancyBboxPatch(DRUM[:2], DRUM[2], DRUM[3],
                          boxstyle='round,pad=0,rounding_size=0.085',
                          fc=th['glass1'], ec=th['stroke'], lw=2, zorder=3)
    ax.add_patch(drum)
    # condensate holdup in the drum, breathing anti-phase to the column sump
    # below: the two inventories read as trading material back and forth.
    dlv = DRUM[1] + DRUM[3]*(0.46 - 0.17*np.sin(2*np.pi*COL_LEVEL_CYCLES
                                                * t/DUR))
    dliq = Rectangle(DRUM[:2], DRUM[2], dlv - DRUM[1], fc=th['peri'],
                     ec='none', alpha=0.55, zorder=3.02)
    ax.add_patch(dliq)
    dliq.set_clip_path(drum)
    ax.add_patch(FancyBboxPatch(DRUM[:2], DRUM[2], DRUM[3],
                                boxstyle='round,pad=0,rounding_size=0.085',
                                fc='none', ec=th['stroke'], lw=2,
                                zorder=3.05))
    # column: plinth + skirt below, glass-gradient shell with tray decks above
    ax.add_patch(Rectangle(COL_PLINTH[:2], COL_PLINTH[2], COL_PLINTH[3],
                           fc=th['stroke'], ec='none', zorder=3.2))
    col = _column_poly(zorder=2)
    ax.add_patch(col)
    gradient_fill(ax, col, (colx, coly - COL_HEAD_B, colx + colw,
                            coly + colh + COL_HEAD_T),
                  th['glass0'], th['glass1'], direction='h', zorder=2)
    # column internals, all under the shell outline (3) and the tray decks:
    # a sump whose level breathes, and vapor ticks rising tray to tray. This
    # is the plant's answer to the reactor's bubble column — the cue that
    # something is happening *inside* the equipment, not just around it.
    lvl = COL_SUMP_Y + 0.030*np.sin(2*np.pi*COL_LEVEL_CYCLES*t/DUR)
    sump = Rectangle((colx, coly - COL_HEAD_B), colw,
                     lvl - (coly - COL_HEAD_B), fc=th['peri'], ec='none',
                     alpha=0.45, zorder=2.45)
    ax.add_patch(sump)
    sump.set_clip_path(col)
    ax.plot([colx + 0.045, colx + colw - 0.045], [lvl, lvl],
            color=th['steel_dark'], lw=1.8, alpha=0.8, zorder=2.5)
    for tx, off in ((cx - 0.11, 0.0), (cx, 0.45), (cx + 0.11, 0.75)):
        for j in (0, 1):
            u = (COL_VAPOR_CYCLES*t/DUR + off + 0.5*j) % 1.0
            ty = COL_TRAY_YS[0] + 0.04 + u*(COL_TRAY_YS[1] - COL_TRAY_YS[0])
            ax.plot([tx - 0.035, tx + 0.035], [ty, ty], color=th['steel_dark'],
                    lw=1.8, solid_capstyle='round',
                    alpha=0.55*np.sin(np.pi*u)**0.6, zorder=2.55)
    ax.add_patch(Rectangle(COL_SKIRT[:2], COL_SKIRT[2], COL_SKIRT[3],
                           fc=th['navy_mid'], ec=th['stroke'], lw=1.8,
                           zorder=3.15))
    ax.add_patch(_column_poly(ec=th['stroke'], lw=3, zorder=3))
    for ty in np.linspace(*COL_TRAY_YS, COL_N_TRAYS):
        ax.plot([colx + 0.055, colx + colw - 0.055], [ty, ty],
                color=th['steel'], lw=2.2, alpha=0.85, zorder=3)
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
    draw_caption(ax, th, MICROBE_C[0], 'Metabolic engineering')
    draw_caption(ax, th, 5.0, 'Kinetics & reactor design')
    draw_caption(ax, th, 8.35, 'Facility-scale economics')
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
# One 10 s loop, one iteration. Ambient motion (impeller spin, bubbles,
# surface wave, vapor plumes, pipe/flowsheet slug trains, column vapor and
# levels) derives from state['t'] and runs in every frame, each with an
# integer number of cycles per loop so the wrap is seamless. On top, three
# comets circulate, each departing its origin stage as that stage's widget
# starts to move and lighting the next stage as it arrives:
#   Comet 1 (0.2-2.0 s)  microbe -> reactor, leaving as the slider starts;
#                        the reactor lights on arrival (~2.0 s).
#   Comet 2 (2.0-3.8 s)  reactor -> plant, leaving as the product curve starts;
#                        the plant lights on arrival (~3.8 s).
#   Comet 3 (3.8-9.8 s)  plant -> microbe on the dashed feedback arc, leaving
#                        as the $ needle starts; slider, curve, and needle ease
#                        home together beneath it, masking the reset, and the
#                        microbe lights again as the next loop opens.


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
    """Scene state at time t (seconds) within the DUR-second loop.

    One clean iteration per loop, carried by three comets that circulate the
    pipeline. Each comet departs a stage at the instant that stage's one
    widget (its last movement) begins to move, and the next stage lights as
    that comet arrives:

    * the microbe's **slider** starts (t 0.2) and sends a comet to the reactor;
    * the reactor's **product curve** starts (t 2.0) and sends one to the plant;
    * the plant's **$ needle** starts (t 3.8) and sends the right-to-left
      feedback comet back to the microbe.

    Each moved widget holds its new value while the downstream stages respond,
    then all three ease home together (``home``) under the feedback comet,
    masking the loop reset.
    """
    s = default_state()
    s['t'] = t
    # Comet 1, microbe -> reactor: departs as the slider starts (0.2), crosses
    # arrow 1 and vanishes into the reactor (pulse_xy returns None past s 0.34),
    # where the reactor then lights.
    if 0.2 <= t < 2.0:
        s['pulse1_s'] = 0.42*(t - 0.2)/1.8
    # Comet 2, reactor -> plant: departs as the curve starts (2.0), emerges from
    # the reactor at arrow 2 (s 0.50) and vanishes into the plant (s > 0.80).
    if 2.0 <= t < 3.8:
        s['pulse2_s'] = 0.50 + 0.40*(t - 2.0)/1.8
    # Comet 3, plant -> microbe: departs as the needle starts (3.8), riding the
    # dashed feedback arc back to the microbe by the loop's end.
    if 3.8 <= t < 9.8:
        s['fb_s'] = (t - 3.8)/6.0
    # Stage glows, each timed to its incoming comet's arrival. Every glow starts
    # a hair before its stage's widget, so the widget reads as the last movement
    # of the stage — the movement the outgoing comet departs with.
    s['glow_microbe'] = bump(0.1, 1.2, t)
    s['glow_reactor'] = bump(1.7, 2.7, t)
    s['glow_plant'] = bump(3.4, 4.4, t)
    # One widget per stage, moving as that stage lights and holding until the
    # shared ease-home under the feedback comet returns all three to base by
    # the loop's end (so frame 0 and the last frame coincide).
    home = smooth(5.2, 9.7, t)
    s['slider_frac'] = 0.35 + 0.40*(smooth(0.2, 1.0, t) - home)
    s['curve_boost'] = smooth(2.0, 2.8, t) - home
    s['needle_frac'] = 0.50 - 0.30*(smooth(3.8, 4.6, t) - home)
    return s


# Theme colors forced into the palette exactly, in priority order. Both are
# large, flat, low-saturation fills, which is precisely where a nearest-match
# miss is most visible: the eye reads a whole uniform panel as the wrong hue.
PINNED_KEYS = ('bg', 'card')


def _pin_exact(pal, rgb, taken):
    """Overwrite the nearest not-yet-pinned palette slot with ``rgb`` exactly.

    Returns the slot index. Editing the *nearest* slot keeps the disturbance
    minimal — that slot moves by at most the pinned color's own quantization
    error — and guarantees the pinned color quantizes to itself afterwards
    (distance 0).
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
    and the light theme's 'card' face otherwise lands on a neighbour whose
    green exceeds its red and blue, flipping the whole control card from
    blue-tinted lavender to mint. Each pin takes its own slot so a later pin
    cannot overwrite an earlier one.

    A ramp of stage-glow stops is pinned the same way. ``_halo`` draws a
    smooth amber-over-background falloff, but MAXCOVERAGE spends its budget
    on the scene's distinct hues and leaves that near-background ramp almost
    nothing: measured on the dark theme, a scan across the plant glow holds
    23 shades in the render and 3 in the GIF, turning the soft falloff back
    into the hard-edged ellipse ``_halo`` exists to avoid. Twelve evenly
    spaced stops up to HALO_PEAK cost ~5% of the palette and restore it.
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
    keyed = []
    for key in PINNED_KEYS:
        rgb = tuple(int(th[key][i:i + 2], 16) for i in (1, 3, 5))
        taken.add(_pin_exact(pal, rgb, taken))
        keyed.append(np.asarray(rgb, float))
    bg, acc = (np.array([int(th[k][i:i + 2], 16) for i in (1, 3, 5)], float)
               for k in ('bg', 'accent'))
    stops = 12
    for i in range(1, stops + 1):
        rgb = np.round(bg + (HALO_PEAK*i/stops)*(acc - bg))
        # Pillow's undithered lookup is cached per 8-level color cube, so a
        # stop sharing a cell with a PINNED_KEYS color captures that cell and
        # silently un-pins it: the faintest light-theme stop lands 5 levels
        # off pure white and stole 'bg'. Skip any stop within one cell of a
        # keyed pin — at that distance it is perceptually redundant anyway.
        if any(np.abs(rgb - p).max() < 8 for p in keyed):
            continue
        taken.add(_pin_exact(pal, tuple(int(v) for v in rgb), taken))
    # The flowsheet panel's face is 'card' at PANEL_ALPHA over the building
    # gradient, so it is *not* the pinned 'card' color and owns no slot. It
    # is also a 240 x 168 px flat region, exactly the case PINNED_KEYS exists
    # for: adding the mini reactor's teal to the panel was enough to pull the
    # light theme's panel from pale lavender to a pale mint. Pinning the
    # composite at the gradient midpoint captures the whole face, which
    # spreads only ~2 levels from left to right. On dark the composite lands
    # within a cell of 'card' itself, so the same guard as above skips it
    # rather than letting it steal that pin.
    card = np.array([int(th['card'][i:i + 2], 16) for i in (1, 3, 5)], float)
    bld = [np.array([int(th[k][i:i + 2], 16) for i in (1, 3, 5)], float)
           for k in ('bldg0', 'bldg1')]
    face = np.round(PANEL_ALPHA*card + (1.0 - PANEL_ALPHA)*0.5*(bld[0] + bld[1]))
    if not any(np.abs(face - p).max() < 8 for p in keyed):
        taken.add(_pin_exact(pal, tuple(int(v) for v in face), taken))
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
    # stills show visible banding in the gradients — on this scene it cost
    # roughly 6x the file size (9.99 / 9.10 MB), over the ~8 MB guardrail.
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
                    help='Floyd-Steinberg dithering (only if gradients band; '
                         'roughly 6x file size on this scene — check the '
                         '~8 MB guardrail)')
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
