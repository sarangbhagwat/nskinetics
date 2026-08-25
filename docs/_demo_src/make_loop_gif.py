# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Render the animated "strain-to-TEA" loop GIFs for the docs landing page.

The scene is a three-stage pipeline drawn left to right — an engineered
microbe with abstract sliders, a stirred-tank reactor, and a minimalistic
chemical plant with a flowsheet nested inside it plus a $ gauge — connected
by forward arrows and a feedback arrow curving back underneath. The
animation sends a pulse down the pipeline, moves a slider, ripples the
change through to the gauge, and pulses the feedback arrow while easing
everything back for a seamless loop.

Run with no arguments to (re)build both theme variants:

    python docs/_demo_src/make_loop_gif.py
    ->  docs/source/_static/images/demo/loop_light.gif
        docs/source/_static/images/demo/loop_dark.gif

Use ``--stills DIR`` to render one PNG per preset state for visual review.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import (Arc, Circle, Ellipse, FancyArrowPatch,
                                FancyBboxPatch, Rectangle)
from PIL import Image

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / 'source' / '_static' / 'images' / 'demo'

DPI = 100    # figsize (10, 3.6) -> 1000 x 360 px (displayed at 200 px height)
DUR = 8.0    # seconds per loop
FPS = 10

# 'bg' must match the docs page background (pydata-sphinx-theme): GIFs are
# opaque, so an off-theme bg shows as a visible rectangle on the page.
THEMES = {
    'light': dict(bg='#ffffff', stroke='#1f2a63', faint='#8a93b8',
                  accent='#f5a623', accent2='#ef7b45'),
    'dark':  dict(bg='#14181e', stroke='#d7dce8', faint='#5a6478',
                  accent='#f5a623', accent2='#ef7b45'),
}

# %% Scene geometry

SLIDER_X0, SLIDER_X1 = 0.65, 2.05
SLIDER_YS = (1.62, 1.34, 1.06)

ARROW1 = ((2.25, 1.95), (4.30, 1.95))
ARROW2 = ((5.70, 1.95), (7.10, 1.95))

# Feedback arrow: quadratic bezier dipping below the pipeline, right to left.
# FB_RAD reproduces the same control point via arc3 for FancyArrowPatch.
FB_P0, FB_P1, FB_P2 = (8.45, 0.92), (4.95, 0.12), (1.45, 0.90)
FB_RAD = -0.113

# Chemical-plant silhouette: building skyline with a stack; a separate
# column (rounded rect) stands to its right, joined by a pipe.
_PLANT_OUTLINE = [
    (7.30, 1.05), (7.30, 2.25), (7.70, 2.25), (7.70, 2.90), (7.90, 2.90),
    (7.90, 2.25), (8.60, 2.25), (8.60, 1.05),
]
# Nested flowsheet (cutaway view inside the building): unit boxes (x, y, w, h)
# and the stream lines connecting their edges.
_FLOWSHEET_BOXES = [(7.48, 1.52, 0.30, 0.24), (7.98, 1.72, 0.30, 0.24),
                    (8.18, 1.22, 0.30, 0.24)]
_FLOWSHEET_LINES = [((7.78, 1.64), (7.98, 1.84)), ((8.13, 1.72), (8.33, 1.46))]

GAUGE_C = (8.25, 3.02)
GAUGE_R = 0.34


def default_state():
    """Scene state for one frame; every animated quantity in one dict."""
    return dict(slider_frac=0.35,   # animated (middle) slider knob, 0..1
                needle_frac=0.50,   # gauge needle, 0 (far left) .. 1 (far right)
                curve_boost=0.0,    # 0 = base product curve, 1 = boosted
                pulse1_s=None,      # forward pulse path params, None = hidden
                pulse2_s=None,
                fb_s=None,          # feedback pulse path param, None = hidden
                glow_microbe=0.0, glow_reactor=0.0, glow_plant=0.0)


# %% Drawing

def _halo(ax, th, xy, rx, ry, strength):
    """Soft amber glow behind a stage while a pulse passes it."""
    if strength <= 0:
        return
    for k, a in ((1.35, 0.10), (1.0, 0.16)):
        ax.add_patch(Ellipse(xy, 2*rx*k, 2*ry*k, fc=th['accent'], ec='none',
                             alpha=a*min(strength, 1.0), zorder=1))


def draw_microbe(ax, th, glow):
    _halo(ax, th, (1.35, 2.45), 0.62, 0.45, glow)
    ax.add_patch(Ellipse((1.35, 2.45), 1.0, 0.6, angle=-16, fc='none',
                         ec=th['stroke'], lw=3, zorder=3))
    for cx, cy, r in ((1.15, 2.52, 0.055), (1.42, 2.33, 0.045),
                      (1.58, 2.52, 0.05)):
        ax.add_patch(Circle((cx, cy), r, fc=th['stroke'], ec='none', zorder=3))
    # flagellum
    xs = np.linspace(0, 0.45, 30)
    ax.plot(0.86 - xs, 2.62 + 0.06*np.sin(xs/0.45*2*np.pi), color=th['stroke'],
            lw=2, solid_capstyle='round', zorder=3)


def draw_sliders(ax, th, frac):
    """Three abstract sliders; the middle (amber) knob is the animated one."""
    knob_fracs = (0.30, frac, 0.72)
    for y, kf, animated in zip(SLIDER_YS, knob_fracs, (False, True, False)):
        ax.plot([SLIDER_X0, SLIDER_X1], [y, y], color=th['faint'], lw=4,
                solid_capstyle='round', zorder=2)
        kx = SLIDER_X0 + kf*(SLIDER_X1 - SLIDER_X0)
        color = th['accent'] if animated else th['stroke']
        ax.add_patch(Circle((kx, y), 0.09, fc=color, ec=th['bg'], lw=1.5,
                            zorder=4))


def draw_arrow(ax, th, p0, p1, rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, connectionstyle=f'arc3,rad={rad}',
                                 arrowstyle='-|>', mutation_scale=22,
                                 lw=2.5, color=th['faint'], zorder=2))


def draw_reactor(ax, th, glow, curve_boost):
    _halo(ax, th, (5.0, 2.05), 0.75, 0.95, glow)
    ax.add_patch(FancyBboxPatch((4.45, 1.3), 1.1, 1.5,
                                boxstyle='round,pad=0,rounding_size=0.14',
                                fc='none', ec=th['stroke'], lw=3, zorder=3))
    # impeller: shaft + blade bar
    ax.plot([5.0, 5.0], [2.8, 2.15], color=th['stroke'], lw=2, zorder=3)
    ax.plot([4.78, 5.22], [2.15, 2.15], color=th['stroke'], lw=2.5,
            solid_capstyle='round', zorder=3)
    # product curve rising inside the vessel wall
    xs = np.linspace(4.62, 5.38, 40)
    amp = 0.30 + 0.28*curve_boost
    ys = 1.48 + amp/(1 + np.exp(-(xs - 4.95)*9))
    ax.plot(xs, ys, color=th['accent'], lw=2.5, solid_capstyle='round',
            zorder=3)


def draw_plant(ax, th, glow, flow_lit):
    # halo kept inside the canvas (rx*1.35 outer ring) so it fades out rather
    # than clipping against the right edge
    _halo(ax, th, (8.30, 1.95), 1.12, 1.05, glow)
    xs, ys = zip(*_PLANT_OUTLINE)
    ax.plot(xs, ys, color=th['stroke'], lw=3, solid_joinstyle='round',
            solid_capstyle='round', zorder=3)
    # column + connecting pipe + ground line
    ax.add_patch(FancyBboxPatch((8.90, 1.05), 0.45, 1.75,
                                boxstyle='round,pad=0,rounding_size=0.16',
                                fc='none', ec=th['stroke'], lw=3, zorder=3))
    ax.plot([8.60, 8.90], [1.90, 1.90], color=th['stroke'], lw=2.5, zorder=3)
    ax.plot([7.15, 9.55], [1.05, 1.05], color=th['stroke'], lw=3,
            solid_capstyle='round', zorder=3)
    # nested flowsheet (cutaway): base pass in 'faint', lit pass in amber
    for color, alpha, lw in ((th['faint'], 1.0, 2.0),
                             (th['accent'], min(flow_lit, 1.0), 2.4)):
        if alpha <= 0:
            continue
        for x, y, w, h in _FLOWSHEET_BOXES:
            ax.add_patch(Rectangle((x, y), w, h, fc='none', ec=color,
                                   lw=lw, alpha=alpha, zorder=4))
        for (x0, y0), (x1, y1) in _FLOWSHEET_LINES:
            ax.plot([x0, x1], [y0, y1], color=color, lw=lw, alpha=alpha,
                    zorder=4)


def draw_gauge(ax, th, needle_frac):
    cx, cy = GAUGE_C
    ax.add_patch(Arc((cx, cy), 2*GAUGE_R, 2*GAUGE_R, theta1=0, theta2=180,
                     ec=th['stroke'], lw=3, zorder=3))
    ax.plot([cx - GAUGE_R, cx + GAUGE_R], [cy, cy], color=th['stroke'], lw=3,
            solid_capstyle='round', zorder=3)
    for ang in (150, 90, 30):
        a = np.deg2rad(ang)
        ax.plot([cx + 0.26*np.cos(a), cx + 0.32*np.cos(a)],
                [cy + 0.26*np.sin(a), cy + 0.32*np.sin(a)],
                color=th['faint'], lw=2, zorder=3)
    a = np.deg2rad(160 - 140*needle_frac)
    ax.plot([cx, cx + 0.24*np.cos(a)], [cy, cy + 0.24*np.sin(a)],
            color=th['accent2'], lw=2.5, solid_capstyle='round', zorder=4)
    ax.add_patch(Circle((cx, cy), 0.045, fc=th['stroke'], ec='none', zorder=5))
    ax.text(cx, cy - 0.28, '$', ha='center', va='center', fontsize=13,
            fontweight='bold', color=th['accent'], zorder=4)


def pulse_xy(s):
    """Forward-pulse position at path parameter s in [0, 1].

    Returns None while the pulse is "inside" a stage (the stage glows
    instead): s 0.35-0.50 crosses the reactor; s > 0.80 is inside the plant.
    """
    if s < 0.35:
        u = s/0.35
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


def draw_pulse(ax, th, xy):
    if xy is None:
        return
    ax.add_patch(Circle(xy, 0.14, fc=th['accent'], ec='none', alpha=0.35,
                        zorder=5))
    ax.add_patch(Circle(xy, 0.075, fc=th['accent'], ec='none', zorder=5))


def draw_scene(ax, th, state):
    draw_microbe(ax, th, state['glow_microbe'])
    draw_sliders(ax, th, state['slider_frac'])
    draw_arrow(ax, th, *ARROW1)
    draw_reactor(ax, th, state['glow_reactor'], state['curve_boost'])
    draw_arrow(ax, th, *ARROW2)
    draw_plant(ax, th, state['glow_plant'], min(state['glow_plant'], 1.0))
    draw_gauge(ax, th, state['needle_frac'])
    draw_arrow(ax, th, FB_P0, FB_P2, rad=FB_RAD)
    for s in (state['pulse1_s'], state['pulse2_s']):
        if s is not None:
            draw_pulse(ax, th, pulse_xy(s))
    if state['fb_s'] is not None:
        draw_pulse(ax, th, fb_xy(state['fb_s']))


def render_frame(th, state):
    """Render one frame to an opaque RGB PIL image (1000 x 360)."""
    fig = plt.figure(figsize=(10, 3.6), dpi=DPI)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.6)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor(th['bg'])
    draw_scene(ax, th, state)
    fig.canvas.draw()
    img = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
    plt.close(fig)
    return img


# %% Animation timeline
#
# One 8 s loop, three acts:
#   Act 1 (0-3 s)    forward pulse microbe -> reactor -> plant; each stage
#                    glows as the pulse passes; the needle settles at a
#                    starting position.
#   Act 2 (3-6.5 s)  a slider knob moves; a second pulse travels the same
#                    path; the product curve redraws higher; the needle
#                    swings to a new (lower-cost) position.
#   Act 3 (6.5-8 s)  the feedback arrow pulses right-to-left while slider,
#                    curve, and needle ease back to their initial values,
#                    masking the reset so the loop is seamless.


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
    if 0.2 <= t < 2.8:
        s['pulse1_s'] = (t - 0.2)/2.6
    if 3.8 <= t < 6.0:
        s['pulse2_s'] = (t - 3.8)/2.2
    if 6.5 <= t < 7.6:
        s['fb_s'] = (t - 6.5)/1.1
    # Glow envelopes are timed to the pulses' stage arrivals (pulse path
    # params 0.35/0.50 = reactor, 0.80 = plant).
    s['glow_microbe'] = bump(0.0, 0.9, t) + bump(3.0, 4.0, t)
    s['glow_reactor'] = bump(0.95, 1.75, t) + bump(4.45, 5.25, t)
    s['glow_plant'] = bump(2.2, 3.1, t) + bump(5.5, 6.4, t)
    s['slider_frac'] = 0.35 + 0.40*smooth(3.1, 3.7, t) - 0.40*smooth(7.0, 7.8, t)
    s['curve_boost'] = smooth(4.5, 5.1, t) - smooth(7.0, 7.8, t)
    s['needle_frac'] = (0.50 + 0.22*smooth(2.5, 3.0, t)
                        - 0.40*smooth(5.7, 6.3, t)
                        + 0.18*smooth(7.1, 7.9, t))
    return s


def _shared_palette(frames, bg_hex):
    """Build the one palette every frame is quantized against.

    Derived from a probe holding *all* frames, not from frame 0: at t = 0 no
    pulse, glow or boosted curve is on screen, so a frame-0 palette has no
    entry for those colors and nearest-match snaps them to wrong hues.
    MAXCOVERAGE (not MEDIANCUT) because the canvas is ~94% background —
    median cut spends its budget subdividing that near-uniform mass and
    starves the sparse accent colors that carry the animation.

    The background itself is then pinned to the exact theme color: the GIF is
    opaque and matted on the docs page, so a near-miss there would show as a
    visibly bounded off-white (or off-dark) rectangle on the page.
    """
    w, h = frames[0].size
    probe = Image.new('RGB', (w, h*len(frames)))
    for k, f in enumerate(frames):
        probe.paste(f, (0, k*h))
    base = probe.quantize(colors=128, method=Image.MAXCOVERAGE)
    bg = tuple(int(bg_hex[i:i + 2], 16) for i in (1, 3, 5))
    idx = Image.new('RGB', (1, 1), bg).quantize(
        palette=base, dither=Image.Dither.NONE).getpixel((0, 0))
    pal = list(base.getpalette())
    pal[3*idx:3*idx + 3] = list(bg)
    out = Image.new('P', (1, 1))
    out.putpalette(pal)
    out.load()   # refresh the core palette so quantize() sees the edit
    return out


def build_gif(theme_name, out_path):
    """Render all frames for one theme and write an infinite-loop GIF."""
    th = THEMES[theme_name]
    n = int(round(DUR*FPS))
    frames = [render_frame(th, timeline(i/FPS)) for i in range(n)]
    # Quantize every frame against one shared palette (flat art, no dither)
    # for a small file with no inter-frame palette flicker.
    base = _shared_palette(frames, th['bg'])
    pal_frames = [f.quantize(palette=base, dither=Image.Dither.NONE)
                  for f in frames]
    pal_frames[0].save(out_path, save_all=True, append_images=pal_frames[1:],
                       duration=int(1000/FPS), loop=0, optimize=True)
    print(f'{out_path} ({out_path.stat().st_size/1024:.0f} KiB, {n} frames)')


# %% CLI

# Hand-picked states for reviewing the scene art without the animation.
_PRESETS = {
    'base': {},
    'pulse': dict(pulse1_s=0.20, glow_microbe=0.4),
    'reactor-glow': dict(glow_reactor=1.0),
    'shifted': dict(slider_frac=0.75, curve_boost=1.0, needle_frac=0.30,
                    glow_plant=1.0),
    'feedback': dict(fb_s=0.5),
}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Build the landing-page loop GIFs (or review stills).')
    ap.add_argument('--stills', metavar='DIR',
                    help='render one PNG per preset state into DIR and exit')
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
    for theme_name in THEMES:
        build_gif(theme_name, OUT_DIR / f'loop_{theme_name}.gif')


if __name__ == '__main__':
    main()
