# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Conceptual diagram of the ``s_cerevisiae_ferm_fb_inhib_mod_ibo``
kinetic model (Antimony model ``bhagwat2026``, extending Lei et al. 2001,
J. Biotechnol. 88:205-21 / BIOMD0000000245 with an engineered isobutanol
pathway, product inhibition, aeration staging, and fed-batch feeding).

The figure shows the reaction network (r1-r11, r13-r17) and the control
structure — exponential product inhibition (r1, r4, r6, r7, r17),
threshold-gated product acceleration of biomass decay (r10), glucose
repression, the
acetaldehyde overflow signal, O2 dependence (``f_O2`` scaling), the AcDH
physiological-state machinery, and the fed-batch glucose-spike loop — sized
to Nature Communications figure specifications (180 mm double-column width,
Arial, 5-7 pt text, vector PDF with editable text plus a 600 dpi PNG).

The network and controls are curated directly from
``s_cerevisiae_ferm_fb_inhib_mod_ibo_antimony.txt``; continuous-mode dilution
terms (``s_glu_in``, ``*_out``; D = 0 in fed-batch use) are not drawn.

Rotary-knob icons mark the strain design variables: the strain-side decision
variables of the isobutanol biorefinery's ``metabolic_split_12d``
kinetic-optimization preset (enzyme capacities k_1l/k_1h/k_1e, k_3, k_6,
k_13, k_14/k_15/k_16, k_17 and the ethanol / isobutanol / acetate
product-inhibition coefficient groups). Its three feeding-policy variables
are process design variables and are not marked.

Run directly to write ``conceptual_diagram.png`` / ``.pdf`` next to this file::

    python conceptual_diagram.py
"""

import os

import numpy as np
import matplotlib

matplotlib.use('Agg') if os.environ.get('MPLBACKEND') is None else None
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib.path import Path
from matplotlib.lines import Line2D

__all__ = ('draw_conceptual_diagram',)

# --- Nature Communications sizing -----------------------------------------
MM = 1 / 25.4                     # mm -> inch
FIG_W_MM = 180.                   # double-column width
Y_TOP_MM = 134.                   # axes top with the process-control row ...
Y_TOP_NO_CONTROLS_MM = 104.       # ... and without it (glucose box tops at 101.5)
LEGEND_ROW_4_MM = 5.2             # depth of the legend strip below y = 0
                                  # (full figure 139.2 mm; Nat. Commun. cap 170 mm)
LEGEND_COL_GAP_MM = 3.5           # clear space between legend columns
LEGEND_PAD_MM = 4.                # legend strip padding left / right of its content

# --- Okabe-Ito palette (colorblind-safe), assigned by control job ----------
C_FLUX = '#3A3A3A'      # mass/reaction flux
C_INHIB = '#D55E00'     # product inhibition (vermillion)
C_REPR = '#CC79A7'      # glucose repression (reddish purple)
C_ACT = '#009E73'       # activation / overflow signal (bluish green)
C_FEED = '#0072B2'      # fed-batch feed & sensing (blue)
C_O2 = '#56B4E9'        # O2-dependence (f_O2 scaling) badge (sky blue)
C_O2_TEXT = '#2E7FB0'   # O2 lettering on white (supply tag, partial badge)
C_LEVER = '#E69F00'     # strain design variable knob (orange)
C_TEXT = '#1A1A1A'
C_MUTED = '#666666'

# Panel tints (light surfaces under white)
TINT_IBO = '#FDF3EA'
TINT_IBO_EDGE = '#E5B896'
TINT_BM = '#EDF6F2'
TINT_BM_EDGE = '#9CCDBB'
TINT_CTRL = '#EFF4FA'
TINT_CTRL_EDGE = '#9DB8D2'
TINT_SUBSTRATE = '#EAF2FA'
TINT_PRODUCT = '#FCEEE3'

FS_SPECIES = 6.5
FS_RXN = 5.0
FS_ENZ = 4.6
FS_TAG = 4.6
FS_NOTE = 5.0
FS_PANEL = 6.8
FS_LEGEND = 5.4

# --- strain design variables (knob positions, mm) --------------------------
# The strain-side decision variables of the isobutanol package's
# metabolic_split_12d kinetic-optimization preset. Capacity variables sit just
# after the enzyme name of their reaction; the three product-tolerance groups
# (inhib_ethanol / inhib_isobutanol / inhib_acetate) sit on every drawn
# product edge -- the inhibition stubs of r1/r4/r6/r7/r17 and the decay edge
# of r10 -- each of which names the effectors acting there.
STRAIN_LEVER_MARKS = {
    'glycolysis (k_1l, k_1h, k_1e) @ r1': (101.8, 87.6),
    'k_3 @ r3': (97.6, 65.4),
    'k_6 @ r6': (97.6, 42.2),
    'k_13 @ r13': (131.6, 79.9),
    'ehrlich_downstream (k_14) @ r14': (166.4, 71.5),
    'ehrlich_downstream (k_15) @ r15': (166.4, 58.5),
    'ehrlich_downstream (k_16) @ r16': (166.4, 45.5),
    'k_17 @ r17': (166.4, 32.5),
    'inhib_* @ r1 stub': (78.6, 90.0),
    'inhib_* @ r4 stub': (62.4, 48.6),
    'inhib_* @ r7 stub': (14.6, 60.4),
    'inhib_acetate, inhib_isobutanol @ r6 stub': (109.8, 36.4),
    'inhib_acetate, inhib_ethanol @ r17 stub': (131.0, 31.1),
    'inhib_* @ r10 decay edge': (31.0, 34.4),
}


def _setup_rcparams():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'pdf.fonttype': 42,   # TrueType -> editable text in the PDF
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
        'text.color': C_TEXT,
    })


# --- drawing helpers (axes coordinates are millimetres) --------------------

def _box(ax, cx, cy, w, h, lines, fc='white', ec=C_FLUX, lw=0.8,
         fs=FS_SPECIES, bold_first=True, rounding=1.6, zorder=4,
         sub_fs=None, text_color=C_TEXT):
    """Rounded species/process box centered at (cx, cy); returns its geometry."""
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=f'round,pad=0,rounding_size={rounding}',
        fc=fc, ec=ec, lw=lw, zorder=zorder))
    sub_fs = FS_ENZ if sub_fs is None else sub_fs
    if len(lines) == 1:
        ax.text(cx, cy, lines[0], ha='center', va='center', fontsize=fs,
                fontweight='bold' if bold_first else 'normal',
                color=text_color, zorder=zorder + 1)
    else:
        ax.text(cx, cy + h * 0.16, lines[0], ha='center', va='center',
                fontsize=fs, fontweight='bold' if bold_first else 'normal',
                color=text_color, zorder=zorder + 1)
        ax.text(cx, cy - h * 0.22, lines[1], ha='center', va='center',
                fontsize=sub_fs, style='italic', color=C_MUTED,
                zorder=zorder + 1)
    return dict(cx=cx, cy=cy, w=w, h=h,
                left=(cx - w / 2, cy), right=(cx + w / 2, cy),
                top=(cx, cy + h / 2), bottom=(cx, cy - h / 2))


def _rounded_path(pts, r=2.5):
    """Polyline through pts with rounded corners (quadratic Beziers)."""
    pts = [np.asarray(p, float) for p in pts]
    verts, codes = [pts[0]], [Path.MOVETO]
    for i in range(1, len(pts) - 1):
        p0, p1, p2 = pts[i - 1], pts[i], pts[i + 1]
        d0, d1 = p1 - p0, p2 - p1
        l0, l1 = np.hypot(*d0), np.hypot(*d1)
        rr = min(r, l0 / 2, l1 / 2)
        a = p1 - d0 / l0 * rr
        b = p1 + d1 / l1 * rr
        verts += [a, p1, b]
        codes += [Path.LINETO, Path.CURVE3, Path.CURVE3]
    verts.append(pts[-1])
    codes.append(Path.LINETO)
    return Path(verts, codes)


def _arrow(ax, pts, color=C_FLUX, lw=1.1, reversible=False, ls='-',
           zorder=3, mutation=6.5, corner_r=3.0):
    style = '<|-|>' if reversible else '-|>'
    if len(pts) == 2:
        path = Path([pts[0], pts[1]], [Path.MOVETO, Path.LINETO])
    else:
        path = _rounded_path(pts, r=corner_r)
    ax.add_patch(FancyArrowPatch(
        path=path, arrowstyle=style, mutation_scale=mutation,
        lw=lw, color=color, linestyle=ls, shrinkA=0, shrinkB=0,
        capstyle='round', joinstyle='round', zorder=zorder,
        fill=True))


def _curve(ax, p0, p1, rad, color, lw=0.9, ls='-', arrow=True, zorder=3,
           mutation=5.5):
    """Curved control edge (arc3) from p0 to p1."""
    ax.add_patch(FancyArrowPatch(
        p0, p1, connectionstyle=f'arc3,rad={rad}',
        arrowstyle='-|>' if arrow else '-', mutation_scale=mutation,
        lw=lw, color=color, linestyle=ls, shrinkA=0, shrinkB=0,
        capstyle='round', zorder=zorder))


def _tbar(ax, tip, tail, color, lw=0.9, bar=1.8, zorder=5, ls='-'):
    """Inhibition stub: line from tail to tip ending in a perpendicular bar."""
    tip, tail = np.asarray(tip, float), np.asarray(tail, float)
    d = tip - tail
    d = d / np.hypot(*d)
    n = np.array([-d[1], d[0]])
    a, b = tip + n * bar / 2, tip - n * bar / 2
    ax.add_line(Line2D([tail[0], tip[0]], [tail[1], tip[1]], color=color,
                       lw=lw, ls=ls, solid_capstyle='round', zorder=zorder))
    ax.add_line(Line2D([a[0], b[0]], [a[1], b[1]], color=color, lw=lw,
                       solid_capstyle='round', zorder=zorder))


def _rxn_marker(ax, x, y, rid, enzyme=None, enz_dxy=(0, -4.0), r=2.1,
                zorder=6):
    ax.add_patch(Circle((x, y), r, fc='white', ec=C_FLUX, lw=0.7,
                        zorder=zorder))
    ax.text(x, y, rid, ha='center', va='center', fontsize=FS_RXN,
            fontweight='bold', zorder=zorder + 1)
    if enzyme:
        ax.text(x + enz_dxy[0], y + enz_dxy[1], enzyme, ha='center',
                va='center', fontsize=FS_ENZ, style='italic', color=C_MUTED,
                zorder=zorder + 1)


def _badge(ax, x, y, text, fc, tc='white', ec='none', fs=4.4, w=6.0, h=3.4,
           zorder=6):
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0,rounding_size=1.4',
        fc=fc, ec=ec, lw=0.6, zorder=zorder))
    ax.text(x, y, text, ha='center', va='center', fontsize=fs,
            fontweight='bold', color=tc, zorder=zorder + 1)


def _tag(ax, x, y, text, color=C_MUTED, ha='center', fs=FS_TAG, zorder=6,
         style='italic'):
    ax.text(x, y, text, ha=ha, va='center', fontsize=fs, color=color,
            style=style, zorder=zorder)


def _knob(ax, x, y, r=1.05, pointer_angle=45., fc=C_LEVER, ec=C_TEXT, lw=0.35,
          zorder=8):
    """Rotary-knob icon centered at (x, y): a filled dial of radius `r` with a
    pointer, inside a 270-degree range arc (open at the bottom) that ends in
    min/max stops. Drawn from vector primitives (no font glyph, so the PDF
    text stays editable); kept smaller than the reaction markers so the two
    circles read apart."""
    sweep = np.deg2rad(np.linspace(-45., 225., 40))
    ax.add_line(Line2D(x + 1.5 * r * np.cos(sweep), y + 1.5 * r * np.sin(sweep),
                       color=ec, lw=1.3 * lw, solid_capstyle='round',
                       zorder=zorder))
    for t in sweep[[0, -1]]:                      # min / max stops
        ax.add_line(Line2D([x + 1.25 * r * np.cos(t), x + 1.75 * r * np.cos(t)],
                           [y + 1.25 * r * np.sin(t), y + 1.75 * r * np.sin(t)],
                           color=ec, lw=1.3 * lw, solid_capstyle='round',
                           zorder=zorder))
    ax.add_patch(Circle((x, y), r, fc=fc, ec=ec, lw=lw, zorder=zorder))
    t = np.deg2rad(pointer_angle)
    ax.add_line(Line2D([x + 0.15 * r * np.cos(t), x + 0.9 * r * np.cos(t)],
                       [y + 0.15 * r * np.sin(t), y + 0.9 * r * np.sin(t)],
                       color=ec, lw=2.0 * lw, solid_capstyle='round',
                       zorder=zorder + 1))


# --- the figure ------------------------------------------------------------

def draw_conceptual_diagram(save_dir=None, formats=('png', 'pdf'),
                            dpi=600, show=False, show_strain_levers=True,
                            show_process_controls=True,
                            filename='conceptual_diagram'):
    """Draw the conceptual reaction-network/controls diagram.

    Parameters
    ----------
    save_dir : str, optional
        Directory for the output files. Defaults to this module's directory.
    formats : tuple of str, optional
        File formats to save (``'png'``, ``'pdf'``, ``'svg'``, ...).
        Defaults to ``('png', 'pdf')``.
    dpi : int, optional
        Raster resolution for PNG output. Defaults to 600 (Nature
        Communications requires >= 300 dpi at final size).
    show : bool, optional
        Call ``plt.show()`` after saving. Defaults to False.
    show_strain_levers : bool, optional
        Mark the strain design variables (``STRAIN_LEVER_MARKS``: the
        strain-side decision variables of the isobutanol package's
        ``metabolic_split_12d`` kinetic-optimization preset) with rotary-knob
        icons. Defaults to True.
    show_process_controls : bool, optional
        Draw the process context: the process-control row above the reactor
        (the fed-batch glucose feeding and two-stage aeration panels), its
        glucose-spike, glucose sensing and O$_2$ supply connectors, the
        "Fed-batch bioreactor" frame they cross, and the feed legend entry.
        False gives the reaction network alone, frameless, on a
        correspondingly shorter figure; the O$_2$ badges stay (they mark
        the reactions whose rate is scaled by ``f_O2`` -- solid for r2, r5
        and r8, outlined for r7, whose ``anaerobic_growth_mult`` share runs
        without O$_2$).
        Defaults to True.
    filename : str, optional
        Output file stem. Defaults to ``'conceptual_diagram'``; pass another
        stem to keep a variant from overwriting the full figure.

    Returns
    -------
    tuple
        ``(fig, ax)`` of the drawn figure.
    """
    _setup_rcparams()
    # the legend's control-edge column always has four entries
    y_floor = -LEGEND_ROW_4_MM
    y_top = Y_TOP_MM if show_process_controls else Y_TOP_NO_CONTROLS_MM
    fig, ax = plt.subplots(figsize=(FIG_W_MM * MM, (y_top - y_floor) * MM))
    ax.set_xlim(0, FIG_W_MM)
    ax.set_ylim(y_floor, y_top)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    # === process-control row (top, outside the reactor) ===================
    if show_process_controls:
        ax.add_patch(FancyBboxPatch((4, 113), 82, 19,
                                    boxstyle='round,pad=0,rounding_size=2',
                                    fc=TINT_CTRL, ec=TINT_CTRL_EDGE, lw=0.9,
                                    zorder=2))
        ax.text(8, 128.2, 'Fed-batch glucose feeding (FeedSpike events)',
                fontsize=FS_PANEL, fontweight='bold', color='#144E7A',
                ha='left', va='center', zorder=3)
        ax.text(8, 119.8,
                'When $s_\\mathrm{glu}$ < threshold (10 g L$^{-1}$): spike '
                'concentrated\nfeed (600 g L$^{-1}$) to restore the target '
                '(100 g L$^{-1}$);\nat most max_n_glu_spikes spikes '
                '(default 5);\neach spike increases broth volume env.',
                fontsize=FS_NOTE, ha='left', va='center', color=C_TEXT,
                linespacing=1.35, zorder=3)

        ax.add_patch(FancyBboxPatch((90, 113), 86, 19,
                                    boxstyle='round,pad=0,rounding_size=2',
                                    fc=TINT_CTRL, ec=TINT_CTRL_EDGE, lw=0.9,
                                    zorder=2))
        ax.text(94, 128.2, 'Two-stage aeration control',
                fontsize=FS_PANEL, fontweight='bold', color='#144E7A',
                ha='left', va='center', zorder=3)
        ax.text(94, 120.2,
                'is_aerobic = 1 while $t$ < stage_1_max_time and '
                '$x$ < stage_1_max_x, then 0.\nf_O2 = is_aerobic × transfer '
                'fraction, capped so respiratory O$_2$ ≤ kLa·C$_{O_2}$*;\nscales '
                'r2, r5, r8 and aerobic growth r7 (anaerobic_growth_mult share '
                'ungated).',
                fontsize=FS_NOTE, ha='left', va='center', color=C_TEXT,
                linespacing=1.35, zorder=3)

    # === bioreactor frame + connectors from the process-control row ========
    # (the frame is the boundary the feed / sensing / O2 connectors cross, so
    # it goes with them: the network-only figure is frameless)
    if show_process_controls:
        ax.add_patch(FancyBboxPatch((2, 16), 176, 94,
                                    boxstyle='round,pad=0,rounding_size=3',
                                    fc='white', ec='#666666', lw=1.1,
                                    zorder=1))
        ax.text(174, 106.3, 'Fed-batch bioreactor',
                fontsize=FS_PANEL, fontweight='bold', color='#444444',
                ha='right', va='center', zorder=3)

        # feed arrow into the reactor + sensing line (the fed-batch loop)
        _arrow(ax, [(45, 113), (45, 104), (78, 104), (78, 101.8)],
               color=C_FEED, lw=1.2, corner_r=3)
        _tag(ax, 60, 106.2, 'glucose spike', color=C_FEED, fs=FS_TAG)
        ax.add_line(Line2D([95, 95, 82, 82], [101.8, 107, 107, 113],
                           color=C_FEED, lw=0.8, ls=(0, (1.2, 1.4)), zorder=3))
        _tag(ax, 89.5, 108.6, 'sense $s_\\mathrm{glu}$', color=C_FEED,
             fs=FS_TAG)

        # O2 supply drop from the aeration panel (left of the reactor caption)
        _arrow(ax, [(120, 113), (120, 103.5)], color=C_O2, lw=1.0,
               ls=(0, (2.4, 1.6)), mutation=5.5)
        _tag(ax, 122.3, 108, 'O$_2$', color=C_O2_TEXT, fs=FS_TAG, ha='left')

    # === species boxes =====================================================
    glu = _box(ax, 86, 98, 24, 7, ['Glucose'], fc=TINT_SUBSTRATE)
    pyr = _box(ax, 86, 76, 24, 7, ['Pyruvate'])
    ald = _box(ax, 86, 54, 26, 7, ['Acetaldehyde'])
    eth = _box(ax, 86, 30, 24, 7, ['Ethanol'], fc=TINT_PRODUCT)
    ace = _box(ax, 48, 54, 20, 7, ['Acetate'])
    tca = _box(ax, 27, 76, 30, 10, ['TCA cycle &', 'respiration'],
               fc='#F4F4F4', ec='#8A8A8A', sub_fs=FS_SPECIES - 0.5)
    # TCA sub-line should not be italic/muted; redraw label cleanly
    # (the helper renders line 2 italic; overwrite with a matching label)
    tca_txt = [t for t in ax.texts if t.get_text() == 'respiration'][-1]
    tca_txt.set_style('normal')
    tca_txt.set_color(C_TEXT)
    tca_txt.set_fontweight('bold')
    tca_txt.set_fontsize(FS_SPECIES)

    # CO2 vent from TCA
    _arrow(ax, [(27, 81.3), (27, 87.5)], color='#8A8A8A', lw=0.8, mutation=5)
    _tag(ax, 27, 89.6, 'CO$_2$', color=C_MUTED, fs=FS_TAG)

    # === engineered isobutanol pathway panel (right) =======================
    # (the panel reaches left of the species column so that r13, the
    # pathway's entry reaction, sits inside it)
    ax.add_patch(FancyBboxPatch((119, 22), 57, 68,
                                boxstyle='round,pad=0,rounding_size=2',
                                fc=TINT_IBO, ec=TINT_IBO_EDGE, lw=0.9,
                                zorder=2))
    ax.text(147.5, 86.0, 'Engineered isobutanol\npathway (r13–r17)',
            fontsize=FS_PANEL - 0.3, fontweight='bold', color='#B04A00',
            ha='center', va='center', linespacing=1.2, zorder=3)

    al = _box(ax, 152, 78, 34, 6.4, ['AL'])
    dhi = _box(ax, 152, 65, 34, 6.4, ['DHIV'])
    kiv = _box(ax, 152, 52, 34, 6.4, ['KIV'])
    iald = _box(ax, 152, 39, 34, 6.4, ['Isobutyraldehyde'],
                fs=FS_SPECIES - 0.8)
    ibo = _box(ax, 152, 26, 34, 6.4, ['Isobutanol'], fc=TINT_PRODUCT)

    # === biomass / physiological-state panel (bottom left) =================
    ax.add_patch(FancyBboxPatch((8, 18), 60, 23.5,
                                boxstyle='round,pad=0,rounding_size=2',
                                fc=TINT_BM, ec=TINT_BM_EDGE, lw=0.9,
                                zorder=2))
    ax.text(38, 38.2, 'Biomass $x$ & physiological state',
            fontsize=FS_PANEL - 0.3, fontweight='bold', color='#1F7A5C',
            ha='center', va='center', zorder=3)

    xa = _box(ax, 24, 29, 21, 6, ['$X_a$ (active)'], fs=FS_SPECIES - 0.5)
    xac = _box(ax, 54, 29, 22, 6, ['$X_\\mathrm{AcDH}$'], fs=FS_SPECIES - 0.5)
    _arrow(ax, [(34.5, 29), (43, 29)], color=C_FLUX, lw=1.1, mutation=5.5)
    _rxn_marker(ax, 38.7, 29, 'r9', enzyme=None)
    # decay r10 / r11 -- angled down-and-outward so the standard marker
    # circles fit on the runs without hiding the arrowheads (the region
    # directly below r9 carries its induction/repression annotations)
    _arrow(ax, [(22, 26), (16.6, 20.2)], color=C_FLUX, lw=1.1)
    _rxn_marker(ax, 19.3, 23.1, 'r10', enzyme=None)
    _tag(ax, 15.3, 19.3, '$\\varnothing$', color=C_MUTED, fs=FS_TAG,
         style='normal')
    _arrow(ax, [(56, 26), (61.4, 20.2)], color=C_FLUX, lw=1.1)
    _rxn_marker(ax, 58.7, 23.1, 'r11', enzyme=None)
    _tag(ax, 62.7, 19.3, '$\\varnothing$', color=C_MUTED, fs=FS_TAG,
         style='normal')

    # === central catabolic backbone ========================================
    # r1 glycolysis
    _arrow(ax, [glu['bottom'], pyr['top']], lw=1.1)
    _rxn_marker(ax, 86, 87.6, 'r1', enzyme='glycolysis', enz_dxy=(9.5, 0))
    _tag(ax, 80.6, 84.3, '+NADH', ha='right')
    # r3 Pdc1
    _arrow(ax, [pyr['bottom'], ald['top']], lw=1.1)
    _rxn_marker(ax, 86, 65.4, 'r3', enzyme='Pdc1', enz_dxy=(6.9, 0))
    _tag(ax, 81.0, 62.2, 'CO$_2$', ha='right')
    # r6 ADH (reversible)
    _arrow(ax, [ald['bottom'], eth['top']], lw=1.1, reversible=True)
    _rxn_marker(ax, 86, 42.2, 'r6', enzyme='Adh1', enz_dxy=(6.9, 0))
    _tag(ax, 81.0, 39.0, '–NADH', ha='right')

    # r2 pyruvate -> TCA (the TCA box shares pyruvate's row, so r2 runs
    # horizontal like r4 and r13)
    _arrow(ax, [pyr['left'], tca['right']], lw=1.1)
    _rxn_marker(ax, 58, 76, 'r2', enzyme=None)
    _badge(ax, 51.5, 79.8, 'O$_2$', C_O2)
    _tag(ax, 58, 72.3, '+NADH', fs=FS_ENZ)

    # r4 acetaldehyde -> acetate (needs AcDH machinery)
    _arrow(ax, [ald['left'], ace['right']], lw=1.1)
    _rxn_marker(ax, 65.5, 54, 'r4', enzyme='Ald6', enz_dxy=(0, 7.5))
    _badge(ax, 65.5, 58.2, 'AcDH', '#FFFFFF', tc='#1F7A5C',
           ec='#1F7A5C', w=8.6)
    _tag(ax, 70, 49.3, '+NADH', fs=FS_ENZ, ha='left')

    # r5 acetate -> TCA
    _arrow(ax, [(44, 57.6), (33.8, 70.9)], lw=1.1)
    _rxn_marker(ax, 38.7, 64.6, 'r5', enzyme='Acs2', enz_dxy=(5.7, 0))
    _badge(ax, 33.0, 64.2, 'O$_2$', C_O2)
    _tag(ax, 42.8, 61.4, '+NADH', fs=FS_ENZ, ha='left')

    # r7 glucose -> biomass (left-margin route)
    _arrow(ax, [(74, 98), (8.5, 98), (8.5, 47), (20, 41.7)], lw=1.1,
           corner_r=4)
    _rxn_marker(ax, 8.5, 58, 'r7', enzyme=None)
    _tag(ax, 11.2, 64.5, 'growth', color=C_MUTED, fs=FS_ENZ, ha='left')
    _tag(ax, 11.2, 53.8, '+NADH, CO$_2$', fs=FS_ENZ, ha='left')
    # only the (1 - anaerobic_growth_mult) share of r7 is scaled by f_O2, so
    # it gets the outlined (partial) badge rather than the solid one
    _badge(ax, 14.3, 49.9, 'O$_2$', '#FFFFFF', tc=C_O2_TEXT, ec=C_O2)

    # r8 acetate -> biomass
    _arrow(ax, [ace['bottom'], (48, 41.7)], lw=1.1)
    _rxn_marker(ax, 48, 46.4, 'r8', enzyme=None)
    _badge(ax, 55.3, 46.4, 'O$_2$', C_O2)
    _tag(ax, 43, 43.5, '+NADH, CO$_2$', fs=FS_ENZ, ha='right')

    # === engineered pathway reactions ======================================
    _arrow(ax, [pyr['right'], (135, 76)], lw=1.1)
    _rxn_marker(ax, 125, 76, 'r13', enzyme='Ilv2+Ilv6', enz_dxy=(0.4, 3.9))
    _tag(ax, 129.5, 73.2, 'CO$_2$', ha='left', fs=FS_ENZ)
    _arrow(ax, [al['bottom'], dhi['top']], lw=1.1)
    _rxn_marker(ax, 152, 71.5, 'r14', enzyme='Ilv5', enz_dxy=(8.6, 0))
    _tag(ax, 144.7, 71.5, '–NADPH', ha='right', fs=FS_ENZ)
    _arrow(ax, [dhi['bottom'], kiv['top']], lw=1.1)
    _rxn_marker(ax, 152, 58.5, 'r15', enzyme='Ilv3', enz_dxy=(8.6, 0))
    # r16 Aro10 decarboxylase (irreversible), r17 Adh6 reductase (reversible,
    # like r6 -- the two alcohol dehydrogenases are drawn alike on purpose)
    _arrow(ax, [kiv['bottom'], iald['top']], lw=1.1)
    _rxn_marker(ax, 152, 45.5, 'r16', enzyme='Aro10', enz_dxy=(8.9, 0))
    _tag(ax, 144.7, 45.5, 'CO$_2$', ha='right', fs=FS_ENZ)
    _arrow(ax, [iald['bottom'], ibo['top']], lw=1.1, reversible=True)
    _rxn_marker(ax, 152, 32.5, 'r17', enzyme='Adh6', enz_dxy=(8.6, 0))
    # (raised off the marker's row to leave room for r17's inhibition stub)
    _tag(ax, 144.7, 34.1, '–NADPH', ha='right', fs=FS_ENZ)

    # === control edges =====================================================
    # acetaldehyde overflow signal activates the high-capacity glycolytic
    # term of r1; routed through the free corridor right of the pyruvate box
    # (crossing only the r13 arrow) so the edge visibly runs acetaldehyde->r1
    _arrow(ax, [(99.2, 55), (105, 61), (105, 79), (88.5, 85.8)],
           color=C_ACT, lw=1.0, ls=(0, (3.2, 1.8)), corner_r=5,
           mutation=5.5)
    _tag(ax, 107, 70, 'overflow\nsignal (+)', color=C_ACT, fs=FS_TAG,
         ha='left')

    # product inhibition stubs (EtOH, acetate, isobutanol -| r1, r4, r7)
    _tbar(ax, (83.3, 87.6), (74.5, 87.6), C_INHIB)
    _tag(ax, 73.5, 87.6, 'ethanol·acetate·isobutanol', color=C_INHIB, ha='right',
         fs=FS_TAG)
    _tbar(ax, (65.5, 51.9), (65.5, 45.5), C_INHIB)
    _tag(ax, 68.5, 43.7, 'ethanol·acetate·isobutanol', color=C_INHIB,
         fs=FS_TAG)
    _tbar(ax, (11.3, 58), (17.5, 58), C_INHIB)
    _tag(ax, 18.3, 58, 'ethanol·acetate·\nisobutanol', color=C_INHIB, ha='left',
         fs=FS_TAG)

    # the two alcohol dehydrogenases carry only the cross-product
    # exponentials (their own product acts through the reversible law):
    # r6 exp(-k_6ia*Ace)*exp(-k_6ii*iBuOH), r17 exp(-k_17ia*Ace)*exp(-k_17ie*EtOH)
    _tbar(ax, (88.1, 40.5), (93.4, 37.0), C_INHIB)
    _tag(ax, 94.2, 36.4, 'acetate·isobutanol', color=C_INHIB, ha='left',
         fs=FS_TAG)
    _tbar(ax, (149.6, 31.1), (145.4, 31.1), C_INHIB)
    _tag(ax, 144.6, 31.1, 'ethanol·acetate', color=C_INHIB, ha='right',
         fs=FS_TAG)

    # product-accelerated decay: above its threshold (P_10e/a/i) each product
    # multiplies r10 by exp(k_10i*(C - P_10)) -- an activating edge in the
    # product color, routed down the panel margin left of the X_a box
    _arrow(ax, [(10.8, 32.7), (10.8, 23.1), (17.0, 23.1)], color=C_INHIB,
           lw=0.9, ls=(0, (3.2, 1.8)), corner_r=2.5, mutation=5, zorder=5)
    _tag(ax, 9.6, 34.3, 'ethanol·acetate·isobutanol', color=C_INHIB, ha='left',
         fs=FS_TAG)

    # glucose repression stubs (glucose -| r2, r5, r8, r9)
    _tbar(ax, (61.5, 77.7), (64.5, 80.9), C_REPR)
    _tag(ax, 65.4, 82.1, 'glucose', color=C_REPR, ha='left', fs=FS_TAG)
    _tbar(ax, (41.4, 66.7), (45.5, 69.3), C_REPR)
    _tag(ax, 46.4, 70.4, 'glucose', color=C_REPR, ha='left', fs=FS_TAG)
    _tbar(ax, (45.3, 46.4), (41.2, 46.4), C_REPR)
    _tag(ax, 40.2, 46.4, 'glucose', color=C_REPR, ha='right', fs=FS_TAG)
    # r9 is glucose-INDUCED (and EtOH-induced) at low glucose and repressed
    # only at high glucose (the 1/(K_9i*s_glu + 1) factor), so it gets a
    # dual annotation rather than the plain repression stub of r2/r5/r8
    _tbar(ax, (40.2, 27.2), (43.2, 23.8), C_REPR)
    _tag(ax, 43.8, 22.7, 'glucose (high)', color=C_REPR, ha='left', fs=FS_TAG)
    _arrow(ax, [(34.2, 23.8), (37.2, 27.2)], color=C_ACT, lw=0.9,
           ls=(0, (3.2, 1.8)), mutation=5, zorder=5)
    _tag(ax, 33.6, 22.7, 'glucose ·\nethanol', color=C_ACT, ha='right', fs=FS_TAG)

    # === strain design variables ===========================================
    if show_strain_levers:
        for x, y in STRAIN_LEVER_MARKS.values():
            _knob(ax, x, y)

    # === legend strip ======================================================
    # three columns: control edges | badges and the knob | the feed entry.
    # Column widths are measured from the rendered labels, so the columns sit
    # LEGEND_COL_GAP_MM apart and the strip is centred under the figure
    # whichever entries are present.
    def _leg_arrow(x, y, color):
        _arrow(ax, [(x, y), (x + 7, y)], color=color, lw=1.0,
               ls=(0, (3, 1.8)), mutation=5.5, zorder=5)

    renderer = fig.canvas.get_renderer()

    def _label_w(label):
        t = ax.text(0, 0, label, fontsize=FS_LEGEND)
        bb = t.get_window_extent(renderer)
        t.remove()
        (x0, _), (x1, _) = ax.transData.inverted().transform(
            [(bb.x0, 0), (bb.x1, 0)])
        return x1 - x0

    edge_labels = (
        'product inhibition',
        'product-accelerated decay',
        'glucose repression',
        'activation',
    )
    o2_label = 'O$_2$-dependent'
    o2_part_label = 'partly O$_2$-dependent'
    acdh_label = 'requires AcDH machinery'
    knob_label = 'strain design variable'
    feed_label = 'fed-batch feed; dotted = sensing'

    # offsets within a column: edge glyphs are 7 mm with the label at +9.5;
    # the badge column is measured from the wide AcDH badge's left edge
    # (badges centred at +4.3, O2 / knob labels at +8.8, AcDH label at +10.3);
    # the feed glyph is 7.5 mm with its label at +9
    widths = [9.5 + max(_label_w(label) for label in edge_labels),
              max(8.8 + _label_w(o2_label), 8.8 + _label_w(o2_part_label),
                  10.3 + _label_w(acdh_label),
                  8.8 + _label_w(knob_label) if show_strain_levers else 0.)]
    if show_process_controls:
        widths.append(9. + _label_w(feed_label))
    total_w = sum(widths) + LEGEND_COL_GAP_MM * (len(widths) - 1)
    xa = (FIG_W_MM - total_w) / 2
    xb = xa + widths[0] + LEGEND_COL_GAP_MM + 4.3
    xc = xb - 4.3 + widths[1] + LEGEND_COL_GAP_MM

    ax.add_patch(FancyBboxPatch((xa - LEGEND_PAD_MM, 1.5 + y_floor),
                                total_w + 2 * LEGEND_PAD_MM, 12 - y_floor,
                                boxstyle='round,pad=0,rounding_size=1.5',
                                fc='#FAFAFA', ec='#CCCCCC', lw=0.7,
                                zorder=1))

    # even 4.2 mm pitch: the badges are 3.4 mm tall, so anything tighter
    # makes the stacked O2 / AcDH badges touch
    y1, y2, y3, y4 = 11.3, 7.1, 2.9, -1.3

    _tbar(ax, (xa + 7, y1), (xa, y1), C_INHIB, lw=1.0)
    _leg_arrow(xa, y2, C_INHIB)
    _tbar(ax, (xa + 7, y3), (xa, y3), C_REPR, lw=1.0)
    _leg_arrow(xa, y4, C_ACT)
    for y, label in zip((y1, y2, y3, y4), edge_labels):
        ax.text(xa + 9.5, y, label, fontsize=FS_LEGEND, va='center',
                zorder=5)

    _badge(ax, xb, y1, 'O$_2$', C_O2)
    ax.text(xb + 4.5, y1, o2_label, fontsize=FS_LEGEND, va='center',
            zorder=5)
    _badge(ax, xb, y2, 'O$_2$', '#FFFFFF', tc=C_O2_TEXT, ec=C_O2)
    ax.text(xb + 4.5, y2, o2_part_label, fontsize=FS_LEGEND, va='center',
            zorder=5)
    _badge(ax, xb, y3, 'AcDH', '#FFFFFF', tc='#1F7A5C', ec='#1F7A5C',
           w=8.6)
    ax.text(xb + 6, y3, acdh_label, fontsize=FS_LEGEND, va='center',
            zorder=5)
    if show_strain_levers:
        _knob(ax, xb, y4)
        ax.text(xb + 4.5, y4, knob_label, fontsize=FS_LEGEND, va='center',
                zorder=5)

    if show_process_controls:
        _arrow(ax, [(xc, y1), (xc + 4.5, y1)], color=C_FEED, lw=1.0,
               mutation=5.5, zorder=5)
        ax.add_line(Line2D([xc + 5.3, xc + 7.5], [y1, y1], color=C_FEED,
                           lw=0.9, ls=(0, (1.2, 1.4)), zorder=5))
        ax.text(xc + 9, y1, feed_label, fontsize=FS_LEGEND, va='center',
                zorder=5)

    # === save ==============================================================
    if save_dir is None:
        save_dir = os.path.dirname(os.path.abspath(__file__))
    for fmt in formats:
        out = os.path.join(save_dir, f'{filename}.{fmt}')
        fig.savefig(out, dpi=dpi, facecolor='white')
        print(f'Saved {out}')
    if show:
        plt.show()
    return fig, ax


if __name__ == '__main__':
    draw_conceptual_diagram()
