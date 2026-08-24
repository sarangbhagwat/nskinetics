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

The figure shows the reaction network (r1-r11, r13-r16) and the control
structure — exponential product inhibition, glucose repression, the
acetaldehyde overflow signal, aerobic gating, the AcDH physiological-state
machinery, and the fed-batch glucose-spike loop — sized to Nature
Communications figure specifications (180 mm double-column width, Arial,
5-7 pt text, vector PDF with editable text plus a 600 dpi PNG).

The network and controls are curated directly from
``s_cerevisiae_ferm_fb_inhib_mod_ibo_antimony.txt``; continuous-mode dilution
terms (``s_glu_in``, ``*_out``; D = 0 in fed-batch use) are noted but not
drawn as edges.

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
FIG_W_MM, FIG_H_MM = 180., 134.   # double-column width; height < 170 mm cap

# --- Okabe-Ito palette (colorblind-safe), assigned by control job ----------
C_FLUX = '#3A3A3A'      # mass/reaction flux
C_INHIB = '#D55E00'     # product inhibition (vermillion)
C_REPR = '#CC79A7'      # glucose repression (reddish purple)
C_ACT = '#009E73'       # activation / overflow signal (bluish green)
C_FEED = '#0072B2'      # fed-batch feed & sensing (blue)
C_O2 = '#56B4E9'        # aerobic gating badge (sky blue)
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


# --- the figure ------------------------------------------------------------

def draw_conceptual_diagram(save_dir=None, formats=('png', 'pdf'),
                            dpi=600, show=False):
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

    Returns
    -------
    tuple
        ``(fig, ax)`` of the drawn figure.
    """
    _setup_rcparams()
    fig, ax = plt.subplots(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
    ax.set_xlim(0, FIG_W_MM)
    ax.set_ylim(0, FIG_H_MM)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    # === process-control row (top, outside the reactor) ===================
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
            '$x$ < stage_1_max_x, then 0.\nGates O$_2$-dependent fluxes '
            '(r2, r5, r8); growth r7 scales by\nanaerobic_growth_mult when '
            'anaerobic.',
            fontsize=FS_NOTE, ha='left', va='center', color=C_TEXT,
            linespacing=1.35, zorder=3)

    # === bioreactor frame ==================================================
    ax.add_patch(FancyBboxPatch((2, 16), 176, 94,
                                boxstyle='round,pad=0,rounding_size=3',
                                fc='white', ec='#666666', lw=1.1, zorder=1))
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
    _tag(ax, 122.3, 108, 'O$_2$', color='#2E7FB0', fs=FS_TAG, ha='left')

    # === species boxes =====================================================
    glu = _box(ax, 86, 98, 24, 7, ['Glucose'], fc=TINT_SUBSTRATE)
    pyr = _box(ax, 86, 76, 24, 7, ['Pyruvate'])
    ald = _box(ax, 86, 54, 26, 7, ['Acetaldehyde'])
    eth = _box(ax, 86, 30, 24, 7, ['Ethanol'], fc=TINT_PRODUCT)
    ace = _box(ax, 48, 54, 20, 7, ['Acetate'])
    tca = _box(ax, 27, 77, 30, 10, ['TCA cycle &', 'respiration'],
               fc='#F4F4F4', ec='#8A8A8A', sub_fs=FS_SPECIES - 0.5)
    # TCA sub-line should not be italic/muted; redraw label cleanly
    # (the helper renders line 2 italic; overwrite with a matching label)
    tca_txt = [t for t in ax.texts if t.get_text() == 'respiration'][-1]
    tca_txt.set_style('normal')
    tca_txt.set_color(C_TEXT)
    tca_txt.set_fontweight('bold')
    tca_txt.set_fontsize(FS_SPECIES)

    # CO2 vent from TCA
    _arrow(ax, [(27, 82.3), (27, 88.5)], color='#8A8A8A', lw=0.8, mutation=5)
    _tag(ax, 27, 90.6, 'CO$_2$', color=C_MUTED, fs=FS_TAG)

    # === engineered isobutanol pathway panel (right) =======================
    ax.add_patch(FancyBboxPatch((128, 22), 48, 68,
                                boxstyle='round,pad=0,rounding_size=2',
                                fc=TINT_IBO, ec=TINT_IBO_EDGE, lw=0.9,
                                zorder=2))
    ax.text(152, 86.0, 'Engineered isobutanol\npathway (r13–r16)',
            fontsize=FS_PANEL - 0.3, fontweight='bold', color='#B04A00',
            ha='center', va='center', linespacing=1.2, zorder=3)

    al = _box(ax, 152, 76, 34, 7, ['AL'])
    dhi = _box(ax, 152, 60.5, 34, 7, ['DHIV'])
    kiv = _box(ax, 152, 45, 34, 7, ['KIV'])
    ibo = _box(ax, 152, 28.5, 34, 7, ['Isobutanol'], fc=TINT_PRODUCT)

    # === biomass / physiological-state panel (bottom left) =================
    ax.add_patch(FancyBboxPatch((8, 18), 60, 23.5,
                                boxstyle='round,pad=0,rounding_size=2',
                                fc=TINT_BM, ec=TINT_BM_EDGE, lw=0.9,
                                zorder=2))
    ax.text(38, 37, 'Biomass $x$ & physiological state',
            fontsize=FS_PANEL - 0.3, fontweight='bold', color='#1F7A5C',
            ha='center', va='center', zorder=3)

    xa = _box(ax, 24, 29, 21, 6, ['$X_a$ (active)'], fs=FS_SPECIES - 0.5)
    xac = _box(ax, 54, 29, 22, 6, ['$X_\\mathrm{AcDH}$'], fs=FS_SPECIES - 0.5)
    _arrow(ax, [(34.5, 29), (43, 29)], color=C_FLUX, lw=0.9, mutation=5.5)
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
    _arrow(ax, [glu['bottom'], pyr['top']], lw=1.3)
    _rxn_marker(ax, 86, 87.6, 'r1', enzyme='glycolysis', enz_dxy=(9.5, 0))
    _tag(ax, 80.6, 84.3, '+NADH', ha='right')
    # r3 PDC
    _arrow(ax, [pyr['bottom'], ald['top']], lw=1.3)
    _rxn_marker(ax, 86, 65.4, 'r3', enzyme='PDC', enz_dxy=(6.6, 0))
    _tag(ax, 81.0, 62.2, 'CO$_2$', ha='right')
    # r6 ADH (reversible)
    _arrow(ax, [ald['bottom'], eth['top']], lw=1.3, reversible=True)
    _rxn_marker(ax, 86, 42.2, 'r6', enzyme='ADH', enz_dxy=(6.6, 0))
    _tag(ax, 81.0, 39.0, '–NADH', ha='right')

    # r2 pyruvate -> TCA
    _arrow(ax, [pyr['left'], tca['right']], lw=1.1)
    _rxn_marker(ax, 58, 76.5, 'r2', enzyme=None)
    _badge(ax, 51.5, 80.3, 'O$_2$', C_O2)
    _tag(ax, 58, 72.8, '+NADH', fs=FS_ENZ)

    # r4 acetaldehyde -> acetate (needs AcDH machinery)
    _arrow(ax, [ald['left'], ace['right']], lw=1.1)
    _rxn_marker(ax, 65.5, 54, 'r4', enzyme=None)
    _badge(ax, 65.5, 58.2, 'AcDH', '#FFFFFF', tc='#1F7A5C',
           ec='#1F7A5C', w=8.6)
    _tag(ax, 70, 49.3, '+NADH', fs=FS_ENZ, ha='left')

    # r5 acetate -> TCA
    _arrow(ax, [(44, 57.6), (33, 71.9)], lw=1.1)
    _rxn_marker(ax, 38.7, 64.6, 'r5', enzyme=None)
    _badge(ax, 31.6, 62.6, 'O$_2$', C_O2)
    _tag(ax, 44.5, 63.5, '+NADH', fs=FS_ENZ, ha='left')

    # r7 glucose -> biomass (left-margin route)
    _arrow(ax, [(74, 98), (8.5, 98), (8.5, 47), (20, 41.7)], lw=1.1,
           corner_r=4)
    _rxn_marker(ax, 8.5, 58, 'r7', enzyme=None)
    _tag(ax, 11.2, 64.5, 'growth', color=C_MUTED, fs=FS_ENZ, ha='left')
    _tag(ax, 11.2, 53.8, '+NADH, CO$_2$', fs=FS_ENZ, ha='left')

    # r8 acetate -> biomass
    _arrow(ax, [ace['bottom'], (48, 41.7)], lw=1.1)
    _rxn_marker(ax, 48, 46.4, 'r8', enzyme=None)
    _badge(ax, 55.3, 46.4, 'O$_2$', C_O2)
    _tag(ax, 43, 43.5, '+NADH, CO$_2$', fs=FS_ENZ, ha='right')

    # === engineered pathway reactions ======================================
    _arrow(ax, [pyr['right'], (135, 76)], lw=1.1)
    _rxn_marker(ax, 113, 76, 'r13', enzyme='ALS', enz_dxy=(0, 3.9))
    _arrow(ax, [al['bottom'], dhi['top']], lw=1.1)
    _rxn_marker(ax, 152, 68.3, 'r14', enzyme='KARI', enz_dxy=(8.6, 0))
    _arrow(ax, [dhi['bottom'], kiv['top']], lw=1.1)
    _rxn_marker(ax, 152, 52.8, 'r15', enzyme='DHAD', enz_dxy=(8.7, 0))
    _arrow(ax, [kiv['bottom'], ibo['top']], lw=1.1, reversible=True)
    _rxn_marker(ax, 152, 36.8, 'r16', enzyme='KDC+ADH', enz_dxy=(10.6, 0))
    _tag(ax, 144.7, 36.8, '–NADH', ha='right', fs=FS_ENZ)

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
    _tag(ax, 73.5, 87.6, 'EtOH·Ace·iBuOH', color=C_INHIB, ha='right',
         fs=FS_TAG)
    _tbar(ax, (65.5, 51.9), (65.5, 45.5), C_INHIB)
    _tag(ax, 66.5, 43.7, 'EtOH·Ace·iBuOH', color=C_INHIB, fs=FS_TAG)
    _tbar(ax, (11.3, 58), (17.5, 58), C_INHIB)
    _tag(ax, 18.3, 58, 'EtOH·Ace·iBuOH', color=C_INHIB, ha='left',
         fs=FS_TAG)

    # glucose repression stubs (glucose -| r2, r5, r8, r9)
    _tbar(ax, (61.5, 78.2), (64.5, 81.4), C_REPR)
    _tag(ax, 65.4, 82.6, 'glc', color=C_REPR, ha='left', fs=FS_TAG)
    _tbar(ax, (41.4, 66.7), (45.5, 69.3), C_REPR)
    _tag(ax, 46.4, 70.4, 'glc', color=C_REPR, ha='left', fs=FS_TAG)
    _tbar(ax, (45.3, 46.4), (41.2, 46.4), C_REPR)
    _tag(ax, 40.2, 46.4, 'glc', color=C_REPR, ha='right', fs=FS_TAG)
    # r9 is glucose-INDUCED (and EtOH-induced) at low glucose and repressed
    # only at high glucose (the 1/(K_9i*s_glu + 1) factor), so it gets a
    # dual annotation rather than the plain repression stub of r2/r5/r8
    _tbar(ax, (40.2, 27.2), (43.2, 23.8), C_REPR)
    _tag(ax, 43.8, 22.7, 'glc (high)', color=C_REPR, ha='left', fs=FS_TAG)
    _arrow(ax, [(34.2, 23.8), (37.2, 27.2)], color=C_ACT, lw=0.9,
           ls=(0, (3.2, 1.8)), mutation=5, zorder=5)
    _tag(ax, 33.6, 22.7, 'glc · EtOH', color=C_ACT, ha='right', fs=FS_TAG)

    # === legend strip ======================================================
    ax.add_patch(FancyBboxPatch((2, 1.5), 176, 12,
                                boxstyle='round,pad=0,rounding_size=1.5',
                                fc='#FAFAFA', ec='#CCCCCC', lw=0.7,
                                zorder=1))

    def _leg_arrow(x, y, color, reversible=False, ls='-', dashed=False):
        _arrow(ax, [(x, y), (x + 7, y)], color=color, lw=1.0,
               reversible=reversible, ls=(0, (3, 1.8)) if dashed else '-',
               mutation=5.5, zorder=5)

    y1, y2, y3 = 11.3, 7.6, 4.3
    _leg_arrow(6, y1, C_FLUX)
    ax.text(15, y1, 'reaction flux (mass basis)', fontsize=FS_LEGEND,
            va='center', zorder=5)
    _leg_arrow(6, y2, C_FLUX, reversible=True)
    ax.text(15, y2, 'reversible (r6, r16)', fontsize=FS_LEGEND, va='center',
            zorder=5)
    _arrow(ax, [(6, y3), (10.5, y3)], color=C_FEED, lw=1.0, mutation=5.5,
           zorder=5)
    ax.add_line(Line2D([11.3, 13.5], [y3, y3], color=C_FEED, lw=0.9,
                       ls=(0, (1.2, 1.4)), zorder=5))
    ax.text(15, y3, 'fed-batch feed; dotted = sensing', fontsize=FS_LEGEND,
            va='center', zorder=5)

    _tbar(ax, (60, y1), (53, y1), C_INHIB, lw=1.0)
    ax.text(62.5, y1, 'product inhibition, $e^{-k_i C}$ (r1, r4, r7)',
            fontsize=FS_LEGEND, va='center', zorder=5)
    _tbar(ax, (60, y2), (53, y2), C_REPR, lw=1.0)
    ax.text(62.5, y2, 'glucose repression (r2, r5, r8; r9 at high glc)',
            fontsize=FS_LEGEND, va='center', zorder=5)
    _leg_arrow(53, y3, C_ACT, dashed=True)
    ax.text(62.5, y3, 'activation (acetaldehyde $\\rightarrow$ r1; '
                      'glc & EtOH $\\rightarrow$ r9)',
            fontsize=FS_LEGEND, va='center', zorder=5)

    _badge(ax, 125, y1, 'O$_2$', C_O2)
    ax.text(129.5, y1, 'aerobic only (is_aerobic gate)',
            fontsize=FS_LEGEND, va='center', zorder=5)
    _badge(ax, 125, y2, 'AcDH', '#FFFFFF', tc='#1F7A5C', ec='#1F7A5C',
           w=8.6)
    ax.text(131, y2, 'requires AcDH machinery ($a\\,X_\\mathrm{AcDH}$)',
            fontsize=FS_LEGEND, va='center', zorder=5)
    ax.text(121.5, y3, 'all rates $\\propto$ $a = X_a x$; conc. in '
                       'g L$^{-1}$; dilution ($D$ = 0) omitted',
            fontsize=FS_LEGEND - 0.4, va='center', ha='left', zorder=5,
            color=C_MUTED)

    # === save ==============================================================
    if save_dir is None:
        save_dir = os.path.dirname(os.path.abspath(__file__))
    for fmt in formats:
        out = os.path.join(save_dir, f'conceptual_diagram.{fmt}')
        fig.savefig(out, dpi=dpi, facecolor='white')
        print(f'Saved {out}')
    if show:
        plt.show()
    return fig, ax


if __name__ == '__main__':
    draw_conceptual_diagram()
