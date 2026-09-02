# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Side-by-side network figure of cumulative flux and inhibition.

Renders one or more :class:`~nskinetics.FluxSummary` objects onto a shared
:class:`FluxMapSpec` node layout: edge width encodes cumulative flux (shared
scale across panels), and each inhibited step carries a compact strip whose
row lengths are the fraction of potential flux lost per inhibitor plus a grey
joint-effect row. Sized for a Nature-family double-column figure (Arial,
5-7 pt, vector PDF with editable text + 600 dpi PNG).
"""

import os
from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

__all__ = ('FluxMapSpec', 'draw_flux_map')

MM = 1 / 25.4
C_FLUX = '#3A3A3A'
C_JOINT = '#666666'
FS_NODE = 6.0
FS_VAL = 4.6
FS_HEAD = 7.0
FS_LEG = 5.2
NODE_HW = 6.0    # node box half-width in mm (see _draw_panel)
NODE_HH = 2.5    # node box half-height in mm
NODE_GAP = 0.8   # clear space in mm between a node box and an arrow

#: rcParams applied *only* while a figure is drawn and saved (Nature-family
#: conventions: Arial, editable vector text). Applied through
#: :func:`matplotlib.pyplot.rc_context` so a caller's later plots keep their
#: own styling; ``savefig`` must run inside that context because
#: ``pdf.fonttype`` is read at save time.
_RC_PARAMS = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'text.color': '#1A1A1A',
}


@dataclass
class FluxMapSpec:
    """Declarative layout for :func:`draw_flux_map`.

    Parameters
    ----------
    nodes : dict
        ``{node_id: (x_mm, y_mm, label)}``.
    edges : dict
        ``{reaction_id: (from_node, to_node)}``.
    inhibitors : dict
        Ordered ``{inhibitor_name: hex_color}``; sets strip row order/colors.
    inhibition_map : dict
        ``{parameter_id: (reaction_id, inhibitor_name)}`` (as for
        :func:`~nskinetics.compute_flux_summary`).
    enhancement_reactions : set
        Reactions whose negative fractions are drawn as enhancement (hatched,
        ``+`` marker) rather than loss.
    products : list
        Species ids whose final titers appear in each panel header.
    strip_offsets : dict, optional
        ``{reaction_id: (dx_mm, dy_mm)}`` offset of the strip from the edge
        midpoint. Missing reactions use ``(4.0, 0.0)``.
    size_mm : tuple
        ``(width_mm, height_mm)`` of one panel.
    """
    nodes: dict
    edges: dict
    inhibitors: dict
    inhibition_map: dict
    enhancement_reactions: set = field(default_factory=set)
    products: list = field(default_factory=list)
    strip_offsets: dict = field(default_factory=dict)
    size_mm: tuple = (88., 120.)


def _width_scale(summaries, spec):
    """Build the flux-to-linewidth mapping shared by every panel.

    Parameters
    ----------
    summaries : list of FluxSummary
        All panels, so that one scale serves the whole figure.
    spec : FluxMapSpec
        Supplies the drawn reactions.

    Returns
    -------
    callable
        ``flux -> linewidth in points`` (``0.0`` for non-positive flux).
    """
    fluxes = [s.cumulative_flux.get(rid, 0.0)
              for s in summaries for rid in spec.edges]
    fmax = max(fluxes) if fluxes else 1.0
    fmax = fmax or 1.0

    def scale(flux):
        # sqrt so area, not width, tracks flux; 0.4-3.0 pt linewidth range
        if flux <= 0:
            return 0.0
        return 0.4 + 2.6 * np.sqrt(flux / fmax)
    return scale


def _trim_to_nodes(p0, p1):
    """Pull an edge's endpoints back to the boundaries of its node boxes.

    ``FancyArrowPatch``'s ``shrinkA``/``shrinkB`` are expressed in points and
    so cannot track the node boxes, which are sized in mm data units: a
    points-based shrink leaves the arrowhead buried under the opaque,
    higher-``zorder`` target box for every near-horizontal edge. Trimming in
    data space instead keeps the head visible at any edge orientation.

    Parameters
    ----------
    p0, p1 : tuple of float
        ``(x_mm, y_mm)`` centers of the source and target nodes.

    Returns
    -------
    tuple
        The trimmed ``(p0, p1)`` endpoints.
    """
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    dist = float(np.hypot(dx, dy))
    if dist == 0.:
        return p0, p1
    ux, uy = dx / dist, dy / dist
    # distance from a node center to its box boundary along (ux, uy)
    reach = min(NODE_HW / abs(ux) if ux else np.inf,
                NODE_HH / abs(uy) if uy else np.inf)
    t = min(reach + NODE_GAP, 0.45 * dist)
    return ((p0[0] + ux * t, p0[1] + uy * t),
            (p1[0] - ux * t, p1[1] - uy * t))


def _draw_edge(ax, p0, p1, w):
    p0, p1 = _trim_to_nodes(p0, p1)
    if w <= 0:
        ax.add_patch(FancyArrowPatch(
            p0, p1, arrowstyle='-|>', mutation_scale=5,
            lw=0.5, ls=(0, (2, 2)), color='#B4B2A9',
            shrinkA=0, shrinkB=0, zorder=3))
        return
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle='-|>', mutation_scale=6, lw=w,
        color=C_FLUX, shrinkA=0, shrinkB=0, zorder=3))


def _draw_bar(ax, x, yy, L, h, frac, color, enh):
    """Draw one strip row: outline, proportional fill, enhancement marker.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Panel to draw on.
    x, yy : float
        Lower-left corner of the row, in mm.
    L, h : float
        Full row length and row height, in mm.
    frac : float
        Fraction of potential flux lost; negative means enhancement. The
        drawn length is clamped to the full row length.
    color : str
        Fill color for this row.
    enh : bool
        Whether this reaction is drawn as enhancement, so that a negative
        ``frac`` is hatched and marked ``+`` instead of read as a loss.
    """
    ax.add_patch(Rectangle((x, yy), L, h, fill=False,
                           ec='#CCCCCC', lw=0.3, zorder=5))
    mag = abs(frac)
    if mag == 0:
        return
    gain = enh and frac < 0
    ax.add_patch(Rectangle((x, yy), L * min(mag, 1.0), h, fc=color, ec='none',
                           hatch='///' if gain else None, zorder=6))
    if gain:
        ax.text(x + L * min(mag, 1.0) + 0.6, yy + h / 2, '+',
                fontsize=FS_VAL, va='center', ha='left', color=color,
                zorder=7)


def _draw_strip(ax, x, y, rid, summary, spec):
    fl = summary.fraction_lost.get(rid, {})
    L = 9.0          # full strip length in mm
    h = 1.1          # row height
    gap = 0.3
    enh = rid in spec.enhancement_reactions
    yy = y
    for inh, color in spec.inhibitors.items():
        _draw_bar(ax, x, yy, L, h, fl.get(inh, 0.0), color, enh)
        yy -= (h + gap)
    # joint row (grey), drawn on the same conventions as the rows above
    _draw_bar(ax, x, yy, L, h, summary.fraction_lost_all.get(rid, 0.0),
              C_JOINT, enh)


def _draw_panel(ax, summary, spec, scale, annotate_values):
    w_mm, h_mm = spec.size_mm
    ax.set_xlim(0, w_mm)
    ax.set_ylim(0, h_mm)
    ax.set_aspect('equal')
    ax.axis('off')
    # header
    header = summary.label or ''
    titer = ', '.join(
        f'{sp} {summary.final_concentrations.get(sp, float("nan")):.0f} g/L'
        for sp in spec.products)
    ax.text(2, h_mm - 3, header, fontsize=FS_HEAD, fontweight='bold',
            ha='left', va='top')
    ax.text(2, h_mm - 8, f'harvest {summary.t_end:.0f} h; {titer}',
            fontsize=FS_LEG, ha='left', va='top', color='#555555')
    # edges
    for rid, (a, b) in spec.edges.items():
        xa, ya, _ = spec.nodes[a]
        xb, yb, _ = spec.nodes[b]
        w = scale(summary.cumulative_flux.get(rid, 0.0))
        _draw_edge(ax, (xa, ya), (xb, yb), w)
        if annotate_values and summary.cumulative_flux.get(rid, 0.0) > 0:
            ax.text((xa + xb) / 2 - 1.5, (ya + yb) / 2,
                    f'{summary.cumulative_flux[rid]:.0f}',
                    fontsize=FS_VAL, ha='right', va='center', color=C_FLUX)
        # strip for inhibited reactions
        if any(v[0] == rid for v in spec.inhibition_map.values()):
            dx, dy = spec.strip_offsets.get(rid, (4.0, 0.0))
            _draw_strip(ax, (xa + xb) / 2 + dx, (ya + yb) / 2 + dy,
                        rid, summary, spec)
    # nodes
    for x, y, label in spec.nodes.values():
        ax.add_patch(FancyBboxPatch((x - NODE_HW, y - NODE_HH),
                     2 * NODE_HW, 2 * NODE_HH,
                     boxstyle='round,pad=0,rounding_size=1.2',
                     fc='white', ec=C_FLUX, lw=0.6, zorder=8))
        ax.text(x, y, label, fontsize=FS_NODE, ha='center', va='center',
                zorder=9)


def draw_flux_map(summaries, spec, save_dir=None, formats=('png', 'pdf'),
                  dpi=600, annotate_values=True, show=False):
    """Draw one panel per summary, side by side, on the shared ``spec`` layout.

    Parameters
    ----------
    summaries : sequence of FluxSummary
        One per panel (typically scenario A then B). All must share the same
        reaction set.
    spec : FluxMapSpec
        Node layout, edge-to-reaction map, inhibitor colors, and coefficient
        mapping.
    save_dir : str, optional
        If given, writes ``flux_map.<fmt>`` for each format.
    formats : sequence of str
        e.g. ``('png', 'pdf')``.
    dpi : int
        Raster resolution.
    annotate_values : bool
        Print the g/L value beside each edge.
    show : bool
        Call ``plt.show()`` after saving.

    Returns
    -------
    (matplotlib.figure.Figure, list of matplotlib.axes.Axes)

    Raises
    ------
    ValueError
        If ``summaries`` is empty, or if the summaries do not share the same
        reaction set.

    Notes
    -----
    Figure styling is confined to a :func:`matplotlib.pyplot.rc_context`, so
    the caller's global ``rcParams`` are unchanged on return. The returned
    figure is not closed; close it with :func:`matplotlib.pyplot.close` when
    drawing many.
    """
    summaries = list(summaries)
    if not summaries:
        raise ValueError('draw_flux_map needs at least one summary.')
    ref = set(summaries[0].reaction_ids)
    for s in summaries[1:]:
        if set(s.reaction_ids) != ref:
            raise ValueError(
                'All summaries must share the same reaction set; '
                f'{s.label!r} differs from {summaries[0].label!r}.')
    w_mm, h_mm = spec.size_mm
    n = len(summaries)
    scale = _width_scale(summaries, spec)
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(1, n, figsize=(w_mm * MM * n, h_mm * MM))
        axes = [axes] if n == 1 else list(axes)
        for ax, s in zip(axes, summaries):
            _draw_panel(ax, s, spec, scale, annotate_values)
        fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99,
                            wspace=0.04)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            for fmt in formats:
                # inside the context: pdf.fonttype is read at save time
                fig.savefig(os.path.join(save_dir, f'flux_map.{fmt}'),
                            dpi=dpi, facecolor='white')
    if show:
        plt.show()
    return fig, axes
