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
joint-effect row. A legend band beneath the panels keys both encodings. Sized
for a Nature-family double-column figure (Arial, 5-8 pt, vector PDF with
editable text + 600 dpi PNG).
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
C_ZERO = '#B4B2A9'   # faint grey: zero-flux edges and their strip outlines
C_ROW = '#CCCCCC'    # strip-row outline for an active reaction
FS_NODE = 6.0
FS_VAL = 4.6
FS_HEAD = 7.0
FS_LETTER = 8.0      # Nature-style bold lowercase panel letter
FS_LEG = 5.2
NODE_HW = 6.0    # node box half-width in mm (see _draw_panel)
NODE_HH = 2.5    # node box half-height in mm
NODE_GAP = 0.8   # clear space in mm between a node box and an arrow
STRIP_LEN = 9.0      # full strip-row length in mm (a row = 0-100 % lost)
STRIP_ROW_H = 1.1    # strip row height in mm
STRIP_ROW_GAP = 0.3  # vertical gap between strip rows in mm
LABEL_OFFSET = 1.6   # edge value label offset from the line, in mm
LEGEND_H = 13.0      # height of the legend band in mm
MAX_MUTATION_SCALE = 9.0   # cap on the flux-scaled arrowhead size

#: Panel letters, in order.
_LETTERS = 'abcdefghijklmnopqrstuvwxyz'

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
    product_labels : dict, optional
        ``{species_id: display name}`` used for those header titers; a species
        with no entry is shown by its raw id.
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
    product_labels: dict = field(default_factory=dict)
    strip_offsets: dict = field(default_factory=dict)
    size_mm: tuple = (88., 120.)

    @property
    def reactions(self):
        """Reactions to summarize for this spec, in drawing order.

        The drawn edges first (in ``edges`` order), then any reaction that is
        named in ``inhibition_map`` but has no edge (in first-seen order),
        so that a mapped-but-undrawn reaction such as biomass decay still gets
        its numbers computed by
        :func:`~nskinetics.compute_flux_summary`.

        Returns
        -------
        list of str
        """
        out = list(self.edges)
        for _p, (rid, _inh) in self.inhibition_map.items():
            if rid not in out:
                out.append(rid)
        return out


def _fmt_value(v):
    """Format a cumulative flux for an edge label or header titer.

    Parameters
    ----------
    v : float
        Value in g/L of final broth.

    Returns
    -------
    str
        ``'0'`` only for an exact zero, one decimal below 10, else an integer.
    """
    if v == 0:
        return '0'
    if abs(v) < 10:
        return f'{v:.1f}'
    return f'{v:.0f}'


def _flux_max(summaries, spec):
    """Largest drawn cumulative flux across every panel (never zero).

    Parameters
    ----------
    summaries : list of FluxSummary
        All panels, so that one scale serves the whole figure.
    spec : FluxMapSpec
        Supplies the drawn reactions.

    Returns
    -------
    float
    """
    fluxes = [s.cumulative_flux.get(rid, 0.0)
              for s in summaries for rid in spec.edges]
    fmax = max(fluxes) if fluxes else 1.0
    return fmax or 1.0


def _width_scale(fmax):
    """Build the flux-to-linewidth mapping shared by every panel.

    Parameters
    ----------
    fmax : float
        Largest drawn cumulative flux (see :func:`_flux_max`).

    Returns
    -------
    callable
        ``flux -> linewidth in points`` (``0.0`` for non-positive flux).
    """
    def scale(flux):
        # sqrt so area, not width, tracks flux; 0.4-3.0 pt linewidth range
        if flux <= 0:
            return 0.0
        return 0.4 + 2.6 * np.sqrt(flux / fmax)
    return scale


def _mutation_scale(w):
    """Arrowhead size for an edge of linewidth ``w`` points."""
    return min(4.0 + 1.2 * w, MAX_MUTATION_SCALE)


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


def _label_placement(p0, p1, offset=LABEL_OFFSET):
    """Place an edge's value label perpendicular to the edge, off the line.

    The label sits ``offset`` mm from the edge midpoint along the unit normal
    ``n = (-uy, ux)``, flipped so that ``n`` points upward (``n_y > 0``) or,
    for an exactly vertical edge, to the right (``n_x > 0``). ``ha``/``va``
    follow the sign of ``n``, so the text is pushed away from the line instead
    of being struck through by it.

    Parameters
    ----------
    p0, p1 : tuple of float
        ``(x_mm, y_mm)`` centers of the source and target nodes.
    offset : float
        Perpendicular clearance in mm.

    Returns
    -------
    tuple
        ``(x_mm, y_mm, ha, va)``.
    """
    mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    dist = float(np.hypot(dx, dy))
    if dist == 0.:
        return mx, my, 'center', 'center'
    ux, uy = dx / dist, dy / dist
    nx, ny = -uy, ux
    tol = 1e-9
    if ny < -tol or (abs(ny) <= tol and nx < 0.):
        nx, ny = -nx, -ny
    ha = 'left' if nx > tol else ('right' if nx < -tol else 'center')
    va = 'bottom' if ny > tol else ('top' if ny < -tol else 'center')
    return mx + nx * offset, my + ny * offset, ha, va


def _draw_edge(ax, p0, p1, w):
    p0, p1 = _trim_to_nodes(p0, p1)
    if w <= 0:
        ax.add_patch(FancyArrowPatch(
            p0, p1, arrowstyle='-|>', mutation_scale=5,
            lw=0.5, ls=(0, (2, 2)), color=C_ZERO,
            shrinkA=0, shrinkB=0, zorder=3))
        return
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle='-|>', mutation_scale=_mutation_scale(w), lw=w,
        color=C_FLUX, shrinkA=0, shrinkB=0, zorder=3))


def _draw_bar(ax, x, yy, L, h, frac, color, enh, ec=C_ROW):
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
    ec : str, optional
        Outline color; ``C_ZERO`` marks a reaction that carried no flux at
        all, so an inactive branch does not read as an uninhibited one.
    """
    ax.add_patch(Rectangle((x, yy), L, h, fill=False,
                           ec=ec, lw=0.3, zorder=5))
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
    enh = rid in spec.enhancement_reactions
    # An inactive branch gets faint outlines, so that an all-empty strip reads
    # as "no flux" rather than as "no inhibition".
    ec = C_ROW if summary.cumulative_flux.get(rid, 0.0) else C_ZERO
    yy = y
    for inh, color in spec.inhibitors.items():
        _draw_bar(ax, x, yy, STRIP_LEN, STRIP_ROW_H, fl.get(inh, 0.0),
                  color, enh, ec=ec)
        yy -= (STRIP_ROW_H + STRIP_ROW_GAP)
    # joint row (grey), drawn on the same conventions as the rows above
    _draw_bar(ax, x, yy, STRIP_LEN, STRIP_ROW_H,
              summary.fraction_lost_all.get(rid, 0.0), C_JOINT, enh, ec=ec)


def _draw_panel(ax, summary, spec, scale, annotate_values, letter):
    w_mm, h_mm = spec.size_mm
    ax.set_xlim(0, w_mm)
    ax.set_ylim(0, h_mm)
    ax.set_aspect('equal')
    ax.axis('off')
    # header: bold lowercase panel letter, then the scenario label on the same
    # baseline, then a sub-header of harvest time and final titers
    ax.text(2, h_mm - 4, letter, fontsize=FS_LETTER, fontweight='bold',
            ha='left', va='baseline')
    ax.text(5.6, h_mm - 4, summary.label or '', fontsize=FS_HEAD,
            fontweight='bold', ha='left', va='baseline')
    parts = [f'harvest {summary.t_end:.0f} h']
    for sp in spec.products:
        name = spec.product_labels.get(sp, sp)
        v = summary.final_concentrations.get(sp, float('nan'))
        parts.append(f'{name} {_fmt_value(v)} g/L')
    ax.text(2, h_mm - 7.8, ' · '.join(parts),
            fontsize=FS_LEG, ha='left', va='baseline', color='#555555')
    # edges
    for rid, (a, b) in spec.edges.items():
        xa, ya, _ = spec.nodes[a]
        xb, yb, _ = spec.nodes[b]
        flux = summary.cumulative_flux.get(rid, 0.0)
        _draw_edge(ax, (xa, ya), (xb, yb), scale(flux))
        if annotate_values and flux > 0:
            lx, ly, ha, va = _label_placement((xa, ya), (xb, yb))
            ax.text(lx, ly, _fmt_value(flux), fontsize=FS_VAL,
                    ha=ha, va=va, color=C_FLUX)
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


def _nice_below(v):
    """Largest "nice" (1-2 significant figure) value not exceeding ``v``."""
    if v <= 0:
        return 0.0
    k = 10.0 ** np.floor(np.log10(v))
    for m in (10., 8., 6., 5., 4., 3., 2.5, 2., 1.5, 1.):
        if m * k <= v * (1 + 1e-9):
            return float(m * k)
    return float(k)


def _legend_flux_values(fmax):
    """Three ascending "nice" flux values for the edge-width key.

    Parameters
    ----------
    fmax : float
        Largest drawn cumulative flux; the top key value is the largest nice
        value at or below it, so no sample is wider than any drawn edge.

    Returns
    -------
    list of float
    """
    vals = []
    for frac in (0.1, 0.4, 1.0):
        v = _nice_below(fmax * frac)
        if v > 0 and v not in vals:
            vals.append(v)
    return vals


def _draw_legend(ax, spec, scale, fmax, width_mm):
    """Draw the legend band: edge-width key (left), strip key (right).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The legend axes, spanning the full figure width.
    spec : FluxMapSpec
        Supplies inhibitor names/colors and the enhancement set.
    scale : callable
        The same flux-to-linewidth mapping the panels use, so a sample line
        is directly comparable to a drawn edge.
    fmax : float
        Largest drawn cumulative flux.
    width_mm : float
        Full figure width in mm (the legend's x data range).
    """
    ax.set_xlim(0, width_mm)
    ax.set_ylim(0, LEGEND_H)
    ax.set_aspect('equal')
    ax.axis('off')
    y_item, y_cap, x0 = 8.0, 2.6, 2.0
    sec2_x = 0.50 * width_mm

    # --- edge-width key ---
    vals = _legend_flux_values(fmax)
    pitch = (sec2_x - 2.0 - x0) / (len(vals) + 1)
    samp = min(6.0, 0.42 * pitch)
    for i, v in enumerate(vals):
        xs = x0 + i * pitch
        ax.plot([xs, xs + samp], [y_item, y_item], color=C_FLUX,
                lw=scale(v), solid_capstyle='butt')
        ax.text(xs + samp + 0.8, y_item, _fmt_value(v), fontsize=FS_LEG,
                ha='left', va='center')
    xs = x0 + len(vals) * pitch
    ax.plot([xs, xs + samp], [y_item, y_item], color=C_ZERO, lw=0.5,
            ls=(0, (2, 2)))
    ax.text(xs + samp + 0.8, y_item, 'pathway inactive', fontsize=FS_LEG,
            ha='left', va='center', color='#555555')
    ax.text(x0, y_cap, 'cumulative flux, g L$^{-1}$ final broth',
            fontsize=FS_LEG, ha='left', va='center', color='#555555')

    # --- strip key ---
    labels = list(spec.inhibitors) + ['all inhibitors']
    colors = list(spec.inhibitors.values()) + [C_JOINT]
    hatches = [None] * len(labels)
    if any(rid in spec.enhancement_reactions for rid in spec.edges):
        labels.append('enhancement')
        colors.append(C_JOINT)
        hatches.append('///')
    pitch = (width_mm - 2.0 - sec2_x) / len(labels)
    sw_w, sw_h = min(4.0, 0.3 * pitch), 2.0
    for i, (lab, col, hat) in enumerate(zip(labels, colors, hatches)):
        x = sec2_x + i * pitch
        ax.add_patch(Rectangle((x, y_item - sw_h / 2), sw_w, sw_h, fc=col,
                               ec=C_ROW, lw=0.3, hatch=hat))
        tx = x + sw_w + 0.8
        if hat:
            ax.text(tx, y_item, '+', fontsize=FS_VAL, ha='left', va='center',
                    color=col)
            tx += 1.4
        ax.text(tx, y_item, lab, fontsize=FS_LEG, ha='left', va='center')
    ax.text(sec2_x, y_cap,
            'fraction of potential flux lost (strip = 0-100 %)',
            fontsize=FS_LEG, ha='left', va='center', color='#555555')


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
        The figure and the PANEL axes only, one per summary. The legend band
        is an additional axes on the figure, reachable as ``fig.axes[-1]``.

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
    fmax = _flux_max(summaries, spec)
    scale = _width_scale(fmax)
    total_w = w_mm * n
    with plt.rc_context(_RC_PARAMS):
        fig = plt.figure(figsize=(total_w * MM, (h_mm + LEGEND_H) * MM))
        gs = fig.add_gridspec(2, n, height_ratios=[h_mm, LEGEND_H],
                              left=0.01, right=0.99, bottom=0.01, top=0.99,
                              wspace=0.04, hspace=0.0)
        axes = [fig.add_subplot(gs[0, i]) for i in range(n)]
        legend_ax = fig.add_subplot(gs[1, :])
        for i, (ax, s) in enumerate(zip(axes, summaries)):
            _draw_panel(ax, s, spec, scale, annotate_values, _LETTERS[i])
        _draw_legend(legend_ax, spec, scale, fmax, total_w)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            for fmt in formats:
                # inside the context: pdf.fonttype is read at save time
                fig.savefig(os.path.join(save_dir, f'flux_map.{fmt}'),
                            dpi=dpi, facecolor='white')
    if show:
        plt.show()
    return fig, axes
