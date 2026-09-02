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
joint-effect row. A legend band beneath the panels keys both encodings, laid
out from measured text widths so items never collide at any figure width.
Sized for a Nature-family double-column figure (Arial, 5-8 pt -- no text below
the 5 pt floor, and 8 pt only for the bold panel letters -- vector PDF with
editable text + 600 dpi PNG).
"""

import os
from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

__all__ = ('FluxMapSpec', 'draw_flux_map')

MM = 1 / 25.4
C_FLUX = '#3A3A3A'
C_JOINT = '#666666'
C_ZERO = '#B4B2A9'   # faint grey: zero-flux edges and their strip outlines
C_ROW = '#CCCCCC'    # strip-row outline for an active reaction
FS_NODE = 6.0
FS_VAL = 5.0         # 5 pt is the Nature-family floor; nothing goes below it
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
LEADER_LW = 0.3      # hairline linking a strip to the edge it belongs to

# --- legend band metrics (mm) ---
LEG_PAD = 1.0        # clear space at the top/bottom of the band
LEG_MARGIN = 2.0     # clear space at the left/right of the band
LEG_GAP = 0.8        # between an item's graphic and its label
LEG_ITEM_GAP = 3.5   # between one item and the next
LEG_ROW_H = 4.0      # maximum row pitch; shrinks when rows wrap
LEG_SAMPLE_W = 5.0   # length of a flux-sample line
LEG_SWATCH_W = 3.5   # width of a strip-color swatch
LEG_SWATCH_H = 2.0   # height of a strip-color swatch

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
        Reactions expected to show enhancement rather than loss; a drawn one
        adds the enhancement key to the legend. It does NOT gate the drawing:
        a negative fraction is hatched and marked ``+`` wherever it occurs.
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


def _panel_letter(i):
    """Panel letter for zero-based panel index ``i``.

    Counts ``a``..``z``, then ``aa``, ``ab``, ... so any number of panels is
    labelled (a 26-entry literal would raise on the 27th).

    Parameters
    ----------
    i : int
        Zero-based panel index.

    Returns
    -------
    str
    """
    out = ''
    i += 1
    while i:
        i, rem = divmod(i - 1, 26)
        out = chr(ord('a') + rem) + out
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
        ``'0'`` only for an exact zero; ``'<0.1'`` for a nonzero value that
        would otherwise round to ``'0.0'``; one decimal below 10, else an
        integer.
    """
    if v == 0:
        return '0'
    if abs(v) < 0.05:
        # a real but tiny flux must not read as nothing at all
        return '<0.1'
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


def _draw_bar(ax, x, yy, L, h, frac, color, ec=C_ROW):
    """Draw one strip row: outline, proportional fill, enhancement marker.

    A negative ``frac`` is an enhancement -- removing the "inhibitor" would
    *lower* the flux -- and is always drawn hatched and marked ``+``, whatever
    the spec's ``enhancement_reactions`` says: that set only decides whether
    the legend carries an enhancement key.

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
    ec : str, optional
        Outline color; ``C_ZERO`` marks a reaction that carried no flux at
        all, so an inactive branch does not read as an uninhibited one.
    """
    ax.add_patch(Rectangle((x, yy), L, h, fill=False,
                           ec=ec, lw=0.3, zorder=5))
    mag = abs(frac)
    if mag == 0:
        return
    gain = frac < 0
    ax.add_patch(Rectangle((x, yy), L * min(mag, 1.0), h, fc=color, ec='none',
                           hatch='///' if gain else None, zorder=6))
    if gain:
        ax.text(x + L * min(mag, 1.0) + 0.6, yy + h / 2, '+',
                fontsize=FS_VAL, va='center', ha='left', color=color,
                zorder=7)


def _strip_bbox(x, y, n_rows):
    """Bounding box ``(x0, y0, x1, y1)`` in mm of a strip anchored at ``x, y``.

    Parameters
    ----------
    x, y : float
        Lower-left corner of the strip's TOP row, in mm; rows hang downward.
    n_rows : int
        Number of rows drawn (inhibitors plus the joint row).

    Returns
    -------
    tuple of float
    """
    y_top = y + STRIP_ROW_H
    y_bot = y - (n_rows - 1) * (STRIP_ROW_H + STRIP_ROW_GAP)
    return x, y_bot, x + STRIP_LEN, y_top


def _draw_leader(ax, bbox, target):
    """Link a strip to the edge it describes with a hairline leader.

    The leader runs from whichever corner of the strip is nearest the edge
    midpoint to that midpoint, so an offset strip is unambiguously attributed
    to its reaction. Nothing is drawn when the midpoint already lies inside
    the strip's bounding box (the strip sits on its own edge).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Panel to draw on.
    bbox : tuple of float
        ``(x0, y0, x1, y1)`` of the strip, in mm.
    target : tuple of float
        ``(x_mm, y_mm)`` edge midpoint to point at.
    """
    x0, y0, x1, y1 = bbox
    tx, ty = target
    if x0 <= tx <= x1 and y0 <= ty <= y1:
        return
    corners = ((x0, y0), (x0, y1), (x1, y0), (x1, y1))
    cx, cy = min(corners, key=lambda c: (c[0] - tx) ** 2 + (c[1] - ty) ** 2)
    ax.add_line(Line2D([cx, tx], [cy, ty], lw=LEADER_LW, color=C_ZERO,
                       solid_capstyle='butt', zorder=2))


def _draw_strip(ax, x, y, rid, summary, spec):
    fl = summary.fraction_lost.get(rid, {})
    # An inactive branch gets faint outlines, so that an all-empty strip reads
    # as "no flux" rather than as "no inhibition".
    ec = C_ROW if summary.cumulative_flux.get(rid, 0.0) else C_ZERO
    yy = y
    for inh, color in spec.inhibitors.items():
        _draw_bar(ax, x, yy, STRIP_LEN, STRIP_ROW_H, fl.get(inh, 0.0),
                  color, ec=ec)
        yy -= (STRIP_ROW_H + STRIP_ROW_GAP)
    # joint row (grey), drawn on the same conventions as the rows above
    _draw_bar(ax, x, yy, STRIP_LEN, STRIP_ROW_H,
              summary.fraction_lost_all.get(rid, 0.0), C_JOINT, ec=ec)


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
        # strip for inhibited reactions, with a hairline back to the edge
        if any(v[0] == rid for v in spec.inhibition_map.values()):
            dx, dy = spec.strip_offsets.get(rid, (4.0, 0.0))
            mx, my = (xa + xb) / 2, (ya + yb) / 2
            sx, sy = mx + dx, my + dy
            _draw_strip(ax, sx, sy, rid, summary, spec)
            _draw_leader(ax, _strip_bbox(sx, sy, len(spec.inhibitors) + 1),
                         (mx, my))
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


def _renderer(fig):
    """A renderer to measure text with, or ``None`` to let matplotlib pick.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure the text lives on.

    Returns
    -------
    matplotlib.backend_bases.RendererBase or None
        ``None`` is a valid argument to
        :meth:`~matplotlib.text.Text.get_window_extent`, which then builds a
        renderer itself; that is the fallback for a canvas (pdf, svg) that
        exposes none.
    """
    get = getattr(fig.canvas, 'get_renderer', None)
    return get() if get is not None else None


def _text_width_mm(ax, artist, renderer):
    """Width of a drawn text artist in the legend's mm data units.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes the text belongs to (its ``transData`` sets the mm scale).
    artist : matplotlib.text.Text
        The text to measure.
    renderer : matplotlib.backend_bases.RendererBase or None
        Renderer to measure with; ``None`` lets matplotlib supply one.

    Returns
    -------
    float
    """
    bb = artist.get_window_extent(renderer)
    inv = ax.transData.inverted()
    (x0, _), (x1, _) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y0)])
    return abs(x1 - x0)


def _legend_entries(spec, fmax):
    """The legend's items and captions, in reading order.

    Parameters
    ----------
    spec : FluxMapSpec
        Supplies inhibitor names/colors and the enhancement set.
    fmax : float
        Largest drawn cumulative flux.

    Returns
    -------
    (list, list)
        Items as ``(kind, payload, label, color)`` -- ``kind`` is ``'line'``
        (flux sample, payload = the flux value), ``'dash'`` (inactive-pathway
        sample) or ``'swatch'`` (payload = hatch or ``None``) -- and the
        caption strings.
    """
    items = [('line', v, _fmt_value(v), C_FLUX) for v in
             _legend_flux_values(fmax)]
    items.append(('dash', None, 'pathway inactive', C_ZERO))
    for inh, color in spec.inhibitors.items():
        items.append(('swatch', None, inh, color))
    items.append(('swatch', None, 'all inhibitors', C_JOINT))
    if any(rid in spec.enhancement_reactions for rid in spec.edges):
        items.append(('swatch', '///', 'enhancement', C_JOINT))
    captions = ["cumulative flux, g of each step's substrate per L "
                'final broth',
                'fraction of potential flux lost (strip = 0-100 %)']
    return items, captions


def _draw_legend(ax, spec, scale, fmax, width_mm):
    """Draw the legend band: edge-width key, then the strip key.

    Items flow left to right and wrap to a further row whenever the next one
    would cross the band's right margin, so the band stays clean at one-panel
    (88 mm) and two-panel (176 mm) widths alike. Each item's width is the
    measured width of its own label (via
    :meth:`~matplotlib.text.Text.get_window_extent`, converted to the band's
    mm data units) plus its graphic, so nothing depends on a guessed pitch.
    The row pitch shrinks if the content needs more rows than the band's
    nominal height allows.

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
    fig = ax.get_figure()
    fig.canvas.draw()                 # so a renderer exists to measure with
    renderer = _renderer(fig)
    items, captions = _legend_entries(spec, fmax)
    usable = width_mm - 2 * LEG_MARGIN

    # --- pass 1: create every text, measure it, and plan the rows ---
    plan = []          # (row, x_mm, entry, artists...)
    row, x = 0, 0.0
    for kind, payload, label, color in items:
        marker = None
        if kind == 'swatch' and payload:
            marker = ax.text(0, 0, '+', fontsize=FS_VAL, ha='left',
                             va='center', color=color)
        text = ax.text(0, 0, label, fontsize=FS_LEG, ha='left', va='center',
                       color='#555555' if kind == 'dash' else None)
        graphic_w = LEG_SAMPLE_W if kind in ('line', 'dash') else LEG_SWATCH_W
        w = graphic_w + LEG_GAP + _text_width_mm(ax, text, renderer)
        if marker is not None:
            w += _text_width_mm(ax, marker, renderer) + LEG_GAP
        if x > 0. and x + w > usable:
            row, x = row + 1, 0.
        plan.append((row, x, kind, payload, color, graphic_w, marker, text))
        x += w + LEG_ITEM_GAP
    row += 1
    cap_artists = []
    x = 0.
    for cap in captions:
        text = ax.text(0, 0, cap, fontsize=FS_LEG, ha='left', va='center',
                       color='#555555')
        w = _text_width_mm(ax, text, renderer)
        if x > 0. and x + w > usable:
            row, x = row + 1, 0.
        cap_artists.append((row, x, text))
        x += w + 2 * LEG_ITEM_GAP
    n_rows = row + 1

    # --- pass 2: place everything on a pitch that fits the band ---
    pitch = min(LEG_ROW_H, (LEGEND_H - 2 * LEG_PAD) / n_rows)
    top = LEGEND_H / 2 + n_rows * pitch / 2

    def _y(r):
        return top - (r + 0.5) * pitch

    for r, x0, kind, payload, color, graphic_w, marker, text in plan:
        x0 += LEG_MARGIN
        y = _y(r)
        if kind == 'line':
            ax.plot([x0, x0 + graphic_w], [y, y], color=color,
                    lw=scale(payload), solid_capstyle='butt')
        elif kind == 'dash':
            ax.plot([x0, x0 + graphic_w], [y, y], color=color, lw=0.5,
                    ls=(0, (2, 2)))
        else:
            ax.add_patch(Rectangle((x0, y - LEG_SWATCH_H / 2), graphic_w,
                                   LEG_SWATCH_H, fc=color, ec=C_ROW, lw=0.3,
                                   hatch=payload))
        tx = x0 + graphic_w + LEG_GAP
        if marker is not None:
            marker.set_position((tx, y))
            tx += _text_width_mm(ax, marker, renderer) + LEG_GAP
        text.set_position((tx, y))
    for r, x0, text in cap_artists:
        text.set_position((x0 + LEG_MARGIN, _y(r)))


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
            _draw_panel(ax, s, spec, scale, annotate_values, _panel_letter(i))
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
