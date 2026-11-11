# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""What the build GIFs of Figure 6 and Extended Data Figure 3 share: drawing
the figure with the unmodified snapshot plot_kin_opt_parameter_sets.py while
recording which artist shows what, and the panel-B (final proteome
allocation) storyboard, which is the same in both figures.

Recording: the plot module is loaded by file path and a few of its drawing
helpers are wrapped at runtime (the file itself is untouched):

    set_colors            the plotted sets and their colours
    outcome_cell          one panel-A cell; while it runs, the cell's scatter,
                          step, axhline and axvspan calls are noted, each
                          artist tagged with what it is (trial cloud,
                          incumbent, baseline reference, preload zone, ...)
                          and, through its colour, the campaign it belongs to
    _draw_relay_preload   marks the scatters drawn inside it as the
                          TRY-informed campaign's preloaded (scout) trials
    _draw_best_irr_marks  marks the scatters drawn inside it as highest-IRR
                          trials
    draw_burden           panel B: its bars (Axes.barh calls, by row) apart
                          from everything else it draws
    _frame_keys           the frame drawn around the panel-A legends

Panel-B storyboard (both figures, after panel A):

    proteome-baseline         the panel title, axes, sector brackets and
                              sector key fade in, then the baseline bar grows
    proteome-uninformed       the uninformed profitability campaign's bar
    proteome-isobutanol       the three isobutanol campaigns' bars
    proteome-ethanol          the three ethanol campaigns' bars
    proteome-try-informed     the TRY-informed profitability campaign's bar

This file is kept byte-identical in Fig6/code and EDFig3/code.

Sim-safe: imports only the standard library, numpy, matplotlib and the
snapshot plot module (which reads CSVs and workbooks only).
"""
import argparse
import dataclasses
import importlib.util
import os
import sys

import matplotlib.axes
from matplotlib.artist import Artist
from matplotlib.colors import to_hex
from matplotlib.legend import Legend
from matplotlib.transforms import Affine2D, Bbox

import gif_engine as ge

HERE = os.path.dirname(os.path.abspath(__file__))
PLOT_SCRIPT = os.path.join(HERE, 'nskinetics', 'biorefineries', 'isobutanol',
                           'plots', 'plot_kin_opt_parameter_sets.py')

#: seconds: a fade, a trajectory being drawn along the trial axis, a bar row
#: growing, the rest before an annotation
FADE_S, TRACE_S, BAR_S, BEAT_S = 0.8, 3.0, 1.6, 0.5
#: white margin kept below panel B's x-axis title (inches); the static figure
#: leaves ~0.9 in there
BOTTOM_MARGIN = 0.25


def load_plot_module():
    spec = importlib.util.spec_from_file_location('_kin_opt_parameter_sets',
                                                  PLOT_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# -----------------------------------------------------------------------------
# campaigns
# -----------------------------------------------------------------------------
def campaign_kind(plot, s):
    """'baseline', 'uninformed', 'try-informed', 'isobutanol' or 'ethanol':
    the role of a plotted set in the story."""
    if s.get('is_baseline'):
        return 'baseline'
    if s.get('is_relay'):
        return 'try-informed'
    if plot._owns_outcome(s.get('objective'), 'IRR'):
        return 'uninformed'
    objective = s.get('objective') or ''
    if objective.startswith('IBO'):
        return 'isobutanol'
    if objective.startswith('EtOH'):
        return 'ethanol'
    raise ValueError(f'set {s["label"]!r}: unexpected objective {objective!r}')


# -----------------------------------------------------------------------------
# recording what the plot script draws
# -----------------------------------------------------------------------------
@dataclasses.dataclass
class Mark:
    kind: str          # cloud | incumbent | incumbent-other | baseline | zone
    label: object      # | preload | preload-mark | best;  campaign label/None
    artist: object


@dataclasses.dataclass
class Cell:
    """One panel-A cell: its axes, outcome column and tagged artists."""
    ax: object
    col: str
    marks: list

    def find(self, kind, labels=None):
        """The artists of one kind, optionally only those of the campaigns
        whose label is in ``labels``."""
        return [m.artist for m in self.marks if m.kind == kind
                and (labels is None or m.label in labels)]

    def find_one(self, kind, label=None):
        found = self.find(kind, None if label is None else {label})
        if len(found) != 1:
            raise LookupError(f'cell {self.col!r}: expected one {kind!r} '
                              f'artist of {label!r}, found {len(found)}')
        return found[0]


@dataclasses.dataclass
class Burden:
    """Panel B: its two axes (the broken x axis), the rows' bars by campaign
    label (top row first) and the figure-level artists it added."""
    axes: list
    rows: dict
    row_y: dict
    figure_artists: list


class Recorder:
    """Wraps the plot module's drawing helpers to tag the artists they draw."""

    def __init__(self, plot):
        self.plot = plot
        self.sets = []
        self.colors = {}
        self.cells = []
        self.burden = None
        self.key_frame = None
        self._context = None
        self._preload_seen = False
        for name in ('set_colors', 'outcome_cell', '_draw_relay_preload',
                     '_draw_best_irr_marks', 'draw_burden', '_frame_keys'):
            setattr(plot, name, getattr(self, '_' + name.lstrip('_'))(
                getattr(plot, name)))

    # -- lookups --------------------------------------------------------------
    def label_of(self, color):
        """The label of the campaign drawn in this colour."""
        hexed = to_hex(color).lower()
        for s in self.sets:
            if to_hex(self.colors[id(s)]).lower() == hexed:
                return s['label']
        raise LookupError(f'no campaign is drawn in {color!r}')

    def labels(self, *kinds):
        """Labels of the sets of these campaign kinds, in set order."""
        return [s['label'] for s in self.sets
                if campaign_kind(self.plot, s) in kinds]

    def set_of(self, label):
        (found,) = [s for s in self.sets if s['label'] == label]
        return found

    # -- wrappers -------------------------------------------------------------
    def _set_colors(self, original):
        def set_colors(sets):
            colors = original(sets)
            self.sets, self.colors = list(sets), colors
            return colors
        return set_colors

    def _outcome_cell(self, original):
        def outcome_cell(ax, sets, colors, col, *args, **kwargs):
            cell = Cell(ax, col, [])
            noted = {}

            def note(method, tagger):
                drawn = getattr(ax, method)

                def call(*a, **k):
                    artist = drawn(*a, **k)
                    cell.marks.extend(tagger(artist, k))
                    return artist
                noted[method] = call
                setattr(ax, method, call)

            note('scatter', self._tag_scatter)
            note('step', lambda lines, k: [
                Mark('incumbent' if k.get('ls') == '-' else 'incumbent-other',
                     self.label_of(k['color']), line) for line in lines])
            note('axhline', lambda line, k: [
                Mark('baseline', self.labels('baseline')[0], line)])
            note('axvspan', lambda patch, k: [Mark('zone', None, patch)])
            try:
                original(ax, sets, colors, col, *args, **kwargs)
            finally:
                for method in noted:
                    delattr(ax, method)      # back to the Axes method
            self.cells.append(cell)
        return outcome_cell

    def _tag_scatter(self, artist, k):
        if self._context == 'best':
            return [Mark('best', self.label_of(k['edgecolors']), artist)]
        if self._context == 'preload':
            # the helper draws the whole preloaded cloud first, then one open
            # circle per marked preloaded trial
            if not self._preload_seen:
                self._preload_seen = True
                return [Mark('preload', None, artist)]
            return [Mark('preload-mark', self.label_of(k['edgecolors']),
                         artist)]
        return [Mark('cloud', self.label_of(k['color']), artist)]

    def _in_context(self, name, original):
        def helper(*args, **kwargs):
            self._context, self._preload_seen = name, False
            try:
                return original(*args, **kwargs)
            finally:
                self._context = None
        return helper

    def _draw_relay_preload(self, original):
        return self._in_context('preload', original)

    def _draw_best_irr_marks(self, original):
        return self._in_context('best', original)

    def _draw_burden(self, original):
        def draw_burden(fig, gs_cell, sets, colors, *args, **kwargs):
            before = _figure_artists(fig)
            axes_before = list(fig.axes)
            bars = []
            barh = matplotlib.axes.Axes.barh

            def noted_barh(ax, y, *a, **k):
                container = barh(ax, y, *a, **k)
                bars.append((y, container))
                return container
            matplotlib.axes.Axes.barh = noted_barh
            try:
                result = original(fig, gs_cell, sets, colors, *args, **kwargs)
            finally:
                matplotlib.axes.Axes.barh = barh
            n = len(sets)
            rows, row_y = {}, {}
            for s_index, s in enumerate(sets):    # row y = n - set index
                row_y[s['label']] = n - s_index
                rows[s['label']] = [patch for y, container in bars
                                    if y == n - s_index for patch in container]
            self.burden = Burden(
                [ax for ax in fig.axes if ax not in axes_before],
                rows, row_y,
                [a for a in _figure_artists(fig) if a not in before])
            return result
        return draw_burden

    def _frame_keys(self, original):
        def frame_keys(fig, *args, **kwargs):
            before = _figure_artists(fig)
            result = original(fig, *args, **kwargs)
            (self.key_frame,) = [a for a in _figure_artists(fig)
                                 if a not in before]
            return result
        return frame_keys


def _figure_artists(fig):
    """The figure-level artists (texts, legends, patches, ...), unique."""
    seen, out = set(), []
    for artist in (fig.texts + fig.legends + fig.patches + fig.lines
                   + fig.artists + fig.images):
        if id(artist) not in seen:
            seen.add(id(artist))
            out.append(artist)
    return out


def draw_figure(variant):
    """Draw one figure variant ('irr' or 'process') with the snapshot plot
    script, as main() would with --variants <variant>, but keep it open
    instead of saving it. Returns (figure, recorder)."""
    font = ge.use_arial_or_metric_twin()
    if font != 'Arial':
        print(f'Arial is not installed; using {font or "the fallback font"}')
    plot = load_plot_module()
    recorder = Recorder(plot)
    # nothing is written (capture_figure switches savefig off); --out-dir only
    # has to pass the path-length check
    fig = ge.capture_figure(plot.main, ['--variants', variant,
                                        '--out-dir', HERE])
    if len(recorder.cells) == 0 or recorder.burden is None:
        raise RuntimeError('the plot script drew no panel-A cells or no '
                           'panel B; has it changed?')
    return fig, recorder


# -----------------------------------------------------------------------------
# figure-level helpers
# -----------------------------------------------------------------------------
def figure_text(fig, content):
    """The figure-level text with exactly this content."""
    (found,) = [t for t in fig.texts if t.get_text() == content]
    return found


def legend_with_entry(fig, label):
    """The figure legend that has an entry with this label."""
    (found,) = [leg for leg in fig.legends
                if any(t.get_text() == label for t in leg.get_texts())]
    return found


def legend_entries(fig, *labels):
    """The handle and text of each legend entry with one of these labels
    (each matched by its full text, or as a prefix for the long relay
    label). Hidden legend children keep their layout space."""
    out = []
    for label in labels:
        found = []
        for leg in fig.legends:
            for handle, text in zip(leg.legend_handles, leg.get_texts()):
                if text.get_text() == label or (
                        text.get_text().startswith(label + ':')):
                    found += [handle, text]
        if len(found) != 2:
            raise LookupError(f'no single legend entry {label!r}')
        out += found
    return out


def campaign_key_frame(fig, recorder):
    """The campaign legend's title and the frame around the panel-A keys:
    shown with the first legend entry."""
    legend = legend_with_entry(fig, recorder.labels('baseline')[0])
    return [legend.get_title(), recorder.key_frame]


def grow_canvas_top(fig, inches):
    """Grow the canvas by ``inches`` at the top (0 is fine), keeping every
    artist where it was, measured from the bottom-left corner in inches.

    The plot script places its axes, figure texts, legends and the key frame
    in figure fractions, which would stretch with the canvas. They are
    re-anchored on transforms built from transFigure, which savefig also
    shifts when it crops (bbox_inches); the legends' own default anchoring
    would not follow a crop."""
    width, height = fig.get_size_inches()
    to_fraction = fig.transFigure.inverted()
    fig.canvas.draw()
    axes = [(ax, ax.get_position(original=True)) for ax in fig.axes
            if ax.get_axes_locator() is None]
    anchors = [(leg, leg.get_bbox_to_anchor().transformed(to_fraction))
               for leg in fig.legends]
    others = [a for a in _figure_artists(fig) if not isinstance(a, Legend)
              and Artist.get_transform(a) is fig.transFigure]
    fig.set_size_inches(width, height + inches)
    scale = height / (height + inches)
    old_fraction = Affine2D().scale(1.0, scale) + fig.transFigure
    for ax, box in axes:
        ax.set_position(Bbox.from_extents(box.x0, box.y0 * scale,
                                          box.x1, box.y1 * scale))
    for leg, box in anchors:
        leg.set_bbox_to_anchor(box, transform=old_fraction)
    for artist in others:
        artist.set_transform(old_fraction)


def inch_transform(fig):
    """Figure inches -> display, for the canvas size at the time of the
    call. Like fig.dpi_scale_trans, but cropped with the rest of the figure
    when savefig is given bbox_inches."""
    width, height = fig.get_size_inches()
    return Affine2D().scale(1.0 / width, 1.0 / height) + fig.transFigure


def bottom_crop(fig, recorder):
    """savefig keyword arguments that drop the empty band below panel B,
    keeping BOTTOM_MARGIN inches under its x-axis title. Call it last,
    once the canvas has its final size."""
    grow_canvas_top(fig, 0.0)    # re-anchor the legends, so they follow the crop
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lowest = min(a.get_window_extent(renderer).y0
                 for a in recorder.burden.figure_artists if a.get_visible())
    bottom = max(0.0, lowest / fig.dpi - BOTTOM_MARGIN)
    width, height = fig.get_size_inches()
    return {'bbox_inches': Bbox.from_extents(0, bottom, width, height)}


# -----------------------------------------------------------------------------
# panel B (identical in both figures)
# -----------------------------------------------------------------------------
def group_panel_b(board, fig, recorder):
    """Register panel B's groups: the scaffold, and one group per bar-row
    block (see the module docstring)."""
    burden = recorder.burden
    board.group('B scaffold', figure_text(fig, 'B'),
                figure_text(fig, 'Final proteome allocation'),
                burden.axes, burden.figure_artists)
    for kind in ('baseline', 'uninformed', 'isobutanol', 'ethanol',
                 'try-informed'):
        board.group(f'B {kind} rows',
                    [burden.rows[label] for label in recorder.labels(kind)])


def _rows_region(recorder, kind):
    """The band of panel B holding the bar rows of one campaign kind, across
    both halves of the broken axis."""
    burden = recorder.burden
    ys = [burden.row_y[label] for label in recorder.labels(kind)]
    return ge.Region(burden.axes, y=(min(ys) - 0.45, max(ys) + 0.45))


def panel_b_steps(board, recorder):
    """The panel-B steps (see the module docstring)."""
    grow = lambda kind, seconds=BAR_S: ge.wipe(  # noqa: E731
        f'B {kind} rows', region=_rows_region(recorder, kind),
        seconds=seconds)
    board.step('proteome-baseline',
               'Panel B shows how each final design allocates the proteome '
               'to sectors, starting with the baseline strain',
               ge.fade('B scaffold', seconds=1.0), ge.pause(BEAT_S),
               grow('baseline', 2.0))
    board.step('proteome-uninformed',
               'The uninformed profitability campaign spends most of its '
               'metabolic budget on glycolysis and ethanol production',
               grow('uninformed'))
    board.step('proteome-isobutanol',
               'The isobutanol campaigns devote a large share of the '
               'proteome to isobutanol production',
               grow('isobutanol', 2.0))
    board.step('proteome-ethanol',
               'The ethanol campaigns put their proteome into glycolysis and '
               'ethanol production',
               grow('ethanol', 2.0))
    board.step('proteome-try-informed',
               'The TRY-informed profitability campaign splits its proteome '
               'between ethanol and isobutanol production',
               grow('try-informed'))


# -----------------------------------------------------------------------------
# command line
# -----------------------------------------------------------------------------
def run(description, variant, stem, build, argv=None):
    """Parse the shared GIF options, draw the figure variant, build the
    storyboard with ``build(fig, recorder)`` and render it into
    <figure folder>/output/gifs (or --out-dir)."""
    parser = argparse.ArgumentParser(description=description)
    ge.add_cli_arguments(parser)
    args = parser.parse_args(argv)
    out_dir = args.out_dir or os.path.join(os.path.dirname(HERE), 'output',
                                           'gifs')
    fig, recorder = draw_figure(variant)
    board = build(fig, recorder)
    board.render(out_dir, stem, **ge.render_kwargs(args))
    return board
