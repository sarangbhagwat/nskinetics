# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Progressive presentation frames of Figure 6 (--variant irr) and Extended
Data Figure 3 (--variant process): optimization trajectories (panel A) and
final proteome allocation (panel B); and of the two panels A side by side
over the proteome allocation (--variant combined).

The same file sits in both figure folders. It loads the folder's snapshot
plots/plot_kin_opt_parameter_sets.py by file path and lets its main() draw
the complete figure once, exactly as make_figure.py does (capture_figures
switches the saving off). Before that, the module's drawing functions are
wrapped so every artist is filed under a group as it is drawn:

  outcome_cell   the cell's scatter / step / axhline / axvspan calls; a
                 trial cloud, incumbent line or highest-IRR circle is
                 assigned to its campaign by the exact colour set_colors gave
                 that campaign (inside _draw_relay_preload: the preload zone)
  draw_burden    every proteome bar segment (Axes.barh) to its row, by the
                 campaign colour of the row's filled / outlined segments;
                 the rest it adds to the figure is panel B's scaffold
  legends        one group per legend row (legend_entry), found by label

Anything these do not account for raises (an artist a cell holds that
none of those calls drew, a legend row or figure-level artist in no group),
so nothing can slip into every frame unnoticed.

Groups (`*` stands for a campaign: the slug of its legend label, e.g.
'isobutanol-yield', 'profitability-uninformed', 'profitability-try-informed',
'baseline-no-optimization'; <cell> is 'irr-left' / 'irr-right' for Figure 6,
'ibo-yield', 'ibo-titer', 'ibo-productivity', 'etoh-yield', 'etoh-titer',
'etoh-productivity' for Extended Data Figure 3):

  A:title                     panel A letter + title
  A:<cell>:axes               a cell's axes (spines, ticks, labels, titles)
  A:<cell>:baseline           its grey dotted baseline line
  A:<cell>:*:cloud            the trial cloud of the cell's own campaign
  A:<cell>:*:incumbent        its thick solid incumbent
  A:<cell>:*:dashed           another campaign's dashed incumbent
  A:<cell>:*:best             a campaign's highest-IRR circle (Figure 6)
  A:irr-right:preload-zone    the grey preload backdrop (Figure 6)
  A:irr-right:preload-trials  the donor-coloured preloaded trials
  A:irr-right:preload-best    the donors' open circles in the preload zone
  legend:*                    one campaign-legend row (the TRY-informed line
                              included)
  key:individual-trials, key:highest-irr-trial, key:incumbent-own,
  key:incumbent-other         the mark-key rows
  B:title                     panel B letter + title
  B:axes                      panel B axes, brackets, dashed lines, labels
  B:sectors                   the sector legend
  B:row:*                     one proteome bar row
  note:uninformed, note:try-traverse, note:scouts, note:try-informed
                              Figure 6's panel-A annotations (box, text,
                              leaders / brace)

The legend frame and title are the background, drawn in every frame. FRAMES
lists the groups each frame adds as fnmatch patterns (each must match at
least one group); re-cut the sequence by editing FRAMES only. Frames are
cumulative and share one canvas, so they overlay exactly; the last frame is
the complete figure, with panel A titled as in the manuscript
(MANUSCRIPT_A_TITLE) and, for Figure 6, the manuscript's annotations above
panel A on a canvas extended upwards by NOTE_SPACE_IN.

--variant combined is one wide figure in two columns (combine): on the
left Figure 6's panel A as drawn above (annotations included) over the
campaign legend; on the right Extended Data Figure 3's panel A (the six
process-metric cells) over the proteome allocation (panel B), each drawn at
the size it has in its own figure. No panel letters, no panel-A title. Its
groups are Figure 6's (the A:title group left out, B:title without the
letter) plus the six process cells' A groups.

Usage (through ../run_plot.py, which redirects the relay manifest path;
../../make_figure_frames.py does this, and for --variant combined
../../make_figure_combined_frames.py, which only the Figure 6 folder has):
    python plot_kin_opt_parameter_sets_frames.py
           --variant irr|process|combined
           [--dpi N] [--formats png [pdf svg ...]]
Frames go to ../../output/frames/.
"""
import argparse
import contextlib
import fnmatch
import importlib.util
import io
import os
import re
import sys

sys.dont_write_bytecode = True   # no __pycache__ under code/ when run directly

import numpy as np
from matplotlib.axes import Axes
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import ArrowStyle, FancyArrowPatch, PathPatch
from matplotlib.path import Path
from matplotlib.transforms import Affine2D, Bbox

HERE = os.path.dirname(os.path.abspath(__file__))
FOLDER = os.path.dirname(os.path.dirname(HERE))
OUT_DIR = os.path.join(FOLDER, 'output', 'frames')
PLOT_SCRIPT = os.path.join(FOLDER, 'code', 'nskinetics', 'biorefineries',
                           'isobutanol', 'plots',
                           'plot_kin_opt_parameter_sets.py')

FRAMES = {
    'irr': [
        ('baseline', ['A:title', 'A:irr-left:axes', 'A:irr-left:baseline',
                      'legend:baseline-no-optimization']),
        ('uninformed-profitability',
         ['A:irr-left:profitability-uninformed:*',
          'legend:profitability-uninformed', 'key:individual-trials',
          'key:highest-irr-trial', 'key:incumbent-own']),
        ('note-uninformed', ['note:uninformed']),
        ('try-campaigns', ['A:irr-left:isobutanol-*', 'A:irr-left:ethanol-*',
                           'legend:isobutanol-*', 'legend:ethanol-*',
                           'key:incumbent-other']),
        ('note-try-traverse', ['note:try-traverse']),
        ('try-scouts', ['A:irr-right:axes', 'A:irr-right:baseline',
                        'A:irr-right:preload-*']),
        ('note-scouts', ['note:scouts']),
        ('try-informed', ['A:irr-right:profitability-try-informed:*',
                          'legend:profitability-try-informed']),
        ('note-try-informed', ['note:try-informed']),
        ('proteome-baseline', ['B:title', 'B:axes', 'B:sectors',
                               'B:row:baseline-no-optimization']),
        ('proteome-uninformed', ['B:row:profitability-uninformed']),
        ('proteome-try', ['B:row:isobutanol-*', 'B:row:ethanol-*']),
        ('proteome-try-informed', ['B:row:profitability-try-informed']),
    ],
    'process': [
        ('baseline', ['A:title', 'A:*:axes', 'A:*:baseline',
                      'legend:baseline-no-optimization']),
        ('isobutanol-campaigns', ['A:*:isobutanol-*:cloud',
                                  'A:*:isobutanol-*:incumbent',
                                  'legend:isobutanol-*',
                                  'key:individual-trials',
                                  'key:incumbent-own']),
        ('ethanol-campaigns', ['A:*:ethanol-*:cloud',
                               'A:*:ethanol-*:incumbent',
                               'legend:ethanol-*']),
        ('cross-campaign-incumbents', ['A:*:isobutanol-*:dashed',
                                       'A:*:ethanol-*:dashed',
                                       'key:incumbent-other']),
        ('uninformed-profitability', ['A:*:profitability-uninformed:dashed',
                                      'legend:profitability-uninformed']),
        ('try-informed-profitability',
         ['A:*:profitability-try-informed:dashed',
          'legend:profitability-try-informed']),
        ('proteome-baseline', ['B:title', 'B:axes', 'B:sectors',
                               'B:row:baseline-no-optimization']),
        ('proteome-try', ['B:row:isobutanol-*', 'B:row:ethanol-*']),
        ('proteome-uninformed', ['B:row:profitability-uninformed']),
        ('proteome-try-informed', ['B:row:profitability-try-informed']),
    ],
    # Figure 6's panel A (irr-left / irr-right cells, notes, legend), Extended
    # Data Figure 3's (the ibo-* / etoh-* cells) and the proteome rows: one
    # frame per step of the story, each ending on its note. A campaign enters
    # every panel in the same frame: its profitability, its incumbent in the
    # six process cells and its proteome row
    'combined': [
        ('baseline', ['A:irr-left:axes', 'A:irr-left:baseline',
                      'A:ibo-*:axes', 'A:etoh-*:axes', 'A:ibo-*:baseline',
                      'A:etoh-*:baseline', 'legend:baseline-no-optimization',
                      'B:title', 'B:axes', 'B:sectors',
                      'B:row:baseline-no-optimization']),
        ('uninformed-profitability',
         ['A:irr-left:profitability-uninformed:*',
          'A:ibo-*:profitability-uninformed:dashed',
          'A:etoh-*:profitability-uninformed:dashed',
          'legend:profitability-uninformed', 'key:individual-trials',
          'key:highest-irr-trial', 'key:incumbent-own', 'key:incumbent-other',
          'note:uninformed', 'B:row:profitability-uninformed']),
        ('try-campaigns',
         ['A:ibo-*:cloud', 'A:etoh-*:cloud', 'A:ibo-*:incumbent',
          'A:etoh-*:incumbent', 'A:ibo-*:isobutanol-*:dashed',
          'A:ibo-*:ethanol-*:dashed', 'A:etoh-*:isobutanol-*:dashed',
          'A:etoh-*:ethanol-*:dashed', 'A:irr-left:isobutanol-*',
          'A:irr-left:ethanol-*', 'legend:isobutanol-*', 'legend:ethanol-*',
          'note:try-traverse', 'B:row:isobutanol-*', 'B:row:ethanol-*']),
        ('try-scouts', ['A:irr-right:axes', 'A:irr-right:baseline',
                        'A:irr-right:preload-*', 'note:scouts']),
        ('try-informed', ['A:irr-right:profitability-try-informed:*',
                          'A:ibo-*:profitability-try-informed:dashed',
                          'A:etoh-*:profitability-try-informed:dashed',
                          'legend:profitability-try-informed',
                          'note:try-informed',
                          'B:row:profitability-try-informed']),
    ],
}
# the plot script's figure each variant starts from
VIEW = {'irr': 'irr', 'process': 'process', 'combined': 'irr'}

# panel A's title as the plot script writes it, and as the manuscript shows
# it in both figures (changed in the manuscript's PowerPoint source)
SCRIPT_A_TITLE = {'irr': 'Profitability optimization trajectories',
                  'process': 'Process-level optimization trajectories'}
MANUSCRIPT_A_TITLE = 'Optimization trajectories'

# the mark-key rows, by their label in the plot script's _style_handles
KEY_ROWS = {'Individual trials': 'key:individual-trials',
            'Highest-IRR trial': 'key:highest-irr-trial',
            'Incumbent (campaign optimizing the shown metric)':
                'key:incumbent-own',
            'Incumbent (other campaign)': 'key:incumbent-other'}

# --- Figure 6's panel-A annotations (as in the manuscript figure) -----------
NOTES = {
    'note:uninformed':
        'Uninformed, the profitability campaign finds a safe ethanol-only\n'
        'region early but cannot find more profitable points elsewhere',
    'note:try-traverse':
        'Campaigns with process-level TRY objectives can\n'
        'traverse more profitable points but cannot recognize them',
    'note:scouts':
        'Campaigns with process-level TRY objectives are\n'
        'used as scouts to chart an initial map of profitability',
    'note:try-informed':
        'The TRY-informed profitability campaign finds a more profitable\n'
        'region co-producing ethanol & isobutanol and optimizes further',
}
# the campaign whose highest-IRR circle note:uninformed points at, and those
# whose circles note:try-traverse points at
UNINFORMED = 'Profitability (uninformed)'
TRAVERSE_TARGETS = ('Isobutanol productivity', 'Isobutanol yield',
                    'Isobutanol titer')
# note:try-informed's first arrow points at the TRY-informed incumbent where
# it has climbed this fraction of the way from its first simulated value to
# its final one (the knee of the line just after the preload zone)
KNEE_FRACTION = 0.8
# the canvas grows by this much at the top to hold the notes (an integer
# number of pixels at any integer dpi, so the original canvas below stays on
# the same pixel grid); panel A's letter and title move up into it
NOTE_SPACE_IN = 1.0
# layout [in]: lower note row above the cells, gap between the two rows,
# title baseline above the upper row, first note from the canvas edge, the
# lower-left note's clearance from the first note's leader, the lower-right
# note's indent from the right cell's left edge, text inset, least gap
# between two notes
CELL_GAP_IN, ROW_GAP_IN, TITLE_GAP_IN = 0.15, 0.13, 0.17
EDGE_IN, LEADER_CLEAR_IN, SCOUTS_INDENT_IN = 0.15, 0.10, 0.05
NOTE_PAD_IN, NOTE_GAP_IN = 0.075, 0.05
# the leader bus of note:try-traverse sits this fraction of the way down from
# the cell top to the highest target; that of note:try-informed this far
# above the cell top, under the scouts brace's height
BUS_DROP_FRACTION, BUS_ABOVE_IN = 0.4, 0.075
BRACE_Y0_IN, BRACE_H_IN, BRACE_R_IN = 0.02, 0.075, 0.04
NOTE_EDGE, NOTE_LW = '0.4', 0.8
LEADER_COLOR, LEADER_LW = '0.3', 0.8
ARROW = ArrowStyle('-|>', head_length=0.5, head_width=0.22)
ARROW_SCALE = 10
TIP_GAP_PT = 1.5        # arrowheads stop this far short of their target

# --- the combined variant ---------------------------------------------------
# what Figure 6 draws that the combined figure leaves out (panel A's title)
# or draws again elsewhere (panel B)
COMBINED_DROPPED = ('A:title', 'B:axes', 'B:sectors', 'B:row:*')
# from the bottom of the process cells to the top of the proteome axes [in]:
# room for the cells' trial-axis titles and the proteome panel's brackets
COMBINED_PANEL_GAP_IN = 0.95
# margin of the combined canvas above and below what is drawn [in]; left and
# right it keeps the figures' own margins
COMBINED_MARGIN_IN = 0.12


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pf = _load('_progressive_frames', os.path.join(HERE, 'progressive_frames.py'))


def _slug(label):
    return re.sub(r'[^0-9a-z]+', '-', label.lower()).strip('-')


def instrument(m, groups):
    """Wrap the plot module's drawing functions so the artists they draw are
    filed into groups ({name: [artists]}) as they are created. Returns a dict
    that receives plot()'s sets ('sets'), the frame drawn around the
    campaign legend and mark key ('key frame'), which is background, per
    cell (name, axes, ids of the children it had before it was drawn on)
    ('cells'), for check_all_filed, and the arguments plot() drew panel B
    with ('burden': sets, colors, args, kwargs), for combine."""
    seen = {'cells': []}
    plot, frame_keys = m.plot, m._frame_keys
    outcome_cell, draw_preload = m.outcome_cell, m._draw_relay_preload
    draw_burden = m.draw_burden

    def add(name, artists):
        groups.setdefault(name, []).extend(artists)

    def plot_(sets, *args, **kwargs):
        seen['sets'] = sets
        return plot(sets, *args, **kwargs)

    def frame_keys_(fig, *args, **kwargs):
        with pf.recording(fig) as new:
            result = frame_keys(fig, *args, **kwargs)
        seen['key frame'] = new
        return result

    def outcome_cell_(ax, sets, colors, col, *args, **kwargs):
        owner = {}
        for s in sets:
            if owner.setdefault(colors[id(s)], s) is not s:
                raise ValueError(f'two campaigns share colour {colors[id(s)]}')
        if col == 'IRR':
            relay = any(s.get('is_relay') for s in sets)
            cell = 'A:irr-right' if relay else 'A:irr-left'
        else:
            cell = 'A:' + _slug(col)
        in_preload = []

        def campaign(colour):
            if colour not in owner:
                raise KeyError(f'{cell}: colour {colour!r} is no campaign of '
                               'this cell')
            return owner[colour]

        def group(method, kw):
            if method == 'axhline':
                return f'{cell}:baseline'
            if method == 'axvspan':
                return f'{cell}:preload-zone'
            if method == 'step':
                s = campaign(kw['color'])
                own = m._owns_outcome(s.get('objective'), col)
                return f'{cell}:{_slug(s["label"])}:' + (
                    'incumbent' if own else 'dashed')
            circle = kw.get('facecolors') == 'white'      # scatter
            if in_preload:
                return f'{cell}:preload-' + ('best' if circle else 'trials')
            s = campaign(kw['edgecolors'] if circle else kw['color'])
            return f'{cell}:{_slug(s["label"])}:' + (
                'best' if circle else 'cloud')

        def tap(method):
            draw = getattr(ax, method)

            def tapped(*a, **kw):
                art = draw(*a, **kw)
                add(group(method, kw), art if isinstance(art, list) else [art])
                return art
            return tapped

        def draw_preload_(*a, **kw):
            in_preload.append(True)
            try:
                return draw_preload(*a, **kw)
            finally:
                in_preload.clear()

        methods = ('scatter', 'step', 'axhline', 'axvspan')
        for method in methods:
            setattr(ax, method, tap(method))
        m._draw_relay_preload = draw_preload_
        seen['cells'].append((cell, ax, {id(a) for a in ax.get_children()}))
        try:
            return outcome_cell(ax, sets, colors, col, *args, **kwargs)
        finally:
            m._draw_relay_preload = draw_preload
            for method in methods:
                delattr(ax, method)
            add(f'{cell}:axes', [ax])

    def draw_burden_(fig, gs_cell, sets, colors, *args, **kwargs):
        seen.setdefault('burden', (sets, colors, args, kwargs))
        bars, barh = [], Axes.barh

        def barh_(self, y, *a, **kw):
            container = barh(self, y, *a, **kw)
            bars.append((y, kw, container))
            return container

        Axes.barh = barh_
        try:
            with pf.recording(fig) as new:
                result = draw_burden(fig, gs_cell, sets, colors, *args,
                                     **kwargs)
        finally:
            Axes.barh = barh
        # a row is the bars at one y; its owner is the set whose colour fills
        # or outlines them (the empty slack segment carries no colour)
        by_colour = {colors[id(s)]: s for s in sets}
        rows = {}
        for y, kw, container in bars:
            rows.setdefault(y, []).append((kw, container))
        owners = set()
        for y, segments in rows.items():
            found = {id(by_colour[c]): by_colour[c] for kw, _ in segments
                     for c in (kw.get('facecolor'), kw.get('edgecolor'))
                     if isinstance(c, str) and c in by_colour}
            if len(found) != 1:
                raise ValueError(f'proteome row at y={y}: {len(found)} owners')
            (s,) = found.values()
            owners.add(id(s))
            add(f'B:row:{_slug(s["label"])}',
                [p for _, container in segments for p in container.patches])
        if owners != {id(s) for s in sets}:
            raise ValueError('not every set has a proteome row')
        for a in new:
            add('B:sectors' if isinstance(a, Legend) else 'B:axes', [a])
        return result

    m.plot, m._frame_keys = plot_, frame_keys_
    m.outcome_cell, m.draw_burden = outcome_cell_, draw_burden_
    return seen


def _one(artists, what):
    if len(artists) != 1:
        raise ValueError(f'expected one {what}, found {len(artists)}')
    return artists[0]


def file_titles_and_legends(m, fig, sets, groups, variant):
    """Group the panel letters / titles and the legend rows (after plot)."""
    def text(s):
        return _one([t for t in fig.texts if t.get_text() == s], repr(s))

    a_title = text(SCRIPT_A_TITLE[variant])
    groups['A:title'] = [text('A'), a_title]
    groups['B:title'] = [text('B'), text('Final proteome allocation')]
    keys = [leg for leg in fig.legends if leg not in groups['B:sectors']]
    a_sets = [s for s in sets if not s.get('is_relay')]

    def row(label):
        return _one([leg for leg in keys if label in
                     [t.get_text() for t in leg.get_texts()]],
                    f'legend with the row {label!r}')

    for s in sets:
        label = (m.relay_legend_label(s, a_sets) if s.get('is_relay')
                 else s['label'])
        groups[f'legend:{_slug(s["label"])}'] = pf.legend_entry(row(label),
                                                                label)
    for leg in keys:
        for t in leg.get_texts():
            if t.get_text() in KEY_ROWS:
                groups[KEY_ROWS[t.get_text()]] = pf.legend_entry(
                    leg, t.get_text())
    return a_title


def check_all_filed(fig, groups, background, cells):
    """Raise if the figure draws anything that is in no group and is not
    `background` (the key frame): it would be shown in every frame. Of a
    campaign legend or mark key only the title may be unfiled. Likewise for
    every panel-A cell: an artist it gained that no tapped call drew would
    be shown with the cell's axes, before its frame."""
    filed = {id(a) for artists in groups.values() for a in artists}
    for cell, ax, scaffold in cells:
        loose = [a for a in ax.get_children()
                 if id(a) not in scaffold and id(a) not in filed]
        if loose:
            raise ValueError(f'{cell}: artists in no group: {loose}')
    for a in fig.get_children():
        if a is fig.patch or id(a) in filed or not a.get_visible():
            continue
        if isinstance(a, Legend):
            loose = [t.get_text() for t in a.get_texts() if id(t) not in filed]
            if loose:
                raise ValueError(f'legend rows in no group: {loose}')
        elif not any(a is b for b in background):
            raise ValueError(f'{a!r} is in no group')


def _inches(fig):
    """Inches from the figure's lower-left corner -> display. Unlike
    fig.dpi_scale_trans it goes through fig.transFigure, so what is placed
    with it stays with the figure when a frame is saved with a bbox that
    does not start at the figure's corner (the combined variant's crop)."""
    W, H = fig.get_size_inches()
    return Affine2D().scale(1 / W, 1 / H) + fig.transFigure


def add_notes(m, fig, sets, groups):
    """Figure 6: grow the canvas by NOTE_SPACE_IN at the top, move panel A's
    letter and title up into it and add the four annotations between the
    title and the cells, every arrow aimed at data computed by the plot
    module. Returns the bbox [in] to save every frame with."""
    W, H = fig.get_size_inches()
    inches = _inches(fig)
    to_in = inches.inverted()
    left = _one(groups['A:irr-left:axes'], 'left IRR cell')
    right = _one(groups['A:irr-right:axes'], 'right IRR cell')
    cell_top = left.get_position().y1 * H
    r_x0, r_x1 = right.get_position().x0 * W, right.get_position().x1 * W

    def at(ax, x, y):                           # data -> figure inches
        return to_in.transform(ax.transData.transform((x, y)))

    best = {s['label']: (s, x, y) for s, x, y in m.best_irr_points(sets)}
    relay = _one([s for s in sets if s.get('is_relay')], 'relay set')
    n_pre = relay['n_preloaded']
    mark_r = np.sqrt(m.BEST_MARK_SIZE) / 2 + m.BEST_MARK_LW / 2   # [pt]

    def arrow(p0, p1, gap_pt):
        return FancyArrowPatch(p0, p1, arrowstyle=ARROW,
                               mutation_scale=ARROW_SCALE, shrinkA=0,
                               shrinkB=gap_pt + TIP_GAP_PT, color=LEADER_COLOR,
                               lw=LEADER_LW, transform=inches)

    def line(xs, ys):
        return Line2D(xs, ys, color=LEADER_COLOR, lw=LEADER_LW,
                      solid_capstyle='butt', transform=inches)

    # the boxes: measure the text, then place each box by its edges
    fs = m.FONTS['legend'] + 1
    pad = NOTE_PAD_IN * 72 / fs                 # boxstyle pad, in font sizes
    texts = {name: fig.text(0, 0, s, fontsize=fs, ha='left', va='bottom',
                            transform=inches, zorder=6,
                            bbox=dict(boxstyle=f'square,pad={pad}',
                                      fc='white', ec=NOTE_EDGE, lw=NOTE_LW))
             for name, s in NOTES.items()}
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    size = {name: to_in.transform(t.get_window_extent(renderer).size)
            - to_in.transform((0, 0)) + 2 * NOTE_PAD_IN
            for name, t in texts.items()}
    h = max(sz[1] for sz in size.values())
    low = cell_top + CELL_GAP_IN                # lower row's box bottom
    up = low + h + ROW_GAP_IN                   # upper row's box bottom
    x_unin = at(left, *best[UNINFORMED][1:])[0]
    lefts = {'note:uninformed': EDGE_IN,
             'note:try-traverse': x_unin + LEADER_CLEAR_IN,
             'note:scouts': r_x0 + SCOUTS_INDENT_IN,
             'note:try-informed': r_x1 - size['note:try-informed'][0]}
    bottoms = {'note:uninformed': up, 'note:try-traverse': low,
               'note:scouts': low, 'note:try-informed': up}
    box = {}
    for name, t in texts.items():
        x0, y0 = lefts[name], bottoms[name]
        box[name] = (x0, y0, x0 + size[name][0], y0 + size[name][1])
        t.set_position((x0 + NOTE_PAD_IN, y0 + NOTE_PAD_IN))
        groups[name] = [t]

    # panel A's letter and title move up above the notes
    title_y = (up + h + TITLE_GAP_IN) / H
    for t in groups['A:title']:
        t.set_y(title_y)

    # note:uninformed -- straight down to the uninformed campaign's circle
    _, xu, yu = best[UNINFORMED]
    p = at(left, xu, yu)
    groups['note:uninformed'].append(arrow((p[0], up), p, mark_r))

    # note:try-traverse -- a trunk down to the leftmost target, a bus across
    # the top of the cell to the others, one drop per target
    pts = sorted((at(left, *best[label][1:]) for label in TRAVERSE_TARGETS),
                 key=lambda q: q[0])
    x_bus = [pts[0][0], pts[-1][0]]
    y_bus = cell_top - BUS_DROP_FRACTION * (cell_top - max(q[1] for q in pts))
    notes = [arrow((pts[0][0], low), pts[0], mark_r), line(x_bus, [y_bus] * 2)]
    notes += [arrow((q[0], y_bus), q, mark_r) for q in pts[1:]]
    groups['note:try-traverse'] += notes

    # note:scouts -- a curly brace over the preload zone (trials 0..N)
    bx0, bx1 = at(right, 0, 0)[0], at(right, n_pre, 0)[0]
    groups['note:scouts'].append(_brace(fig, bx0, bx1, cell_top + BRACE_Y0_IN))

    # note:try-informed -- a trunk down to the TRY-informed campaign's circle
    # and, along a bus just above the cell, a drop to its incumbent's knee
    _, xr, yr = best[m.DEFAULT_RELAY_LABEL]
    tx = np.asarray(relay['traj_x'], dtype=float)
    ty = np.asarray(relay['traj']['IRR'], dtype=float)
    tx, ty = tx[tx >= n_pre], ty[tx >= n_pre]
    if not (np.all(np.isfinite(ty)) and ty[-1] > ty[0]):
        raise ValueError('the TRY-informed incumbent does not climb from its '
                         'seed: no knee for note:try-informed to point at')
    k = int(np.argmax(ty >= ty[0] + KNEE_FRACTION * (ty[-1] - ty[0])))
    pr, pk = at(right, xr, yr), at(right, tx[k], ty[k])
    y_knee_bus = cell_top + BUS_ABOVE_IN
    incumbent = _one(
        groups[f'A:irr-right:{_slug(relay["label"])}:incumbent'],
        'TRY-informed incumbent line')
    groups['note:try-informed'] += [
        arrow((pr[0], up), pr, mark_r),
        line([pk[0], pr[0]], [y_knee_bus] * 2),
        arrow((pk[0], y_knee_bus), pk, incumbent.get_linewidth() / 2)]

    for a in [a for name in NOTES for a in groups[name][1:]]:
        fig.add_artist(a)
    _check_note_layout(box, W, H, x_unin, pts[0][0], pr[0], bx1, pk[0],
                       title_y, groups['A:title'], fig)
    return Bbox([[0, 0], [W, H + NOTE_SPACE_IN]])


def combine(m, fig, sets, groups, seen):
    """The combined variant, from the finished Figure 6 (notes added): a
    second column one figure width to the right, i.e. each panel at the x it
    has in its own figure, holding

      - Extended Data Figure 3's panel A (the six process-metric cells,
        m.draw_process_outcomes) in a box the size of Figure 6's panel A,
        its top level with the notes';
      - under it panel B, drawn again by m.draw_burden with the arguments
        plot() gave it, its axes the size they have in Figure 6. Its sector
        legend, which draw_burden hangs at the figure's left edge, and its
        title (without the letter) move over with it.

    Panel A's title and the first drawing of panel B (COMBINED_DROPPED) are
    hidden for good; the new artists are filed like any other (instrument).
    Returns the bbox [in] to save every frame with: both columns, with
    COMBINED_MARGIN_IN above and below. Raises if anything drawn falls
    outside it."""
    W, H = fig.get_size_inches()
    to_fig = fig.transFigure.inverted()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pair = Bbox.union([_one(groups[f'A:{cell}:axes'], cell).get_position()
                       for cell in ('irr-left', 'irr-right')])
    bars = Bbox.union([a.get_position() for a in groups['B:axes']
                       if isinstance(a, Axes)])
    letter, title = groups['B:title']
    title_rise = title.get_position()[1] - bars.y1
    for name in [g for g in groups for p in COMBINED_DROPPED
                 if fnmatch.fnmatchcase(g, p)]:
        for a in groups.pop(name):
            a.set_visible(False)

    top = max(to_fig.transform(
        groups[name][0].get_bbox_patch().get_window_extent(renderer))[1, 1]
        for name in NOTES)
    cells = fig.add_gridspec(1, 1, left=pair.x0 + 1, right=pair.x1 + 1,
                             top=top, bottom=top - pair.height)
    m.draw_process_outcomes(fig, cells[0], sets, m.set_colors(sets),
                            m._trial_xmax(sets))
    bars_top = top - pair.height - COMBINED_PANEL_GAP_IN / H
    burden = fig.add_gridspec(1, 1, left=bars.x0 + 1, right=bars.x1 + 1,
                              top=bars_top, bottom=bars_top - bars.height)
    b_sets, b_colors, b_args, b_kwargs = seen['burden']
    m.draw_burden(fig, burden[0], b_sets, b_colors, *b_args, **b_kwargs)
    sectors = _one(groups['B:sectors'], 'sector legend')
    anchor = sectors.get_bbox_to_anchor().transformed(to_fig)
    sectors.set_bbox_to_anchor((anchor.x0 + 1, anchor.y0),
                               transform=fig.transFigure)
    letter.set_visible(False)
    title.set_position((letter.get_position()[0] + 1, bars_top + title_rise))
    groups['B:title'] = [title]

    fig.canvas.draw()
    drawn = fig.get_tightbbox(renderer)
    # whole hundredths of an inch: whole pixels at 100, 200, 300, 600 dpi
    crop = Bbox([[0, np.floor((drawn.y0 - COMBINED_MARGIN_IN) * 100) / 100],
                 [2 * W, np.ceil((drawn.y1 + COMBINED_MARGIN_IN) * 100) / 100]])
    if not (crop.x0 <= drawn.x0 and drawn.x1 <= crop.x1):
        raise ValueError(f'the combined figure {drawn.extents} [in] does not '
                         f'fit its canvas {crop.extents}')
    return crop


def _brace(fig, x0, x1, y0):
    """A horizontal curly brace from x0 to x1 [in], ends at y0, tip up."""
    xm, ym, yt, r = (x0 + x1) / 2, y0 + BRACE_H_IN / 2, y0 + BRACE_H_IN, \
        BRACE_R_IN
    verts = [(x0, y0), (x0, ym), (x0 + r, ym), (xm - r, ym), (xm, ym),
             (xm, yt), (xm, ym), (xm + r, ym), (x1 - r, ym), (x1, ym),
             (x1, y0)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.LINETO,
             Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.LINETO,
             Path.CURVE3, Path.CURVE3]
    return PathPatch(Path(verts, codes), fill=False, ec=LEADER_COLOR,
                     lw=LEADER_LW, transform=_inches(fig))


def _check_note_layout(box, W, H, x_unin, x_traverse, x_relay, x_brace,
                       x_knee, title_y, title, fig):
    """Raise if the notes overlap each other or the title, leave the canvas
    or have a leader run through another note (data could move them)."""
    names = list(box)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            (ax0, ay0, ax1, ay1), (bx0, by0, bx1, by1) = box[a], box[b]
            g = NOTE_GAP_IN
            if ax0 < bx1 + g and bx0 < ax1 + g and ay0 < by1 + g \
                    and by0 < ay1 + g:
                raise ValueError(f'{a} and {b} overlap or touch')
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    top = max(t.get_window_extent(r).y1 for t in title) / fig.dpi
    if top > H + NOTE_SPACE_IN - 0.05 or \
            title_y * H < max(b[3] for b in box.values()):
        raise ValueError('panel A title does not fit above the notes')
    if min(b[0] for b in box.values()) < 0.05 or \
            max(b[2] for b in box.values()) > W - 0.05:
        raise ValueError('a note leaves the canvas')
    inside = [(x_unin, 'note:uninformed'), (x_traverse, 'note:try-traverse'),
              (x_relay, 'note:try-informed')]
    for x, name in inside:
        if not box[name][0] + 0.05 < x < box[name][2] - 0.05:
            raise ValueError(f'the leader of {name} misses its box')
    if not (x_unin < box['note:try-traverse'][0]
            and x_relay > box['note:scouts'][2] + NOTE_GAP_IN
            and x_knee > x_brace + NOTE_GAP_IN):
        raise ValueError('a leader runs through another note')


def expand(frames, groups):
    """FRAMES with every fnmatch pattern replaced by the groups it matches."""
    out = []
    for name, patterns in frames:
        names = []
        for p in patterns:
            hits = [g for g in groups if fnmatch.fnmatchcase(g, p)]
            if not hits:
                raise KeyError(f'frame {name!r}: no group matches {p!r}')
            names += [g for g in hits if g not in names]
        out.append((name, names))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--variant', choices=sorted(FRAMES), required=True,
                    help='the figure (set by the make_figure_*frames.py '
                         'wrapper): irr = Figure 6, process = Extended Data '
                         'Figure 3, combined = their two panels A side by '
                         'side over the proteome allocation')
    ap.add_argument('--dpi', type=int, default=300,
                    help='resolution of the frames (default: %(default)s, '
                         "as make_figure.py's figure)")
    ap.add_argument('--formats', nargs='+', default=['png'],
                    help='file formats to write (default: png)')
    args = ap.parse_args()

    view = VIEW[args.variant]
    m = _load('_plot_kin_opt_parameter_sets', PLOT_SCRIPT)
    groups = {}
    seen = instrument(m, groups)
    # main()'s console report and "wrote ..." lines are dropped: with saving
    # switched off it writes nothing (its --out-dir only names the stems)
    with contextlib.redirect_stdout(io.StringIO()):
        _, figures = pf.capture_figures(
            m.main, ['--variants', view, '--dpi', str(args.dpi),
                     '--out-dir', OUT_DIR])
    fig = _one(figures, 'figure')
    a_title = file_titles_and_legends(m, fig, seen['sets'], groups, view)
    a_title.set_text(MANUSCRIPT_A_TITLE)
    bbox = add_notes(m, fig, seen['sets'], groups) if view == 'irr' else None
    if args.variant == 'combined':
        bbox = combine(m, fig, seen['sets'], groups, seen)
    check_all_filed(fig, groups, seen['key frame'], seen['cells'])
    pf.save_frames(fig, groups, expand(FRAMES[args.variant], groups),
                   OUT_DIR, f'kinBO_parameter_sets_{args.variant}', args.dpi,
                   formats=args.formats, bbox_inches=bbox)


if __name__ == '__main__':
    main()
