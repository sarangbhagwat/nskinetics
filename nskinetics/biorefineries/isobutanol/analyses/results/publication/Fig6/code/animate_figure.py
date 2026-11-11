# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Build GIFs of Figure 6 (panel A: IRR optimization trajectories; panel B:
final proteome allocation), one GIF per step towards the full figure, for
presentations.

The figure is drawn by the unmodified snapshot
nskinetics/biorefineries/isobutanol/plots/plot_kin_opt_parameter_sets.py
(--variants irr), recorded and grouped with kin_opt_animation.py and
animated with gif_engine.Storyboard. On top of the static figure, the GIFs
carry the four panel-A annotation boxes of the manuscript figure (added
there by hand), so the canvas gets headroom above panel A for them, the
panel title reads "Optimization trajectories" and sits at the top, and the
empty band below panel B is trimmed. Everything else is drawn as
make_figure.py draws it.

Storyboard:

    start   the panel-A title and the empty left IRR cell
    1  baseline                 the baseline IRR reference line
    2  profitability-uninformed the uninformed profitability campaign's trials
                                and incumbent, drawn along the trial axis,
                                its highest-IRR trial, annotation 1
    3  process-level-incumbents the six process-level campaigns' incumbents
    4  process-level-best-trials their highest-IRR trials, annotation 2
    5  scout-trials             the right cell: the 1000 preloaded trials of
                                the process-level campaigns (scouts) and
                                their marked trials, annotation 4 (bracket)
    6  profitability-try-informed  the TRY-informed campaign's own trials and
                                incumbent from trial 1000, its highest-IRR
                                trial, annotation 3
    7-11                        panel B, as in kin_opt_animation.py

Legend entries appear with what they describe. The annotation anchors (the
highest-IRR markers, the start of the TRY-informed incumbent, the preload
zone) are taken from the drawn artists. Run through ../make_gifs.py; see
gif_engine.py for the files written and the GIF options.

Sim-safe: reads the campaign CSVs, two workbooks and the antimony text only.
"""
import os
import sys

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, PathPatch, Rectangle
from matplotlib.path import Path
from matplotlib.transforms import Bbox

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import gif_engine as ge  # noqa: E402
import kin_opt_animation as ka  # noqa: E402

STEM = 'kinBO_parameter_sets_irr'
TITLE = 'Optimization trajectories'

ANNOTATION_1 = ('Uninformed, the profitability campaign finds a safe '
                'ethanol-only\nregion early but cannot find more profitable '
                'points elsewhere')
ANNOTATION_2 = ('Campaigns with process-level TRY objectives can\ntraverse '
                'more profitable points but cannot recognize them')
ANNOTATION_3 = ('The TRY-informed profitability campaign finds a more '
                'profitable\nregion co-producing ethanol & isobutanol and '
                'optimizes further')
ANNOTATION_4 = ('Campaigns with process-level TRY objectives are\nused as '
                'scouts to chart an initial map of profitability')

# Layout of the annotations, in inches, as in the manuscript figure. Heights
# are measured up from the top edge of the IRR cells.
BOX_PAD_X, BOX_PAD_Y = 0.07, 0.075     # text to box edge
LOWER_ROW_GAP = 0.155                  # cells -> bottom of boxes 2 and 4
ROW_GAP = 0.11                         # boxes 2/4 -> boxes 1/3
TITLE_GAP = 0.17                       # boxes 1/3 -> panel title baseline
TOP_MARGIN = 0.24                      # panel title baseline -> canvas top
BOX_1_LEFT = 0.17                      # from the canvas's left edge
STEM_TO_BOX = 0.08                     # annotation 1's stem -> box 2
BRACE_Y, BRACE_DROP = 0.072, 0.04      # annotation 4's brace above the cells
BRACE_GAP = 0.02                       # brace tip -> box 4
BRANCH_DROP = 0.26                     # annotation 2's branch below the cells
TRY_RISE = 0.9                         # annotation 3 points where the
#                                        TRY-informed incumbent has made this
#                                        share of its first climb
LEADER = dict(color='0.35', lw=0.8)
BOX_EDGE = dict(edgecolor='0.35', lw=0.8)
ANNOTATION_ZORDER = 10


# -----------------------------------------------------------------------------
# geometry
# -----------------------------------------------------------------------------
def to_inches(fig, ax, xy):
    """Data point(s) of ``ax`` -> figure inches."""
    return ax.transData.transform(xy) / fig.dpi


def extent_inches(fig, artist):
    fig.canvas.draw()
    box = artist.get_window_extent(fig.canvas.get_renderer())
    return Bbox(box.get_points() / fig.dpi)


# -----------------------------------------------------------------------------
# annotation pieces (all in figure inches)
# -----------------------------------------------------------------------------
class Annotation:
    """A text box and its leader lines, collected for one storyboard group."""

    def __init__(self, fig, text, fontsize):
        self.fig = fig
        self.artists = []
        self.text = fig.text(0, 0, text, ha='left', va='bottom',
                             fontsize=fontsize, zorder=ANNOTATION_ZORDER + 1)
        self.artists.append(self.text)
        self.box = None
        self.inches = None

    def size(self):
        """(width, height) of the box, in inches."""
        ext = extent_inches(self.fig, self.text)
        return ext.width + 2 * BOX_PAD_X, ext.height + 2 * BOX_PAD_Y

    def place(self, x0, y0):
        """Put the box's bottom-left corner at (x0, y0), once the canvas
        has its final size."""
        width, height = self.size()
        self.inches = ka.inch_transform(self.fig)
        self.text.set_transform(self.inches)
        self.text.set_position((x0 + BOX_PAD_X, y0 + BOX_PAD_Y))
        self.box = Bbox.from_bounds(x0, y0, width, height)
        frame = Rectangle((x0, y0), width, height, facecolor='white',
                          transform=self.inches,
                          zorder=ANNOTATION_ZORDER, **BOX_EDGE)
        self.fig.add_artist(frame)
        self.artists.append(frame)

    def line(self, *points):
        xs, ys = zip(*points)
        line = Line2D(xs, ys, transform=self.inches,
                      solid_capstyle='butt', zorder=ANNOTATION_ZORDER,
                      **LEADER)
        self.fig.add_artist(line)
        self.artists.append(line)

    def arrow(self, start, tip, shrink=0.0):
        """An arrow from start to tip, stopping ``shrink`` points short."""
        arrow = FancyArrowPatch(
            start, tip, transform=self.inches,
            arrowstyle='-|>,head_length=0.5,head_width=0.22',
            mutation_scale=10, shrinkA=0, shrinkB=shrink,
            zorder=ANNOTATION_ZORDER, **LEADER)
        arrow.set_facecolor(LEADER['color'])
        self.fig.add_artist(arrow)
        self.artists.append(arrow)

    def brace(self, x0, x1, y, drop, peak):
        """A flat curly brace from x0 to x1 at height y: its ends curl down
        by ``drop``, its middle rises to a point at height ``peak``."""
        xm = 0.5 * (x0 + x1)
        r = drop
        vertices = [(x0, y - drop), (x0, y), (x0 + r, y),
                    (xm - r, y), (xm, y), (xm, peak),
                    (xm, y), (xm + r, y),
                    (x1 - r, y), (x1, y), (x1, y - drop)]
        codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3,
                 Path.LINETO, Path.CURVE3, Path.CURVE3,
                 Path.CURVE3, Path.CURVE3,
                 Path.LINETO, Path.CURVE3, Path.CURVE3]
        patch = PathPatch(Path(vertices, codes), fill=False,
                          transform=self.inches,
                          capstyle='butt', joinstyle='round',
                          zorder=ANNOTATION_ZORDER,
                          edgecolor=LEADER['color'], lw=LEADER['lw'])
        self.fig.add_artist(patch)
        self.artists.append(patch)


def marker_radius_points(collection):
    """Radius of a scatter's circle marker in points, its edge included."""
    return 0.5 * np.sqrt(collection.get_sizes()[0]) \
        + 0.5 * collection.get_linewidths()[0]


def marker_xy(fig, cell, collection):
    return to_inches(fig, cell.ax, collection.get_offsets()[0])


def try_informed_climb_point(fig, cell, line):
    """Where the TRY-informed incumbent, which starts after its preloaded
    trials, has made TRY_RISE of its first climb (its rise up to the first
    trial that comes within 1 % of the final value): (point in inches, the
    line's half width in points)."""
    x, y = np.asarray(line.get_xdata(), float), np.asarray(line.get_ydata(),
                                                           float)
    final = y[-1]
    plateau = int(np.argmax(y >= final - 0.01 * abs(final)))
    target = y[0] + TRY_RISE * (y[plateau] - y[0])
    i = int(np.argmax(y >= target))
    return to_inches(fig, cell.ax, (x[i], y[i])), 0.5 * line.get_linewidth()


def add_annotations(fig, recorder, left, right):
    """Make room above panel A, retitle it and draw the four annotations.
    Returns {group name: artists}."""
    # the size of the campaign legend's entries
    fontsize = recorder.plot.FONTS['legend'] + 1
    notes = {name: Annotation(fig, text, fontsize) for name, text in (
        ('annotation 1', ANNOTATION_1), ('annotation 2', ANNOTATION_2),
        ('annotation 3', ANNOTATION_3), ('annotation 4', ANNOTATION_4))}
    upper = max(notes[n].size()[1] for n in ('annotation 1', 'annotation 3'))
    lower = max(notes[n].size()[1] for n in ('annotation 2', 'annotation 4'))

    width, height = fig.get_size_inches()
    top = left.ax.get_position().y1 * height
    needed = top + LOWER_ROW_GAP + lower + ROW_GAP + upper + TITLE_GAP \
        + TOP_MARGIN
    ka.grow_canvas_top(fig, needed - height)

    # the panel title and its letter, at the top
    title_y = needed - TOP_MARGIN
    for content in ('A', 'Profitability optimization trajectories'):
        text = ka.figure_text(fig, content)
        x = text.get_position()[0] * width     # a figure fraction before
        text.set_transform(ka.inch_transform(fig))
        text.set_position((x, title_y))
    ka.figure_text(fig, 'Profitability optimization trajectories') \
        .set_text(TITLE)

    fig.canvas.draw()
    cells_top = left.ax.get_position().y1 * fig.get_size_inches()[1]
    lower_y = cells_top + LOWER_ROW_GAP
    upper_y = lower_y + lower + ROW_GAP
    right_box = extent_inches(fig, right.ax.patch)

    # 1: the uninformed campaign's highest-IRR trial
    uninformed = recorder.labels('uninformed')[0]
    best_1 = left.find_one('best', uninformed)
    m1 = marker_xy(fig, left, best_1)
    note = notes['annotation 1']
    note.place(BOX_1_LEFT, upper_y)
    note.arrow((m1[0], note.box.y0), m1, marker_radius_points(best_1))

    # 2: the isobutanol campaigns' highest-IRR trials, one stem branching
    note = notes['annotation 2']
    note.place(m1[0] + STEM_TO_BOX, lower_y)
    marks = sorted(((marker_xy(fig, left, c), marker_radius_points(c))
                    for c in left.find('best', recorder.labels('isobutanol'))),
                   key=lambda mark: mark[0][0])
    branch_y = cells_top - BRANCH_DROP
    (stem_xy, stem_r), others = marks[0], marks[1:]
    note.arrow((stem_xy[0], note.box.y0), stem_xy, stem_r)
    note.line((stem_xy[0], branch_y), (marks[-1][0][0], branch_y))
    for xy, r in others:
        note.arrow((xy[0], branch_y), xy, r)

    # 3: the start of the TRY-informed incumbent and its highest-IRR trial
    try_informed = recorder.labels('try-informed')[0]
    best_3 = right.find_one('best', try_informed)
    m3 = marker_xy(fig, right, best_3)
    note = notes['annotation 3']
    width_3 = note.size()[0]
    note.place(right_box.x1 - width_3, upper_y)
    note.arrow((m3[0], note.box.y0), m3, marker_radius_points(best_3))
    climb, half_width = try_informed_climb_point(
        fig, right, right.find_one('incumbent', try_informed))
    branch_3 = cells_top + BRACE_Y
    note.line((m3[0], branch_3), (climb[0], branch_3))
    note.arrow((climb[0], branch_3), climb, half_width)

    # 4: a brace over the preload zone (the scouts' trials)
    zone = extent_inches(fig, right.find_one('zone'))
    note = notes['annotation 4']
    note.place(right_box.x0 + 0.02, lower_y)
    note.brace(zone.x0, zone.x1, cells_top + BRACE_Y, BRACE_DROP,
               note.box.y0 - BRACE_GAP)

    check_layout(notes, m1, marks, m3, climb, zone, branch_y)
    return {name: note.artists for name, note in notes.items()}


def check_layout(notes, m1, marks, m3, climb, zone, branch_y):
    """Fail loudly if new data would make the annotations collide."""
    box = {name: note.box for name, note in notes.items()}
    gap = 0.04
    problems = []
    if box['annotation 1'].x1 > box['annotation 3'].x0 - gap:
        problems.append('boxes 1 and 3 overlap')
    if box['annotation 2'].x1 > box['annotation 4'].x0 - gap:
        problems.append('boxes 2 and 4 overlap')
    if not box['annotation 1'].x0 < m1[0] < box['annotation 2'].x0:
        problems.append('stem 1 does not run between box 1 and box 2')
    if not box['annotation 2'].x0 < marks[0][0][0] < box['annotation 2'].x1:
        problems.append('stem 2 misses box 2')
    if not box['annotation 4'].x1 + gap < m3[0] < box['annotation 3'].x1:
        problems.append('stem 3 crosses box 4 or misses box 3')
    if not zone.x1 < climb[0] < m3[0]:
        problems.append('branch 3 runs into the brace')
    if max(xy[1] for xy, _ in marks) > branch_y - 0.1:
        problems.append('a marker of annotation 2 sits above its branch')
    if problems:
        raise RuntimeError('annotation layout: ' + '; '.join(problems))


# -----------------------------------------------------------------------------
# storyboard
# -----------------------------------------------------------------------------
def build_storyboard(fig, recorder):
    left, right = recorder.cells          # the uninformed cell, then the
    if right.find('zone') == []:          # TRY-informed cell with its zone
        raise RuntimeError('expected the TRY-informed cell second')
    notes = add_annotations(fig, recorder, left, right)

    baseline = recorder.labels('baseline')[0]
    uninformed = recorder.labels('uninformed')[0]
    process = recorder.labels('isobutanol', 'ethanol')
    try_informed = recorder.labels('try-informed')[0]
    entry = lambda *labels: ka.legend_entries(fig, *labels)  # noqa: E731

    board = ge.Storyboard(fig, ka.bottom_crop(fig, recorder))
    board.group('baseline', left.find('baseline'),
                ka.campaign_key_frame(fig, recorder), entry(baseline))
    board.group('uninformed', left.find('cloud', {uninformed}),
                left.find('incumbent', {uninformed}),
                entry(uninformed, 'Individual trials',
                      'Incumbent (campaign optimizing the shown metric)'))
    board.group('uninformed best', left.find('best', {uninformed}),
                entry('Highest-IRR trial'))
    board.group('process incumbents', left.find('incumbent-other'),
                entry(*process), entry('Incumbent (other campaign)'))
    board.group('process best', left.find('best', set(process)))
    board.group('right cell', right.ax, right.find('baseline'),
                right.find('zone'))
    board.group('scout trials', right.find('preload'))
    board.group('scout marks', right.find('preload-mark'))
    board.group('try-informed', right.find('cloud', {try_informed}),
                right.find('incumbent', {try_informed}), entry(try_informed))
    board.group('try-informed best', right.find('best', {try_informed}))
    for name, artists in notes.items():
        board.group(name, artists)
    ka.group_panel_b(board, fig, recorder)
    ka.require_cells_grouped(board, recorder)

    n_preloaded = recorder.set_of(try_informed)['n_preloaded']
    left_cell = ge.Region(left.ax)
    scouts = ge.Region(right.ax, x=(0, n_preloaded))
    own_trials = ge.Region(right.ax, x=(n_preloaded, right.ax.get_xlim()[1]))
    fade, wipe, pause = ge.fade, ge.wipe, ge.pause
    S, T, B = ka.FADE_S, ka.TRACE_S, ka.BEAT_S

    board.step('baseline',
               'The dotted grey line is the IRR of the baseline strain, '
               'before any optimization',
               fade('baseline', seconds=1.0))
    board.step('profitability-uninformed',
               'The uninformed profitability campaign climbs to a safe '
               'ethanol-only optimum early and stays there',
               wipe('uninformed', region=left_cell, seconds=T),
               fade('uninformed best', seconds=S), pause(B),
               fade('annotation 1', seconds=S))
    board.step('process-level-incumbents',
               'Plotted as IRR, the incumbents of the six campaigns with '
               'process-level TRY objectives all end below the uninformed one',
               wipe('process incumbents', region=left_cell, seconds=T))
    board.step('process-level-best-trials',
               'Yet several of them tried designs far more profitable than '
               'the uninformed optimum',
               fade('process best', seconds=S), pause(B),
               fade('annotation 2', seconds=S))
    board.step('scout-trials',
               'The TRY-informed campaign starts from 1000 trials of the '
               'process-level campaigns, used as scouts',
               fade('right cell', seconds=S),
               wipe('scout trials', region=scouts, seconds=2.0),
               fade('scout marks', seconds=S), pause(B),
               fade('annotation 4', seconds=S))
    board.step('profitability-try-informed',
               'From there it finds a more profitable region co-producing '
               'ethanol and isobutanol and optimizes further',
               wipe('try-informed', region=own_trials, seconds=T),
               fade('try-informed best', seconds=S), pause(B),
               fade('annotation 3', seconds=S))
    ka.panel_b_steps(board, recorder)
    return board


def main(argv=None):
    ka.run(__doc__.split('\n\n')[0], 'irr', STEM, build_storyboard, argv)


if __name__ == '__main__':
    main()
