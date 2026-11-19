# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Regenerate the kinetic-model figures of Kim et al., Appl. Catal. B: Environ.
Energy 386 (2026) 126403, from the NSKinetics models of this package.

Figures
-------
``fig9``
    Fig. 9: normalized concentration profiles with (a) H2O2 and (b) O2 as
    the oxidant, next to the reaction lists. Circles are the experimental
    points and solid lines are the NSKinetics model; the black dotted line
    is the model curve the authors deposited (a polyline through their
    solver's steps), drawn for comparison.
``fig10``
    Fig. 10: (a-c) normalized concentration profiles with O2 at 25, 75 and
    100 C: experimental points, the reversible O2 variant that was inferred
    from the deposited model curves (solid lines) and those deposited curves
    (black dotted lines); (d) the Arrhenius plot: the deposited points
    (taken to be rate constants of the simplified reversible model, see
    ``NOTES.md``) and the three dashed lines, recomputed by least squares
    through the points each line was fitted to.
``fig10_table3_o2``
    The experiments of Fig. 10a-c compared with the Table 3 O2 model
    extended with the Arrhenius relation (as the text of the paper
    describes for Fig. 11), with the per-temperature constant sets as
    dashed lines; and the Arrhenius plot with the lines that the activation
    energies quoted in the text imply next to the regression lines.
``figS12``
    SI Fig. S12: the same profiles with the predicted (a) H2O2 and (b) OPA
    as dashed lines and in the insets. The insets also show the published
    inset curve, digitized from the SI figure, and in (a) the H2O2 trace
    the authors deposited.
``fig11``
    Fig. 11: heat maps of (a) MA "selectivity" and (b) MA "yield" over
    temperature and time, computed with the authors' simplified reversible
    model as their heat-map script runs it, each panel with its own
    constant set. Contour lines at 0.1, 0.3 and 0.5, the red dashed line
    through the row maxima, a red star at the largest cell and a purple
    star at the maximum of the 50 C row; the insets are zooms of the same
    matrix. Two things differ from the published figure on purpose: the
    largest cell of each map is annotated with its value, and the
    colour-bar label of panel a is spelled "MA Selectivity" (the published
    panel prints "MA Selectvity").
``fig11_verification``
    The reproduced Fig. 11 maps with the contour lines, red dashed line
    and stars digitized from the published figure drawn on top as fine
    black dots, and the offset of every digitized contour sample from the
    reproduced contour.
``fig11_table3_o2``
    The same two maps computed with the Table 3 O2 model extended with the
    Arrhenius relation, which is what the text of the paper describes for
    Fig. 11, with the published (digitized) contour lines drawn on top.

The layout, the series colors, the colors of the reaction lists and the
ticks (inward with minor ticks in Figs. 9 and 10, outward in SI Fig. S12)
follow the published panels; the fonts are Arial where it is installed.
Black dotted deposited curves, digitized overlays and the measured R2
values are additions, named in the entries above. Run as a script to write
every figure (PNG and PDF) to a folder::

    python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR

Add ``--figures fig10 fig10_table3_o2`` to write only the figures named,
``--formats png`` to choose the formats and ``--dpi 100`` to set the
resolution of the PNG files (default 600).
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.offsetbox import (
    AnchoredOffsetbox, DrawingArea, HPacker, TextArea, VPacker)
from matplotlib.ticker import MultipleLocator

from .heat_maps import (
    FIG11_C_PER_PX, FIG11_LEVELS, FIG11_MIN_PER_PX, contour_offsets,
    fig11_panel_heat_map, heat_map_maximum, heat_map_ridge,
    table3_o2_heat_maps)
from .pd_c_3hp_oxidation import (
    NORMALIZED, load_h2o2_model, load_o2_model, load_o2_reversible_model,
    simulate_at_times)
from .reference import (
    read_fig9a_h2o2, read_fig9b_o2, read_fig10_o2,
    read_fig11a_selectivity_digitized, read_fig11b_yield_digitized,
    read_figS12a_h2o2_inset_digitized, read_figS12b_opa_inset_digitized)
from .temperature_series import (
    FIG10_TEMPERATURES_C, FIG10D_LINE_ENDS_C, fig10d_points,
    fig10d_regression_lines, fig10d_text_lines, pooled_r2)

__all__ = ('draw_fig9a', 'draw_fig9b', 'draw_fig10abc', 'draw_fig10d',
           'draw_figS12a', 'draw_figS12b', 'draw_fig11a', 'draw_fig11b',
           'plot_fig9', 'plot_fig10', 'plot_fig10_table3_o2', 'plot_figS12',
           'plot_fig11', 'plot_fig11_verification', 'plot_fig11_table3_o2',
           'save_figures',)

#%% Style
# Arial is the font of the published figure; the others are fallbacks for
# machines without it. Applied through rc_context only, never globally.
_RC_PARAMS = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans',
                        'DejaVu Sans'],
    'font.size': 12,
    'mathtext.default': 'regular',
    'axes.linewidth': 1.5,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    }

_LABELS = {'X_HP3': '3HP', 'X_MA': 'MA', 'X_AAc': 'AAc', 'X_OXA': 'OXA',
           'X_D': 'Others'}
_COLORS = {'X_HP3': '#4d4d4d', 'X_MA': '#a66bd6', 'X_AAc': '#1f6fdb',
           'X_OXA': '#f26a1b', 'X_D': '#5b9bd5'}
# Fig. 10 draws "Others" in a grayer blue than Fig. 9, in all four panels.
_C_FIG10_OTHERS = '#6f9bc4'
_FIG10_COLORS = dict(_COLORS, X_D=_C_FIG10_OTHERS)
_C_H2O2 = '#f28e8e'
_C_OPA = '#2aa35a'
_C_INK = '#1a1a1a'
_C_MUTED = '#5f5f5f'
# Reaction lists of Fig. 9: "Rxn n" labels and the k over each arrow.
_C_RXN = '#3653a0'
_C_RATE_CONSTANT = '#ee2c32'
# Oxidant tag in the corner of each panel: text, face color, edge color.
_TAGS = {'H2O2': ('H$_2$O$_2$', '#fbb4b4', '#e98686'),
         'O2': ('O$_2$', '#c3c3dc', '#8f8fb8')}

# Temperature tag in the corner of each panel of Fig. 10a-c: face color and
# edge color.
_TEMPERATURE_TAGS = {25: ('#b9d3e6', '#2f8fd0'), 75: ('#ffc9bb', '#ff7a3d'),
                     100: ('#ff8585', '#e5262b')}
# Fig. 10d: the three plotted constants, named as the constants of the
# simplified reversible model they are taken to be, with the label and the
# color of the published series.
_FIG10D_SERIES = {'k_1': ('3-HP $\\rightarrow$ MA', _COLORS['X_MA']),
                  'k_2': ('MA $\\rightarrow$ AAc', _COLORS['X_AAc']),
                  'k_5': ('MA $\\rightarrow$ Others', _C_FIG10_OTHERS)}

_T_END = 1440.                       # right edge of the published time axis
_T_DENSE = np.linspace(0., _T_END, 2881)

# Fig. 11: color map and contour colors of the authors' heat-map script.
_HEAT_MAP_CMAP = LinearSegmentedColormap.from_list('kim2026_heat_map', [
    (0.0, (0.737, 0.675, 0.89)), (0.15, (0.659, 0.788, 1.)),
    (0.3, (0.69, 0.89, 0.776)), (0.5, (1., 0.922, 0.588)),
    (0.7, (1., 0.584, 0.392)), (1.0, (1., 0.353, 0.353))])
_CONTOUR_COLORS = {0.1: (0.30, 0.50, 0.82), 0.3: (0.35, 0.57, 0.40),
                  0.5: (0.82, 0.56, 0.09)}
_C_RIDGE = '#e41a1c'
_C_RED_STAR = '#f0525a'
_C_PURPLE_STAR = '#7b6fc4'
# The two panels: mapped quantity, its index in the tuple returned by
# simulate_heat_maps, and the reader of the lines digitized from the
# published panel.
_FIG11_PANELS = {
    'a': ('MA Selectivity', 2, read_fig11a_selectivity_digitized),
    'b': ('MA Yield', 3, read_fig11b_yield_digitized),
    }
# Where the contour labels sit in the published panels and insets, as
# (font size, [(time in min, temperature in C), ...]) groups; each label
# goes on the nearest line. The last group of each inset holds the labels
# of the steep lines before the red dashed line.
_FIG11_LABELS = {
    'a': [(9, [(350., 76.5), (233., 71.), (156., 68.8), (25., 50.),
               (58., 48.7), (227., 42.)])],
    'b': [(9, [(683., 66.4), (522., 59.3), (359., 54.2), (236., 49.),
               (1467., 54.9), (1298., 45.6), (1082., 36.4), (968., 32.),
               (107., 40.), (282., 37.)])],
    }
_FIG11_INSET_LABELS = {
    'a': [(7, [(430., 53.7), (800., 53.7), (1480., 55.)]),
          (4, [(25., 49.6), (64., 50.4), (141., 50.)])],
    'b': [(7, [(330., 55.4), (480., 49.5), (700., 55.3), (985., 49.5),
               (1440., 55.5)]),
          (7, [(65., 47.8), (131., 48.8), (262., 48.2)])],
    }

# Line widths of the computed maps in the overlay figures: wide, so that
# the fine black dots of the digitized published lines show inside them.
_OVERLAY_LINES = dict(contour_lw=2.6, ridge_lw=2.6)

# Reaction lists of Fig. 9: (reactants, products) of each reaction.
_REACTIONS = {
    'H2O2': (('3-HP + H$_2$O$_2$', 'OPA + 2H$_2$O'),
             ('OPA + H$_2$O$_2$', 'MA + H$_2$O'),
             ('MA', 'AAc + CO$_2$'),
             ('MA + H$_2$O$_2$', 'OXA + 2H$_2$ + CO$_2$'),
             ('3-HP', 'Others'),
             ('MA', 'Others'),
             ('H$_2$O$_2$', 'H$_2$O + 0.5 O$_2$')),
    'O2': (('3-HP + 0.5 O$_2$', 'OPA + H$_2$O'),
           ('OPA + 0.5 O$_2$', 'MA'),
           ('MA', 'AAc + CO$_2$'),
           ('MA + O$_2$', 'OXA + H$_2$ + CO$_2$'),
           ('3-HP', 'Others'),
           ('MA', 'Others')),
    }


#%% Shared drawing helpers
def _species(ax, sim, experiment, ms=6.5, lw=1.8, colors=_COLORS):
    """Model lines and experimental circles of the five measured species."""
    # Points on the frame are drawn whole (clip_on=False), so the samples
    # taken after the end of the published time axis are left out here.
    experiment = experiment[experiment['time'] <= _T_END]
    for name in NORMALIZED:
        ax.plot(sim['time'], sim[name], color=colors[name], lw=lw,
                zorder=2)
        ax.plot(experiment['time'], experiment[name], ls='none', marker='o',
                ms=ms, mfc='none', mec=colors[name], mew=lw*0.8, zorder=3,
                clip_on=False)


def _frame(ax, ylim=(-0.05, 1.0), labelsize=15, inward_ticks=True):
    ax.set_xlim(-30., _T_END)
    ax.set_ylim(*ylim)
    ax.set_xticks([0, 500, 1000])
    ax.set_yticks(np.arange(0., 1.01, 0.2))
    ax.tick_params(which='major', length=5)
    if inward_ticks:
        # Inward ticks with minor ticks half-way, as in the published
        # Figs. 9 and 10 (SI Fig. S12 has outward major ticks only).
        ax.xaxis.set_minor_locator(MultipleLocator(250))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params(which='both', direction='in')
        ax.tick_params(which='minor', length=3, width=1.2)
    ax.set_xlabel('Time (min)', fontsize=labelsize)
    ax.set_ylabel('Normalized\nConcentration', fontsize=labelsize,
                  linespacing=0.9)


def _tag(ax, oxidant, xy=(0.965, 0.955), fontsize=16):
    text, face, edge = _TAGS[oxidant]
    ax.text(*xy, text, transform=ax.transAxes, ha='right', va='top',
            fontsize=fontsize, color=_C_INK, zorder=6,
            bbox=dict(boxstyle='square,pad=0.25', fc=face, ec=edge, lw=1.5))


def _line_handle(color, label, **kwargs):
    return Line2D([], [], color=color, label=label, **kwargs)


def _species_legend(ax, anchor, loc='upper left', colors=_COLORS):
    """Two key columns as published (model | experiment), rows aligned.

    The key is built cell by cell and not with ``ax.legend``: a legend
    packs each of its columns on its own, and the rows of the two columns
    then drift apart because a label with the superscript "mod" is taller
    than one with "exp". Here every cell has the same height, whatever the
    font. ``anchor`` is the corner named by ``loc``, in axes fractions.
    """
    size = 8.5                              # font size in points
    ascent, descent = 1.02*size, 0.24*size  # of a cell, about the baseline
    columns = []
    for superscript, width, handle in (
            ('mod', 2.4*size, dict(lw=1.8)),
            ('exp', 1.0*size, dict(ls='none', marker='o', ms=5, mfc='none',
                                   mew=1.3))):
        # A line across the handle, or one marker in its middle.
        xdata = [0., width] if superscript == 'mod' else [0.5*width]
        cells = []
        for name in NORMALIZED:
            box = DrawingArea(width, ascent + descent, 0., descent)
            box.add_artist(Line2D(xdata, [0.33*size]*len(xdata),
                                  color=colors[name], **handle))
            label = TextArea(_LABELS[name] + '$^{%s}$' % superscript,
                             textprops=dict(size=size))
            cells.append(HPacker(children=[box, label], align='baseline',
                                 pad=0., sep=0.5*size))
        columns.append(VPacker(children=cells, align='left', pad=0., sep=0.))
    key = AnchoredOffsetbox(
        loc, child=HPacker(children=columns, align='top', pad=0.,
                           sep=1.2*size),
        pad=0., borderpad=0., frameon=False, bbox_to_anchor=anchor,
        bbox_transform=ax.transAxes)
    key.set_gid('species_key')
    ax.add_artist(key)
    return key


def _draw_profiles(ax, km, reference, legend_anchor, deposited_anchor,
                   legend_loc='upper left', colors=_COLORS):
    """Experiment, model and deposited model curve of one profile panel."""
    sim = simulate_at_times(km, _T_DENSE)
    _species(ax, sim, reference['experiment'], colors=colors)
    deposited = reference['model']
    for name in NORMALIZED:
        ax.plot(deposited['time'], deposited[name], color='black', lw=0.9,
                ls=(0, (1, 2)), zorder=4)
    _frame(ax)
    # A separate key for the deposited curve, placed where the panel is
    # empty.
    _species_legend(ax, legend_anchor, legend_loc, colors)
    ax.legend(handles=[_line_handle('black', 'Deposited model curve',
                                    lw=0.9, ls=(0, (1, 2)))],
              loc='upper right', bbox_to_anchor=deposited_anchor,
              frameon=False, fontsize=8.5, handlelength=2.4,
              handletextpad=0.5, borderaxespad=0.)


def _draw_fig9_panel(ax, km, reference, oxidant, legend_anchor,
                     deposited_anchor):
    _draw_profiles(ax, km, reference, legend_anchor, deposited_anchor)
    _tag(ax, oxidant)


def _temperature_tag(ax, temperature_C, r2_text, y=0.775):
    """Temperature box and, below it, a box with ``r2_text`` (the pooled R2
    against the experiment)."""
    face, edge = _TEMPERATURE_TAGS[temperature_C]
    ax.text(0.945, y, '%d \N{DEGREE SIGN}C' % temperature_C,
            transform=ax.transAxes, ha='right', va='center', fontsize=17,
            color=_C_INK, zorder=6,
            bbox=dict(boxstyle='square,pad=0.22', fc=face, ec=edge, lw=1.5))
    ax.text(0.945, y - 0.08, r2_text, transform=ax.transAxes,
            ha='right', va='top', fontsize=11, color=_C_INK, zorder=6,
            linespacing=1.25,
            bbox=dict(boxstyle='round,pad=0.3', fc='#e6e6e6', ec='none'))


def _inset(ax, bounds):
    inset = ax.inset_axes(bounds, zorder=5)
    inset.patch.set_alpha(0.75)       # the main panel shows faintly through
    inset.tick_params(labelsize=8.5, length=3.5, pad=2)
    inset.set_xlabel('Time (min)', fontsize=9, labelpad=0.5)
    inset.set_ylabel('Normalized\nConcentration', fontsize=9,
                     linespacing=0.9, labelpad=2)
    return inset


#%% Panels
def draw_fig9a(ax):
    """Draw Fig. 9a (H2O2): experiment, model and deposited curve.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    """
    _draw_fig9_panel(ax, load_h2o2_model(), read_fig9a_h2o2(), 'H2O2',
                     legend_anchor=(0.15, 0.97), deposited_anchor=(0.97, 0.66))


def draw_fig9b(ax):
    """Draw Fig. 9b (O2): experiment, model and deposited curve.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    """
    _draw_fig9_panel(ax, load_o2_model(), read_fig9b_o2(), 'O2',
                     legend_anchor=(0.07, 0.97), deposited_anchor=(0.97, 0.50))


# Where the species key (as published) and the key of the deposited curve
# sit in the panels of Fig. 10a-c. The species key is anchored by its right
# edge, so that a wider fallback font makes it grow away from the frame.
_FIG10_LEGEND_ANCHOR = (0.975, 0.555)
_FIG10_KEY_ANCHOR = (0.975, 0.215)


def draw_fig10abc(ax, temperature_C):
    """Draw Fig. 10a, b or c: O2 at 25, 75 or 100 C.

    Circles are the experimental points, solid lines the reversible O2
    variant with the constant set of that temperature, and the black dotted
    line the model curve the authors deposited (a polyline through their
    solver's steps), from which that variant was inferred. The box below
    the temperature gives the pooled R2 of the solid lines against all
    experimental points of the panel (the published panels print 0.99).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    temperature_C : {25, 75, 100}
        Temperature of the experiment in degrees Celsius.
    """
    reference = read_fig10_o2(temperature_C)
    km = load_o2_reversible_model('%dC' % temperature_C)
    _draw_profiles(ax, km, reference, _FIG10_LEGEND_ANCHOR,
                   _FIG10_KEY_ANCHOR, legend_loc='upper right',
                   colors=_FIG10_COLORS)
    experiment = reference['experiment']
    _temperature_tag(ax, temperature_C, 'R$^2$ = %.3f' % pooled_r2(
        experiment, simulate_at_times(km, experiment['time'])))


def _arrhenius_line(ax, line, ends_C, color, **kwargs):
    inv_T = 1./(np.asarray(ends_C) + 273.15)
    ax.plot(1000.*inv_T, line['slope']*inv_T + line['intercept'],
            color=color, **kwargs)


def _draw_arrhenius_points(ax):
    points = fig10d_points()
    for name, (_, color) in _FIG10D_SERIES.items():
        ax.plot(1000.*points['inv_T'], points['ln_' + name], ls='none',
                marker='o', ms=8, mfc='none', mec=color, mew=1.8, zorder=3)


def _arrhenius_frame(ax, labelsize=15):
    ax.set_xlim(2.6, 3.45)
    ax.set_ylim(-13., -3.)
    ax.set_xticks(np.arange(2.6, 3.41, 0.2))
    ax.set_yticks(np.arange(-12, -3, 2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(which='both', direction='in')
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', length=3, width=1.2)
    ax.set_xlabel('1000/T (K$^{-1}$)', fontsize=labelsize)
    ax.set_ylabel('ln k', fontsize=labelsize)


def draw_fig10d(ax):
    """Draw Fig. 10d: the Arrhenius plot.

    The circles are the deposited points: ``ln k`` (k in 1/min) of three
    rate constants at 25, 50, 75 and 100 C, taken to be 3-HP to MA
    (``k_1``), MA to AAc (``k_2``) and MA to others (``k_5``) of the
    simplified reversible model (see ``NOTES.md``; the paper does not name
    the model). Each dashed line is the least-squares line through three
    of the four points, recomputed here and drawn over the range of the
    published line: ``k_1`` through 25, 75 and 100 C, ``k_2`` and ``k_5``
    through 50, 75 and 100 C.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    """
    lines = fig10d_regression_lines()
    for name, (label, color) in _FIG10D_SERIES.items():
        _arrhenius_line(ax, lines[name], FIG10D_LINE_ENDS_C[name], color,
                        lw=1.4, ls=(0, (3, 2.5)), alpha=0.55, zorder=2)
    _draw_arrhenius_points(ax)
    _arrhenius_frame(ax)
    ax.legend(handles=[
        _line_handle(color, label, lw=1.4, ls=(0, (3, 2.5)), marker='o',
                     ms=8, mfc='none', mew=1.8)
        for label, color in _FIG10D_SERIES.values()],
        loc='upper right', frameon=False, fontsize=12, handlelength=3.2,
        handletextpad=0.5, labelspacing=0.25, borderaxespad=0.6)


def draw_figS12a(ax):
    """Draw SI Fig. S12a: H2O2 model with the predicted H2O2.

    The predicted H2O2 (equivalents of the initial 3-HP) is the dashed line
    of the main axes, where only its tail below 1 is in view, and of the
    0-15 inset. The inset adds the published curve (digitized) and the
    trace the authors deposited; both lie up to 0.38 equivalents below the
    model line and are not reproduced by the Table 3 constants.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    """
    reference = read_fig9a_h2o2()
    sim = simulate_at_times(load_h2o2_model(), _T_DENSE,
                            list(NORMALIZED) + ['eq_H2O2'])
    _species(ax, sim, reference['experiment'])
    ax.plot(sim['time'], sim['eq_H2O2'], color=_C_H2O2, lw=1.8, ls='--',
            zorder=2)
    _frame(ax, inward_ticks=False)

    inset = _inset(ax, [0.50, 0.52, 0.455, 0.44])
    _species(inset, sim, reference['experiment'], ms=2.5, lw=1.0)
    published = read_figS12a_h2o2_inset_digitized()
    trace = reference['h2o2_trace']
    inset.plot(published['time'], published['eq_H2O2'], color=_C_INK, lw=0.8,
               zorder=4)
    inset.plot(trace['time'], trace['eq_H2O2'], ls='none', marker='s',
               ms=2.8, mfc='none', mec=_C_INK, mew=0.6, zorder=4)
    inset.plot(sim['time'], sim['eq_H2O2'], color=_C_H2O2, lw=1.6, ls='--',
               zorder=5)
    inset.set_xlim(-30., _T_END)
    inset.set_ylim(-0.3, 15.)
    inset.set_xticks([0, 500, 1000])
    inset.set_yticks([0, 5, 10, 15])
    _tag(inset, 'H2O2', xy=(0.95, 0.42), fontsize=12)
    inset.legend(handles=[
        _line_handle(_C_H2O2, 'H$_2$O$_2$$^{predicted}$', lw=1.6, ls='--'),
        _line_handle(_C_INK, 'Published (digitized)', lw=0.8),
        _line_handle(_C_INK, 'Deposited trace', ls='none', marker='s',
                     ms=2.8, mfc='none', mew=0.6)],
        loc='upper right', frameon=False, fontsize=7.5, handlelength=2.2,
        labelspacing=0.2, borderaxespad=0.3)


def draw_figS12b(ax):
    """Draw SI Fig. S12b: O2 model with the predicted OPA.

    The predicted OPA liquid mole fraction is the dashed line of the main
    axes (indistinguishable from zero there) and of the inset, whose axis
    is in units of 1e-5. The inset adds the published curve, a polyline
    whose vertices were digitized from the SI figure.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    """
    reference = read_fig9b_o2()
    sim = simulate_at_times(load_o2_model(), _T_DENSE,
                            list(NORMALIZED) + ['x_OPA'])
    _species(ax, sim, reference['experiment'])
    ax.plot(sim['time'], sim['x_OPA'], color=_C_OPA, lw=1.8, ls='--',
            zorder=2)
    _frame(ax, ylim=(-0.05, 1.05), inward_ticks=False)
    handles = [_line_handle(_COLORS[n], _LABELS[n].replace('Others', 'D'),
                            lw=1.8, marker='o', ms=5.5, mfc='white', mew=1.4)
               for n in NORMALIZED]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.10, 0.99),
              frameon=False, fontsize=10, handlelength=2.6, labelspacing=0.25)

    inset = _inset(ax, [0.50, 0.52, 0.455, 0.42])
    published = read_figS12b_opa_inset_digitized()
    inset.plot(published['time'], published['x_OPA']*1e5, color=_C_INK,
               lw=0.8, marker='D', ms=2.8, mfc='none', mew=0.6, zorder=4)
    inset.plot(sim['time'], sim['x_OPA']*1e5, color=_C_OPA, lw=1.6, ls='--',
               zorder=5)
    inset.set_xlim(-50., 600.)
    inset.set_ylim(-0.3, 2.0)
    inset.set_xticks([0, 200, 400, 600])
    inset.set_yticks([0., 0.5, 1., 1.5, 2.])
    inset.text(0.0, 1.01, r'$\times10^{-5}$', transform=inset.transAxes,
               fontsize=8.5, color=_C_INK, va='bottom')
    _tag(inset, 'O2', xy=(0.95, 0.50), fontsize=12)
    inset.legend(handles=[
        _line_handle(_C_OPA, 'OPA$^{predicted}$', lw=1.6, ls='--'),
        _line_handle(_C_INK, 'Published (digitized)', lw=0.8, marker='D',
                     ms=2.8, mfc='none', mew=0.6)],
        loc='upper right', frameon=False, fontsize=7.5, handlelength=2.6,
        labelspacing=0.2, borderaxespad=0.3)


def _draw_reactions(ax, oxidant):
    """Reaction list drawn to the right of a Fig. 9 panel.

    Each reaction is drawn in pieces that share a baseline: reactants,
    arrow and products follow one another, each placed at the right end of
    the piece before it, and the rate constant sits over the arrow.
    """
    ax.axis('off')
    reactions = _REACTIONS[oxidant]
    size = 12.5
    for i, (reactants, products) in enumerate(reactions, start=1):
        y = 1. - (i - 0.5)/len(reactions) - 0.025
        ax.text(0.0, y, 'Rxn %d' % i, transform=ax.transAxes, va='baseline',
                fontsize=size, color=_C_RXN)
        left = ax.text(0.19, y, reactants, transform=ax.transAxes,
                       va='baseline', fontsize=size, color=_C_INK)
        arrow = ax.annotate(
            r'$\longrightarrow$', (1., y), xycoords=(left, 'axes fraction'),
            xytext=(4, 0), textcoords='offset points', va='baseline',
            fontsize=size, color=_C_INK)
        ax.annotate(
            'k$_%d$' % i, (0.45, y), xycoords=(arrow, 'axes fraction'),
            xytext=(0, 0.62*size), textcoords='offset points', ha='center',
            va='baseline', fontsize=0.8*size, color=_C_RATE_CONSTANT)
        ax.annotate(
            products, (1., y), xycoords=(arrow, 'axes fraction'),
            xytext=(4, 0), textcoords='offset points', va='baseline',
            fontsize=size, color=_C_INK)


#%% Heat-map panels (Fig. 11)
def _star(ax, t, T, color, ms=13, **kwargs):
    ax.plot([t], [T], ls='none', marker='*', ms=ms, mfc=color, mec=color,
            mew=0.6, alpha=0.9, zorder=6, **kwargs)


def _heat_map(ax, T_C, times, matrix, contour_lw=1.0, labels=(),
              ridge_lw=1.6):
    """Flat 1 C cells, contour lines, red dashed line and the two stars."""
    mesh = ax.pcolormesh(times, T_C, matrix, shading='nearest',
                         cmap=_HEAT_MAP_CMAP, vmin=0., vmax=1.,
                         rasterized=True)
    contours = ax.contour(times, T_C, matrix, levels=list(FIG11_LEVELS),
                          colors=[_CONTOUR_COLORS[level]
                                  for level in FIG11_LEVELS],
                          linewidths=contour_lw)
    for label_size, positions in labels:
        ax.clabel(contours, fmt='%g', fontsize=label_size, manual=positions)
    ax.plot(heat_map_ridge(times, matrix), T_C, color=_C_RIDGE, ls='--',
            lw=ridge_lw, zorder=4)
    _, T_max, t_max = heat_map_maximum(T_C, times, matrix)
    _star(ax, t_max, T_max, _C_RED_STAR)
    row_50 = matrix[np.argmin(np.abs(T_C - 50.))]
    _star(ax, times[np.argmax(row_50)], 50., _C_PURPLE_STAR)
    return mesh


def _heat_map_frame(ax, labelsize=15):
    ax.set_xlim(7.2, 1682.8)            # half a cell beyond the grid
    ax.set_ylim(29.5, 100.5)
    ax.set_xticks(np.arange(200, 1601, 200))
    ax.set_yticks(np.arange(30, 101, 10))
    ax.tick_params(length=4, direction='in')
    ax.set_xlabel('Time (min)', fontsize=labelsize)
    ax.set_ylabel('Temperature (\N{DEGREE SIGN}C)', fontsize=labelsize)


def _colorbar(ax, mesh, label, labelsize=15):
    # Drawn beside the axes without shrinking them.
    bar = ax.figure.colorbar(mesh, cax=ax.inset_axes([1.025, 0., 0.05, 1.]),
                             ticks=np.arange(0., 1.01, 0.1), format='%g')
    bar.ax.tick_params(length=3)
    bar.set_label(label, fontsize=labelsize, rotation=270, labelpad=20)
    return bar


def _annotate_maximum(ax, T_C, times, matrix, xytext, prefix='Max'):
    value, T_max, t_max = heat_map_maximum(T_C, times, matrix)
    ax.annotate(
        '%s %.3f\n%g \N{DEGREE SIGN}C, %.0f min' % (prefix, value, T_max,
                                                    t_max),
        (t_max, T_max), xytext=xytext, fontsize=9, color=_C_INK, zorder=7,
        ha='left', va='center', linespacing=1.1,
        arrowprops=dict(arrowstyle='-', color=_C_INK, lw=0.7, shrinkA=1,
                        shrinkB=6))


def _draw_fig11_panel(ax, panel, annotation_xytext):
    label = _FIG11_PANELS[panel][0]
    T_C, times, matrix = fig11_panel_heat_map(panel)
    mesh = _heat_map(ax, T_C, times, matrix, labels=_FIG11_LABELS[panel])
    _heat_map_frame(ax)
    _annotate_maximum(ax, T_C, times, matrix, annotation_xytext)
    _colorbar(ax, mesh, label)

    # Inset as published: the 47 to 57 C rows over the whole time axis.
    # The limits are set first: the contour labels are rotated and the
    # lines opened for them as the axes are scaled when they are placed.
    inset = ax.inset_axes([0.535, 0.525, 0.435, 0.435], zorder=5)
    inset.set_xlim(7.2, 1682.8)
    inset.set_ylim(46.5, 57.5)
    _heat_map(inset, T_C, times, matrix, contour_lw=0.9,
              labels=_FIG11_INSET_LABELS[panel], ridge_lw=1.3)
    inset.set_xticks([200, 600, 1000, 1400])
    inset.set_yticks([48, 50, 52, 54, 56])
    # Inward ticks on all four sides, every 200 min and every 1 C, as
    # published.
    inset.xaxis.set_minor_locator(MultipleLocator(200))
    inset.yaxis.set_minor_locator(MultipleLocator(1))
    inset.tick_params(which='both', direction='in', top=True, right=True,
                      length=3)
    inset.tick_params(labelsize=8.5, pad=2)
    inset.set_xlabel('Time (min)', fontsize=9, labelpad=0.5)
    inset.set_ylabel('Temperature (\N{DEGREE SIGN}C)', fontsize=9, labelpad=1)


def draw_fig11a(ax):
    """Draw Fig. 11a: heat map of MA "selectivity" (``X_MA``).

    Computed with the simplified reversible model and its ``'fig11a'``
    constant set, as the authors' heat-map script runs it. The red star is
    the largest cell of the map and the purple star the maximum of its
    50 C row. The inset is a zoom of the same matrix; in the published
    figure the inset lines sit 7 to 10 min later than those of the main
    panel, which a zoom cannot reproduce.

    Two things are not as in the published panel: the largest cell is
    annotated with its value, temperature and plotted time, and the
    colour-bar label reads "MA Selectivity" (published: "MA Selectvity").

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on. The colour bar is drawn beside it.
    """
    _draw_fig11_panel(ax, 'a', annotation_xytext=(130., 96.5))


def draw_fig11b(ax):
    """Draw Fig. 11b: heat map of MA "yield" (``(1 - X_HP3)*X_MA``).

    Computed with the simplified reversible model and its ``'fig11b'``
    constant set; see :func:`draw_fig11a`. In the published figure the
    inset lines are up to 15 min off those of the main panel (41 min for
    the 0.3 contour on the late side). The largest cell is annotated,
    which the published panel does not do; the stars sit on the computed
    maxima, whereas the published main panel draws them 7 and 13 min
    later (see ``NOTES.md``).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on. The colour bar is drawn beside it.
    """
    _draw_fig11_panel(ax, 'b', annotation_xytext=(240., 88.))


def _draw_digitized(ax, digitized, ms=1.5):
    """Lines and stars digitized from a published Fig. 11 panel, in black.

    The line samples are drawn as fine dots, narrower than the lines of the
    computed map they are compared with, so that both stay visible where
    they coincide.
    """
    contours, ridge = digitized['contours'], digitized['ridge']
    for samples in (contours, ridge):
        ax.plot(samples['time'], samples['T_C'], ls='none', marker='.',
                ms=ms, mew=0., color='black', zorder=5)
    stars = digitized['stars']
    ax.plot(stars['time'], stars['T_C'], ls='none', marker='*', ms=12,
            mfc='none', mec='black', mew=0.8, zorder=7)


def _overlay_legend(ax, model_label):
    """Key of an overlay panel: computed map against the published lines."""
    model_contours = tuple(_line_handle(_CONTOUR_COLORS[level], None, lw=2.6)
                           for level in FIG11_LEVELS)
    model_stars = tuple(_line_handle(color, None, ls='none', marker='*',
                                     ms=10)
                        for color in (_C_RED_STAR, _C_PURPLE_STAR))
    handles = [
        model_contours,
        _line_handle(_C_RIDGE, None, lw=2.6, ls='--'),
        model_stars,
        _line_handle('black', None, lw=1.1, ls=(0, (0.1, 2.4)),
                     dash_capstyle='round'),
        _line_handle('black', None, ls='none', marker='*', ms=9, mfc='none',
                     mew=0.8)]
    labels = [model_label + ': contours 0.1, 0.3, 0.5',
              model_label + ': row maxima',
              model_label + ': maximum, 50 \N{DEGREE SIGN}C maximum',
              'Published: lines (digitized)',
              'Published: stars (digitized)']
    ax.legend(handles, labels, loc='upper right', fontsize=8, frameon=True,
              framealpha=0.85, edgecolor='none', handlelength=3.4,
              labelspacing=0.25, borderaxespad=0.6,
              handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})


def _draw_offsets(ax, T_C, times, matrix, digitized, ylim):
    """Offset of each digitized contour sample from the map's contour."""
    contours = digitized['contours']
    rms = []
    for (level, branch), samples in contours.groupby(['level', 'branch']):
        inside = samples[samples['time'].between(times[0], times[-1])
                         & samples['T_C'].between(T_C[0], T_C[-1])]
        offsets = contour_offsets(T_C, times, matrix, level, inside['time'],
                                  inside['T_C'])
        rms.append(offsets)
        early = branch == 'early'
        ax.plot(inside['T_C'], offsets, ls='none',
                marker='o' if early else 's', ms=2.6,
                mfc='none' if early else _CONTOUR_COLORS[level],
                mec=_CONTOUR_COLORS[level], mew=0.5,
                label='%g, %s' % (level, branch))
    rms = np.sqrt(np.mean(np.concatenate(rms)**2))
    ax.axhline(0., color=_C_INK, lw=0.8)
    ax.set_xlim(29.5, 100.5)
    ax.set_ylim(*ylim)
    ax.set_xticks(np.arange(30, 101, 10))
    ax.tick_params(length=4, direction='in')
    ax.set_xlabel('Temperature of the sample (\N{DEGREE SIGN}C)', fontsize=13)
    ax.set_ylabel('Offset from the\ncomputed contour (px)', fontsize=13,
                  linespacing=0.95)
    ax.text(0.02, 0.05, 'rms %.2f px (1 px = %.2f min or %.3f '
            '\N{DEGREE SIGN}C)' % (rms, FIG11_MIN_PER_PX, FIG11_C_PER_PX),
            transform=ax.transAxes, va='bottom', fontsize=10, color=_C_INK)
    ax.legend(title='Contour level, side of the red dashed line', ncol=3,
              loc='upper right', fontsize=8, title_fontsize=8, frameon=False,
              handletextpad=0.2, columnspacing=1.0, labelspacing=0.25,
              markerscale=1.6)
    return rms


#%% Figures
def plot_fig9():
    """Regenerate Fig. 9 (panel a: H2O2, panel b: O2, with reaction lists).

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(
            2, 2, figsize=(9.6, 7.4),
            gridspec_kw=dict(width_ratios=[1., 1.05], wspace=0.07,
                             hspace=0.32))
        for row, oxidant, draw in ((0, 'H2O2', draw_fig9a),
                                   (1, 'O2', draw_fig9b)):
            draw(axes[row, 0])
            _draw_reactions(axes[row, 1], oxidant)
            axes[row, 0].text(-0.26, 1.0, 'ab'[row], fontsize=24,
                              transform=axes[row, 0].transAxes, va='top',
                              color=_C_INK)
    return fig


def _panel_letters(axes, x=-0.25):
    for ax, letter in zip(axes, 'abcd'):
        ax.text(x, 1.0, letter, transform=ax.transAxes, fontsize=24,
                va='top', color=_C_INK)


def plot_fig10():
    """Regenerate Fig. 10 (a-c: O2 at 25, 75 and 100 C; d: Arrhenius plot).

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(2, 2, figsize=(11.6, 9.0),
                                 gridspec_kw=dict(wspace=0.36, hspace=0.30))
        axes = axes.ravel()
        for ax, temperature_C in zip(axes, FIG10_TEMPERATURES_C):
            draw_fig10abc(ax, temperature_C)
        draw_fig10d(axes[3])
        _panel_letters(axes)
    return fig


# Height of the temperature tag in panels a to c of the comparison figure,
# chosen so that the tag and the R2 box below it lie clear of the curves.
_FIG10_TABLE3_TAG_Y = {25: 0.80, 75: 0.80, 100: 0.735}


def plot_fig10_table3_o2():
    """Compare the experiments of Fig. 10 with the temperature-extended
    Table 3 O2 model.

    Panels a to c: the experimental points of Fig. 10a-c with the Table 3
    O2 model at 25, 75 and 100 C (solid lines; the model of Fig. 9b
    extended with the Arrhenius relation, as the text of the paper
    describes for Fig. 11; the paper does not show this comparison) and
    the curves of the per-temperature constant sets (dashed lines; the
    reversible O2 variant, which reproduces the published curves). The box
    gives the pooled R2 of each against all experimental points of the
    panel: first the Table 3 model, then the per-temperature set.

    Panel d: the points of Fig. 10d with the published regression lines
    (dashed) and the lines that the activation energies quoted in the text
    imply (solid; 59.8, 91.6 and 59.5 kJ/mol through the 50 C points, as in
    the authors' heat-map script), over 25 to 100 C.

    The curves of the Table 3 model cross the whole of panels a to c, so
    the keys are placed above the panels and not inside them.

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(2, 2, figsize=(11.6, 9.0),
                                 gridspec_kw=dict(wspace=0.36, hspace=0.30))
        axes = axes.ravel()
        for ax, temperature_C in zip(axes, FIG10_TEMPERATURES_C):
            experiment = read_fig10_o2(temperature_C)['experiment']
            table3 = load_o2_model(T=temperature_C + 273.15)
            _species(ax, simulate_at_times(table3, _T_DENSE), experiment,
                     colors=_FIG10_COLORS)
            reversible = load_o2_reversible_model('%dC' % temperature_C)
            fit = simulate_at_times(reversible, _T_DENSE)
            for name in NORMALIZED:
                ax.plot(fit['time'], fit[name], color=_FIG10_COLORS[name],
                        lw=1.0, ls=(0, (4, 2)), zorder=2)
            _frame(ax)
            r2 = [pooled_r2(experiment,
                            simulate_at_times(km, experiment['time']))
                  for km in (table3, reversible)]
            _temperature_tag(
                ax, temperature_C,
                'R$^2$ = %.3f (Table 3)\n'
                'R$^2$ = %.3f (reversible variant)' % tuple(r2),
                y=_FIG10_TABLE3_TAG_Y[temperature_C])
        for y, handles in (
                (0.930, [_line_handle(_FIG10_COLORS[n], _LABELS[n], lw=1.8,
                                      marker='o', ms=5, mfc='none', mew=1.3)
                         for n in NORMALIZED]),
                (0.900, [_line_handle(_C_MUTED, 'Experiment', ls='none',
                                      marker='o', ms=5, mfc='none', mew=1.3),
                         _line_handle(_C_MUTED, 'Table 3 O$_2$ model with '
                                      'Arrhenius temperature dependence',
                                      lw=1.8),
                         _line_handle(_C_MUTED, 'Reversible variant, one '
                                      'constant set per temperature',
                                      lw=1.0, ls=(0, (4, 2)))])):
            fig.legend(handles=handles, ncol=len(handles),
                       loc='lower center', bbox_to_anchor=(0.5, y),
                       bbox_transform=fig.transFigure, frameon=False,
                       fontsize=10, handlelength=2.6, handletextpad=0.5,
                       columnspacing=1.6, borderaxespad=0.)

        ax = axes[3]
        regression, text = fig10d_regression_lines(), fig10d_text_lines()
        for name, (_, color) in _FIG10D_SERIES.items():
            _arrhenius_line(ax, regression[name], FIG10D_LINE_ENDS_C[name],
                            color, lw=1.4, ls=(0, (3, 2.5)), alpha=0.55,
                            zorder=2)
            _arrhenius_line(ax, text[name], (100., 25.), color, lw=1.6,
                            zorder=2)
        _draw_arrhenius_points(ax)
        _arrhenius_frame(ax)
        handles = [_line_handle(color, label, ls='none', marker='o', ms=8,
                                mfc='none', mew=1.8)
                   for label, color in _FIG10D_SERIES.values()]
        handles += [_line_handle(_C_MUTED, 'Regression line (published)',
                                 lw=1.4, ls=(0, (3, 2.5)), alpha=0.55),
                    _line_handle(_C_MUTED, 'E$_a$ of the text, from 50 '
                                 '\N{DEGREE SIGN}C', lw=1.6)]
        ax.legend(handles=handles, loc='lower left', frameon=False,
                  fontsize=10, handlelength=3.0, handletextpad=0.5,
                  labelspacing=0.25, borderaxespad=0.5)
        _panel_letters(axes)
    return fig


def plot_figS12():
    """Regenerate SI Fig. S12 (a: predicted H2O2, b: predicted OPA).

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.3),
                                 gridspec_kw=dict(wspace=0.36))
        for ax, letter, draw in zip(axes, 'ab', (draw_figS12a, draw_figS12b)):
            draw(ax)
            ax.text(-0.25, 1.0, letter, transform=ax.transAxes, fontsize=24,
                    va='top', color=_C_INK)
    return fig


def plot_fig11():
    """Regenerate Fig. 11 (a: MA "selectivity", b: MA "yield" heat maps).

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.6),
                                 gridspec_kw=dict(wspace=0.48))
        for ax, letter, draw in zip(axes, 'ab', (draw_fig11a, draw_fig11b)):
            draw(ax)
            ax.text(-0.21, 1.03, letter, transform=ax.transAxes, fontsize=24,
                    va='top', color=_C_INK)
    return fig


def plot_fig11_verification():
    """Compare the reproduced Fig. 11 maps with the published figure.

    Top row: the reproduced maps of :func:`plot_fig11` with the contour
    lines, red dashed line and stars digitized from the published panels
    drawn on top in black. Bottom row: the offset of every digitized
    contour sample from the reproduced contour line, in pixels of the
    published raster, against the temperature of the sample.

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(
            2, 2, figsize=(13.2, 8.0),
            gridspec_kw=dict(wspace=0.48, hspace=0.30,
                             height_ratios=[1.6, 1.]))
        for column, panel in enumerate('ab'):
            label, _, read = _FIG11_PANELS[panel]
            digitized = read()
            T_C, times, matrix = fig11_panel_heat_map(panel)
            top, bottom = axes[0, column], axes[1, column]
            mesh = _heat_map(top, T_C, times, matrix, **_OVERLAY_LINES)
            _draw_digitized(top, digitized)
            _heat_map_frame(top)
            _colorbar(top, mesh, label)
            _overlay_legend(top, 'Reproduced')
            top.text(-0.21, 1.03, panel, transform=top.transAxes,
                     fontsize=24, va='top', color=_C_INK)
            _draw_offsets(bottom, T_C, times, matrix, digitized,
                          ylim=(-1.5, 1.5))
    return fig


def plot_fig11_table3_o2():
    """Draw the Fig. 11 maps as the Table 3 O2 model predicts them.

    The maps are computed with :func:`~.heat_maps.table3_o2_heat_maps`: the
    six-reaction O2 model of Table 3 and Fig. 9b with its Arrhenius
    temperature dependence, started from the feed at time 0. This is the
    calculation the text of the paper describes; the authors' heat-map
    script computes the maps with the simplified reversible model instead.
    The published contour lines, red dashed line and stars (digitized) are
    drawn on top in black for comparison, and the maximum of each computed
    map is annotated.

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.6),
                                 gridspec_kw=dict(wspace=0.48))
        maps = table3_o2_heat_maps()
        T_C, times = maps[:2]
        for ax, panel, xytext in zip(axes, 'ab', ((420., 80.), (560., 76.))):
            label, index, read = _FIG11_PANELS[panel]
            matrix = maps[index]
            mesh = _heat_map(ax, T_C, times, matrix, **_OVERLAY_LINES)
            _draw_digitized(ax, read())
            _heat_map_frame(ax)
            _annotate_maximum(ax, T_C, times, matrix, xytext)
            _colorbar(ax, mesh, label)
            _overlay_legend(ax, 'Table 3 O$_2$ model')
            ax.text(-0.21, 1.03, panel, transform=ax.transAxes, fontsize=24,
                    va='top', color=_C_INK)
    return fig


# Name and plot function of every figure of this module.
_FIGURES = {'fig9': plot_fig9, 'fig10': plot_fig10,
            'fig10_table3_o2': plot_fig10_table3_o2, 'figS12': plot_figS12,
            'fig11': plot_fig11,
            'fig11_verification': plot_fig11_verification,
            'fig11_table3_o2': plot_fig11_table3_o2}


def save_figures(output_dir, formats=('png', 'pdf'), figures=None, dpi=600):
    """Write the figures of this module to ``output_dir``.

    Parameters
    ----------
    output_dir : str
        Folder to write to; created if missing.
    formats : tuple of str, optional
        File formats, by default PNG and vector PDF.
    figures : sequence of str, optional
        Names of the figures to write (see the module docstring), by
        default all of them.
    dpi : float, optional
        Resolution of the raster formats, by default 600. The images of
        the documentation page are written with ``formats=('png',)`` and
        ``dpi=100``.

    Returns
    -------
    list of str
        Paths of the files written, named
        ``pd_c_3hp_oxidation_<figure>.<format>``.
    """
    if figures is None:
        figures = list(_FIGURES)
    unknown = [name for name in figures if name not in _FIGURES]
    if unknown:
        raise ValueError(f'unknown figures {unknown}; '
                         f'choose from {list(_FIGURES)}')
    os.makedirs(output_dir, exist_ok=True)
    paths = []
    for name in figures:
        with plt.rc_context(_RC_PARAMS):
            fig = _FIGURES[name]()
            for fmt in formats:
                path = os.path.join(output_dir,
                                    f'pd_c_3hp_oxidation_{name}.{fmt}')
                fig.savefig(path, dpi=dpi, bbox_inches='tight',
                            facecolor='white')
                paths.append(path)
            plt.close(fig)
    return paths


def main():
    parser = argparse.ArgumentParser(
        prog='python -m nskinetics.models.pd_c_3hp_oxidation.figures',
        description='Regenerate Fig. 9, Fig. 10, SI Fig. S12 and Fig. 11 '
                    'of Kim et al. (2026) from the NSKinetics models.')
    parser.add_argument('output_dir', help='folder to write the figures to')
    parser.add_argument('--formats', nargs='+', default=['png', 'pdf'],
                        help='file formats (default: png pdf)')
    parser.add_argument('--figures', nargs='+', default=None,
                        choices=list(_FIGURES), metavar='NAME',
                        help='figures to write (default: all): '
                             + ', '.join(_FIGURES))
    parser.add_argument('--dpi', type=float, default=600.,
                        help='resolution of the PNG files (default: 600)')
    args = parser.parse_args()
    if os.environ.get('MPLBACKEND') is None:
        plt.switch_backend('Agg')    # no display needed to write files
    for path in save_figures(args.output_dir, tuple(args.formats),
                             args.figures, args.dpi):
        print(f'Saved {path}')


if __name__ == '__main__':
    main()
