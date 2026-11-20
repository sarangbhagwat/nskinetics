# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Regenerate Fig. 9 and SI Fig. S12 of Kim et al., Appl. Catal. B: Environ.
Energy 386 (2026) 126403, from the NSKinetics models of this package, and
draw a comparison figure and a verification figure.

Figures
-------
``fig9``
    Fig. 9: normalized concentration profiles with (a) H2O2 and (b) O2 as
    the oxidant, next to the reaction lists. Circles are the experimental
    points and solid lines are the NSKinetics model; the black dotted line
    is the model curve the authors deposited (a polyline through their
    solver's steps), drawn for comparison.
``fig10_table3_o2``
    The experiments of Fig. 10a-c compared with the Table 3 O2 model
    extended with the Arrhenius relation, with the pooled R2 of the model
    against the experimental points in a box in each panel; and the
    Arrhenius plot (Fig. 10d) with the lines that the activation energies
    quoted in the text imply next to the regression lines.
``figS12``
    SI Fig. S12: the same profiles with the predicted (a) H2O2 and (b) OPA
    as dashed lines and in the insets. The insets also show the published
    inset curve, digitized from the SI figure, and in (a) the H2O2 trace
    the authors deposited.

``solver_emulation``
    Verification figure, not a figure of the paper: the three model curves
    the authors deposited (Fig. 9a, Fig. 9b, the H2O2 trace of SI
    Fig. S12a) as open circles at their nodes, the steps of the repeated
    solver run (``solver_emulation``) as crosses on top, and the accurate
    solution of the same model as lines. Below each panel, on a
    logarithmic axis, how far the emulated steps and the accurate solution
    are from the deposited nodes.

The layout of Fig. 9 and SI Fig. S12, the series colors, the colors of the
reaction lists and the ticks (inward with minor ticks in Figs. 9 and 10,
outward in SI Fig. S12) follow the published panels; the fonts are Arial
where it is installed. Black dotted deposited curves, digitized overlays
and the measured R2 values are additions, named in the entries above. Run
as a script to write every figure (PNG and PDF) to a folder::

    python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR

Add ``--figures fig9 fig10_table3_o2`` to write only the figures named,
``--formats png`` to choose the formats and ``--dpi 100`` to set the
resolution of the PNG files (default 600).
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.offsetbox import (
    AnchoredOffsetbox, DrawingArea, HPacker, TextArea, VPacker)
from matplotlib.ticker import MultipleLocator

from .pd_c_3hp_oxidation import (
    NORMALIZED, load_h2o2_model, load_o2_model, simulate_at_times)
from .reference import (
    read_fig9a_h2o2, read_fig9b_o2, read_fig10_o2,
    read_figS12a_h2o2_inset_digitized, read_figS12b_opa_inset_digitized)
from .solver_emulation import emulate_authors_run, node_deviations
from .temperature_series import (
    FIG10_TEMPERATURES_C, FIG10D_LINE_ENDS_C, fig10d_points,
    fig10d_regression_lines, fig10d_text_lines, pooled_r2)

__all__ = ('draw_fig9a', 'draw_fig9b', 'draw_figS12a', 'draw_figS12b',
           'plot_fig9', 'plot_fig10_table3_o2', 'plot_figS12',
           'plot_solver_emulation', 'save_figures',)

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
# Fig. 10d: the three plotted constants, named as read_fig10d_arrhenius
# names them, with the label and the color of the published series.
_FIG10D_SERIES = {'k_1': ('3-HP $\\rightarrow$ MA', _COLORS['X_MA']),
                  'k_2': ('MA $\\rightarrow$ AAc', _COLORS['X_AAc']),
                  'k_5': ('MA $\\rightarrow$ Others', _C_FIG10_OTHERS)}

_T_END = 1440.                       # right edge of the published time axis
_T_DENSE = np.linspace(0., _T_END, 2881)

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


def _frame(ax, ylim=(-0.05, 1.0), inward_ticks=True):
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
    ax.set_xlabel('Time (min)', fontsize=15)
    ax.set_ylabel('Normalized\nConcentration', fontsize=15, linespacing=0.9)


def _tag(ax, oxidant, xy=(0.965, 0.955), fontsize=16):
    text, face, edge = _TAGS[oxidant]
    ax.text(*xy, text, transform=ax.transAxes, ha='right', va='top',
            fontsize=fontsize, color=_C_INK, zorder=6,
            bbox=dict(boxstyle='square,pad=0.25', fc=face, ec=edge, lw=1.5))


def _line_handle(color, label, **kwargs):
    return Line2D([], [], color=color, label=label, **kwargs)


def _species_legend(ax, anchor):
    """Two key columns as published (model | experiment), rows aligned.

    The key is built cell by cell and not with ``ax.legend``: a legend
    packs each of its columns on its own, and the rows of the two columns
    then drift apart because a label with the superscript "mod" is taller
    than one with "exp". Here every cell has the same height, whatever the
    font. ``anchor`` is the upper left corner, in axes fractions.
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
                                  color=_COLORS[name], **handle))
            label = TextArea(_LABELS[name] + '$^{%s}$' % superscript,
                             textprops=dict(size=size))
            cells.append(HPacker(children=[box, label], align='baseline',
                                 pad=0., sep=0.5*size))
        columns.append(VPacker(children=cells, align='left', pad=0., sep=0.))
    key = AnchoredOffsetbox(
        'upper left', child=HPacker(children=columns, align='top', pad=0.,
                           sep=1.2*size),
        pad=0., borderpad=0., frameon=False, bbox_to_anchor=anchor,
        bbox_transform=ax.transAxes)
    key.set_gid('species_key')
    ax.add_artist(key)
    return key


def _draw_fig9_panel(ax, km, reference, oxidant, legend_anchor,
                     deposited_anchor):
    """Experiment, model and deposited model curve of one panel of Fig. 9."""
    sim = simulate_at_times(km, _T_DENSE)
    _species(ax, sim, reference['experiment'])
    deposited = reference['model']
    for name in NORMALIZED:
        ax.plot(deposited['time'], deposited[name], color='black', lw=0.9,
                ls=(0, (1, 2)), zorder=4)
    _frame(ax)
    # A separate key for the deposited curve, placed where the panel is
    # empty.
    _species_legend(ax, legend_anchor)
    ax.legend(handles=[_line_handle('black', 'Deposited model curve',
                                    lw=0.9, ls=(0, (1, 2)))],
              loc='upper right', bbox_to_anchor=deposited_anchor,
              frameon=False, fontsize=8.5, handlelength=2.4,
              handletextpad=0.5, borderaxespad=0.)
    _tag(ax, oxidant)


def _temperature_tag(ax, temperature_C, r2_text, y):
    """Temperature box at height ``y`` (axes fraction) and, below it, a box
    with ``r2_text`` (the pooled R2 against the experiment)."""
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


def _panel_letters(axes, x=-0.25):
    for ax, letter in zip(axes, 'abcd'):
        ax.text(x, 1.0, letter, transform=ax.transAxes, fontsize=24,
                va='top', color=_C_INK)


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


def _arrhenius_line(ax, line, ends_C, color, **kwargs):
    inv_T = 1./(np.asarray(ends_C) + 273.15)
    ax.plot(1000.*inv_T, line['slope']*inv_T + line['intercept'],
            color=color, **kwargs)


def _draw_arrhenius_points(ax):
    points = fig10d_points()
    for name, (_, color) in _FIG10D_SERIES.items():
        ax.plot(1000.*points['inv_T'], points['ln_' + name], ls='none',
                marker='o', ms=8, mfc='none', mec=color, mew=1.8, zorder=3)


def _arrhenius_frame(ax):
    ax.set_xlim(2.6, 3.45)
    ax.set_ylim(-13., -3.)
    ax.set_xticks(np.arange(2.6, 3.41, 0.2))
    ax.set_yticks(np.arange(-12, -3, 2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(which='both', direction='in')
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', length=3, width=1.2)
    ax.set_xlabel('1000/T (K$^{-1}$)', fontsize=15)
    ax.set_ylabel('ln k', fontsize=15)


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
        _panel_letters(axes[:, 0], x=-0.26)
    return fig


# Height of the temperature tag in panels a to c of the comparison figure,
# chosen so that the tag and the R2 box below it lie clear of the curves.
_FIG10_TABLE3_TAG_Y = {25: 0.80, 75: 0.80, 100: 0.735}


def plot_fig10_table3_o2():
    """Compare the experiments of Fig. 10 with the temperature-extended
    Table 3 O2 model.

    Panels a to c: the experimental points of Fig. 10a-c with the Table 3
    O2 model at 25, 75 and 100 C (solid lines; the model of Fig. 9b
    extended with the Arrhenius relation; the model curves of the
    published panels are not this extension). The box gives the pooled R2
    of the model against all deposited experimental rows of that
    temperature. The time axis ends at 1440 min, as in the published
    panels; the samples after it (1740 and 1890 min at 25 C, 1560 min at
    100 C) are not drawn but are counted in R2.

    Panel d: the points of Fig. 10d with the published regression lines
    (dashed) and the lines that the activation energies quoted in the text
    imply (solid; 59.8, 91.6 and 59.5 kJ/mol through the 50 C points), over
    25 to 100 C. The model has these activation energies for R1 and R2, R3
    and R6, so each solid line is the plotted 50 C constant multiplied by
    the model's Arrhenius factor ``k(T)/k_ref``. Whether the solid lines
    are the model's own constants depends on the rate law behind the
    plotted constants, which the paper does not state (see
    ``fig10d_text_lines`` and ``NOTES.md``).

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
            _frame(ax)
            _temperature_tag(
                ax, temperature_C, 'R$^2$ = %.3f' % pooled_r2(
                    experiment, simulate_at_times(table3,
                                                  experiment['time'])),
                y=_FIG10_TABLE3_TAG_Y[temperature_C])
        for y, handles in (
                (0.930, [_line_handle(_FIG10_COLORS[n], _LABELS[n], lw=1.8,
                                      marker='o', ms=5, mfc='none', mew=1.3)
                         for n in NORMALIZED]),
                (0.900, [_line_handle(_C_MUTED, 'Experiment', ls='none',
                                      marker='o', ms=5, mfc='none', mew=1.3),
                         _line_handle(_C_MUTED, 'Table 3 O$_2$ model with '
                                      'Arrhenius temperature dependence',
                                      lw=1.8)])):
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
        for ax, draw in zip(axes, (draw_figS12a, draw_figS12b)):
            draw(ax)
        _panel_letters(axes)
    return fig


def _emulation_cases():
    """Title, model, deposited tables, name of the deposited nodes, compared
    columns and colors of the three panels of the solver-emulation figure."""
    fig9a = read_fig9a_h2o2()
    return [
        ('Fig. 9a (H$_2$O$_2$, 50 \N{DEGREE SIGN}C)',
         load_h2o2_model('fig9a_nodes'), fig9a, 'model', NORMALIZED, _COLORS),
        ('Fig. 9b (O$_2$, 50 \N{DEGREE SIGN}C)',
         load_o2_model(parameter_set='fig9b_nodes'), read_fig9b_o2(),
         'model', NORMALIZED, _COLORS),
        ('SI Fig. S12a, H$_2$O$_2$ (equivalents)',
         load_h2o2_model('figS12a_h2o2_nodes'), fig9a, 'h2o2_trace',
         ('eq_H2O2',), {'eq_H2O2': _C_H2O2})]


_NODE_MARKER = dict(ls='none', marker='o', ms=6, mfc='none', mew=1.2)
_STEP_MARKER = dict(ls='none', marker='+', ms=6, mew=1.0)


def plot_solver_emulation():
    """Draw the deposited model nodes with the emulated solver steps on top.

    A verification figure, not a figure of the paper. One pair of panels
    for each deposited model curve (Fig. 9a, Fig. 9b, the H2O2 trace of SI
    Fig. S12a). Upper panel: the deposited nodes (open circles), the steps
    of the authors' solver run repeated with the shipped model
    (:func:`~.solver_emulation.emulate_authors_run`; black crosses) and the
    accurate solution of the same model (lines). Lower panel, logarithmic:
    at each deposited node, the largest absolute difference over the
    plotted quantities between the node and the emulated step (crosses) and
    between the node and the accurate solution (circles), with a dotted
    line at 1e-5, one unit of the last decimal the deposited files print.
    Differences below 1e-8 (the initial state) are not drawn. The models
    carry the constant sets with which the emulation returns the nodes
    (the sets inferred from the nodes, not Table 3 as printed). The title
    of each pair gives the number of steps and the largest difference of
    the emulated steps from the nodes; it is kept short enough for the
    widest fallback font (DejaVu Sans).

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(_RC_PARAMS):
        fig = plt.figure(figsize=(13.2, 5.2))
        grid = fig.add_gridspec(2, 3, height_ratios=[2.2, 1.3], hspace=0.10,
                                wspace=0.30, bottom=0.245, top=0.86)
        for column, (title, km, tables, nodes_name, columns,
                     colors) in enumerate(_emulation_cases()):
            top = fig.add_subplot(grid[0, column])
            bottom = fig.add_subplot(grid[1, column], sharex=top)
            columns = list(columns)
            nodes = tables[nodes_name]
            t_end = 1.1*tables['experiment']['time'].iloc[-1]
            emulated = emulate_authors_run(km, t_end)
            dense = simulate_at_times(km, np.linspace(0., t_end, 1201),
                                      columns)
            accurate = simulate_at_times(km, nodes['time'], columns)
            n = len(nodes)
            for name in columns:
                color = colors[name]
                top.plot(dense['time'], dense[name], color=color, lw=1.2,
                         zorder=2)
                top.plot(nodes['time'], nodes[name], color=color, zorder=3,
                         **_NODE_MARKER)
                top.plot(emulated['time'], emulated[name], color=_C_INK,
                         zorder=4, **_STEP_MARKER)
            for values, marker in ((accurate, _NODE_MARKER),
                                   (emulated.iloc[:n], _STEP_MARKER)):
                difference = np.abs(
                    values[columns].to_numpy()
                    - nodes[columns].to_numpy()).max(axis=1)
                shown = difference >= 1e-8      # the floor of the axis
                bottom.plot(nodes['time'][shown], difference[shown],
                            color=_C_INK, zorder=3, **marker)
            bottom.axhline(1e-5, color=_C_MUTED, lw=1.0, ls=(0, (1, 2)),
                           zorder=1)
            deviations = node_deviations(emulated, nodes, columns)
            top.set_title('%s\n%d steps; within %.1e of the nodes'
                          % (title, deviations['steps'], deviations['max']),
                          fontsize=11.5, color=_C_INK, linespacing=1.3)
            top.set_xlim(-0.02*t_end, 1.02*t_end)
            top.tick_params(labelbottom=False, length=4)
            if nodes_name == 'h2o2_trace':
                top.set_ylim(-0.75, 15.75)
                top.set_yticks([0, 5, 10, 15])
            else:
                top.set_ylim(-0.05, 1.05)
            bottom.set_yscale('log')
            bottom.set_ylim(1e-8, 1e-1)
            bottom.set_yticks([1e-7, 1e-5, 1e-3, 1e-1])
            bottom.tick_params(length=4)
            bottom.tick_params(which='minor', left=False)
            bottom.set_xlabel('Time (min)', fontsize=12)
            if column == 0:
                top.set_ylabel('Normalized\nConcentration', fontsize=12,
                               linespacing=0.9)
                bottom.set_ylabel('Largest |difference|\nfrom the node',
                                  fontsize=12, linespacing=0.9)
        handles = [_line_handle(_COLORS[name], _LABELS[name], lw=2.2)
                   for name in NORMALIZED]
        handles += [
            _line_handle(_C_H2O2, 'H$_2$O$_2$', lw=2.2),
            _line_handle(_C_INK, 'Shipped model, accurate solution', lw=1.2),
            _line_handle(_C_INK, 'Deposited node (below: accurate solution '
                         'against the node)', **_NODE_MARKER),
            _line_handle(_C_INK, 'Emulated solver step (below: step against '
                         'the node)', **_STEP_MARKER),
            _line_handle(_C_MUTED, 'Last printed decimal (1e-5)', lw=1.0,
                         ls=(0, (1, 2)))]
        fig.legend(handles=handles, loc='lower center', ncol=5,
                   frameon=False, fontsize=10.5, handlelength=2.0,
                   columnspacing=1.6, bbox_to_anchor=(0.5, -0.005))
    return fig


# Name and plot function of every figure of this module.
_FIGURES = {'fig9': plot_fig9, 'fig10_table3_o2': plot_fig10_table3_o2,
            'figS12': plot_figS12, 'solver_emulation': plot_solver_emulation}


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
        description='Regenerate Fig. 9 and SI Fig. S12 of Kim et al. '
                    '(2026) from the NSKinetics models, and draw the '
                    'comparison and verification figures.')
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
