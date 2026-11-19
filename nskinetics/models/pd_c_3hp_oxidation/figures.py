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
``figS12``
    SI Fig. S12: the same profiles with the predicted (a) H2O2 and (b) OPA
    as dashed lines and in the insets. The insets also show the published
    inset curve, digitized from the SI figure, and in (a) the H2O2 trace
    the authors deposited.

Colors and layout follow the published panels. Run as a script to write
every figure (PNG and PDF) to a folder::

    python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .pd_c_3hp_oxidation import (
    NORMALIZED, load_h2o2_model, load_o2_model, simulate_at_times)
from .reference import (
    read_fig9a_h2o2, read_fig9b_o2, read_figS12a_h2o2_inset_digitized,
    read_figS12b_opa_inset_digitized)

__all__ = ('draw_fig9a', 'draw_fig9b', 'draw_figS12a', 'draw_figS12b',
           'plot_fig9', 'plot_figS12', 'save_figures',)

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

LABELS = {'X_HP3': '3HP', 'X_MA': 'MA', 'X_AAc': 'AAc', 'X_OXA': 'OXA',
          'X_D': 'Others'}
COLORS = {'X_HP3': '#4d4d4d', 'X_MA': '#a66bd6', 'X_AAc': '#1f6fdb',
          'X_OXA': '#f26a1b', 'X_D': '#5b9bd5'}
C_H2O2 = '#f28e8e'
C_OPA = '#2aa35a'
C_INK = '#1a1a1a'
C_MUTED = '#5f5f5f'
# Oxidant tag in the corner of each panel: text, face color, edge color.
_TAGS = {'H2O2': ('H$_2$O$_2$', '#fbb4b4', '#e98686'),
         'O2': ('O$_2$', '#c3c3dc', '#8f8fb8')}

T_END = 1440.                       # right edge of the published time axis
_T_DENSE = np.linspace(0., T_END, 2881)

REACTIONS = {
    'H2O2': ('3-HP + H$_2$O$_2$ {} OPA + 2H$_2$O',
             'OPA + H$_2$O$_2$ {} MA + H$_2$O',
             'MA {} AAc + CO$_2$',
             'MA + H$_2$O$_2$ {} OXA + 2H$_2$ + CO$_2$',
             '3-HP {} Others',
             'MA {} Others',
             'H$_2$O$_2$ {} H$_2$O + 0.5 O$_2$'),
    'O2': ('3-HP + 0.5 O$_2$ {} OPA + H$_2$O',
           'OPA + 0.5 O$_2$ {} MA',
           'MA {} AAc + CO$_2$',
           'MA + O$_2$ {} OXA + H$_2$ + CO$_2$',
           '3-HP {} Others',
           'MA {} Others'),
    }


#%% Shared drawing helpers
def _species(ax, sim, experiment, ms=6.5, lw=1.8):
    """Model lines and experimental circles of the five measured species."""
    # Points on the frame are drawn whole (clip_on=False), so the samples
    # taken after the end of the published time axis are left out here.
    experiment = experiment[experiment['time'] <= T_END]
    for name in NORMALIZED:
        ax.plot(sim['time'], sim[name], color=COLORS[name], lw=lw,
                zorder=2)
        ax.plot(experiment['time'], experiment[name], ls='none', marker='o',
                ms=ms, mfc='none', mec=COLORS[name], mew=lw*0.8, zorder=3,
                clip_on=False)


def _frame(ax, ylim=(-0.05, 1.0), labelsize=15):
    ax.set_xlim(-30., T_END)
    ax.set_ylim(*ylim)
    ax.set_xticks([0, 500, 1000])
    ax.set_yticks(np.arange(0., 1.01, 0.2))
    ax.tick_params(length=5)
    ax.set_xlabel('Time (min)', fontsize=labelsize)
    ax.set_ylabel('Normalized\nConcentration', fontsize=labelsize,
                  linespacing=0.9)


def _tag(ax, oxidant, xy=(0.965, 0.955), fontsize=16):
    text, face, edge = _TAGS[oxidant]
    ax.text(*xy, text, transform=ax.transAxes, ha='right', va='top',
            fontsize=fontsize, color=C_INK, zorder=6,
            bbox=dict(boxstyle='square,pad=0.25', fc=face, ec=edge, lw=1.5))


def _line_handle(color, label, **kwargs):
    return Line2D([], [], color=color, label=label, **kwargs)


def _draw_fig9_panel(ax, km, reference, oxidant, legend_anchor,
                     deposited_anchor):
    sim = simulate_at_times(km, _T_DENSE)
    _species(ax, sim, reference['experiment'])
    deposited = reference['model']
    for name in NORMALIZED:
        ax.plot(deposited['time'], deposited[name], color='black', lw=0.9,
                ls=(0, (1, 2)), zorder=4)
    _frame(ax)
    _tag(ax, oxidant)
    # Two legend columns as published (model | experiment), and a separate
    # key for the deposited curve, placed where the panel is empty.
    handles = [_line_handle(COLORS[n], LABELS[n] + '$^{mod}$', lw=1.8)
               for n in NORMALIZED]
    handles += [_line_handle(COLORS[n], LABELS[n] + '$^{exp}$', ls='none',
                             marker='o', ms=5, mfc='none', mew=1.3)
                for n in NORMALIZED]
    species_legend = ax.legend(
        handles=handles, ncol=2, loc='upper left',
        bbox_to_anchor=legend_anchor, frameon=False, fontsize=8.5,
        handlelength=2.4, handletextpad=0.5, columnspacing=1.2,
        labelspacing=0.1, borderaxespad=0.)
    ax.add_artist(species_legend)
    ax.legend(handles=[_line_handle('black', 'Deposited model curve',
                                    lw=0.9, ls=(0, (1, 2)))],
              loc='upper right', bbox_to_anchor=deposited_anchor,
              frameon=False, fontsize=8.5, handlelength=2.4,
              handletextpad=0.5, borderaxespad=0.)


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
    """Draw Fig. 9a (H2O2) on ``ax``: experiment, model, deposited curve."""
    _draw_fig9_panel(ax, load_h2o2_model(), read_fig9a_h2o2(), 'H2O2',
                     legend_anchor=(0.15, 0.97), deposited_anchor=(0.97, 0.66))


def draw_fig9b(ax):
    """Draw Fig. 9b (O2) on ``ax``: experiment, model, deposited curve."""
    _draw_fig9_panel(ax, load_o2_model(), read_fig9b_o2(), 'O2',
                     legend_anchor=(0.07, 0.97), deposited_anchor=(0.97, 0.50))


def draw_figS12a(ax):
    """Draw SI Fig. S12a on ``ax``: H2O2 model with the predicted H2O2.

    The predicted H2O2 (equivalents of the initial 3-HP) is the dashed line
    of the main axes, where only its tail below 1 is in view, and of the
    0-15 inset. The inset adds the published curve (digitized) and the
    trace the authors deposited; both come from a parameter set other than
    Table 3 and lie up to 0.38 equivalents below the model line.
    """
    reference = read_fig9a_h2o2()
    sim = simulate_at_times(load_h2o2_model(), _T_DENSE,
                            list(NORMALIZED) + ['eq_H2O2'])
    _species(ax, sim, reference['experiment'])
    ax.plot(sim['time'], sim['eq_H2O2'], color=C_H2O2, lw=1.8, ls='--',
            zorder=2)
    _frame(ax)

    inset = _inset(ax, [0.50, 0.52, 0.455, 0.44])
    _species(inset, sim, reference['experiment'], ms=2.5, lw=1.0)
    published = read_figS12a_h2o2_inset_digitized()
    trace = reference['h2o2_trace']
    inset.plot(published['time'], published['eq_H2O2'], color=C_INK, lw=0.8,
               zorder=4)
    inset.plot(trace['time'], trace['eq_H2O2'], ls='none', marker='s',
               ms=2.8, mfc='none', mec=C_INK, mew=0.6, zorder=4)
    inset.plot(sim['time'], sim['eq_H2O2'], color=C_H2O2, lw=1.6, ls='--',
               zorder=5)
    inset.set_xlim(-30., T_END)
    inset.set_ylim(-0.3, 15.)
    inset.set_xticks([0, 500, 1000])
    inset.set_yticks([0, 5, 10, 15])
    _tag(inset, 'H2O2', xy=(0.95, 0.42), fontsize=12)
    inset.legend(handles=[
        _line_handle(C_H2O2, 'H$_2$O$_2$$^{predicted}$', lw=1.6, ls='--'),
        _line_handle(C_INK, 'Published (digitized)', lw=0.8),
        _line_handle(C_INK, 'Deposited trace', ls='none', marker='s',
                     ms=2.8, mfc='none', mew=0.6)],
        loc='upper right', frameon=False, fontsize=7.5, handlelength=2.2,
        labelspacing=0.2, borderaxespad=0.3)


def draw_figS12b(ax):
    """Draw SI Fig. S12b on ``ax``: O2 model with the predicted OPA.

    The predicted OPA liquid mole fraction is the dashed line of the main
    axes (indistinguishable from zero there) and of the inset, whose axis
    is in units of 1e-5. The inset adds the published curve, a polyline
    whose vertices were digitized from the SI figure.
    """
    reference = read_fig9b_o2()
    sim = simulate_at_times(load_o2_model(), _T_DENSE,
                            list(NORMALIZED) + ['x_OPA'])
    _species(ax, sim, reference['experiment'])
    ax.plot(sim['time'], sim['x_OPA'], color=C_OPA, lw=1.8, ls='--',
            zorder=2)
    _frame(ax, ylim=(-0.05, 1.05))
    handles = [_line_handle(COLORS[n], LABELS[n].replace('Others', 'D'),
                            lw=1.8, marker='o', ms=5.5, mfc='white', mew=1.4)
               for n in NORMALIZED]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.10, 0.99),
              frameon=False, fontsize=10, handlelength=2.6, labelspacing=0.25)

    inset = _inset(ax, [0.50, 0.52, 0.455, 0.42])
    published = read_figS12b_opa_inset_digitized()
    inset.plot(published['time'], published['x_OPA']*1e5, color=C_INK,
               lw=0.8, marker='D', ms=2.8, mfc='none', mew=0.6, zorder=4)
    inset.plot(sim['time'], sim['x_OPA']*1e5, color=C_OPA, lw=1.6, ls='--',
               zorder=5)
    inset.set_xlim(-50., 600.)
    inset.set_ylim(-0.3, 2.0)
    inset.set_xticks([0, 200, 400, 600])
    inset.set_yticks([0., 0.5, 1., 1.5, 2.])
    inset.text(0.0, 1.01, r'$\times10^{-5}$', transform=inset.transAxes,
               fontsize=8.5, color=C_INK, va='bottom')
    _tag(inset, 'O2', xy=(0.95, 0.50), fontsize=12)
    inset.legend(handles=[
        _line_handle(C_OPA, 'OPA$^{predicted}$', lw=1.6, ls='--'),
        _line_handle(C_INK, 'Published (digitized)', lw=0.8, marker='D',
                     ms=2.8, mfc='none', mew=0.6)],
        loc='upper right', frameon=False, fontsize=7.5, handlelength=2.6,
        labelspacing=0.2, borderaxespad=0.3)


def _draw_reactions(ax, oxidant):
    """Reaction list drawn to the right of a Fig. 9 panel."""
    ax.axis('off')
    reactions = REACTIONS[oxidant]
    for i, reaction in enumerate(reactions, start=1):
        y = 1. - (i - 0.5)/len(reactions)
        arrow = r'$\overset{k_%d}{\longrightarrow}$' % i
        ax.text(0.0, y, 'Rxn %d' % i, transform=ax.transAxes, va='center',
                fontsize=12.5, color=C_MUTED)
        ax.text(0.19, y, reaction.format(arrow), transform=ax.transAxes,
                va='center', fontsize=12.5, color=C_INK)


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
                              color=C_INK)
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
                    va='top', color=C_INK)
    return fig


def save_figures(output_dir, formats=('png', 'pdf')):
    """Write every figure of this module to ``output_dir``.

    Parameters
    ----------
    output_dir : str
        Folder to write to; created if missing.
    formats : tuple of str, optional
        File formats, by default PNG (600 dpi) and vector PDF.

    Returns
    -------
    list of str
        Paths of the files written.
    """
    os.makedirs(output_dir, exist_ok=True)
    paths = []
    for name, plot in (('fig9', plot_fig9), ('figS12', plot_figS12)):
        with plt.rc_context(_RC_PARAMS):
            fig = plot()
            for fmt in formats:
                path = os.path.join(output_dir, f'pd_c_3hp_oxidation_{name}.{fmt}')
                fig.savefig(path, dpi=600, bbox_inches='tight',
                            facecolor='white')
                paths.append(path)
            plt.close(fig)
    return paths


def main():
    parser = argparse.ArgumentParser(
        description='Regenerate Fig. 9 and SI Fig. S12 of Kim et al. (2026) '
                    'from the NSKinetics models.')
    parser.add_argument('output_dir', help='folder to write the figures to')
    parser.add_argument('--formats', nargs='+', default=['png', 'pdf'],
                        help='file formats (default: png pdf)')
    args = parser.parse_args()
    if os.environ.get('MPLBACKEND') is None:
        plt.switch_backend('Agg')    # no display needed to write files
    for path in save_figures(args.output_dir, tuple(args.formats)):
        print(f'Saved {path}')


if __name__ == '__main__':
    main()
