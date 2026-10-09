# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Figures of the isobutanol biorefinery analyses.

Importing this package imports nothing heavy: the two plotting libraries
(:mod:`.spearman_matrix` and :mod:`.trajectory_contourplot`, which need the
optional ``contourplots`` package) are imported lazily, on first access of
one of their public names or of the submodule itself, e.g.::

    from nskinetics.biorefineries.isobutanol.plots.spearman_matrix import plot_spearman_matrix

Every other file in this folder is a stand-alone figure script. Run it as a
file (``python <path>/plot_<name>.py``): the scripts load their siblings,
``paths.py``, ``kinetic_optimization.py`` and ``enzyme_burden.py`` BY FILE
PATH, so they never import ``nskinetics``, BioSTEAM or the biorefinery and are
safe to run beside a simulation. Default figure outputs go to
``paths.publication_dir(<theme>)``.
"""
import importlib

#: public name -> submodule that defines it
_LAZY_NAMES = {
    'plot_spearman_matrix': 'spearman_matrix',
    'plot_metric_with_trajectories': 'trajectory_contourplot',
    'JBEI_UCB_colormap': 'trajectory_contourplot',
}

#: submodules importable as attributes of this package
_LAZY_SUBMODULES = ('spearman_matrix', 'trajectory_contourplot')

__all__ = tuple(_LAZY_NAMES)


def __getattr__(name):
    if name in _LAZY_SUBMODULES:
        return importlib.import_module(f'{__name__}.{name}')
    if name in _LAZY_NAMES:
        module = importlib.import_module(f'{__name__}.{_LAZY_NAMES[name]}')
        return getattr(module, name)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')


def __dir__():
    return sorted(set(globals()) | set(_LAZY_NAMES) | set(_LAZY_SUBMODULES))
