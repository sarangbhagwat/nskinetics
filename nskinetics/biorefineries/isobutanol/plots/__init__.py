# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Figures of the isobutanol biorefinery analyses.

The two plotting libraries (:mod:`.spearman_matrix` and
:mod:`.trajectory_contourplot`, which need the optional ``contourplots``
package) are imported lazily, on first access of one of their public names or
of the submodule itself, so the biorefinery package never pulls in
``contourplots``. Importing anything under this package is NOT light,
though: Python first imports the parent ``nskinetics.biorefineries.isobutanol``
package, whose ``system`` module imports BioSTEAM, ``hensmith`` and the
upstream ``biorefineries`` package (building and simulating nothing), e.g.::

    from nskinetics.biorefineries.isobutanol.plots.spearman_matrix import plot_spearman_matrix

Beside a running simulation, load the library file by path instead
(``paths.load_by_path('spearman_matrix', <this folder>/spearman_matrix.py)``).

Every other file in this folder is a stand-alone figure script. Run it as a
file (``python <path>/plot_<name>.py``): the scripts load their siblings,
``paths.py``, ``kinetic_optimization.py`` and ``enzyme_burden.py`` BY FILE
PATH, so they do not import ``nskinetics``, BioSTEAM or the biorefinery and
are safe to run beside a simulation. The exception is
``plot_feeding_strategy_baseline_A.py``, whose stage 1 builds and simulates
scenario A when its trajectory cache is missing (or with ``--resimulate``).
Default figure outputs go to ``paths.publication_dir(<theme>)``.
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
