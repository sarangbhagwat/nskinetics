# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Filesystem locations used by the isobutanol biorefinery.

This module is **stdlib-only and has no package-relative imports**, so it can
be loaded by file path (see :func:`load_by_path`) from code that must stay
"sim-safe" -- ``kinetic_optimization.py``, the plotting scripts and the offline
tests load their dependencies by path precisely so that ``nskinetics``,
``biosteam`` and the upstream ``biorefineries`` package are never imported.

Every analysis output goes under one flat results root,
``<package>/analyses/results`` by default. Set the environment variable
``NSKINETICS_IBO_RESULTS_DIR`` to put it elsewhere, e.g. to reuse the results
of an earlier run. Keep that root short: kinetic-optimization file names are
long and Windows paths are limited to 260 characters.
"""
import importlib.util
import os

__all__ = (
    'PACKAGE_DIR',
    'ANALYSES_DIR',
    'DATA_DIR',
    'PLOTS_DIR',
    'RESULTS_DIR_ENV',
    'DEFAULT_RESULTS_DIR',
    'results_dir',
    'publication_dir',
    'analyses_dir',
    'workbook_dir',
    'workbook_path',
    'require_input',
    'load_by_path',
)

#: Folder of the ``nskinetics.biorefineries.isobutanol`` package.
PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

#: Folder holding the analysis scripts (thematic subfolders).
ANALYSES_DIR = os.path.join(PACKAGE_DIR, 'analyses')

#: Folder holding package data read at run time (parameter workbooks, bands).
DATA_DIR = os.path.join(PACKAGE_DIR, 'data')

#: Folder holding the plotting modules and figure scripts.
PLOTS_DIR = os.path.join(PACKAGE_DIR, 'plots')

#: Environment variable that overrides the results root.
RESULTS_DIR_ENV = 'NSKINETICS_IBO_RESULTS_DIR'

#: Results root used when ``NSKINETICS_IBO_RESULTS_DIR`` is not set.
DEFAULT_RESULTS_DIR = os.path.join(ANALYSES_DIR, 'results')


def results_dir(*parts, create=False):
    """Return a path under the results root.

    Parameters
    ----------
    *parts : str
        Path components appended to the results root.
    create : bool, optional
        If True, create the resulting directory (and parents). Only pass
        ``create=True`` when ``parts`` name a directory.

    Returns
    -------
    str
        ``<results root>/<parts...>``, where the root is
        ``$NSKINETICS_IBO_RESULTS_DIR`` if set, else :data:`DEFAULT_RESULTS_DIR`.
    """
    root = os.environ.get(RESULTS_DIR_ENV) or DEFAULT_RESULTS_DIR
    path = os.path.join(root, *parts)
    if create:
        os.makedirs(path, exist_ok=True)
    return path


def publication_dir(*parts, create=False):
    """Return a path under ``<results root>/publication`` (see :func:`results_dir`)."""
    return results_dir('publication', *parts, create=create)


def analyses_dir(*parts):
    """Return a path under the analyses folder, e.g. ``analyses_dir('sensitivity')``."""
    return os.path.join(ANALYSES_DIR, *parts)


def workbook_dir():
    """Return the folder holding the parameter-distribution workbooks."""
    return os.path.join(DATA_DIR, 'parameter_distributions')


def workbook_path(scenario):
    """Return the parameter-distribution workbook path for a scenario.

    Parameters
    ----------
    scenario : str
        Scenario key, e.g. ``'A'``, ``'B'`` or ``'opt_PI_TRY_informed'``.

    Returns
    -------
    str
        ``<data>/parameter_distributions/parameter-distributions_corn_IBO_EtOH_<scenario>.xlsx``
    """
    return os.path.join(workbook_dir(),
                        f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx')


def require_input(path, produced_by):
    """Return ``path`` if it exists, else raise an informative FileNotFoundError.

    Only ``analyses/results/publication`` ships with the repository; every
    other analysis output must be regenerated (or ``NSKINETICS_IBO_RESULTS_DIR``
    pointed at an existing results folder) before scripts that read it can run.

    Parameters
    ----------
    path : str
        File or directory the caller needs.
    produced_by : str
        What generates it, e.g. ``'analyses/sensitivity/evaluate_sobol_split12d.py'``.

    Returns
    -------
    str
        ``path``, unchanged.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f'missing input {path!r}: it is produced by {produced_by}. Run that '
            f'analysis first, or set {RESULTS_DIR_ENV} to a results folder that '
            f'already contains it.')
    return path


def load_by_path(name, path, register=False):
    """Import a module from a file path without importing its parent package.

    Parameters
    ----------
    name : str
        Name to give the module object.
    path : str
        Path of the ``.py`` file.
    register : bool, optional
        If True, add the module to ``sys.modules`` under ``name`` before
        executing it (needed e.g. when the module defines dataclasses).

    Returns
    -------
    module
    """
    import sys
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    if register:
        sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
