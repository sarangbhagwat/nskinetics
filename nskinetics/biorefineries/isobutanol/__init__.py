# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Corn-to-ethanol + isobutanol biorefinery built around NSKinetics' fed-batch
fermentation model.

Importing this package is **build-free**: it imports ``system``, ``models``,
``scenarios``, ``utils`` and ``process_settings`` (and, through ``system``,
``separations`` and ``enzyme_burden``) but builds and simulates nothing; the
submodules marked "not imported here" below must be imported explicitly.
Call :func:`load` once per process to build the
flowsheet, simulate it to its baseline state and build the uncertainty
``Model``::

    from nskinetics.biorefineries import isobutanol
    isobutanol.load()                       # both separation trains (default)
    bundle = isobutanol.load_scenario('B')  # workbook kinetics + feeding strategy
    isobutanol.solve_TEA()

Submodules
----------
system
    Flowsheet build (:func:`system.load`), TEA solves, the convergence loop
    and its recovery scaffolding.
models
    The BioSTEAM ``Model`` (metrics + the workbook ``exec`` namespace).
scenarios
    The scenario registry (workbooks, feeding strategies, smoke-test pins).
separations, process_settings, utils
    Separation-train factories, cost/price-year settings, workbook generator.
enzyme_burden
    Proteome-burden constraint (imported through ``system``; also loadable
    by file path without importing this package).
kinetic_optimization, sensitivity_analysis
    Kinetic Bayesian optimization and Sobol' sensitivity analysis; not
    imported here (``from nskinetics.biorefineries.isobutanol import
    kinetic_optimization``, or load them by file path without importing
    this package).
smoke
    The 16 pinned smoke-test cases (``python -m
    nskinetics.biorefineries.isobutanol.smoke N``); not imported here.
plots
    Plotting library and figure scripts; not imported here (import e.g.
    ``nskinetics.biorefineries.isobutanol.plots.spearman_matrix`` explicitly).
"""
import types as _types

from . import (
    models,
    scenarios,
    system,
    utils,
    process_settings
    )

__all__ = [
            *models.__all__,
            *scenarios.__all__,
            *system.__all__,
            *utils.__all__,
            *process_settings.__all__,
]

from .models import *
from .scenarios import *
from .system import *
from .utils import *
from .process_settings import *

# NOTE: `load` must be defined AFTER `from .system import *`, which binds
# `load` to system.load (that one does not build the uncertainty Model).
def load(simulate_baseline=True,
         separation_processes=('IBO_EtOH', 'ethanol')):
    """Build and baseline-simulate the biorefinery and its uncertainty Model.

    Runs :func:`system.load` (flowsheet build + baseline simulation) and
    :func:`models.models_EtOH_IBO_corn.create_model`, then publishes every
    built name (``corn_EtOH_IBO_sys``, ``fbs_spec``, ``V406``, ``model``,
    ``namespace_dict``, ...) into this package's namespace, into
    ``nskinetics.biorefineries.isobutanol.models`` and into the defining
    modules. Names of submodules (``system``, ``models``, ...) are never
    overwritten by a published object of the same name.

    Parameters
    ----------
    simulate_baseline : bool, optional
        If True (default), run the baseline ``model_specification`` at the
        end of the build.
    separation_processes : Iterable[str], optional
        Non-empty subset of ``('IBO_EtOH', 'ethanol')``: which separation
        train(s) to build. Default both (the gated-parallel baseline
        configuration; see :func:`system.load`).

    Returns
    -------
    dict
        All published names, the ``create_model`` names winning ties.

    Notes
    -----
    One configuration per process: a repeat call is a no-op that ignores a
    different ``separation_processes`` (the build mutates global thermo,
    price and flowsheet state). Run a different configuration in a fresh
    process.
    """
    published = system.load(simulate_baseline=simulate_baseline,
                            separation_processes=separation_processes)
    published_models = models.models_EtOH_IBO_corn.create_model()
    models.__dict__.update(published_models)
    merged = {**published, **published_models}
    # Protect every module attribute (submodules, including any added
    # later) from being shadowed by a same-named build object.
    protected = {k for k, v in globals().items()
                 if isinstance(v, _types.ModuleType)}
    globals().update({k: v for k, v in merged.items()
                      if k not in protected})
    return merged
