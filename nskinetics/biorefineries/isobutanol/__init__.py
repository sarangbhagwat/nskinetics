# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from . import (
    models,
    plots,
    scenarios,
    system,
    units,
    utils,
    process_settings
    )

__all__ = [
            *models.__all__,
            *plots.__all__,
            *scenarios.__all__,
            *system.__all__,
            *units.__all__,
            *utils.__all__,
            *process_settings.__all__,
]

from .models import *
from .plots import *
from .scenarios import *
from .system import *
from .units import *
from .utils import *
from .process_settings import *

def load(simulate_baseline=True,
         separation_processes=('IBO_EtOH', 'ethanol')):
    """Build and baseline-simulate the biorefinery (see system.load;
    `separation_processes` selects which separation train(s) are built --
    non-empty subset of ('IBO_EtOH', 'ethanol'), one configuration per
    kernel, a repeat call ignores different arguments), build the
    uncertainty Model (see models.models_EtOH_IBO_corn.create_model), and
    publish every built name here, in nskinetics.biorefineries.isobutanol.models, and
    in their defining modules -- reproducing the namespaces the former
    import-time build created. Submodule names (system, models, ...) are
    never overwritten. Idempotent; returns the merged published dict."""
    published = system.load(simulate_baseline=simulate_baseline,
                            separation_processes=separation_processes)
    published_models = models.models_EtOH_IBO_corn.create_model()
    models.__dict__.update(published_models)
    merged = {**published, **published_models}
    _submodules = {'system', 'models', 'plots', 'units', 'utils',
                   'process_settings', 'separations', 'scenarios'}
    globals().update({k: v for k, v in merged.items()
                      if k not in _submodules})
    return merged
