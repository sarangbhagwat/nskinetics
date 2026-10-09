# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Uncertainty/sensitivity ``Model`` of the isobutanol biorefinery.

Build-free at import. After ``nskinetics.biorefineries.isobutanol.load()``
this subpackage's namespace also holds every name published by
:func:`models_EtOH_IBO_corn.create_model` (``namespace_dict``, ``fbs_spec``,
``plot_kinetic_results``, the ``optimize_*`` feeding-strategy optimizers,
...); the package ``load()`` copies them here.
"""
from . import models_EtOH_IBO_corn

__all__ = [*models_EtOH_IBO_corn.__all__]

from .models_EtOH_IBO_corn import *
