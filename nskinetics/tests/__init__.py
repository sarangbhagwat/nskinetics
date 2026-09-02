# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

# NOTE: `nskinetics/__init__.py` imports this package, so every module
# registered below is imported by a plain `import nskinetics`. Test modules
# that import pytest, biosteam or the shipped kinetic model at top level, or
# that otherwise have import-time side effects, are therefore deliberately NOT
# registered here (currently `test_processes`, `test_flux_map` and
# `test_flux_map_render`): registering them would make `import nskinetics`
# fail on an install without those extras, would pollute `sys.modules` for the
# import-guard tests, and -- in the case of `test_flux_map_render`, which calls
# `matplotlib.use('Agg')` at import time -- would silently force a headless
# matplotlib backend on every user who imports the package. pytest still
# collects them by path.

from . import test_events
from . import test_kinetic_model_reset
from . import test_flux_analysis

__all__ = (
     'test_events',
     'test_kinetic_model_reset',
     'test_flux_analysis',
     )
