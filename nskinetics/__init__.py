# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

__version__ = '0.6.0'

# %% Initialize NSKinetics 

from . import engine
from .engine import *

from . import exceptions
from .exceptions import *

from . import utils

from . import tests

from . import units

from . import models

# An empty shell: the biorefinery case studies (e.g.
# ``nskinetics.biorefineries.isobutanol``) need optional dependencies and are
# imported explicitly, so this pulls in nothing heavy.
from . import biorefineries

from . import processes

from . import visualization

__all__ = (
     *exceptions.__all__,
     *engine.__all__,
     'utils',
     'tests',
     'units',
     'models',
     'biorefineries',
     'processes',
     'visualization'
     )
