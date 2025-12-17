# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

__version__ = '0.2.4'

# %% Initialize NSKinetics 

from . import species
from .species import *

from . import reactions
from .reactions import *

from . import reaction_systems 
from .reaction_systems import *

from . import doe
# from .doe import *

from . import utils

from . import tests

from . import units

__all__ = (
     *species.__all__,
     *reactions.__all__,
     *reaction_systems.__all__,
     'utils',
     'doe',
     'tests',
     'units'
     )
