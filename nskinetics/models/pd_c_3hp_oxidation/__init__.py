# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

# Kinetic models of 3-HP oxidation to malonic acid over Pd/C (Kim et al.,
# Appl. Catal. B: Environ. Energy 386 (2026) 126403). Importing this
# subpackage loads no model and runs no simulation: the loaders build a model
# when called. The figure module (``figures``) is not imported here; import
# it explicitly, or run it with ``python -m`` to write the figures.

from . import pd_c_3hp_oxidation
from .pd_c_3hp_oxidation import *

from . import reference
from .reference import *

__all__ = (
     *pd_c_3hp_oxidation.__all__,
     *reference.__all__,
     )
