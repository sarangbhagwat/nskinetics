# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

# Kinetic models of 3-HP oxidation to malonic acid over Pd/C (Kim et al.,
# Appl. Catal. B: Environ. Energy 386 (2026) 126403). Four models, each
# with a loader in ``pd_c_3hp_oxidation``:
#
# * H2O2 model (``load_h2o2_model``; Table 3, H2O2 column): reproduces
#   Fig. 9a. Its predicted H2O2 is compared with SI Fig. S12a, whose H2O2
#   curve lies up to 0.38 equivalents lower;
# * O2 model (``load_o2_model``; Table 3, O2 column): reproduces Fig. 9b
#   and SI Fig. S12b at 50 C. Its Arrhenius temperature dependence is an
#   NSKinetics addition and reproduces no published item;
# * reversible O2 variant (``load_o2_reversible_model``), inferred from the
#   deposited model curves: reproduces the curves of Fig. 10a-c (25, 75
#   and 100 C), one constant set per curve;
# * simplified reversible O2 model (``load_o2_simplified_reversible_model``)
#   of the authors' heat-map script: reproduces Fig. 11; the 50 C point of
#   Fig. 10d supplies three of its rate constants.
#
# Time is in min, concentrations are in M and temperatures are in K where an
# argument is named ``T`` and in degrees Celsius where its name ends in
# ``_C``.
#
# ``NOTES.md`` in this folder records what is printed in the paper, what was
# inferred from the authors' data deposit and where the two disagree; the
# documentation page "Pd/C oxidation of 3-HP to malonic acid"
# (docs/source/case_studies/pd_c_3hp_oxidation.rst) summarizes it.
#
# Importing this subpackage loads no model and runs no simulation: the
# loaders build a model when called. The figure module (``figures``) is not
# imported here; import it explicitly, or run it with ``python -m`` to
# write the figures.

from . import pd_c_3hp_oxidation
from .pd_c_3hp_oxidation import *

from . import reference
from .reference import *

from . import heat_maps
from .heat_maps import *

from . import temperature_series
from .temperature_series import *

__all__ = (
     *pd_c_3hp_oxidation.__all__,
     *reference.__all__,
     *heat_maps.__all__,
     *temperature_series.__all__,
     )
