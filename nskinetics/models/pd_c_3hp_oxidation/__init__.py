# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

# Kinetic models of 3-HP oxidation to malonic acid over Pd/C (Kim et al.,
# Appl. Catal. B: Environ. Energy 386 (2026) 126403). Three models, each
# with a loader in ``pd_c_3hp_oxidation``:
#
# * H2O2 model (``load_h2o2_model``; Table 3, H2O2 column): reproduces
#   Fig. 9a. Its predicted H2O2 is compared with SI Fig. S12a, whose H2O2
#   curve lies up to 0.38 equivalents lower;
# * O2 model (``load_o2_model``; Table 3, O2 column): reproduces Fig. 9b
#   and SI Fig. S12b. It is the published model at 50 C and has no
#   temperature parameter;
# * simplified reversible O2 model (``load_o2_simplified_reversible_model``;
#   the name is this package's): the model of the authors' deposited
#   heat-map script (Fig. 11), with its temperature dependence. Every
#   constant the deposit gives is used as given (where it gives one twice,
#   the full-precision value): k_1, k_2 and k_5 from the deposited Arrhenius
#   table, k_4 = 0.0023 and K_6 to K_10 from the vector k0 written in the
#   script. The one constant it does not give, k_3 (0.0000 in k0), is
#   fitted to the 50 C data of the script
#   (``fit_simplified_reversible_constants`` in ``fitting``).
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
# ``solver_emulation`` repeats the solver run behind the model curves the
# authors deposited for Fig. 9 (a port of MATLAB ``ode23s`` driven by the
# rates of the shipped models); it is used for verification.
#
# Importing this subpackage loads no model, runs no simulation and fits
# nothing: the loaders build a model when called. The figure module
# (``figures``) is not imported here; import it explicitly, or run it with
# ``python -m`` to write the figures.

from . import pd_c_3hp_oxidation
from .pd_c_3hp_oxidation import *

from . import reference
from .reference import *

from . import heat_maps
from .heat_maps import *

from . import temperature_series
from .temperature_series import *

from . import solver_emulation
from .solver_emulation import *

from . import fitting
from .fitting import *

__all__ = (
     *pd_c_3hp_oxidation.__all__,
     *reference.__all__,
     *heat_maps.__all__,
     *temperature_series.__all__,
     *solver_emulation.__all__,
     *fitting.__all__,
     )
