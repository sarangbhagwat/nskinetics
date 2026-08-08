# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

# Re-export the model module's public names so the historical import path
# ``from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r``
# keeps resolving now that the model and its Antimony/SBML data files live in
# this subpackage. Importing this subpackage builds and configures the model
# (it is heavy on purpose — downstream consumers import it for ``te_r``).

from . import s_cerevisiae_ferm_fb_inhib_mod_ibo
from .s_cerevisiae_ferm_fb_inhib_mod_ibo import *

__all__ = (
     *s_cerevisiae_ferm_fb_inhib_mod_ibo.__all__,
     )
