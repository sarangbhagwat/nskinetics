# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from . import (
    spearman_matrix,
    trajectory_contourplot,
    # system,
    # units,
    # utils,
    )

__all__ = [
            *spearman_matrix.__all__,
            *trajectory_contourplot.__all__,
            # *system.__all__,
            # *units.__all__,
            # *utils.__all__,
]

from .spearman_matrix import *
from .trajectory_contourplot import *
# from .system import *
# from .units import *
# from .utils import *

# def load(*args, **kwargs):
#     br = Biorefinery(*args, **kwargs)
#     globals().update(br.__dict__)
#     globals().update({
#         'biorefinery': br,
#         'system': br.system,
#         'tea': br.TEA,
#     })
