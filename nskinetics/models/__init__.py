# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

# This file intentionally imports nothing. Its sole purpose is to make
# ``nskinetics.models`` a real (non-namespace) subpackage so that
# ``setuptools.find_packages()`` includes it in the built distribution and
# ``from . import models`` (in nskinetics/__init__.py) resolves on an
# installed copy of NSKinetics — not just on an editable clone.
#
# The individual model scripts are NOT imported here on purpose: each runs a
# full simulation and imports heavy optional dependencies (biosteam, tellurium,
# simplesbml) at module level, and importing them during package initialization
# would be slow and would create a circular import back into nskinetics.
