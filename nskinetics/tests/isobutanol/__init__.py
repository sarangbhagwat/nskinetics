# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

# Tests of the isobutanol biorefinery (nskinetics.biorefineries.isobutanol).
#
# Deliberately NOT registered in nskinetics/tests/__init__.py (which a plain
# `import nskinetics` executes), and no module here imports the biorefinery,
# biosteam or the upstream `biorefineries` package at module level: simulating
# tests run each script in a fresh subprocess (see _support.run_isolated), and
# offline tests load the sim-safe modules by file path.
