# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Complete biorefinery case studies built on NSKinetics.

Each biorefinery is a subpackage that is imported explicitly, e.g.::

    from nskinetics.biorefineries import isobutanol
    isobutanol.load()

Nothing is imported here on purpose: the biorefineries need optional
dependencies (see the ``isobutanol`` extra) and ``import nskinetics`` must not
pull them in.
"""
