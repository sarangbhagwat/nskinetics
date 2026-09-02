# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Flux-map tests exercising a fully simulated ``NSKBatchReactor``.

Marked ``slow``: the module-scoped fixture builds and simulates the shipped
sugar-prep + fermentation system, which imports the shipped ``te_r`` kinetic
model and the corn-biorefinery chemical set. Every heavy import therefore
lives inside the fixture/test bodies -- ``nskinetics.tests`` is imported by
``nskinetics/__init__.py``, so a module-level ``biosteam``/``biorefineries``
import here would pollute ``sys.modules`` and break the import-guard tests.
"""

import pytest

pytestmark = pytest.mark.slow


@pytest.fixture(scope='module')
def simulated_V406():
    pytest.importorskip('biorefineries')
    import biosteam as bst
    import thermosteam as tmo
    from biorefineries import corn
    chems = [c for c in corn.chemicals.create_chemicals()]
    chems.append(tmo.Chemical('Isobutanol'))
    chems.append(tmo.Chemical('AceticAcid'))
    tmo.settings.set_thermo(chems)
    bst.main_flowsheet.set_flowsheet('test_flux_map')
    from nskinetics.processes import create_sugar_prep_and_fermentation_system
    sys = create_sugar_prep_and_fermentation_system()
    sys.simulate()
    return sys.flowsheet.unit.V406


def test_reactor_records_state_selection_columns(simulated_V406):
    V406 = simulated_V406
    km = V406.nsk_kinetic_model
    cols = set(V406.nsk_results_df.columns)
    for sel in km.state_selections():
        assert sel in cols, f'missing state column {sel}'
