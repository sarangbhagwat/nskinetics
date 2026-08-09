# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Fast contract test pinning the public interface of the
``nskinetics.processes`` sugar-prep + fermentation system factory: its factory
type, the IDs and order of its declared ``ins``/``outs``, and its fixed inlet
and outlet sizes. This is the interface a drop-in replacement of the isobutanol
biorefinery's inline system relies on, so it is pinned separately from the
``slow`` structure tests in ``test_processes``.

Deliberately needs only ``biosteam``: no chemical set, no thermo, no
simulation, and no ``te_r`` kinetic model (the factory only imports it when
*called*, which this module never does). It is therefore *not* ``slow``-marked
and runs in the default ``-m "not slow"`` suite and in CI.
"""

import sys
import inspect

import biosteam as bst

from nskinetics.processes import create_sugar_prep_and_fermentation_system

EXPECTED_INS_IDS = ['saccharified_slurry', 'seed']

EXPECTED_OUTS_IDS = ['fermentation_vent',
                     'fermentation_effluent',
                     'initial_feed_condensate',
                     'spike_feed_condensate']


def test_factory_interface():
    """The factory is a biosteam SystemFactory declaring fixed-size ins/outs
    with exactly these IDs, in this order."""
    f = create_sugar_prep_and_fermentation_system
    assert isinstance(f, bst.SystemFactory)
    assert [i['ID'] for i in f.ins] == EXPECTED_INS_IDS
    assert [o['ID'] for o in f.outs] == EXPECTED_OUTS_IDS
    assert f.fixed_ins_size is True
    assert f.fixed_outs_size is True


def test_import_is_lightweight():
    """Importing the factory must not pull in the shipped kinetic model or the
    biorefineries package: both are call-time (and slow-test) dependencies."""
    assert not [m for m in sys.modules
                if 's_cerevisiae_ferm_fb_inhib_mod_ibo' in m]
    assert 'biorefineries' not in sys.modules


EXPECTED_FBS_KWARG_DEFAULTS = {
    'target_conc': 220.0,
    'threshold_conc': 210.0,
    'spike_conc': 600.0,
    'spike_control_variables': None,
    'baseline_specifications': None,
    'fbs_spec_kwargs': None,
}


def test_fed_batch_strategy_kwargs_in_signature():
    """The factory's wrapped function accepts the fed-batch-strategy kwargs
    with defaults matching the isobutanol biorefinery's inline
    FedBatchStrategySpecification construction. Signature-only on purpose:
    introspecting ``.f`` keeps this module free of call-time dependencies
    (kinetic model, chemical set)."""
    params = inspect.signature(
        create_sugar_prep_and_fermentation_system.f).parameters
    for name, default in EXPECTED_FBS_KWARG_DEFAULTS.items():
        assert name in params, f'missing factory kwarg {name!r}'
        assert params[name].default == default, (
            f'{name!r} default is {params[name].default!r}, '
            f'expected {default!r}')


def test_set_thermo_flag_contract():
    """The factory accepts an opt-in set_thermo keyword, defaulting to False,
    handled by a SystemFactory subclass BEFORE stream creation (a body-level
    kwarg would come too late: SystemFactory.__call__ creates the declared
    ins/outs streams before the wrapped function runs)."""
    from nskinetics.processes.sugar_prep_and_fermentation import (
        _SetThermoSystemFactory)
    f = create_sugar_prep_and_fermentation_system
    assert isinstance(f, _SetThermoSystemFactory)
    assert isinstance(f, bst.SystemFactory)  # existing contract preserved
    params = inspect.signature(_SetThermoSystemFactory.__call__).parameters
    assert params['set_thermo'].default is False
