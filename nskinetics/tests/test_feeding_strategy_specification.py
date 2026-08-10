# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests for the chemical/model-agnostic FedBatchStrategySpecification."""

import numpy as np
import pytest


def test_feeding_strategy_error_is_nsk_error():
    from nskinetics.exceptions import FeedingStrategyError, NSKError, __all__
    assert issubclass(FeedingStrategyError, NSKError)
    assert 'FeedingStrategyError' in __all__


class _StubUnit:
    """Bare attribute bag standing in for a biosteam unit."""
    def __init__(self, **attrs):
        self.ID = attrs.pop('ID', 'U0')
        for k, v in attrs.items():
            setattr(self, k, v)


class _StubSpikeRetry:
    """Stands in for nskinetics.units.SpikeReduceRetry."""
    max_count_var = 'max_n_glu_spikes'


class _StubKineticModel:
    """Records set_value writes the way KineticModel forwards them to _te."""
    def __init__(self):
        self.values = {}

    def set_value(self, selection, value):
        self.values[selection] = value


class _StubReactor:
    volume_var = 'curr_env'
    feed_volume_added_var = 'curr_tot_vol_added'
    spike_feed_index = 2
    spike_retry = _StubSpikeRetry()

    def __init__(self):
        self.nsk_kinetic_model = _StubKineticModel()


def test_concentration_actuator_set_get_repr():
    from nskinetics.units import ConcentrationActuator
    u = _StubUnit(ID='F301', V=0.0)
    act = ConcentrationActuator(u, 'V', 0.0, 0.8)
    act.set(0.5)
    assert u.V == 0.5
    assert act.get() == 0.5
    assert (act.lb, act.ub) == (0.0, 0.8)
    assert 'F301' in repr(act) and 'V' in repr(act)


def test_spike_control_variables_explicit_columns_win():
    from nskinetics.units import SpikeControlVariables
    cv = SpikeControlVariables(
        spike_conc_var='a', target_conc_var='b', threshold_conc_var='c',
        volume_col='my_vol', feed_volume_added_col='my_added')
    r = _StubReactor()
    assert cv.resolve_volume_col(r) == 'my_vol'
    assert cv.resolve_feed_volume_added_col(r) == 'my_added'


def test_spike_control_variables_default_from_reactor():
    from nskinetics.units import SpikeControlVariables
    cv = SpikeControlVariables(
        spike_conc_var='a', target_conc_var='b', threshold_conc_var='c')
    r = _StubReactor()
    assert cv.resolve_volume_col(r) == 'curr_env'
    assert cv.resolve_feed_volume_added_col(r) == 'curr_tot_vol_added'


def test_spike_control_variables_unresolvable_raises():
    from nskinetics.units import SpikeControlVariables
    cv = SpikeControlVariables(
        spike_conc_var='a', target_conc_var='b', threshold_conc_var='c')
    bare = _StubUnit()  # no volume_var / feed_volume_added_var
    with pytest.raises(ValueError, match='volume'):
        cv.resolve_volume_col(bare)
    with pytest.raises(ValueError, match='volume'):
        cv.resolve_feed_volume_added_col(bare)


def test_spike_control_variables_explicit_max_n_spikes_names_win():
    from nskinetics.units import SpikeControlVariables
    cv = SpikeControlVariables(
        spike_conc_var='a', target_conc_var='b', threshold_conc_var='c',
        max_n_spikes_var='my_cap', default_max_n_spikes_attr='my_default_cap')
    r = _StubReactor()
    assert cv.resolve_max_n_spikes_var(r) == 'my_cap'
    assert cv.resolve_default_max_n_spikes_attr(r) == 'my_default_cap'


def test_spike_control_variables_max_n_spikes_defaults_from_reactor():
    from nskinetics.units import SpikeControlVariables
    cv = SpikeControlVariables(
        spike_conc_var='a', target_conc_var='b', threshold_conc_var='c')
    r = _StubReactor()
    assert cv.resolve_max_n_spikes_var(r) == 'max_n_glu_spikes'
    assert cv.resolve_default_max_n_spikes_attr(r) == 'default_max_n_glu_spikes'


def test_spike_control_variables_max_n_spikes_default_attr_follows_explicit_var():
    from nskinetics.units import SpikeControlVariables
    cv = SpikeControlVariables(
        spike_conc_var='a', target_conc_var='b', threshold_conc_var='c',
        max_n_spikes_var='cap_x')
    assert cv.resolve_default_max_n_spikes_attr(_StubReactor()) == 'default_cap_x'


def test_spike_control_variables_max_n_spikes_unresolvable_raises():
    from nskinetics.units import SpikeControlVariables
    cv = SpikeControlVariables(
        spike_conc_var='a', target_conc_var='b', threshold_conc_var='c')
    bare = _StubUnit()  # no spike_retry
    with pytest.raises(ValueError, match='spike'):
        cv.resolve_max_n_spikes_var(bare)
    with pytest.raises(ValueError, match='spike'):
        cv.resolve_default_max_n_spikes_attr(bare)


class _Indexer:
    def __init__(self, data):
        self._data = data
    def __getitem__(self, key):
        if isinstance(key, (list, tuple)):
            return np.array([self._data[k] for k in key])
        return self._data[key]


class _StubStream:
    def __init__(self, imass, ivol):
        self.imass = _Indexer(imass)
        self.ivol = _Indexer(ivol)


def _make_spec(target_conc=220.0, threshold_conc=210.0, spike_conc=600.0,
               **overrides):
    from nskinetics.units import (FedBatchStrategySpecification,
                                  SpikeControlVariables, ConcentrationActuator)
    cv = SpikeControlVariables(
        spike_conc_var='conc_spike', target_conc_var='target_conc_v',
        threshold_conc_var='threshold_conc_v')
    kwargs = dict(
        target_conc=target_conc, threshold_conc=threshold_conc,
        spike_conc=spike_conc, tau_max=72.0,
        fermentation_reactor=_StubReactor(), splitter=_StubUnit(ID='S301'),
        control_variables=cv,
        feed_concentrator=ConcentrationActuator(_StubUnit(ID='F301', V=0.0), 'V', 0.0, 0.8),
        feed_diluter=ConcentrationActuator(_StubUnit(ID='M301', dil=0.0), 'dil', 0.0, 1e5),
        spike_concentrator=ConcentrationActuator(_StubUnit(ID='F302', V=0.0), 'V', 0.0, 0.8),
        spike_diluter=ConcentrationActuator(_StubUnit(ID='M302', dil=0.0), 'dil', 0.0, 1e5),
        species_IDs=['Glucose'],
    )
    kwargs.update(overrides)
    return FedBatchStrategySpecification(**kwargs)


def test_spec_constructor_and_current_specifications():
    spec = _make_spec()
    assert spec.current_specifications == {
        'target_conc': 220.0, 'threshold_conc': 210.0,
        'spike_conc': 600.0, 'tau_max': 72.0, 'max_n_spikes': None}


def test_spec_validation_rejects_bad_ordering():
    with pytest.raises(ValueError):
        _make_spec(threshold_conc=230.0)          # threshold > target
    with pytest.raises(ValueError):
        _make_spec(spike_conc=100.0)              # target > spike
    with pytest.raises(ValueError):
        _make_spec(target_conc=-1.0)              # non-positive


def test_load_specifications_rejects_bad_ordering_before_simulating():
    spec = _make_spec()
    with pytest.raises(ValueError):
        spec.load_specifications(threshold_conc=250.0)


def test_get_conc_respects_species_and_solvent():
    spec = _make_spec(species_IDs=['Glucose', 'Xylose'], solvent_ID='Ethanol')
    stream = _StubStream(imass={'Glucose': 30.0, 'Xylose': 10.0},
                         ivol={'Ethanol': 2.0, 'Water': 100.0})
    assert spec.get_conc(stream) == pytest.approx(20.0)


def test_spec_stores_max_n_spikes():
    spec = _make_spec(max_n_spikes=12)
    assert spec.max_n_spikes == 12
    assert spec.current_specifications['max_n_spikes'] == 12


def test_spec_rejects_negative_max_n_spikes():
    with pytest.raises(ValueError, match='max_n_spikes'):
        _make_spec(max_n_spikes=-1)


def test_load_specifications_rejects_negative_max_n_spikes():
    spec = _make_spec()
    with pytest.raises(ValueError, match='max_n_spikes'):
        spec.load_specifications(max_n_spikes=-3)


def test_load_max_n_spikes_none_writes_nothing():
    spec = _make_spec()
    km = spec.fermentation_reactor.nsk_kinetic_model
    spec.load_max_n_spikes(None)
    assert km.values == {}
    assert not hasattr(km, 'default_max_n_glu_spikes')


def test_load_max_n_spikes_writes_both_targets():
    spec = _make_spec()
    km = spec.fermentation_reactor.nsk_kinetic_model
    spec.load_max_n_spikes(12)
    assert km.values['max_n_glu_spikes'] == 12
    assert km.default_max_n_glu_spikes == 12


def test_load_specifications_imposes_cap_before_the_deriving_run():
    """load_specifications must impose the cap, and impose it BEFORE the run
    that derives the splitter split.

    load_desired_concs and load_threshold_conc_and_tau_max are replaced on the
    instance with recorders (they need a simulating reactor the stub cannot
    provide); load_max_n_spikes is left REAL, so this pins the actual cap
    writes and their real ordering rather than asserting a mock was called.
    """
    spec = _make_spec()
    km = spec.fermentation_reactor.nsk_kinetic_model
    seen = {}

    def _snapshot():
        return (dict(km.values), getattr(km, 'default_max_n_glu_spikes', None))

    def _record_desired_concs(target_conc, spike_conc):
        seen['at_actuators'] = _snapshot()

    def _record_threshold_and_tau(threshold_conc, tau_max):
        seen['at_deriving_run'] = _snapshot()

    spec.load_desired_concs = _record_desired_concs
    spec.load_threshold_conc_and_tau_max = _record_threshold_and_tau

    spec.load_specifications(max_n_spikes=12)

    # The cap was imposed at all, in both places it must live.
    assert km.values['max_n_glu_spikes'] == 12
    assert km.default_max_n_glu_spikes == 12

    # ...and was already imposed by the time the deriving run would have gone,
    # and before the concentration actuators too.
    values_at_run, default_at_run = seen['at_deriving_run']
    assert values_at_run.get('max_n_glu_spikes') == 12
    assert default_at_run == 12

    values_at_actuators, default_at_actuators = seen['at_actuators']
    assert values_at_actuators.get('max_n_glu_spikes') == 12
    assert default_at_actuators == 12


def test_load_max_n_spikes_accepts_zero_and_floats():
    spec = _make_spec()
    km = spec.fermentation_reactor.nsk_kinetic_model
    spec.load_max_n_spikes(0)
    assert km.values['max_n_glu_spikes'] == 0
    spec.load_max_n_spikes(12.0)
    assert km.values['max_n_glu_spikes'] == 12.0
