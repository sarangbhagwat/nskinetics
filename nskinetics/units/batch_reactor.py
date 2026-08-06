# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np

from ..utils import get_index_nearest_element_from_sorted_array

__all__ = ('NSKBatchReactor', 'AerationSpec', 'SpikeReduceRetry',
           'select_tau_index')


def select_tau_index(results, col_names, tau, policy, n_decimals=2):
    """Select the results-array row index for the reported ``tau``.

    Parameters
    ----------
    results : numpy.ndarray
        2-D array of simulation results (rows = timesteps).
    col_names : list of str
        Column names for ``results``.
    tau : float or None
        Reaction time; used only when ``policy is None``.
    policy : None or tuple
        ``None`` -> nearest-time row. ``('max'|'min', var)`` -> first row where
        ``var`` (rounded to ``n_decimals``) equals its rounded max/min.
        ``('equals', var, value)`` -> first row where ``var`` (rounded) equals
        ``value`` (rounded).
    n_decimals : int
        Rounding used by the ``max``/``min``/``equals`` policies.

    Returns
    -------
    (index, success) : tuple[int, bool]
        ``success`` is ``False`` (and ``index`` is ``-1``) only when an
        ``'equals'`` policy finds no match.
    """
    if policy is None:
        idx = get_index_nearest_element_from_sorted_array(
            results[:, col_names.index('time')], tau)
        return idx, True

    if policy[0] in ('max', 'min'):
        col = col_names.index(policy[1])
        column = results[:, col]
        target = getattr(column, policy[0])()
        idx = np.where(np.round(column, n_decimals) == np.round(target, n_decimals))[0][0]
        return idx, True

    if policy[0] == 'equals':
        col = col_names.index(policy[1])
        column = results[:, col]
        matches = np.where(np.round(column, n_decimals) == np.round(policy[2], n_decimals))[0]
        if len(matches):
            return int(matches[0]), True
        return -1, False

    raise ValueError(f'Unknown tau_update_policy {policy!r}.')


class SpikeReduceRetry:
    """Reduce-then-retry loop for fed-batch spike counts.

    Parameters
    ----------
    max_count_var : str
        Model attribute capping the number of feed spikes (the same parameter a
        :class:`FeedSpike` points ``max_count`` at).
    stop_when : callable
        ``stop_when(model) -> bool``; the loop stops when this returns ``True``
        or the cap reaches 0.
    """
    def __init__(self, *, max_count_var='max_n_glu_spikes', stop_when=lambda r: True):
        self.max_count_var = max_count_var
        self.stop_when = stop_when

    def run(self, model, simulate_once):
        """Run ``simulate_once`` while decrementing the cap on ``model`` until
        ``stop_when(model)`` holds or the cap hits 0."""
        var = self.max_count_var
        while (not self.stop_when(model)) and (getattr(model, var) >= 0):
            simulate_once()
            setattr(model, var, getattr(model, var) - 1)


class AerationSpec:
    """Opt-in aeration / O2 model for :class:`NSKBatchReactor`.

    Parameters
    ----------
    qO2_var, is_aerobic_var, biomass_var : str
        Result column names for specific O2 uptake, the aerobic flag, and
        biomass concentration.
    volume_var : str
        Result column name for the working volume (default ``'curr_env'``).
    biomass_chemical : str
        biosteam chemical ID for biomass (default ``'Yeast'``); used to convert
        model biomass grams to biosteam moles when sizing the air stream.
    safety_factor : float
        Multiplier applied to the stoichiometric air requirement.
    air_index : int
        Index of the compressed-air inlet stream on the unit's ``ins``.
    stop_when_cell_density_plateaus : bool
        If ``True``, stop aeration at the biomass-growth plateau instead of the
        aerobic-flag cutoff.
    factor_for_cell_density_plateau : float
        Fraction-of-peak-growth-rate threshold used for plateau detection.
    stage_1_max_time, stage_1_max_x : float
        Stage-1 aerobic cutoffs (mirrored onto the model by the unit).
    """
    def __init__(self, *, qO2_var='qO2', is_aerobic_var='is_aerobic',
                 biomass_var='[x]', volume_var='curr_env', biomass_chemical='Yeast',
                 safety_factor=2.0, air_index=3,
                 stop_when_cell_density_plateaus=False,
                 factor_for_cell_density_plateau=0.5,
                 stage_1_max_time=np.inf, stage_1_max_x=np.inf):
        self.qO2_var = qO2_var
        self.is_aerobic_var = is_aerobic_var
        self.biomass_var = biomass_var
        self.volume_var = volume_var
        self.biomass_chemical = biomass_chemical
        self.safety_factor = safety_factor
        self.air_index = air_index
        self.stop_when_cell_density_plateaus = stop_when_cell_density_plateaus
        self.factor_for_cell_density_plateau = factor_for_cell_density_plateau
        self.stage_1_max_time = stage_1_max_time
        self.stage_1_max_x = stage_1_max_x

    def _tau_index_cell_density_plateau(self, unit):
        x_arr = unit.nsk_results_dict[self.biomass_var]
        growth_rate = np.diff(x_arr)
        max_growth = growth_rate.max()
        i_max = np.where(growth_rate == max_growth)[0][0]
        factor = self.factor_for_cell_density_plateau
        for i in range(i_max, len(growth_rate)):
            if growth_rate[i] < factor * max_growth:
                return i - 1
        return len(x_arr) - 1

    def resolve_stop_index(self, unit, tau_index):
        d = unit.nsk_results_dict
        if self.stop_when_cell_density_plateaus:
            plateau = self._tau_index_cell_density_plateau(unit)
            unit.tau_index_cell_density_plateau = plateau
            unit.tau_cell_density_plateau = d['time'][plateau]
            return plateau
        aerobic_off = np.where(d[self.is_aerobic_var] == 0.0)[0]
        if len(aerobic_off):
            return min(tau_index, aerobic_off[0])
        return tau_index

    def compute_cumulative_O2(self, unit, tau_index):
        d = unit.nsk_results_dict
        stop_index = self.resolve_stop_index(unit, tau_index)
        unit.tau_index_stop_aeration = stop_index
        unit.tau_stop_aeration = d['time'][stop_index]
        stepwise = d[self.qO2_var][:-1] * d[self.biomass_var][:-1] * np.diff(d['time'])
        unit._stepwise_O2 = stepwise
        unit.cumulative_O2 = stepwise[:stop_index].sum()
        return unit.cumulative_O2

    def set_air_stream(self, unit, effluent, vent):
        compressed_air = unit.compressed_air
        compressed_air.empty()
        xO2_air = 0.21
        compressed_air.imol['O2'] = xO2_air
        compressed_air.imol['N2'] = 1.0 - xO2_air
        d = unit.nsk_results_specific_tau_dict
        cumulative_gram_x = d[self.biomass_var] * d[self.volume_var]
        material_factor = effluent.imass[self.biomass_chemical] * 1e3 / cumulative_gram_x
        mmol_to_kmol = 1e-6
        stoich_O2 = mmol_to_kmol * unit.cumulative_O2 * material_factor
        unit.stoich_O2_flow_req = stoich_O2
        unit.stoich_air_flow_req = stoich_air = stoich_O2 / xO2_air
        compressed_air.F_mol = stoich_air * self.safety_factor
        vent.imol['O2'] += max(0.0, compressed_air.imol['O2'] - stoich_O2)
        vent.imol['N2'] += compressed_air.imol['N2']


class NSKBatchReactor:  # replaced in Task 9
    def __init__(self, *a, **k):
        raise NotImplementedError
