# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np

from biosteam.units import BatchBioreactor

from ..utils import get_index_nearest_element_from_sorted_array
from ..engine.kinetic_model import KineticModel
from ..exceptions import KineticSimulationError, MassBalanceError

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
        :class:`FeedSpike` points ``max_count`` at). Required.
    count_var : str, optional
        Result-column name reporting the current spike count. Read by
        :class:`NSKBatchReactor` to seed the cap before the retry loop; not used
        by :meth:`run` itself.
    stop_when : callable
        ``stop_when(model) -> bool``; the loop stops when this returns ``True``
        or the cap reaches 0.
    """
    def __init__(self, *, max_count_var, count_var=None, stop_when=lambda r: True):
        self.max_count_var = max_count_var
        self.count_var = count_var
        self.stop_when = stop_when

    def run(self, model, simulate_once):
        """Run ``simulate_once`` while decrementing the cap on ``model`` until
        ``stop_when(model)`` holds or the cap hits 0."""
        var = self.max_count_var
        while (not self.stop_when(model)) and (getattr(model, var) > 0):
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
    biomass_chemical : str, optional
        biosteam chemical ID for biomass; used to convert model biomass grams to
        biosteam moles when sizing the air stream. Required whenever the air
        stream is sized (:meth:`set_air_stream`); defaults to ``None``.
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
                 biomass_var='[x]', volume_var='curr_env', biomass_chemical=None,
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

    def result_columns(self):
        """Result columns this aeration model reads from the simulation."""
        return (self.qO2_var, self.is_aerobic_var, self.biomass_var,
                self.volume_var)

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


class NSKBatchReactor(BatchBioreactor):
    """
    General batch bioreactor driven by a Tellurium kinetic model.

    Owns the reusable machinery (kinetic simulate loop, ``tau`` selection,
    species-to-chemical mapping, effluent construction) and four opt-in
    features: aeration, fed-batch spike-count retry, feed pre-reactions, and
    post-simulation validators. Chemistry-agnostic: subclasses configure the
    features and may override :meth:`_finalize_effluent`.

    The core kinetic integration is delegated to the kinetic model's
    :meth:`~nskinetics.KineticModel.simulate`, which is the source of
    truth for full-trajectory results; the ``nsk_results``,
    ``nsk_results_dict``, ``nsk_results_col_names``, and ``nsk_results_df``
    attributes are read-only properties delegating to it (the full-trajectory
    DataFrame is exposed as ``nsk_results_df`` so it does not shadow biosteam's
    ``Unit.results()`` design/cost method). The ``tau``-selected results
    (``nsk_results_specific_tau`` / ``nsk_results_specific_tau_dict``) remain
    owned by the reactor. :meth:`plot_simulation_results` (aliases
    :meth:`plot_time_course`, :meth:`plot_trajectory`) plots the run with the
    fermentation- and aeration-end times marked.

    Parameters
    ----------
    ins :
        Inlet fluids to be mixed into the reactor.
    outs :
        * [0] Vent
        * [1] Effluent
    nsk_kinetic_model : KineticModel
        Kinetic model driving the reactor.
    map_species_to_chemicals : dict
        ``{model var (or '[conc]' selection): biosteam chemical ID}``.
    track_vars : sequence of str
        Extra model selections to record as result columns.
    tau : float
        Reaction time [h].
    tau_max : float
        Maximum simulated time [h].
    n_simulation_steps : int
        Number of integration output steps.
    tau_update_policy : None or tuple
        Passed to :func:`select_tau_index`.
    n_decimal_places_for_tau_update_policy : int
        Rounding for the ``max``/``min``/``equals`` policies.
    volume_var : str, optional
        Result column name for the working volume used to scale effluent
        ``F_vol`` (e.g. ``'curr_env'``).
    feed_volume_added_var : str, optional
        Result column name for the cumulative feed volume added during the run.
        When given (with ``volume_var``), the effluent ``F_vol`` is corrected for
        volume added by fed-batch feeding. Chemistry-agnostic; a fermentation
        subclass maps this to its own bookkeeping variable.
    aeration : AerationSpec, optional
    spike_retry : SpikeReduceRetry, optional
    pre_reactions : sequence of thermosteam.Reaction
        Applied to feed (and spike-feed) before kinetics.
    validators : sequence of callable
        ``validator(model) -> None``; each raises on failure.
    spike_feed_index : int, optional
        Index of the spike-feed inlet stream (excluded from the initial mix).
    N, V, T, P, Nmin, Nmax :
        Standard ``BatchBioreactor`` sizing parameters.
    """
    line = 'NSKBatchReactor'
    _ins_size_is_fixed = False
    autoselect_N = True
    # Class-level placeholder so `hasattr(NSKBatchReactor, 'simulate_kinetics')`
    # holds even before instantiation; `_init` overrides this per-instance with
    # the bound `_nsk_te_simulate_kinetics` method.
    simulate_kinetics = None

    def _init(self, *, nsk_kinetic_model,
              tau, tau_max,
              map_species_to_chemicals={},
              track_vars=(),
              n_simulation_steps=1000,
              tau_update_policy=None,
              n_decimal_places_for_tau_update_policy=2,
              volume_var=None,
              feed_volume_added_var=None,
              aeration=None, spike_retry=None,
              pre_reactions=(), validators=(),
              spike_feed_index=None,
              N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36):
        BatchBioreactor._init(self, tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax)
        self._load_components()

        nkm = nsk_kinetic_model
        if not isinstance(nkm, KineticModel):
            raise NotImplementedError(
                'NSKBatchReactor currently supports only KineticModel '
                f'instances; got {type(nkm).__name__}.')
        self.nsk_kinetic_model = nkm
        nkm.validate_units()

        self.map_species_to_chemicals = dict(map_species_to_chemicals)
        self.track_vars = list(track_vars)
        self.n_simulation_steps = n_simulation_steps
        self.tau_max = tau_max
        self.tau_update_policy = tau_update_policy
        self.n_decimal_places_for_tau_update_policy = n_decimal_places_for_tau_update_policy
        self.volume_var = volume_var
        self.feed_volume_added_var = feed_volume_added_var
        self.aeration = aeration
        self.spike_retry = spike_retry
        self.pre_reactions = list(pre_reactions)
        self.validators = list(validators)
        self.spike_feed_index = spike_feed_index

        self.run_type = 'simulate kinetics'
        # Full-trajectory results live on the kinetic model; the
        # reactor exposes them via delegating properties (see below).

        self.simulate_kinetics = self._nsk_te_simulate_kinetics
        self._validate_model_selections()

    # --- guard rails --------------------------------------------------------
    def _validate_model_selections(self):
        nkm = self.nsk_kinetic_model
        for selection in list(self.map_species_to_chemicals) + list(self.track_vars):
            try:
                nkm.get_value(selection)
            except Exception as e:
                raise KineticSimulationError(
                    f'Model selection {selection!r} is not a valid model '
                    f'variable/selection: {e}') from e

    # --- alias for back-compat ---------------------------------------------
    @property
    def map_chemicals_nsk_to_bst(self):
        return self.map_species_to_chemicals

    @map_chemicals_nsk_to_bst.setter
    def map_chemicals_nsk_to_bst(self, value):
        self.map_species_to_chemicals = dict(value)

    # --- full-trajectory results (delegated to the kinetic model) ---------
    # The kinetic model is the source of truth for full-trajectory
    # results; these read-only properties expose it under the reactor's
    # `nsk_results*` names. They return None until `nsk_kinetic_model`
    # is assigned during `_init`. The full-trajectory DataFrame is exposed as
    # `nsk_results_df` (not `results`) so it does not shadow biosteam's
    # `Unit.results()` design/cost method.
    @property
    def nsk_results_df(self):
        nkm = getattr(self, 'nsk_kinetic_model', None)
        return nkm.results_df if nkm is not None else None

    @property
    def nsk_results(self):
        nkm = getattr(self, 'nsk_kinetic_model', None)
        return nkm.results_array if nkm is not None else None

    @property
    def nsk_results_dict(self):
        nkm = getattr(self, 'nsk_kinetic_model', None)
        return nkm.results_dict if nkm is not None else None

    @property
    def nsk_results_col_names(self):
        nkm = getattr(self, 'nsk_kinetic_model', None)
        return nkm.results_col_names if nkm is not None else None

    # --- plotting -----------------------------------------------------------
    def plot_simulation_results(self, *args, show_fermentation_end=True,
                                show_aeration_end=True, **kwargs):
        """Plot the kinetic model's most recent simulation, marking
        the fermentation- and aeration-end times.

        Delegates to
        :meth:`nskinetics.KineticModel.plot_simulation_results`,
        injecting dashed vertical event lines at the fermentation-end time
        (``tau``) and, when aeration is configured, the aeration-end time
        (``tau_stop_aeration``). All other positional and keyword arguments are
        forwarded unchanged; keyword arguments override the injected lines when
        their names collide.

        Parameters
        ----------
        show_fermentation_end : bool
            Draw the fermentation-end (``tau``) line. Defaults to ``True``.
        show_aeration_end : bool
            Draw the aeration-end (``tau_stop_aeration``) line when aeration is
            configured. Defaults to ``True``.
        *args, **kwargs
            Forwarded to the kinetic model's
            :meth:`~nskinetics.KineticModel.plot_simulation_results`.
        """
        events = {}
        if show_fermentation_end:
            tau = getattr(self, 'tau', None)
            if tau is not None:
                events['fermentation_end'] = tau
        if show_aeration_end and self.aeration is not None:
            events['aeration_end'] = getattr(self, 'tau_stop_aeration', None)
        # User kwargs win on name collisions and may add further event lines.
        events.update(kwargs)
        return self.nsk_kinetic_model.plot_simulation_results(
            *args, **events)

    plot_time_course = plot_simulation_results
    plot_trajectory = plot_simulation_results

    # --- kinetic simulation -------------------------------------------------
    def _result_columns(self):
        """Assemble the ordered, de-duplicated result columns to request.

        Built from the configured features (working volume, spike-count retry,
        feed-volume tracking, aeration) plus ``track_vars`` and the
        species->chemical map. Chemistry-agnostic: no feature-specific column
        name is hardcoded here.
        """
        cols = []

        def add(name):
            if name is not None and name not in cols:
                cols.append(name)

        add('time')
        add(self.volume_var)
        if self.spike_retry is not None:
            add(self.spike_retry.count_var)
        add(self.feed_volume_added_var)
        if self.aeration is not None:
            for name in self.aeration.result_columns():
                add(name)
        for name in self.track_vars:
            add(name)
        for name in self.map_species_to_chemicals:
            add(name)
        return cols

    def _reset_and_simulate(self, feed, reset_spike_cap=False):
        nkm = self.nsk_kinetic_model
        nkm.reset(reset_spike_cap=reset_spike_cap)
        volume = getattr(feed, nkm.material_indexer.replace('imol', 'ivol')
                         .replace('imass', 'ivol'))['Water']
        indexer = getattr(feed, nkm.material_indexer)
        for c_nsk, c_bst in self.map_species_to_chemicals.items():
            nkm.set_concentration_from_stream(c_nsk, indexer[c_bst], volume)
        cols = self._result_columns()
        # Delegate the core kinetic integration to the kinetic model, which
        # stores the full-trajectory results; the reactor reads them back via
        # its nsk_results* / results properties.
        nkm.simulate(0, self.tau_max * nkm.time_factor,
                     self.n_simulation_steps, cols)

    def _nsk_te_simulate_kinetics(self, feed, tau):
        """Bound implementation of ``simulate_kinetics``: run the kinetic model
        on ``feed`` for reaction time ``tau`` and return the effluent stream."""
        nkm = self.nsk_kinetic_model
        # initial full simulation
        self._reset_and_simulate(feed, reset_spike_cap=True)
        # optional spike-count retry
        if self.spike_retry is not None:
            model = nkm._te
            sr = self.spike_retry
            # Seed the cap from the count reached by the full run, then let the
            # retry loop honor the SAME cap attribute it decrements.
            setattr(model, sr.max_count_var,
                    list(self.nsk_results_dict[sr.count_var])[-1] - 1)
            sr.run(
                model, lambda: self._reset_and_simulate(feed, reset_spike_cap=False))
        # validators
        for validate in self.validators:
            validate(nkm._te)
        # tau selection
        tau_index, ok = select_tau_index(
            self.nsk_results, self.nsk_results_col_names, tau,
            self.tau_update_policy, self.n_decimal_places_for_tau_update_policy)
        self._tau_update_success = ok
        self._load_specific_tau(tau_index)
        if self.aeration is not None:
            self.aeration.compute_cumulative_O2(self, tau_index)
        return self._build_effluent(feed)

    def _load_specific_tau(self, tau_index):
        cols = self.nsk_results_col_names
        row = self.nsk_results[tau_index]
        self.tau_index = tau_index
        self.nsk_results_specific_tau = row
        self.nsk_results_specific_tau_dict = {cols[i]: row[i] for i in range(len(cols))}
        self.tau = row[cols.index('time')]

    def _select_tau_index_from_saved(self, tau):
        tau_index, ok = select_tau_index(
            self.nsk_results, self.nsk_results_col_names, tau, None,
            self.n_decimal_places_for_tau_update_policy)
        self._load_specific_tau(tau_index)
        if self.aeration is not None:
            self.aeration.compute_cumulative_O2(self, tau_index)

    def _build_effluent(self, minimal_feed):
        effluent = minimal_feed.copy()
        nkm = self.nsk_kinetic_model
        cols = self.nsk_results_col_names
        row = self.nsk_results_specific_tau
        indexer = getattr(effluent, nkm.material_indexer)
        volume = effluent.ivol['Water']
        for c_nsk, c_bst in self.map_species_to_chemicals.items():
            indexer[c_bst] = row[cols.index(c_nsk)] * volume
        if self.volume_var is not None:
            d = self.nsk_results_specific_tau_dict
            curr_vol = d[self.volume_var]
            added = (d[self.feed_volume_added_var]
                     if self.feed_volume_added_var is not None else 0.0)
            effluent.F_vol *= curr_vol / (curr_vol - added)
        return effluent

    # --- run flow -----------------------------------------------------------
    def _run(self):
        vent, effluent = self.outs
        ins = self.ins
        vent.empty()
        effluent.empty()

        excluded = []
        spike_feed = None
        if self.spike_feed_index is not None:
            spike_feed = ins[self.spike_feed_index]
            excluded.append(spike_feed)
        if self.aeration is not None:
            self.compressed_air = ins[self.aeration.air_index]
            excluded.append(self.compressed_air)

        effluent.mix_from(i for i in ins if i not in excluded)

        for rxn in self.pre_reactions:
            rxn.force_reaction(effluent)
            if spike_feed is not None:
                rxn.force_reaction(spike_feed)

        keep = set(self.map_species_to_chemicals.values()) | {'Water'}
        minimal_feed = effluent.copy()
        for chem in minimal_feed.chemicals:
            if chem.ID not in keep:
                minimal_feed.imol[chem.ID] = 0.0

        run_type = self.run_type
        if run_type == 'simulate kinetics':
            minimal_effluent = self.simulate_kinetics(minimal_feed, self._tau)
        elif run_type == 'index saved nsk_results by tau':
            self._select_tau_index_from_saved(self.tau)
            minimal_effluent = self._build_effluent(minimal_feed)
        else:
            raise KineticSimulationError(f'Unknown run_type {run_type!r}.')

        for chem in minimal_effluent.chemicals:
            if chem.ID not in keep:
                minimal_effluent.imol[chem.ID] += effluent.imol[chem.ID]
                if spike_feed is not None:
                    minimal_effluent.imol[chem.ID] += spike_feed.imol[chem.ID]

        effluent.copy_like(minimal_effluent)
        self._finalize_effluent(effluent, vent, minimal_feed)
        if self.aeration is not None:
            self.aeration.set_air_stream(self, effluent, vent)

    def _finalize_effluent(self, effluent, vent, feed):
        """Subclass hook for chemistry-specific effluent finishing. Default:
        drop negative flows and route the vent."""
        effluent.empty_negative_flows()
        vent.receive_vent(effluent, energy_balance=False)

    def set_tolerances_kinetic_simulation(self, atol, rtol):
        """Set the kinetic integrator's absolute/relative tolerances.

        Parameters
        ----------
        atol, rtol : float
            Absolute and relative tolerances applied to the underlying
            RoadRunner integrator before simulation.
        """
        nkm = self.nsk_kinetic_model
        integrator = nkm._te.getIntegrator()
        integrator.absolute_tolerance = atol
        integrator.relative_tolerance = rtol
