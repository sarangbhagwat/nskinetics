# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
References
----------
.. [1] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden
    National. Renewable Energy Laboratory Golden, Colorado. P. Schoen,
    J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group
    Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics
    for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid
    Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical
    Report NREL/TP-5100-47764
"""

import numpy as np
from thermosteam import Reaction

from .batch_reactor import NSKBatchReactor, AerationSpec, SpikeReduceRetry

__all__ = ('NSKFermentation',)


def _negative_concentration_validator(model):
    from ..exceptions import MassBalanceError
    if np.any(np.round([model.x, model.s_glu, model.s_EtOH,
                        model.s_IBO, model.s_acetate], 2) < 0):
        raise MassBalanceError(
            'Negative concentrations in final kinetic simulation.')


def _yield_over_theoretical_validator(model):
    from ..exceptions import MassBalanceError
    if np.any([model.y_EtOH_glu_added > 0.511, model.y_IBO_glu_added > 0.410]):
        raise MassBalanceError(
            'Yield over theoretical maximum in final kinetic simulation.')


class NSKFermentation(NSKBatchReactor):
    """
    S. cerevisiae fed-batch fermentation for 1st-generation ethanol (and
    isobutanol) production. Thin :class:`NSKBatchReactor` subclass that
    configures aeration, glucose fed-batch spike-count retry, sucrose
    hydrolysis, and yield/mass-balance validators.

    Parameters
    ----------
    ins :
        Inlet fluids to be mixed into the fermentor. The first inlets are:
        [0] initial feed, [1] seed culture, [2] spike feed (fed-batch),
        [3] compressed air (aeration).
    outs :
        * [0] Vent
        * [1] Effluent
    tau : float
        Reaction time [h].
    kinetic_reaction_system : TelluriumReactionSystem
        Kinetic model driving the fermentation.
    map_chemicals_nsk_to_bst : dict
        ``{model var: biosteam chemical ID}``. Preferred alias:
        ``map_species_to_chemicals``.
    perform_hydrolysis : bool
        If ``True``, apply ``Sucrose + Water -> 2 Glucose`` to feed and spike-feed.
    try_fewer_n_spikes_until : callable
        ``stop_when(model) -> bool`` for the spike-count retry loop.
    aeration_safety_factor, stage_1_max_time, stage_1_max_x, stop_aeration_when_cell_density_plateaus, factor_for_cell_density_plateau :
        Aeration configuration (see :class:`AerationSpec`).
    sugar_IDs : tuple
        Sugar chemical IDs (retained for back-compat).

    Notes
    -----
    Either ``N`` or ``V`` must be given.
    """
    line = 'NSKFermentation'

    def _init(self, tau, kinetic_reaction_system,
              map_chemicals_nsk_to_bst=None,
              map_species_to_chemicals=None,
              track_vars=None,
              n_simulation_steps=1000,
              N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36,
              sugar_IDs=('Sucrose', 'Glucose', 'Xylose'),
              tau_max=24. * 7.,
              tau_update_policy=None,
              n_decimal_places_for_tau_update_policy=2,
              try_fewer_n_spikes_until=lambda r: True,
              perform_hydrolysis=True,
              aeration_safety_factor=2.0,
              stage_1_max_time=np.inf,
              stage_1_max_x=np.inf,
              stop_aeration_when_cell_density_plateaus=False,
              factor_for_cell_density_plateau=0.5):

        # Resolve the species->chemical map (new name preferred, old name accepted).
        if map_species_to_chemicals is None:
            map_species_to_chemicals = map_chemicals_nsk_to_bst or {}
        if track_vars is None:
            track_vars = []

        aeration = AerationSpec(
            qO2_var='qO2', is_aerobic_var='is_aerobic', biomass_var='[x]',
            volume_var='curr_env', biomass_chemical='Yeast',
            safety_factor=aeration_safety_factor, air_index=3,
            stop_when_cell_density_plateaus=stop_aeration_when_cell_density_plateaus,
            factor_for_cell_density_plateau=factor_for_cell_density_plateau,
            stage_1_max_time=stage_1_max_time, stage_1_max_x=stage_1_max_x)

        spike_retry = SpikeReduceRetry(
            max_count_var='max_n_glu_spikes', count_var='curr_n_glu_spikes',
            stop_when=try_fewer_n_spikes_until)

        NSKBatchReactor._init(
            self, kinetic_reaction_system=kinetic_reaction_system,
            tau=tau, tau_max=tau_max,
            map_species_to_chemicals=map_species_to_chemicals,
            track_vars=track_vars, n_simulation_steps=n_simulation_steps,
            tau_update_policy=tau_update_policy,
            n_decimal_places_for_tau_update_policy=n_decimal_places_for_tau_update_policy,
            volume_var='curr_env',
            feed_volume_added_var='curr_tot_vol_glu_feed_added',
            aeration=aeration, spike_retry=spike_retry,
            pre_reactions=(), validators=(
                _negative_concentration_validator,
                _yield_over_theoretical_validator),
            spike_feed_index=2,
            N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax)

        chemicals = self.chemicals
        self.perform_hydrolysis = perform_hydrolysis
        if perform_hydrolysis:
            self.hydrolysis_reaction = Reaction(
                'Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
            self.pre_reactions = [self.hydrolysis_reaction]

        self.sugar_IDs = sugar_IDs
        # These assignments go through property setters (below) that mirror the
        # value onto self.aeration, so post-construction mutation takes effect at
        # simulate time instead of silently no-op'ing.
        self.stop_aeration_when_cell_density_plateaus = \
            stop_aeration_when_cell_density_plateaus
        self.factor_for_cell_density_plateau = factor_for_cell_density_plateau
        self.aeration_safety_factor = aeration_safety_factor
        # mirror stage-1 cutoffs onto the model (property setters below)
        self.stage_1_max_time = stage_1_max_time
        self.stage_1_max_x = stage_1_max_x

    # --- aeration config mirrors onto the AerationSpec ----------------------
    @property
    def aeration_safety_factor(self):
        return self._aeration_safety_factor

    @aeration_safety_factor.setter
    def aeration_safety_factor(self, val):
        self._aeration_safety_factor = val
        if self.aeration is not None:
            self.aeration.safety_factor = val

    @property
    def stop_aeration_when_cell_density_plateaus(self):
        return self._stop_aeration_when_cell_density_plateaus

    @stop_aeration_when_cell_density_plateaus.setter
    def stop_aeration_when_cell_density_plateaus(self, val):
        self._stop_aeration_when_cell_density_plateaus = val
        if self.aeration is not None:
            self.aeration.stop_when_cell_density_plateaus = val

    @property
    def factor_for_cell_density_plateau(self):
        return self._factor_for_cell_density_plateau

    @factor_for_cell_density_plateau.setter
    def factor_for_cell_density_plateau(self, val):
        self._factor_for_cell_density_plateau = val
        if self.aeration is not None:
            self.aeration.factor_for_cell_density_plateau = val

    # --- stage-1 cutoffs mirror onto the kinetic model ---------------------
    @property
    def stage_1_max_time(self):
        return self._stage_1_max_time

    @stage_1_max_time.setter
    def stage_1_max_time(self, val):
        self._stage_1_max_time = val
        self.kinetic_reaction_system._te.stage_1_max_time = val
        if self.aeration is not None:
            self.aeration.stage_1_max_time = val

    @property
    def stage_1_max_x(self):
        return self._stage_1_max_x

    @stage_1_max_x.setter
    def stage_1_max_x(self, val):
        self._stage_1_max_x = val
        self.kinetic_reaction_system._te.stage_1_max_x = val
        if self.aeration is not None:
            self.aeration.stage_1_max_x = val

    # --- isobutanol-specific effluent finishing ----------------------------
    def _finalize_effluent(self, effluent, vent, feed):
        effluent.imol['NH3'] = 0.  # NH3 in ins must be based on final Yeast mass
        effluent.imol['CO2'] = (sum(i.get_atomic_flow('C') for i in self.ins)
                                - sum(i.get_atomic_flow('C') for i in self.outs))
        effluent.empty_negative_flows()
        vent.receive_vent(effluent, energy_balance=False)
        effluent.imol['Ethanol'] += max(0.0, vent.imol['Ethanol'])
        vent.imol['Ethanol'] = 0.0
