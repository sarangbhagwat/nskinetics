# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
System factory for sugar-solution preparation and fed-batch fermentation, as
configured in the isobutanol biorefinery: a splitter feeding parallel
initial-feed and spike-feed conditioning trains (multi-effect evaporator,
pumps, dilution-water mixer, heat exchanger), a
:class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol`
fermentor, and a compressed-air (compressor + valve) aeration loop.
"""

import biosteam as bst

from ..units import FermentationSaccharomycesEthanolIsobutanol

__all__ = ('create_sugar_prep_and_fermentation_system',)

_DEFAULT_MAP_CHEMICALS_NSK_TO_BST = {
    '[s_glu]': 'Glucose',
    '[x]': 'Yeast',
    '[s_EtOH]': 'Ethanol',
    '[s_IBO]': 'Isobutanol',
    '[s_acetate]': 'AceticAcid',
}

_DEFAULT_TRACK_VARS = (
    'y_EtOH_glu_added',
    'y_EtOH_glu_consumed',
    'y_IBO_glu_added',
    'y_IBO_glu_consumed',
    'y_EtOH_IBO_glu_added',
    'curr_n_glu_spikes',
    'curr_a',
    'prod_EtOH',
    'curr_tot_vol_glu_feed_added',
    'curr_env',
)


@bst.SystemFactory(
    ID='sugar_prep_and_fermentation_sys',
    ins=[dict(ID='saccharified_slurry'),
         dict(ID='seed')],
    outs=[dict(ID='fermentation_vent'),
          dict(ID='fermentation_effluent'),
          dict(ID='initial_feed_condensate'),
          dict(ID='spike_feed_condensate')],
)
def create_sugar_prep_and_fermentation_system(
        ins, outs,
        nsk_kinetic_model=None,
        map_chemicals_nsk_to_bst=None,
        track_vars=None,
        tau=3 * 24,
        tau_max=3 * 24,
        tau_update_policy=('max', '[s_EtOH]'),
        n_decimal_places_for_tau_update_policy=0,
        try_fewer_n_spikes_until=None,
        perform_hydrolysis=False,
        sugar_IDs=('Glucose',),
        stage_1_max_x=5.0,
        stage_1_max_time=25.0,
        n_simulation_steps=1000,
        split=0.8,
        evaporator_pressures=(101325, 73581, 50892, 32777, 20000),
        evaporator_V=0.1,
        water_to_sugar_mol_ratio=100.,
        feed_T=32 + 273.15,
        air_pressure=3e7,
        compressor_eta=0.6,
        fermentor_kwargs=None,
    ):
    """
    Create the sugar-solution preparation and fermentation system of the
    isobutanol biorefinery: saccharified sugar slurry in, fermentation
    effluent out. Defaults reproduce the inline configuration of the
    isobutanol biorefinery's ``system.py`` exactly.

    The splitter/evaporator/mixer settings (``split``, ``evaporator_V``,
    ``water_to_sugar_mol_ratio``) are initial values, intended to be tuned
    after construction by a
    :class:`~nskinetics.units.FedBatchStrategySpecification`.

    Deliberately NOT built here (add caller-side after calling the factory):
    any fermentor specification coupling to upstream flowsheet streams (e.g.
    yeast/enzyme/ammonia feed-flow corrections), and the docking of the vent
    and effluent outlets to downstream units. Note that the inline system's
    fermentor specification ends by re-simulating the aeration loop
    (``V406.simulate(); K330.simulate(); V330.simulate()``); a caller re-adding
    feed-flow corrections must re-simulate ``K330`` and ``V330`` likewise, or
    the compressed-air streams go stale with respect to the fermentor.

    Three internal feed streams are created unconnected and therefore surface
    as system inlets beyond the declared ``ins``: ``atmospheric_air`` (feeding
    ``K330``) and one dilution-water stream per mixer. Faithful to the inline
    system, both mixers' water streams share the ID ``'dilution_water'``, so a
    flowsheet lookup by that ID resolves only to M302's stream; reach M301's
    through ``M301.ins[1]`` instead.

    Parameters
    ----------
    ins :
        * [0] Saccharified sugar slurry.
        * [1] Seed culture (docked to the fermentor's inlet 1).
    outs :
        * [0] Fermentation vent.
        * [1] Fermentation effluent.
        * [2] Initial-feed evaporator condensate.
        * [3] Spike-feed evaporator condensate.
    nsk_kinetic_model : KineticModel, optional
        Kinetic model driving the fermentor. Defaults to the shipped ``te_r``
        S. cerevisiae ethanol/isobutanol model (imported lazily — the model
        module configures and simulates itself on first import).
    map_chemicals_nsk_to_bst : dict, optional
        ``{model selection: biosteam chemical ID}`` for the fermentor.
    track_vars : list of str, optional
        Extra model selections recorded as result columns.
    tau, tau_max : float
        Reaction time and maximum simulated time [h].
    tau_update_policy : None or tuple
        Passed to the fermentor (see
        :func:`~nskinetics.units.select_tau_index`).
    n_decimal_places_for_tau_update_policy : int
        Rounding for the ``max``/``min``/``equals`` policies.
    try_fewer_n_spikes_until : callable, optional
        ``stop_when(model) -> bool`` for the spike-count retry loop. Defaults
        to stopping once residual glucose (rounded to 2 decimals) is zero.
    perform_hydrolysis : bool
        If ``True``, apply sucrose hydrolysis to feed and spike-feed.
    sugar_IDs : tuple of str
        Sugar chemical IDs; also drive the dilution-water mixers.
    stage_1_max_x, stage_1_max_time : float
        Stage-1 aerobic cutoffs (see :class:`~nskinetics.units.AerationSpec`).
    n_simulation_steps : int
        Number of kinetic integration output steps.
    split : float
        Initial S301 split fraction to the initial-feed train.
    evaporator_pressures : tuple of float
        Effect pressures [Pa] for both multi-effect evaporators.
    evaporator_V : float
        Initial molar vapor fraction for both evaporators.
    water_to_sugar_mol_ratio : float
        Initial dilution-water-to-sugar molar ratio for both mixers.
    feed_T : float
        Feed temperature [K] set by both heat exchangers.
    air_pressure : float
        Compressed-air pressure [Pa] (K330).
    compressor_eta : float
        Isentropic efficiency of the air compressor (K330).
    fermentor_kwargs : dict, optional
        Extra keyword arguments forwarded verbatim to
        :class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol`
        (e.g. ``stop_aeration_when_cell_density_plateaus``).
    """
    saccharified_slurry, seed = ins
    (fermentation_vent, fermentation_effluent,
     initial_feed_condensate, spike_feed_condensate) = outs

    if nsk_kinetic_model is None:
        # Heavy on purpose: importing the model module builds, configures,
        # and simulates te_r. Deferred to call time so that importing
        # nskinetics.processes stays cheap.
        from ..models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r as nsk_kinetic_model
    if map_chemicals_nsk_to_bst is None:
        map_chemicals_nsk_to_bst = dict(_DEFAULT_MAP_CHEMICALS_NSK_TO_BST)
    if track_vars is None:
        track_vars = list(_DEFAULT_TRACK_VARS)
    if try_fewer_n_spikes_until is None:
        try_fewer_n_spikes_until = lambda r_te: round(r_te.s_glu, 2) == 0.0

    ## Splitter
    S301 = bst.Splitter('S301', ins=saccharified_slurry,
                        outs=('fermentation_initial_feed', 'fermentation_spike'),
                        split=split, # initial value, tuned by FedBatchStrategySpecification
                        )

    ## Initial feed evaporator, pumps, mixer, and hx
    F301 = bst.MultiEffectEvaporator('F301', ins=S301-0, outs=('F301_l', 'F301_g'),
                                     P=evaporator_pressures, V=evaporator_V,
                                     flash=False)
    F301.V = evaporator_V # initial value, tuned by FedBatchStrategySpecification
    F301_design = F301._design
    F301_cost = F301._cost

    @F301.add_specification(run=False)
    def F301_spec():
        feed = F301.ins[0]
        if feed.F_mol:
            F301._run()
            F301._design = F301_design
            F301._cost = F301_cost
        else:
            F301.outs[1].empty()
            F301.outs[0].copy_like(feed)
            F301._design = lambda: 0
            F301._cost = lambda: 0

    F301_P0 = bst.units.Pump('F301_P0', ins=F301-0, outs='', P=101325.)
    F301_P1 = bst.units.Pump('F301_P1', ins=F301-1, outs=initial_feed_condensate,
                             P=101325.)

    M301 = bst.units.Mixer('M301', ins=(F301_P0-0, 'dilution_water'))
    M301.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio # initial value

    @M301.add_specification(run=False)
    def adjust_M301_water():
        M301_ins_1 = M301.ins[1]
        M301_ins_1.imol['Water'] = (M301.water_to_sugar_mol_ratio
                                    * M301.ins[0].imol[V406.sugar_IDs].sum())
        M301._run()

    H301 = bst.units.HXutility('H301', ins=M301-0, outs=('glucose_initial_feed',),
                               T=feed_T, rigorous=True)

    @H301.add_specification(run=False)
    def H301_spec():
        H301._run()
        H301.outs[0].phase = 'l'

    ## Spike evaporator, pumps, mixer, and hx
    F302 = bst.MultiEffectEvaporator('F302', ins=S301-1, outs=('F302_l', 'F302_g'),
                                     P=evaporator_pressures, V=evaporator_V,
                                     flash=False)
    F302.V = evaporator_V # initial value, tuned by FedBatchStrategySpecification
    F302_design = F302._design
    F302_cost = F302._cost

    @F302.add_specification(run=False)
    def F302_spec():
        feed = F302.ins[0]
        if feed.F_mol:
            F302._run()
            F302._design = F302_design
            F302._cost = F302_cost
        else:
            F302.outs[1].empty()
            F302.outs[0].copy_like(feed)
            F302._design = lambda: 0
            F302._cost = lambda: 0

    F302_P0 = bst.units.Pump('F302_P0', ins=F302-0, outs='', P=101325.)
    F302_P1 = bst.units.Pump('F302_P1', ins=F302-1, outs=spike_feed_condensate,
                             P=101325.)

    M302 = bst.units.Mixer('M302', ins=(F302_P0-0, 'dilution_water'))
    M302.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio # initial value

    @M302.add_specification(run=False)
    def adjust_M302_water():
        M302_ins_1 = M302.ins[1]
        M302_ins_1.imol['Water'] = (M302.water_to_sugar_mol_ratio
                                    * M302.ins[0].imol[V406.sugar_IDs].sum())
        M302._run()

    H302 = bst.units.HXutility('H302', ins=M302-0, outs=('glucose_spike_feed',),
                               T=feed_T, rigorous=True)

    @H302.add_specification(run=False)
    def H302_spec():
        H302._run()
        H302.outs[0].phase = 'l'

    ## Fermentor
    V406 = FermentationSaccharomycesEthanolIsobutanol(
        'V406',
        ins=(H301-0, seed, H302-0),
        outs=(fermentation_vent, fermentation_effluent),
        nsk_kinetic_model=nsk_kinetic_model,
        n_simulation_steps=n_simulation_steps,
        map_chemicals_nsk_to_bst=map_chemicals_nsk_to_bst,
        track_vars=track_vars,
        tau=tau,
        tau_max=tau_max,
        sugar_IDs=sugar_IDs,
        tau_update_policy=tau_update_policy,
        n_decimal_places_for_tau_update_policy=n_decimal_places_for_tau_update_policy,
        try_fewer_n_spikes_until=try_fewer_n_spikes_until,
        perform_hydrolysis=perform_hydrolysis,
        stage_1_max_x=stage_1_max_x,
        stage_1_max_time=stage_1_max_time,
        **(fermentor_kwargs or {}),
    )

    ## Compressed air system
    K330 = bst.units.IsothermalCompressor('K330', ins='atmospheric_air',
                                          outs=('pressurized_air'),
                                          P=air_pressure,
                                          eta=compressor_eta,
                                          driver='Electric motor',
                                          )

    @K330.add_specification(run=False)
    def K330_spec():
        K330_ins_0 = K330.ins[0]
        K330_ins_0.T = V406.T
        K330_ins_0.phase = 'g'
        K330_ins_0.mol[:] = K330.outs[0].mol[:]
        K330._run()

    V330 = bst.units.IsenthalpicValve('V330', ins=K330-0,
                                      P=101325.,
                                      vle=False,
                                      )
    V330.line = 'Valve'

    @V330.add_specification(run=False)
    def V330_spec():
        V330.ins[0].mol[:] = V330.outs[0].mol[:]
        V330._run()

    V330-0-3-V406
