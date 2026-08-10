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

The factory's dedicated chemical set is built by
:func:`create_sugar_prep_and_fermentation_chemicals`.
"""

import biosteam as bst
import thermosteam as tmo
from thermosteam import functional as fn

from ..units import (FermentationSaccharomycesEthanolIsobutanol,
                     FedBatchStrategySpecification,
                     SpikeControlVariables,
                     ConcentrationActuator)

__all__ = ('create_sugar_prep_and_fermentation_system',
           'create_sugar_prep_and_fermentation_chemicals')


def create_sugar_prep_and_fermentation_chemicals():
    """
    Create the chemical set for
    :func:`create_sugar_prep_and_fermentation_system` — this factory's own,
    dedicated set, activated by the factory's ``set_thermo=True`` option or
    usable directly via ``bst.settings.set_thermo(...)``. It is *not* a
    package-wide set: future factories in :mod:`nskinetics.processes` ship
    their own ``create_<factory>_chemicals`` beside them.

    Contains exactly what this factory's units touch: ``Water``, ``Ethanol``,
    ``Isobutanol``, ``AceticAcid``, ``Glucose``/``Sucrose`` (liquid-phase
    solutes; ``Sucrose`` covers ``perform_hydrolysis=True``), ``CO2``/``O2``/
    ``N2`` (aeration and vent), ``NH3`` (zeroed in the fermentor's effluent),
    and a ``Yeast`` pseudochemical replicating the cane/corn-biorefinery
    definition this package is validated against (solid, ``CH1.61O0.56``,
    rho = 1540 kg/m^3, glucose-tied Cp and Hf, synonym ``DryYeast``).

    Returns
    -------
    thermosteam.Chemicals
        Compiled chemical set with synonyms ``H2O`` (Water) and ``DryYeast``
        (Yeast).
    """
    chemicals = tmo.Chemicals([
        'Water', 'Ethanol', 'Isobutanol', 'AceticAcid',
        tmo.Chemical('Glucose', phase='l'),
        tmo.Chemical('Sucrose', phase='l'),
        tmo.Chemical('CO2', phase='g'),
        tmo.Chemical('O2', phase='g'),
        tmo.Chemical('N2', phase='g'),
        'NH3',
    ])
    Glucose = chemicals.Glucose
    Sucrose = chemicals.Sucrose
    Glucose.N_solutes = 1
    Sucrose.N_solutes = 2
    Yeast = tmo.Chemical('Yeast', phase='s', phase_ref='s', search_db=False,
                         formula='CH1.61O0.56', rho=1540,
                         Cp=Glucose.Cp(298.15), default=True)
    # Hf tied to glucose per unit mass, ignoring heats related to growth
    # (same convention as the cane/corn biorefinery chemical sets).
    Yeast.Hf = Glucose.Hf / Glucose.MW * Yeast.MW
    chemicals.append(Yeast)
    Yeast.V.add_model(fn.rho_to_V(rho=1540, MW=Yeast.MW), top_priority=True)
    for sugar in (Glucose, Sucrose):
        # Dissolved solutes occupy negligible volume.
        sugar.V.add_model(fn.rho_to_V(rho=1e5, MW=sugar.MW), top_priority=True)
    chemicals.NH3.at_state('l')
    for chemical in chemicals:
        chemical.default()
    chemicals.compile()
    chemicals.set_synonym('Water', 'H2O')
    chemicals.set_synonym('Yeast', 'DryYeast')
    return chemicals


class _SetThermoSystemFactory(bst.SystemFactory):
    """A :class:`biosteam.SystemFactory` accepting an opt-in ``set_thermo``
    keyword at call time. When ``True``, the factory's own chemical set
    (:func:`create_sugar_prep_and_fermentation_chemicals`) is activated via
    ``bst.settings.set_thermo`` *before* the system's inlet/outlet streams
    are created — a plain factory-function keyword would come too late, as
    ``SystemFactory.__call__`` creates the declared streams before the
    wrapped function body runs. Defaults to ``False`` (no thermo side
    effects; any existing thermo is used unchanged).

    When the factory returns a real :class:`biosteam.System` (i.e. not
    ``mockup=True``), a system specification is attached that imposes the
    fed-batch strategy (``fbs_spec.load_specifications()``) before each
    ``system.simulate()``; see
    :func:`_attach_fed_batch_strategy_system_specification`."""

    def __call__(self, *args, set_thermo=False, **kwargs):
        if set_thermo:
            bst.settings.set_thermo(
                create_sugar_prep_and_fermentation_chemicals())
        result = super().__call__(*args, **kwargs)
        # With udct=True the factory returns (system, unit_dct); with
        # mockup=True the "system" is a MockSystem, which cannot carry
        # specifications (callers absorb its units into their own system).
        system = result[0] if isinstance(result, tuple) else result
        if isinstance(system, bst.System):
            _attach_fed_batch_strategy_system_specification(system)
        return result


def _attach_fed_batch_strategy_system_specification(system):
    """Attach a system specification to ``system`` that runs
    ``fbs_spec.load_specifications`` before every ``system.simulate()``,
    re-imposing the fed-batch strategy (kinetic-model concentrations,
    evaporator/dilution actuators, and splitter ratio) so the simulated
    flowsheet always matches the spec's current values.

    The specification function accepts the optional keyword arguments of
    :meth:`~nskinetics.units.FedBatchStrategySpecification.load_specifications`
    (``target_conc``, ``threshold_conc``, ``spike_conc``, ``tau_max``) for
    direct calls — e.g. ``system.specifications[0].f(target_conc=240.)`` —
    mirroring the isobutanol biorefinery's model specification; when invoked
    by ``system.simulate()`` it runs with the spec's stored values.
    """
    fermentor, = [u for u in system.units
                  if isinstance(u, FermentationSaccharomycesEthanolIsobutanol)]

    def fed_batch_strategy_specification(**kwargs):
        fermentor.fbs_spec.load_specifications(**kwargs)

    system.add_specification(fed_batch_strategy_specification, simulate=True)


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

# The isobutanol biorefinery's baseline_specifications dict (scenario baselines
# differ from the constructor's initial values on purpose; see its system.py).
_DEFAULT_BASELINE_SPECIFICATIONS = {
    'target_conc': 221.25,
    'threshold_conc': 217.125,
    'spike_conc': 600.0,
    'tau_max': 120.0,
    'max_n_spikes': 16,
}


@_SetThermoSystemFactory(
    ID='sugar_prep_and_fermentation_sys',
    # Default compositions: the isobutanol biorefinery's scenario-A baseline
    # values of these two streams (measured 2026-08-09 at MPSP 0.7451371545),
    # restricted to the chemicals the factory's shipped chemical set can
    # represent (the slurry's corn-mash solids — Fiber, proteins, TriOlein,
    # Ash, H2SO4, CaO — are dropped; they drive neither the sugar trains nor
    # the kinetics). Applied only when the factory creates the streams itself;
    # caller-passed streams are used as-is. SystemFactory drops flow keys for
    # chemicals absent from the active thermo, so these are safe under any
    # caller thermo.
    ins=[dict(ID='saccharified_slurry', units='kg/hr', T=305.35,
              Water=162685.7, Glucose=38624.7, NH3=209.1),
         dict(ID='seed', units='kg/hr',
              Water=11.41, Yeast=0.548)],
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
        target_conc=220.0,
        threshold_conc=210.0,
        spike_conc=600.0,
        max_n_spikes=None,
        spike_control_variables=None,
        baseline_specifications=None,
        fbs_spec_kwargs=None,
    ):
    """
    Create the sugar-solution preparation and fermentation system of the
    isobutanol biorefinery: saccharified sugar slurry in, fermentation
    effluent out. Defaults match the configuration the isobutanol
    biorefinery uses — historically written inline in its ``system.py``,
    which now builds this section by calling this factory.

    The splitter/evaporator/mixer settings (``split``, ``evaporator_V``,
    ``water_to_sugar_mol_ratio``) are initial values, tuned by the
    :class:`~nskinetics.units.FedBatchStrategySpecification` this factory
    builds from its own units and attaches to the fermentor
    (``V406.fed_batch_strategy_specification``, short alias
    ``V406.fbs_spec``). Building is side-effect-free; the strategy is
    imposed by the spec's ``load_specifications``. When the factory returns
    a real :class:`biosteam.System` (not ``mockup=True``), that call is also
    attached as a system specification, so ``system.simulate()`` first
    re-imposes the strategy at the spec's current values (see
    :func:`_attach_fed_batch_strategy_system_specification`); with
    ``mockup=True`` the caller keeps full control and must invoke
    ``load_specifications`` itself.

    Deliberately NOT built here (add caller-side after calling the factory):
    any fermentor specification coupling to upstream flowsheet streams (e.g.
    yeast/enzyme/ammonia feed-flow corrections), and the docking of the vent
    and effluent outlets to downstream units. Aeration-loop convergence is
    built into the fermentor itself for ALL builds (including
    ``mockup=True``): via :class:`~nskinetics.units.NSKBatchReactor`'s
    ``converge_air_supply`` behavior, every ``V406`` run ends by
    re-simulating ``V330`` then ``K330`` at the freshly written air demand
    (the aeration model writes the demand into ``V330``'s outlet, ``V330``'s
    specification pulls it onto its inlet — which is ``K330``'s outlet — and
    only then can ``K330`` pull the fresh demand onto its own inlet), so
    caller-side specifications must NOT re-simulate the aeration loop
    themselves. Avoid adding more than one self-simulating specification to
    the fermentor (each specification's nested ``simulate`` runs the other,
    doubling the kinetic simulation per fermentor run). Any feed set from
    the fermentor's effluent (e.g. an ammonia-per-yeast correction) must be
    computed *after* ``V406.simulate()``, not before. The
    biorefinery also attaches ``fbs_spec.product_stream`` and
    ``fbs_spec.n_tea_solves`` to the spec after construction (plain
    attributes read at MPSP-evaluation time), so a caller driving a TEA
    from this system must add them caller-side likewise.

    Three internal feed streams are created unconnected and therefore surface
    as system inlets beyond the declared ``ins``: ``atmospheric_air`` (feeding
    ``K330``) and one dilution-water stream per mixer. Faithful to the
    biorefinery's original inline construction, both mixers' water streams
    share the ID ``'dilution_water'``, so a
    flowsheet lookup by that ID resolves only to M302's stream; reach M301's
    through ``M301.ins[1]`` instead.

    Parameters
    ----------
    ins :
        * [0] Saccharified sugar slurry. When created by the factory,
          defaults to the isobutanol biorefinery's scenario-A baseline
          composition (162,685.7 kg/hr Water, 38,624.7 kg/hr Glucose,
          209.1 kg/hr NH3 at 305.35 K — the corn-mash solids the shipped
          chemical set cannot represent are dropped).
        * [1] Seed culture (docked to the fermentor's inlet 1). When created
          by the factory, defaults to the scenario-A baseline composition
          (11.41 kg/hr Water, 0.548 kg/hr Yeast).
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
    target_conc, threshold_conc, spike_conc : float
        Initial-feed target, spike-trigger threshold, and spike-feed
        concentrations of the attached
        :class:`~nskinetics.units.FedBatchStrategySpecification`. Its
        ``tau_max`` reuses this factory's ``tau_max`` parameter.
    max_n_spikes : float, optional
        Cap on the number of fed-batch spikes, held by the same
        specification. ``None`` (default) leaves the kinetic model's own cap
        untouched; the default ``baseline_specifications`` still carry the
        isobutanol biorefinery's baseline cap of 16.
    spike_control_variables : SpikeControlVariables, optional
        Kinetic-model variable names through which the strategy is imposed.
        Defaults to the shipped ibo model's names (``conc_glu_feed_spike``,
        ``target_conc_glu_spike``, ``threshold_conc_glu_spike``).
    baseline_specifications : dict, optional
        Baseline values of the five specifications, keyed by
        ``'target_conc'``, ``'threshold_conc'``, ``'spike_conc'``,
        ``'tau_max'``, ``'max_n_spikes'``. Defaults to a per-call copy of
        the isobutanol biorefinery's baseline dict.
    fbs_spec_kwargs : dict, optional
        Overrides merged over the factory-computed
        :class:`~nskinetics.units.FedBatchStrategySpecification`
        constructor kwargs (last-write-wins). Overrides that need
        references to factory-internal units (e.g. actuator bounds) are
        easier applied after construction on the spec's plain attributes,
        e.g. ``V406.fbs_spec.feed_concentrator.ub = 0.9``.
    set_thermo : bool
        Call-time keyword (handled by the factory object itself, before the
        system's streams are created — it is therefore not part of this
        function's signature). If ``True``, activate this factory's own
        chemical set, :func:`create_sugar_prep_and_fermentation_chemicals`,
        via ``bst.settings.set_thermo`` — overwriting any thermo currently
        set. Defaults to ``False``: the caller's pre-set thermo (e.g. the
        isobutanol biorefinery's corn chemicals) is used unchanged.
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

    ## Fed-batch strategy specification (built, NOT imposed here: it is
    ## imposed by load_specifications, which simulates upstream units —
    ## called by the system specification attached post-build for real
    ## System returns, or caller-side for mockup=True)
    if spike_control_variables is None:
        spike_control_variables = SpikeControlVariables(
            spike_conc_var='conc_glu_feed_spike',
            target_conc_var='target_conc_glu_spike',
            threshold_conc_var='threshold_conc_glu_spike',
            )
    if baseline_specifications is None:
        baseline_specifications = dict(_DEFAULT_BASELINE_SPECIFICATIONS)

    fbs_kwargs = dict(
        target_conc=target_conc,
        threshold_conc=threshold_conc,
        spike_conc=spike_conc,
        max_n_spikes=max_n_spikes,
        tau_max=tau_max,
        fermentation_reactor=V406,
        splitter=S301,
        control_variables=spike_control_variables,
        feed_concentrator=ConcentrationActuator(F301, 'V', 0.0, 0.8),
        feed_diluter=ConcentrationActuator(
            M301, 'water_to_sugar_mol_ratio', 0.0, 100_000),
        spike_concentrator=ConcentrationActuator(F302, 'V', 0.0, 0.8),
        spike_diluter=ConcentrationActuator(
            M302, 'water_to_sugar_mol_ratio', 0.0, 100_000),
        feed_units_sequential=[F301, F301_P0, F301_P1, M301, H301],
        spike_units_sequential=[F302, F302_P0, F302_P1, M302, H302],
        species_IDs=list(sugar_IDs),
        solvent_ID='Water',
        baseline_specifications=baseline_specifications,
        )
    fbs_kwargs.update(fbs_spec_kwargs or {})
    V406.fed_batch_strategy_specification = V406.fbs_spec = \
        FedBatchStrategySpecification(**fbs_kwargs)
