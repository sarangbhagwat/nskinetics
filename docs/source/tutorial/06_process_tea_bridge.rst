Coupling to a process and TEA
================================

:doc:`05_real_model` ran a kinetic model on its own — species concentrations
in, species concentrations out, no notion of feed streams, equipment sizing,
or cost. :class:`~nskinetics.units.NSKBatchReactor` is the bridge: a
biosteam ``BatchBioreactor`` subclass that owns a
:class:`~nskinetics.TelluriumReactionSystem`, runs it to pick a reaction
time, and turns the result into a biosteam effluent stream — so the same
kinetic model that produced the numbers in the last two pages can also drive
equipment design, installed cost, and (through biosteam's TEA machinery)
minimum product selling price.

A generic bridge model
------------------------

This section builds and simulates a small, self-contained
``NSKBatchReactor`` on a *generic* two-species kinetic model (substrate
``S`` consumed to product ``P``) with one :class:`~nskinetics.FeedSpike` —
deliberately not the isobutanol model from :doc:`05_real_model`, to keep the
bridge mechanics separate from any one chemistry. It mirrors the package's
own ``nskinetics.examples.generic_batch_reactor_demo`` module inline, so it
is fully runnable on its own:

.. code-block:: python

   import warnings
   import tellurium as te
   import biosteam as bst
   import nskinetics as nsk

   _MODEL = """
   model generic_batch()
     compartment env; species S in env, P in env;
     S = 50; P = 0; env = 1; k = 0.5;
     threshold = 5; target = 50; feed_conc = 500;
     n_spk = 0; n_max = 3; last_vol = 0; tot_vol = 0;
     n_spk has dimensionless; n_max has dimensionless;
     J: S => P; k*S*env;
     curr_env := env;
   end
   """

   r = te.loadAntimonyModel(_MODEL)
   te_r = nsk.TelluriumReactionSystem(r, units={'time': 'h', 'conc': 'g/L'})

   spike = nsk.FeedSpike(
       species='S', when='S <= threshold', target='target',
       feed_conc='feed_conc', volume_var='env', max_count='n_max',
       count_var='n_spk', last_vol_var='last_vol', tot_vol_var='tot_vol',
       priority=5, name='S_spike')
   for event in spike.expand():
       te_r.add_event(event)
   te_r.compile_events()

   bst.settings.set_thermo(['Water', 'Glucose', 'Ethanol'], cache=True)
   feed = bst.Stream('feed', Water=1000., Glucose=50., units='kg/hr')
   seed = bst.Stream('seed', Water=10., units='kg/hr')
   spike_feed = bst.Stream('spike_feed', Glucose=500., units='kg/hr')

   reactor = nsk.units.NSKBatchReactor(
       'R_demo', ins=(feed, seed, spike_feed), outs=('vent', 'effluent'),
       kinetic_reaction_system=te_r,
       map_species_to_chemicals={'S': 'Glucose', 'P': 'Ethanol'},
       track_vars=['[S]', '[P]'],
       tau=48., tau_max=120., volume_var='curr_env',
       feed_volume_added_var='tot_vol',
       spike_feed_index=2, V=100.)

The kinetic model itself is the same shape as :doc:`04_fed_batch_feeding`'s
toy example — one rate law (``J: S => P``) plus a ``FeedSpike`` on ``S`` — it
is the reactor construction around it that is new. Note ``curr_env := env``:
an assignment rule that mirrors the compartment volume ``env`` into its own
variable, so it can be handed to the reactor as a plain result column (see
``volume_var`` below).

The bridge-specific keyword arguments
----------------------------------------

- ``map_species_to_chemicals`` — ``{model variable: biosteam chemical ID}``.
  The kinetic model is simulated once over ``[0, tau_max]``, a single row is
  then picked out at the selected ``tau`` (via ``select_tau_index``, using
  ``tau_update_policy`` if one is given), and *only that one row's* value of
  ``S`` and ``P`` is converted into the matching chemical's mass in the
  effluent stream (``Glucose`` and ``Ethanol`` here) — not every simulated
  timestep. Any model variable or bracketed ``'[conc]'`` selection can be
  mapped this way.
- ``tau`` / ``tau_max`` — ``tau`` is the reaction time requested (48 h here);
  ``tau_max`` (120 h) is how far the kinetic model is actually integrated,
  so that a ``tau_update_policy`` (not used in this minimal example) could
  pick a different, data-driven tau from within the simulated window without
  re-running the kinetics. With no policy given, the reactor simply reads
  the kinetic trajectory at ``tau``.
- ``volume_var`` / ``feed_volume_added_var`` — together they correct the
  effluent's volumetric flow for fed-batch feeding: ``volume_var``
  (``'curr_env'``) is the result column holding the working volume, and
  ``feed_volume_added_var`` (``'tot_vol'``) is the cumulative volume added
  by ``FeedSpike`` events during the run. Without this pair the reactor
  would size the effluent off the *initial* mixed volume only, ignoring any
  fed-batch spikes.
- ``spike_feed_index`` — the index into ``ins`` (here, ``2`` →
  ``spike_feed``) of the stream that supplies spike-feed composition/cost,
  excluded from the initial reactor charge so it is not double-counted.
- ``track_vars`` — extra model selections to record as result columns beyond
  those the reactor needs internally. The ``map_species_to_chemicals`` keys
  (``'S'``, ``'P'``) record species *amounts* (what the effluent mass is built
  from); adding the bracketed ``'[S]'`` / ``'[P]'`` selections here also records
  the corresponding *concentrations*, so the plot below can show concentration
  rather than amount. The distinction matters in a fed-batch run: each
  ``FeedSpike`` enlarges the working volume, so an amount trace keeps climbing
  even though the spike resets the *concentration* back to the same target.

Simulating and reading results
----------------------------------

Standalone ``simulate()`` runs biosteam's design/cost logic over several
internal passes without clearing results between them; the kinetic ``_run``
adds no costs itself, so the only effect is a benign, expected
``RuntimeWarning`` on later passes, silenced here:

.. code-block:: python

   with warnings.catch_warnings():
       warnings.filterwarnings(
           'ignore', message=r'.*method added unit results.*',
           category=RuntimeWarning)
       reactor.simulate()

   reactor.show(N=100)
   print(reactor.results())

Running the full snippet above end-to-end in the ``HP_2024`` environment
prints:

.. code-block:: text

   NSKBatchReactor: R_demo
   ins...
   [0] feed
       phase: 'l', T: 298.15 K, P: 101325 Pa
       flow (kmol/hr): Water    55.5
                       Glucose  0.278
   [1] seed
       phase: 'l', T: 298.15 K, P: 101325 Pa
       flow (kmol/hr): Water  0.555
   [2] spike_feed
       phase: 'l', T: 298.15 K, P: 101325 Pa
       flow (kmol/hr): Glucose  2.78
   outs...
   [0] vent
       phase: 'g', T: 305.15 K, P: 101325 Pa
       flow: 0
   [1] effluent
       phase: 'l', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water    74.6
                       Glucose  1.81e-08
                       Ethanol  6.29
   NSKBatchReactor                                    Units    R_demo
   Electricity         Power                             kW      5.47
                       Cost                          USD/hr     0.428
   Chilled water       Duty                           kJ/hr -9.48e+06
                       Flow                         kmol/hr  6.35e+03
                       Cost                          USD/hr      47.4
   Design              Reactor volume                    m3      97.7
                       Batch time                        hr       102
                       Loading time                      hr        51
                       Number of reactors                           2
                       Recirculation flow rate        m3/hr     0.861
                       Reactor duty                   kJ/hr -3.16e+06
                       Cleaning and unloading time       hr         3
                       Working volume fraction                    0.9
   Purchase cost       Heat exchangers (x2)             USD  1.38e+04
                       Reactors (x2)                    USD  2.95e+05
                       Agitators (x2)                   USD  1.83e+04
                       Cleaning in place                USD   5.1e+04
                       Recirculation pumps (x2)         USD  2.81e+03
   Total purchase cost                                  USD  3.81e+05
   Utility cost                                      USD/hr      47.8

``reactor.show(N=100)`` prints the unit with its inlet and outlet streams —
the three feeds (initial ``feed``, ``seed``, and fed-batch ``spike_feed``)
and the two outlets (``vent`` and ``effluent``). The ``effluent`` carries the
kinetic result mapped back to chemicals: ~6.29 kmol/hr of ``Ethanol`` (product
``P``) and essentially no ``Glucose`` (substrate ``S`` fully consumed).
``N=100`` raises the per-stream cap on how many chemicals are listed, so every
flow is shown.

``reactor.results()`` is biosteam's standard ``Unit.results()`` design-and-cost
table: reactor volume, batch/loading/cleaning times, number of reactors, the
purchase-cost breakdown (total ~$381,000), and utility duties and costs.
``NSKBatchReactor`` exposes its own full kinetic trajectory as ``nsk_results_df``
precisely so it does not shadow this method. biosteam's
``reactor.installed_cost`` applies installation factors on top of the purchase
cost shown here — the figure a downstream TEA carries forward into capital cost
and minimum product selling price.

The kinetic detail the earlier pages focused on is still one attribute away:
``reactor.tau`` (48 h) for the selected reaction time,
``reactor.nsk_results_specific_tau_dict`` for the model state at ``tau``
(substrate ``[S]`` ~1.82e-06 g/L, product ``[P]`` ~161 g/L), and
``te_r.get_value('n_spk')`` for the fed-batch spike count (3, its ``n_max``
cap over the 48 h window). ``FeedSpike`` uses the same
threshold/target/feed_conc mechanics from :doc:`04_fed_batch_feeding`, just
driving a real biosteam effluent this time.

The same run can be visualized straight off the reactor, which marks the
fermentation-end (``tau``) time on the kinetic trajectory:

.. code-block:: python

   reactor.plot_simulation_results(variables=['[S]', '[P]'], labels=['S', 'P'])

.. figure:: /_static/images/examples/tutorial_06_process_tea_bridge.png
   :width: 400

   Substrate ``S`` is consumed to product ``P`` (concentrations, ``[S]`` and
   ``[P]``); each ``FeedSpike`` resets ``[S]`` to the same target, and the
   dashed line marks the selected fermentation-end time ``tau`` (48 h). When an
   :class:`~nskinetics.units.AerationSpec` is configured (as in
   :class:`~nskinetics.units.NSKFermentation`), the aeration-end time is
   marked too.

From here to a real biorefinery: ``NSKFermentation``
--------------------------------------------------------

The generic reactor above demonstrates the bridge mechanics; a real
fermentation needs more. :class:`~nskinetics.units.NSKFermentation` is
:class:`~nskinetics.units.NSKBatchReactor`'s fed-batch fermentation
subclass, adding: aeration bookkeeping through
:class:`~nskinetics.units.AerationSpec` (tracking cumulative oxygen uptake
and the aerobic/anaerobic stage switch seen in :doc:`05_real_model`,
including an option to stop aeration once cell density plateaus); optional
sucrose hydrolysis (``Sucrose + Water -> 2 Glucose``) applied to the feed and
spike-feed before kinetics; a fed-batch spike-count retry
(:class:`~nskinetics.units.SpikeReduceRetry`) that re-simulates with a lower
spike cap if the requested count cannot be reached; and post-simulation
validators that raise on negative concentrations or yields exceeding the
theoretical maximum. ``NSKFermentation`` is what actually drives the
*S. cerevisiae* fed-batch model from :doc:`05_real_model` inside a process
flowsheet: the downstream ``isobutanol`` biorefinery (a real biorefinery
model built on `BioSTEAM <https://biosteam.readthedocs.io/>`_) instantiates
one ``NSKFermentation`` unit wired to that exact kinetic model to drive its
corn-to-ethanol(-and-isobutanol) fermentation step, feeding sized equipment
and installed cost onward into a full techno-economic analysis — the same
``reactor.installed_cost`` read above, at biorefinery scale.
