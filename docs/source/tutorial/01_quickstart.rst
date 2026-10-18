Quickstart
==========

One factory call builds a complete, industrially configured process around a
real kinetic model: ``nskinetics.processes.create_sugar_prep_and_fermentation_system``
assembles the sugar-solution preparation and fed-batch fermentation section of
an actual biorefinery model — a splitter feeding parallel initial-feed and
spike-feed conditioning trains (multi-effect evaporator, pumps, dilution-water
mixer, heat exchanger), a
:class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol` fermentor
driven by the shipped *S. cerevisiae* kinetic model, and a compressed-air
aeration loop. The complete biorefinery built around this section ships with
NSKinetics as the :doc:`isobutanol biorefinery case study
<../case_studies/isobutanol/index>`. This page runs the section end to end:
build, simulate, change the fed-batch strategy, and inspect the reactor.
Total simulation time is a few seconds once the imports are loaded; every
number and figure below is the output of the code shown.

Build the process
-----------------

The factory needs only a chemical set — ``set_thermo=True`` activates its
own, shipped set. One call builds the whole section on biosteam's main
flowsheet, with its two inlets (the saccharified sugar slurry and the seed
culture) created with default compositions matching the baseline simulation
of the :doc:`isobutanol biorefinery <../case_studies/isobutanol/index>`, ready
to use as-is or overwrite:

.. code-block:: python

   import biosteam as bst
   import nskinetics as nsk

   sugar_ferm_sys = nsk.processes.create_sugar_prep_and_fermentation_system(
       set_thermo=True)
   f = bst.main_flowsheet
   sugar_ferm_sys.diagram()

.. figure:: /_static/images/examples/tutorial_01_quickstart_flowsheet.png
   :width: 700

   The factory's flowsheet: ``S301`` splits the slurry between the
   initial-feed train (``F301`` → ``M301`` → ``H301``) and the spike-feed
   train (``F302`` → ``M302`` → ``H302``), both feeding the ``V406``
   fermentor; ``K330``/``V330`` supply compressed air for aeration.

Simulate
--------

One call runs everything. The factory builds a fed-batch strategy — a
:class:`~nskinetics.units.FedBatchStrategySpecification` wired to its own
units, attached to the fermentor as ``V406.fbs_spec`` — and
``sugar_ferm_sys.simulate()`` begins by imposing it: the initial feed is
conditioned to the spec's ``target_conc`` (220 g/L glucose), a fed-batch
spike fires whenever the reactor's glucose falls to ``threshold_conc``
(210 g/L), each spike draws from a concentrated ``spike_conc`` (600 g/L)
feed, and the evaporators, dilution water, and splitter are solved so the
physical streams actually deliver those concentrations. Then the flowsheet
runs: ``V406`` hands its mixed feed to the kinetic model, integrates the
fed-batch fermentation over the full ``tau_max`` = 72 h window — spikes,
aeration switching, and a spike-count retry included — picks the harvest
time ``tau`` where ethanol peaks, converts that one time point into its
effluent stream, and finishes by re-simulating the aeration loop at the
freshly computed air demand. Plotting the reactor's kinetic trajectory shows
the fermentation behind the process result (the ``simulate()`` call emits a
handful of biosteam ``CostWarning``/``DesignWarning`` messages — transient
solver states outside cost-correlation validity ranges — that are benign
here):

.. code-block:: python

   sugar_ferm_sys.simulate()
   f.V406.plot_simulation_results()

.. figure:: /_static/images/examples/tutorial_01_quickstart_kinetics.png
   :width: 500

   The full kinetic trajectory. Early on, 3 spikes (between ~10 and
   ~12 h) hold glucose (``[s_glu]``) in the 210–220 g/L band while aerobic
   growth builds cells (``[x]``); aeration ends at ~12.4 h, when the cell
   density passes its stage-1 cutoff. Glucose is then drawn down to ~0 as
   ethanol (``[s_EtOH]``) climbs to ~114.5 g/L, and the dashed
   ``fermentation end`` line marks the selected harvest time ``tau`` ≈
   61.8 h.

Change the fed-batch strategy
-----------------------------

The strategy imposed above is held by the specification the factory attached
to the fermentor, ``f.V406.fbs_spec`` (alias
``fed_batch_strategy_specification``). Calling its ``load_specifications()``
with new values re-imposes the strategy immediately — writing the
concentrations into the kinetic model and re-solving the evaporators,
dilution water, and splitter so the physical streams deliver them — and the
values persist on the spec, so every later ``simulate()`` keeps the new
strategy. Dropping the initial-feed target to 200 g/L and letting glucose
fall to 180 g/L between spikes:

.. code-block:: python

   f.V406.fbs_spec.load_specifications(target_conc=200, threshold_conc=180,
                                       spike_conc=600)
   sugar_ferm_sys.simulate()
   f.V406.plot_simulation_results()

.. figure:: /_static/images/examples/tutorial_01_quickstart_kinetics_retuned.png
   :width: 500

   The same fermentation under the retuned strategy. The batch now starts at
   200 g/L glucose (the re-solved initial-feed evaporator delivers
   ~200.1 g/L), and the wider 180–200 g/L band means fewer, larger spikes —
   2 instead of 3 — so less total glucose is fed: it is exhausted by ~41 h,
   ethanol peaks earlier and lower (~109.5 g/L), and the selected harvest
   time drops to ``tau`` ≈ 43.0 h.

Inspect the reactor
-------------------

The same simulation drives a real process unit. ``show()`` lists the
fermentor's inlet and outlet streams — the effluent carries the kinetic
result at the harvest time, mapped onto biosteam chemicals — and
``results()`` is biosteam's standard design-and-cost table, sized and costed
from the kinetic ``tau``. Both reflect the retuned strategy simulated above:

.. code-block:: python

   f.V406.show(N=100)

.. code-block:: text

   FermentationSaccharomycesEthanolIsobutanol: V406
   ins...
   [0] glucose_initial_feed  from  HXutility-H301
       phase: 'l', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water    8.16e+03
                       Glucose  164
                       NH3      9.39
   [1] seed
       phase: 'l', T: 298.15 K, P: 101325 Pa
       flow (kmol/hr): Water  0.633
                       Yeast  0.0243
   [2] glucose_spike_feed  from  HXutility-H302
       phase: 'l', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water    836
                       Glucose  50.4
                       NH3      2.89
   [3] s5  from  IsenthalpicValve-V330
       phase: 'g', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): O2  14.3
                       N2  54
   outs...
   [0] fermentation_vent
       phase: 'g', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water       17.7
                       Ethanol     7.47
                       AceticAcid  0.0122
                       CO2         414
                       O2          7.17
                       N2          54
   [1] fermentation_effluent
       phase: 'l', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water       8.98e+03
                       Ethanol     380
                       AceticAcid  3.32
                       Glucose     9.55e-14
                       Yeast       91.9

.. code-block:: python

   f.V406.results()

.. code-block:: text

   FermentationSaccharomycesEthanolIsobutanol              Units      V406
   Electricity              Power                             kW       148
                            Cost                          USD/hr      11.5
   Chilled water            Duty                           kJ/hr -1.61e+07
                            Flow                         kmol/hr  1.08e+04
                            Cost                          USD/hr      80.6
   Design                   Reactor volume                    m3  3.17e+03
                            Batch time                        hr      61.3
                            Loading time                      hr      15.3
                            Number of reactors                           4
                            Recirculation flow rate        m3/hr      46.6
                            Reactor duty                   kJ/hr -1.61e+07
                            Cleaning and unloading time       hr         3
                            Working volume fraction                    0.9
   Purchase cost            Heat exchangers (x4)             USD  8.66e+04
                            Reactors (x4)                    USD  3.36e+06
                            Agitators (x4)                   USD  2.09e+05
                            Cleaning in place                USD  4.12e+05
                            Recirculation pumps (x4)         USD  1.37e+05
   Total purchase cost                                       USD  4.21e+06
   Installed equipment cost                                  USD   6.6e+06
   Utility cost                                           USD/hr      92.1

Where to next
-------------

The rest of the tutorial builds this picture up one layer at a time:
:doc:`02_build_and_simulate` writes and simulates a kinetic model from
scratch; :doc:`03_events` and :doc:`04_fed_batch_feeding` add the event API
and ``FeedSpike`` feeding used inside the fermentor here; :doc:`05_real_model`
tours the shipped *S. cerevisiae* model itself;
:doc:`06_process_tea_bridge` unpacks the kinetics-to-biosteam bridge that
made ``V406`` possible; and :doc:`07_sensitivity_identifiability` and
:doc:`08_calibration` ask which parameters matter and which ones
measurements can determine, and fit them to fed-batch data.
:doc:`../concepts` covers the ideas in more depth.
