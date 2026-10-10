Feeding strategies
==================

.. currentmodule:: nskinetics

How should the fermentor be fed? In NSKinetics a fed-batch policy comes down
to three numbers. The first is the glucose concentration the fermentation
starts at and is topped back up to (``target_conc``). The second is the
concentration that triggers a top-up (``threshold_conc``). The third is how
many top-ups ("spikes") are allowed (``max_n_spikes``). The kinetic model, not
the user, then decides how many spikes actually fire and when the
fermentation ends. This page asks which policy minimizes the minimum ethanol
selling price (MESP) of the scenario-A corn biorefinery (the wild-type ethanol
strain, with the isobutanol pathway off). It also asks whether fed-batch
operation beats a plain batch. To answer, it sweeps a 40 × 40 threshold ×
target grid and optimizes the spike cap for MESP at every point.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: +4.6 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      MESP penalty for a sugar target 10 g/L too high (230 → 240 g/L),
      ten times the +0.46 % penalty for one 10 g/L too low

      +++
      Threshold × target sweep, build 2b5b27d-equivalent

   .. grid-item-card:: +0.41 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      MESP of a plain batch at 230 g/L, compared with the best fed-batch
      policy in the grid

      +++
      Threshold × target sweep, build 2b5b27d-equivalent

   .. grid-item-card:: +31.6 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      MESP of the titer-maximizing policy, even though its capital and
      operating costs both fall

      +++
      Threshold × target sweep, build 2b5b27d-equivalent

   .. grid-item-card:: 291 of 819
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      feasible policies that lose money (IRR < 0); MESP spans
      3.82–5.03 $/GGE across them

      +++
      Threshold × target sweep, build 2b5b27d-equivalent

.. admonition:: Key insights
   :class: tip

   - **Sugar loading, not feed mode, sets cost.** What matters is where the
     target sits relative to the concentration the saccharified slurry
     delivers without evaporation (between 230 and 240 g/L of water here).
     Batch versus fed-batch matters much less.
   - **Err on the dilute side.** Overshooting that boundary by 10 g/L sends
     the initial feed through multi-effect evaporation and adds 11.6 MM$ of
     capital. Undershooting by the same margin costs a tenth as much.
   - **Chasing one fermentation metric is costly, and each fails through a
     different cost channel.** Maximizing titer raises MESP by +31.6 %,
     productivity by +6.2 %, cell density by +4.8 % and yield by +3.7 %.
   - **Give productivity up first.** The cost optimum keeps 97.6 % of the
     attainable titer and 99.6 % of the attainable yield, but only 64 % of
     the attainable productivity.
   - **The baseline policy is near-optimal.** Scenario A's 217/221 g/L
     policy (7 spikes) lies 0.54 % above the grid optimum (226/230 g/L,
     4 spikes), on the safe, dilute side of the evaporation cliff.

Fed-batch feeding in NSKinetics terms
-------------------------------------

Inside the fermentor ``V406``, the shipped *S. cerevisiae* model carries a
glucose :class:`FeedSpike` (see :doc:`/tutorial/04_fed_batch_feeding` for
the mechanics and :doc:`/tutorial/05_real_model` for this model). Each time
``[s_glu]`` falls to ``threshold_conc_glu_spike``, the spike adds enough
600 g/L spike feed to restore ``target_conc_glu_spike``, which grows the
working volume. It stops once ``max_n_glu_spikes`` spikes have fired. At
the process level, a
:class:`~nskinetics.units.FedBatchStrategySpecification` owns those model
variables. The sugar-preparation and fermentation factory
(:obj:`~nskinetics.processes.create_sugar_prep_and_fermentation_system`)
attaches it to the fermentor as ``V406.fbs_spec``. Its five specifications
are:

.. list-table::
   :header-rows: 1
   :widths: 21 43 16 20

   * - Specification
     - Meaning
     - Scenario A
     - In this sweep
   * - ``target_conc``
     - Glucose concentration (g per L of water) of the initial feed. The
       fermentation starts here, and every spike restores it.
     - 221.25 g/L
     - 10–400 g/L, 40 steps
   * - ``threshold_conc``
     - Reactor glucose concentration that triggers a spike. It has no
       effect when no spike may fire.
     - 217.125 g/L
     - 1–400 g/L, 40 steps (only threshold < target is valid)
   * - ``spike_conc``
     - Glucose concentration of the spike feed.
     - 600 g/L
     - 600 g/L
   * - ``max_n_spikes``
     - Upper bound on the number of spikes; ``0`` means batch.
     - 16
     - optimized for MESP at each point, over 0–40
   * - ``tau_max``
     - Longest allowed fermentation.
     - 120 h
     - 120 h

Three behaviours are decided by the model and the flowsheet, not by the
policy:

- **The number of spikes that fire is emergent.** ``max_n_spikes`` is only
  an upper bound. The fermentor's
  :class:`~nskinetics.units.SpikeReduceRetry` re-simulates with the cap
  lowered by one until the final glucose rounds to 0.00 g/L or the cap
  reaches 0. The last attempt runs at a cap of 1, so a fermentation that
  never exhausts its glucose ends with at most one spike fired (with
  scenario A's threshold just below its target, exactly one: the "stalled"
  signature on :doc:`kinetic_sweeps`). A sweep therefore reports the number
  of spikes that actually fired, not the cap it was given.
- **The harvest time τ is read off the trajectory.** It is the first time,
  at or after the last spike, at which glucose (rounded to whole g/L)
  reaches its minimum (below 0.5 g/L). It is capped at ``tau_max``. The
  biorefinery's fermentor uses ``tau_update_policy=['min', '[s_glu]']``.
  Productivity is titer / τ.
- **The upstream trains deliver the requested concentrations.**
  :meth:`~nskinetics.units.FedBatchStrategySpecification.load_specifications`
  tunes two :class:`~nskinetics.units.ConcentrationActuator` knobs per feed
  train. For the initial feed these are the vapour fraction ``V`` of
  evaporator ``F301`` (the concentrator) and the water-to-sugar ratio of
  dilution mixer ``M301`` (the diluter). The spike feed uses ``F302`` and
  ``M302``. If the slurry is already at least as concentrated as
  ``target_conc``, ``F301`` stays at ``V = 0`` and ``M301`` dilutes. Only a
  target above the slurry's own concentration makes the solver open the
  evaporator. This switch is where the evaporation cliff in the results
  comes from.

Two guards keep a policy physical. ``load_specifications`` rejects any
policy that does not satisfy ``threshold_conc < target_conc < spike_conc``.
The biorefinery also refuses a run whose working volume would grow more than
20-fold (``max_fed_batch_volume_ratio``). Each spike multiplies the volume by
(``spike_conc`` − ``threshold_conc``) / (``spike_conc`` − ``target_conc``),
a factor of 1.0106 at the grid optimum.

To try a policy yourself, set it through the model specification (the top of
the specification's precedence order) and read the fermentor and the TEA:

.. code-block:: python

   from nskinetics.biorefineries import isobutanol
   # build the flowsheet; simulate the baseline once
   isobutanol.load()
   from nskinetics.biorefineries.isobutanol import scenarios

   # A's kinetics + policy: cap 16, 217.125 / 221.25 g/L
   h = scenarios.load_scenario('A')
   V406 = h['V406']
   print(h['fbs_spec'].current_specifications)

   # A plain batch at 230 g/L: with max_n_spikes=0 the threshold has no effect
   h['model_specification'](threshold_conc=226.1, target_conc=230.,
                            max_n_spikes=0)
   r = V406.nsk_results_specific_tau_dict
   # τ [h], spikes fired, titer [g/L water]
   print(V406.tau, r['curr_n_glu_spikes'], r['[s_EtOH]'])
   # {'IRR': ..., 'MPSPs': {'ethanol': $/kg, 'isobutanol': nan}}
   print(h['solve_TEA']())

Run it on the current model build, and the economics will differ slightly
from this page's numbers, which use an earlier build (see
`Caveats and limitations`_).

Method
------

- **Strain and process.** The analysis uses scenario-A kinetics, applied
  from the scenario's parameter workbook (wild-type ethanol strain,
  isobutanol pathway off). The script sets the workbook and feeding policy
  directly rather than through ``load_scenario``, so it runs without the
  enzyme-burden constraint (:doc:`enzyme_burden`), which would be inert at
  this strain anyway. The rest of the flowsheet and the TEA basis are as in
  :doc:`process_and_tea`.
- **Grid.** The threshold takes 40 values (``linspace(1, 400, 40)``, step
  10.23 g/L) and the target takes 40 values (``linspace(10, 400, 40)``,
  step 10 g/L), giving 1,600 cells. Only the 819 cells with threshold <
  target are valid policies, and all 819 returned finite results. The other
  781 are left empty by construction.
- **Spike cap, optimized per cell.** The biorefinery helper
  ``optimize_max_n_glu_spikes_for_MPSP(bounds=(0, 40), Ns=41)`` tries caps
  0, 1, 2, … and stops once a higher cap adds no spike that actually fires.
  It then re-simulates the best cap. Between 0 and 30 spikes fired at the
  chosen caps, so the cap bound of 40 was never reached. Batch (0 spikes)
  won at 631 cells and fed-batch at 188. Individual cap evaluations were
  not logged; the early stop implies between about 3,400 and 34,398 of
  them, each a converged flowsheet simulation.
- **Fixed settings.** Spike feed 600 g/L, ``tau_max`` 120 h, and the
  glucose-depletion harvest rule above.
- **Metrics.** The sweep records 19 metrics per cell. They include the
  purity-adjusted ethanol MPSP at a 15 % IRR, the IRR at default prices
  (ethanol 0.835 $/kg), TCI, AOC, titer, productivity, yield, the number of
  spikes that fired, τ, the heating duty of the two sugar-solution
  evaporators, cell density and aeration.

.. note::

   **Conventions.** MESP is the purity-adjusted ethanol MPSP at a 15 % IRR,
   converted from $/kg to $ per gallon of gasoline equivalent:
   $/GGE = $/kg × 0.789 kg/L × 3.785 L/gal ÷ 0.67 GGE/gal = $/kg × 4.45775.
   Titer, productivity and glucose concentrations are per litre of
   **water**, not broth. Yield is g ethanol per g glucose **added**. The
   sweep's MESP and IRR values come from code equivalent to
   heat-exchanger-network build ``2b5b27d`` (2b5b27d-equivalent). On that
   build the scenario-A baseline MESP is 3.844 $/GGE.

Results
-------

.. figure:: /_static/images/isobutanol/feeding/feeding_strategy_mesp_map_and_baseline_time_course.png
   :target: ../../_static/images/isobutanol/feeding/feeding_strategy_mesp_map_and_baseline_time_course.png
   :alt: Panel A is a heat map of MESP over threshold glucose concentration
         (x, 0 to 400 g/L) and target glucose concentration (y, 0 to 400 g/L).
         Only the triangle above the diagonal, where threshold is below target,
         is coloured. MESP is lowest (yellow-orange, about 3.8 $/GGE) in a band
         at target 200 to 230 g/L and rises steeply above 240 g/L to dark blue
         (about 5 $/GGE) at 400 g/L. A cyan star marks the MESP optimum near
         226/230 g/L. White markers label the titer, cell-density,
         productivity and yield optima, and a cyan pentagon labels the
         TCI minimum. A white staircase line separates batch policies (above)
         from fed-batch policies (below). Panel B shows the concentrations of
         glucose, ethanol, cell density, active cell density, acetate and
         isobutanol over 80 hours for the baseline policy. Glucose holds near
         220 g/L with small spikes at 9 to 12 hours, then falls to zero by
         about 60 hours while ethanol rises to about 115 g/L. Dashed lines mark
         the end of aeration and the end of fermentation.
   :class: pub-figure
   :width: 100%

   **Feeding-strategy map (A) and the baseline fermentation (B).** MESP is
   lowest at targets just below the colour jump between 230 and 240 g/L, the
   evaporation cliff, and depends far less on the threshold. Panel B is the
   scenario-A baseline policy, *not* the panel-A optimum. Panel A uses build
   ``2b5b27d``-equivalent code; panel B's concentrations come from build
   ``6d4776f``, which does not affect them. Click the figure to open it at
   full resolution.

.. dropdown:: How to read this figure

   - **(A)** MESP over the 819 valid threshold × target policies, each with
     its spike cap optimized for MESP. The empty lower-right triangle holds
     the invalid policies (threshold ≥ target).
   - **Markers in (A).** The cyan star is the grid optimum: threshold
     226.1 / target 230 g/L, with 4 spikes and an MESP of 3.8234 $/GGE. The
     white markers are the first grid points that reach the best titer, cell
     density, productivity and yield, and the cyan pentagon is the lowest
     TCI. The titer and cell-density markers sit on broad plateaus (248 and
     213 cells within 0.1 % of their maxima), so read them as regions rather
     than points.
   - **The staircase in (A).** The white staircase line is the lowest target
     at which the optimized cap is zero, with batch above it and fed-batch
     below.
   - **(B)** The scenario-A **baseline** policy (217.125 / 221.25 g/L,
     7 spikes). All 7 spikes fall between 9.0 and 11.5 h, before aeration
     ends at 12.13 h, which shows up as the small ripple on the glucose
     trace. The fermentation ends at τ = 61.62 h, when glucose first falls
     below 0.5 g/L. Ethanol then stands at 114.54 g/L (water basis) and is
     still creeping up (114.95 g/L at 120 h).

Sugar loading, not feed mode, sets cost
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The table below takes the best policy at each target over all thresholds.
MESP falls gently as the target rises toward 230 g/L, then jumps:

.. list-table:: Best MESP at each target (envelope over thresholds; build ``2b5b27d``-equivalent)
   :header-rows: 1
   :stub-columns: 1
   :widths: auto
   :class: ibo-num

   * - Target (g/L)
     - 160
     - 200
     - 220
     - **230**
     - **240**
     - 250
     - 300
     - 400
   * - ΔMESP vs optimum
     - +2.72 %
     - +1.29 %
     - +0.46 %
     - **0**
     - **+4.63 %**
     - +4.82 %
     - +12.5 %
     - +31.6 %
   * - Sugar-evaporation duty, both trains (GJ/h)\ :sup:`*`
     - 41.0
     - 23.0
     - 14.0
     - 9.1
     - **62.2**
     - 61.9
     - 56.7
     - 73.5
   * - TCI (MM$)
     - 144.4
     - 139.8
     - 137.6
     - 136.4
     - **148.0**
     - 149.5
     - 135.1
     - 125.0

:sup:`*` Heating duty of the initial-feed evaporator ``F301`` plus the
spike-feed evaporator ``F302``. The best policies at targets 160–230 g/L are
fed-batch (1–6 spikes), so below 240 g/L this duty is the spike-feed train
concentrating its feed to 600 g/L; ``F301`` stays idle up to 230 g/L and
switches on at 240 g/L, which is the cliff.

- **The evaporation cliff.** Batch policies need no sugar-evaporation duty
  at any target up to 230 g/L, but they need 60.8–73.5 GJ/h at 240 g/L and
  above. The slurry's own glucose concentration therefore lies in
  (230, 240] g/L of water. Going from a 230 to a 240 g/L batch adds
  11.6 MM$ of capital (TCI 135.03 → 146.64 MM$) and 3.5 MM$/y of operating
  cost, but raises titer only from 105.2 to 110.1 g/L.
- **Why the cliff is a step.** The duty jumps rather than growing with the
  small amount of extra water removed. The actuator logic above explains
  when the evaporator turns on. That ``F301`` then processes the *whole*
  initial feed is inferred from the code, not read from unit duty tables. A
  unit-level check at targets 230 and 240 g/L would settle whether the step
  is water removal or preheating of the whole stream. That check has not
  been run.
- **The tolerance ceiling.** Titer saturates at 115.0–115.5 g/L for every
  target of 250 g/L and above. In batch, yield is exactly titer / target,
  so sugar above roughly 250 g/L is never converted. The cost optimum sits
  at the lower of the two ceilings: the no-evaporation ceiling (about
  230–240 g/L) comes before the biological tolerance ceiling (about
  250 g/L).

In this flowsheet the solids loading upstream is fixed, so the dilute side of
the boundary is the safe side. The *size* of the penalty (+4.6 % and
11.6 MM$ per 10 g/L) depends on the flowsheet. The rule does not: never
thermally concentrate the main feed. Raising the solids loading upstream is
the lever that would move the ceiling.

Batch versus fed-batch
~~~~~~~~~~~~~~~~~~~~~~

The best fed-batch policy beats the best batch by only 0.016 $/GGE (0.41 %).
The near-optimal set (within 1 % of the optimum) spans thresholds of
1–226 g/L and 0–8 spikes, but **only targets 210–230 g/L**. The optimum is
sharp in the target and flat in the threshold and the spike count.

.. list-table:: Three policies compared (MESP, IRR and TCI on build ``2b5b27d``-equivalent)
   :header-rows: 1
   :stub-columns: 1

   * -
     - Grid optimum
     - Best batch
     - Scenario-A baseline
   * - Threshold / target (g/L)
     - 226.1 / 230
     - no effect / 230
     - 217.125 / 221.25
   * - Spikes fired
     - 4 (small top-ups)
     - 0
     - 7
   * - MESP ($/GGE)
     - 3.8234
     - 3.8390 (+0.41 %)
     - 3.8440 (+0.54 %)
   * - IRR
     - 12.98 %
     - 12.64 %
     - 12.60 %
   * - TCI (MM$)
     - 136.41
     - 135.03
     - 139.3
   * - Titer (g/L water)
     - 112.73
     - 105.16
     - 114.54
   * - Productivity (g/L/h)
     - 2.31
     - 2.98
     - 1.86
   * - Yield (g/g)
     - 0.4596
     - 0.4572
     - 0.4602
   * - τ (h)
     - 48.89
     - 35.3
     - 61.62

The baseline's titer, productivity, yield and τ come from its ``6d4776f``
time course. Fermentation metrics do not depend on the build, and they match
the ``2b5b27d`` values. The baseline is in effect a front-loaded top-up of
the aerobic growth phase, and the rest of its run is batch. Its titer
(114.5 g/L) and yield (0.460 g/g) lie inside reported U.S. high-gravity
ranges (100–142 g/L and 0.42–0.48 g/g). Its productivity (1.86 g/L/h) lies
**below** the reported 2.0–4.4 g/L/h. The grid optimum (2.31) and the
230 g/L batch (2.98) lie inside it. The model reports per litre of water,
whereas industrial figures are usually per litre of broth.

Single-metric optima are costly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each fermentation metric, maximized on its own, costs money, and each fails
through a different cost channel:

.. list-table:: Grid optima of single metrics (first argmax on the grid; build ``2b5b27d``-equivalent)
   :header-rows: 1
   :widths: 14 20 12 18 36

   * - Metric maximized
     - Policy (threshold / target, spikes)
     - ΔMESP
     - ΔTCI vs cost optimum
     - What drives the penalty
   * - Titer
     - 72.6 / 400 g/L, batch
     - **+31.6 %**
     - −11.4 MM$ (AOC −1.5 MM$/y)
     - Lost product. Yield collapses to 0.289 g/g because sugar above the
       tolerance ceiling is never converted, even though capital and
       operating costs fall.
   * - Productivity
     - 41.9 / 180 g/L, batch
     - +6.2 %
     - +7.2 MM$ (AOC +3.1 MM$/y)
     - Dilute beer (80.45 g/L at 3.60 g/L/h). Attributing the cost to
       downstream separation is an inference; the costs were not broken
       down by unit.
   * - Cell density
     - 185.2 / 250 g/L, batch
     - +4.8 %
     - +13.1 MM$
     - The target sits past the evaporation cliff (61.9 vs 9.1 GJ/h of
       sugar-evaporation duty).
   * - Yield
     - 11.2 / 140 g/L, 1 spike
     - +3.7 %
     - +13.2 MM$
     - A large, evaporated spike stream (49 vs 9 GJ/h).

The cost optimum is a compromise. It reaches 97.6 % of the best titer, 99.6 %
of the best yield and 99.5 % of the best cell density, but only **64 %** of
the best productivity. Of the three fermentation metrics (titer, rate and
yield, or TRY), productivity is the one to give up first.

Titer runs the wrong way
~~~~~~~~~~~~~~~~~~~~~~~~

Across the 819 policies, higher titer goes with *higher* MESP: Spearman
ρ(MESP, titer) = +0.87. Of the policies, 63 % reach at least 99 % of the
best titer, and their MESP spans 3.88–5.03 $/GGE, so titer cannot rank
feeding designs. Yield is the best single proxy (ρ = −0.94), but even its
strict maximum costs +3.7 %. The reason is the tolerance ceiling. At this
fixed strain, titer saturates near 115 g/L, so the only way to raise it is to
over-supply sugar that is never converted. This is a statement about a fixed,
wild-type strain. :doc:`kinetic_sweeps` turns to the strain-side levers.

Caveats and limitations
-----------------------

- **Model builds.** The sweep (panel A and every table here) ran on 2026-09-20
  on code equivalent to the heat-exchanger-network build ``2b5b27d``. The
  current scenario-A pin comes from build ``6d4776f``: MPSP 0.86021 $/kg,
  or 3.835 $/GGE. The two builds' baselines differ by about 0.25 %. The
  0.54 % baseline-to-optimum gap is computed within one build (``2b5b27d``).
  Do not set this page's 3.82 $/GGE next to the uncertainty analysis's
  baseline (:doc:`uncertainty`) as though they were one basis.
- **Grid resolution and the cliff edge.** The target step is 10 g/L. The
  true optimum probably lies between 230 g/L and the slurry's own
  concentration, right at the cliff edge. It is not robust to feedstock or
  solids variation. That margin is why the scenario-A baseline keeps the
  dilute-side 217/221 policy instead of moving to 226/230. In the 6,000-run
  uncertainty analysis (build ``6d4776f``), sugar-preparation heating
  followed the spike count (ρ = 0.976) rather than the starch content
  (ρ = 0.019). No sample went above 35.07 GJ/h, so none crossed the cliff.
- **Where the baseline came from.** The 217/221 g/L baseline, and scenario
  B's batch policy (34.25 / 140 g/L, with the threshold having no effect),
  were the optima of an older 25 × 25 sweep (2026-09-02) on a superseded
  model. They are **not** the optimum of the grid shown here.
- **Fixed strain.** Every result is for scenario-A kinetics with isobutanol
  off. The optimum moves with the strain, so the strain-design campaigns
  re-optimize feeding together with kinetics (:doc:`kinetic_optimization`).
- **Concentration basis.** Titers and concentrations are per litre of water.
  The water-to-broth conversion has been verified for scenario B only.
- **Evaporation mechanism.** It is inferred from the code (see above); the
  unit-level check that would confirm it has not been run.
- **Reproducibility margin.** Within a sweep, numerical path dependence is
  about 1e-4 relative, so the 0.4 % differences above are resolvable. Treat
  differences under about 0.05 % as ties, for example 3.8234 versus
  3.8242 $/GGE at target 230 g/L with 4 versus 1 spikes.
- **Markers are argmax points.** The panel-A markers mark the first grid
  point that reaches each extreme. The titer and cell-density markers in
  particular stand for plateaus, not unique optima.

Reproduce it
------------

No analysis results ship with NSKinetics: the sweep CSVs and the trajectory
cache behind the figure must be regenerated, or ``NSKINETICS_IBO_RESULTS_DIR``
set to a results folder that already holds them. Long simulations are
opt-in. Run one simulation process at a time (:doc:`running_analyses`).
Script paths are relative to the package folder.

.. list-table::
   :header-rows: 1
   :widths: 30 42 28

   * - Script
     - What it does and writes
     - Runtime
   * - ``analyses/feeding_strategies/evaluate_feeding_strategies.py``
     - The threshold × target sweep, with the cap optimized per cell.
       Writes one CSV per metric to the results root,
       ``ibo_(40, 40, 1)_Thres_Targe_Max n_A__<metric>.csv`` (19 metrics),
       plus contour plots. It defaults to scenario **B** on a 25 × 25 grid,
       so set ``IBO_SWEEP_SCENARIO=A`` and ``IBO_SWEEP_STEPS=40`` to
       reproduce this page. ``IBO_SWEEP_REPLOT_FROM_CSV=1`` redraws the plots
       from saved CSVs without simulating. The sweep is **not** checkpointed,
       so a crash restarts it from the beginning.
     - The longest grid sweep: thousands of converged flowsheet simulations
       (see `Method`_). At 15–30 s each that is at least ~15–30 h; plan
       for more than a day.
   * - ``plots/plot_feeding_strategy_baseline_A.py``
     - The two-panel figure above. Stage 1 runs one ``load()`` and
       simulates the scenario-A baseline, caching its time course to
       ``feeding_strategy_baseline_A_trajectory.{csv,json}``. It runs only
       if that cache is missing, or with ``--resimulate``. Stage 2 reads
       the sweep CSVs and the cache without importing the package (it is
       safe to run beside a running simulation). It writes
       ``publication/Feed-strat/feeding_strategy_baseline_A.{png,pdf}``.
     - Stage 1: one baseline load and simulation. Stage 2: seconds.

.. tab-set::

   .. tab-item:: Bash

      .. code-block:: bash

         # from the package folder (see running_analyses to find it)
         export MPLBACKEND=Agg
         IBO_SWEEP_SCENARIO=A IBO_SWEEP_STEPS=40 \
             python analyses/feeding_strategies/evaluate_feeding_strategies.py
         python plots/plot_feeding_strategy_baseline_A.py

   .. tab-item:: PowerShell

      .. code-block:: powershell

         # from the package folder (see running_analyses to find it)
         $env:MPLBACKEND = "Agg"
         $env:IBO_SWEEP_SCENARIO = "A"; $env:IBO_SWEEP_STEPS = "40"
         python analyses\feeding_strategies\evaluate_feeding_strategies.py
         python plots\plot_feeding_strategy_baseline_A.py

Next: :doc:`uncertainty` holds this feeding policy fixed and asks how much
of the plant's cost risk comes from the strain rather than from prices and
plant context; :doc:`kinetic_sweeps` then asks which strain changes pay.
