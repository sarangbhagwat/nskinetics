.. _isobutanol-case-study:

Isobutanol biorefinery
======================

.. toctree::
   :maxdepth: 1
   :hidden:

   process_and_tea
   feeding_strategies
   uncertainty
   kinetic_sweeps
   kinetic_optimization
   objective_landscape
   enzyme_burden
   sensitivity
   running_analyses

``nskinetics.biorefineries.isobutanol`` is a complete corn dry-grind
biorefinery, built in `BioSTEAM <https://biosteam.readthedocs.io/>`__, that
turns corn into fuel ethanol, isobutanol, distillers' dried grains with
solubles (DDGS) and corn oil. Its fermentor, ``V406``, is not a fixed-yield
black box: every simulation integrates an NSKinetics fed-batch kinetic model
of *Saccharomyces cerevisiae* (the :ref:`Lei et al. 2001 <ibo-ref-lei2001>`
core plus an engineered
Ehrlich isobutanol pathway, reactions r13–r17), and the titer, yield, rate and
fermentation time it returns size the separations, utilities and capital
downstream. That makes it a test bed for one question: **when strain kinetics
and the feeding policy set titer, rate and yield (TRY) together, which strain
and process decisions actually make the plant profitable — and are they the
ones a TRY target would pick?** The package ships the flowsheet, its
scenarios and every analysis script, so each result below can be re-run
(:doc:`running_analyses`).

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 12.6 % → 27.3 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      IRR of the wild-type ethanol plant vs the best co-producing strain
      found in silico, at essentially the same capital (TCI 139.3 vs
      137.7 MM$). A prospective, bounds-limited target: 39.6 g/L
      isobutanol is 15–20× the best reported *S. cerevisiae* titers, at one
      isobutanol price, from two random seeds (:ref:`scope <ibo-opt-scope>`).

      +++
      Campaign record, build 2b5b27d

   .. grid-item-card:: ~92 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      of the baseline's ethanol-cost (MESP) variance comes from context —
      corn and DDGS prices, uptime, capacity, starch. All 40 kinetic
      parameters together explain ~5 %.

      +++
      6,000-sample uncertainty run, build 6d4776f

   .. grid-item-card:: 11 of 12
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      titer-, rate- or yield-optimal strains missed the 15 % IRR hurdle;
      none exceeded 15.5 % IRR.

      +++
      Six objectives × two seeds, build 2b5b27d

   .. grid-item-card:: +21 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      MESP from fully removing ethanol inhibition, even though
      productivity rose 3.7× — IRR fell from 12.6 % to −9.4 %.

      +++
      Glycolysis × ethanol-inhibition plane, build 2b5b27d

Each card names the model version its numbers were computed on. The
versions (``2b5b27d``, ``6d4776f``, ``95e9362``) differ only in the
heat-exchanger-network synthesis package (hensmith), which shifts TEA results
by up to a few percent, so numbers from different versions are never compared
directly. See :ref:`isobutanol-conventions`.


What the analyses found
-----------------------

**The story, in the order of the pages.** At the wild-type ethanol plant,
cost is set by sugar loading and commodity context, not by the strain:
about 92 % of the MESP variance is prices, uptime and scale. The classic
strain levers have little headroom: the wild type sits within 1.5 % of the
best design on the glycolysis × ethanol-inhibition plane, and full ethanol
tolerance raises MESP by 21 %. A step change needs a higher-value second
product, yet strains optimized for titer, rate or yield miss it: 11 of 12
TRY optima fall short of the 15 % IRR hurdle. A search that maximizes
profitability directly can find a co-producer near 27 % IRR, but one of two
random seeds stalled on an ethanol-only design at 15.8 %. Preloading that
search with the TRY campaigns' trials made it fast and repeatable and
returned the flagship design: 27.3 % IRR at unchanged capital (a prospective
target; see its :ref:`scope <ibo-opt-scope>`). The landscape, burden and
sensitivity pages explain why: profitable designs form small, separated
islands that need ALS, the downstream Ehrlich enzymes and Adh6 raised
together within the cell's proteome budget.

Each page states the question, the method, the results with their figures,
and the caveats; the build behind each headline is in the card footer.

.. grid:: 1 1 2 2
   :gutter: 3

   .. grid-item-card:: Process and techno-economics
      :link: process_and_tea
      :link-type: doc
      :shadow: sm

      **Co-production pays on the revenue side.** The flagship co-producer
      lifts sales 19 % at 1.2 % lower capital; fermentation is only ~10 %
      of installed cost, so strain changes pay through what the plant sells.

      +++
      Flowsheet, TEA basis, scenarios A/B · builds 6d4776f and 2b5b27d

   .. grid-item-card:: Feeding strategies
      :link: feeding_strategies
      :link-type: doc
      :shadow: sm

      **Sugar loading, not feed mode, sets cost.** Overshooting the target
      glucose concentration by 10 g/L costs +4.6 % MESP, about ten times
      undershooting it; a pure batch at a 230 g/L charge is only +0.41 % off
      the best of 819 feasible policies.

      +++
      Threshold × target sweep · build 2b5b27d-equivalent

   .. grid-item-card:: Uncertainty analysis
      :link: uncertainty
      :link-type: doc
      :shadow: sm

      **Context, not kinetics, sets cost risk.** ~92 % of MESP variance is
      context, ~5 % all 40 kinetic parameters — yet per unit change, yield
      moves MESP almost as much as the corn price (elasticity −0.54 vs
      +0.65).

      +++
      6,000 Monte Carlo simulations, Spearman · build 6d4776f

   .. grid-item-card:: Kinetic sweeps
      :link: kinetic_sweeps
      :link-type: doc
      :shadow: sm

      **More ethanol tolerance can cost more.** Removing ethanol inhibition
      entirely makes the strain 3.7× more productive yet raises MESP by
      21 %, because once the feeding policy stops delivering sugar the
      extra tolerance feeds growth. The wild-type strain already sits within
      1.5 % of the best design on the glycolysis × ethanol-inhibition plane,
      and no point on it reaches the 15 % hurdle.

      +++
      40 × 40 kinetic planes · build 2b5b27d

   .. grid-item-card:: Kinetic optimization
      :link: kinetic_optimization
      :link-type: doc
      :shadow: sm

      **TRY optima miss the profit optimum.** Eleven of twelve titer-,
      rate- or yield-optimal strains miss the hurdle; profit-directed search
      reached 27.3 % IRR, but whether it got there depended on the random
      seed — preloading TRY designs made it fast and repeatable.

      +++
      Bayesian-optimization campaigns · build 2b5b27d

   .. grid-item-card:: Objective landscape
      :link: objective_landscape
      :link-type: doc
      :shadow: sm

      **Profitable designs form small, separated islands** — an
      ethanol-only shelf, a Pdc-null isobutanol island and the
      co-production island. Profitability looks needle-like only in
      optimizer-sampled data; under unbiased sampling it is as smooth as
      the TRY objectives.

      +++
      Pooled campaign trials and 1-D slices · build 2b5b27d

   .. grid-item-card:: Enzyme burden
      :link: enzyme_burden
      :link-type: doc
      :shadow: sm

      **Yield optima overload the proteome.** The isobutanol- and
      ethanol-yield optima run at 157 % and 130 % of the penalty-free
      enzyme budget, throttling growth; the flagship profit optimum uses
      97 % with no derating.

      +++
      Proteome-allocation constraint · build 2b5b27d

   .. grid-item-card:: Global sensitivity analysis
      :link: sensitivity
      :link-type: doc
      :shadow: sm

      **Titer's driver is not profit's driver.** The branch-point enzyme
      ALS carries 36 % of ethanol-titer variance but only 3.5 % of the
      variance in profitability (Shapley effects, screening measure).

      +++
      Sobol'/Shapley analysis of 12 decisions · build 2b5b27d

   .. grid-item-card:: Running the analyses
      :link: running_analyses
      :link-type: doc
      :shadow: sm
      :columns: 12

      **Reproduce any figure.** Which script produces which result, what it
      writes, how long it runs, and how to point the scripts at existing
      results with ``NSKINETICS_IBO_RESULTS_DIR``.

      +++
      Scripts, runtimes and the results folder


Why it is an NSKinetics showcase
--------------------------------

Most biorefinery TEAs take titer, yield and productivity as *inputs*. Here
they are *outputs* of one kinetic parameter vector plus a feeding policy, so
a change to an enzyme capacity, a tolerance or a glucose-spike setpoint
propagates all the way to the selling price. One simulation runs this loop:

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 1 · Kinetics
      :shadow: sm

      ``V406`` integrates the *S. cerevisiae* ODE model with
      :class:`~nskinetics.FeedSpike` glucose spikes and a two-stage
      (aerobic → anaerobic) aeration switch.

   .. grid-item-card:: 2 · Fermentation outcome
      :shadow: sm

      The run is cut at the fermentation time τ; titer, yield, productivity,
      cell mass and the spike-feed volume become biosteam streams.

   .. grid-item-card:: 3 · Process
      :shadow: sm

      Feed evaporation, aeration, distillation, DDGS drying, wastewater
      treatment and the boiler are re-simulated, sized and costed around
      that broth.

   .. grid-item-card:: 4 · Economics
      :shadow: sm

      The TEA solves the minimum product selling price (MPSP) at a 15 %
      IRR, the IRR at market prices, and the profitability index.

The NSKinetics pieces doing the work are the
:class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol` bridge
unit (``V406``), the
:obj:`~nskinetics.processes.create_sugar_prep_and_fermentation_system`
factory that builds the sugar-preparation trains, fermentor and aeration loop
around it, and the shipped kinetic model
``nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo`` (walked through in
:doc:`/tutorial/05_real_model`). Two safeguards keep the results physical.
On every simulation the fermentor closes its spike-feed balance against what
the upstream feed train delivers (``V406.spike_feed_residual``, held to
1e-3). Wherever an enzyme-burden model is installed, a proteome-allocation
constraint derates growth for strains whose enzyme load would not fit in the
cell; it is on for scenario A (where it is inert), the ``opt_*`` scenarios,
the kinetic-sweep planes reported here, the optimization campaigns and the
Sobol' analysis, and off for scenario B, the feeding sweep and the
uncertainty run (:doc:`enzyme_burden`).

.. tab-set::

   .. tab-item:: Kinetic model inside V406

      .. figure:: /_static/images/isobutanol/overview/kinetic_model_conceptual_diagram.png
         :target: ../../_static/images/isobutanol/overview/kinetic_model_conceptual_diagram.png
         :class: pub-figure
         :width: 100%
         :alt: Reaction network of the S. cerevisiae fed-batch model: glucose to pyruvate (r1), pyruvate to acetaldehyde (r3) and ethanol (r6), acetate and respiration branches (r2, r4, r5, r8), growth and biomass states (r7, r9 to r11), and the engineered isobutanol pathway from pyruvate through alpha-acetolactate, 2,3-dihydroxyisovalerate, alpha-ketoisovalerate and isobutyraldehyde to isobutanol (r13 to r17), with product-inhibition, glucose-repression and oxygen-dependence markers, plus boxes for fed-batch glucose feeding and two-stage aeration control.

         The kinetic model behind ``V406``. Read it top-down: glucose enters
         glycolysis (r1) and pyruvate splits three ways — to ethanol via
         Pdc1 (r3) and Adh1 (r6), to respiration (r2), and into the
         engineered isobutanol pathway on the right (Ilv2+Ilv6, Ilv5, Ilv3,
         Aro10, Adh6; r13–r17). Orange bars mark product inhibition by
         ethanol, acetate and isobutanol; dial icons mark the enzyme
         capacities and tolerances the strain-design analyses vary. The
         feeding and aeration settings printed in the two top boxes are the
         model file's own defaults; the biorefinery scenarios override them
         (scenario A spikes at 217.125 → 221.25 g/L with a 16-spike cap).
         Drawn by ``nskinetics/models/s_cerevisiae_ferm_fb_inhib_mod_ibo/conceptual_diagram.py``.

   .. tab-item:: Process around it

      .. figure:: /_static/images/isobutanol/overview/block_flow_diagram.png
         :target: ../../_static/images/isobutanol/overview/block_flow_diagram.png
         :class: pub-figure
         :width: 100%
         :alt: Block flow diagram. Panel A, fermentation process: saccharified corn slurry is split into two multi-effect evaporation and dilution trains giving an initial feed and a spike feed to the fed-batch fermentor, which receives compressed air and yeast; vented gases are scrubbed and the recovered alcohols rejoin the broth. Panel B, separation process: beer column, ethanol rectifier, molecular sieves and denaturant mixing give ethanol; stripper, decanter and drying column give isobutanol above 99.9 wt%; stillage goes through centrifugation, evaporation, oil separation and drying to DDGS and corn oil.

         Block flow of the biorefinery around the fermentor. **A** — the
         saccharified slurry is split between an initial-feed and a
         spike-feed train, fermented with compressed air in the fed-batch
         (or batch) fermentor, and the vent is scrubbed so stripped alcohols
         rejoin the beer sent to distillation. **B** — a solvent-free distillation train
         recovers ethanol and isobutanol; dashed units run only when
         isobutanol is produced; the stillage becomes DDGS and corn oil.
         Wastewater treatment, the boiler-turbogenerator, utilities and the
         idle ethanol-primary train are left out. Details:
         :doc:`process_and_tea`. Drawn by
         ``plots/plot_block_flow_diagram.py``.


Scenarios at a glance
---------------------

A *scenario* bundles a parameter workbook (baseline kinetics plus uncertainty
distributions), a feeding strategy and the expected selling prices that the
package's smoke tests pin. All of them live in
``nskinetics.biorefineries.isobutanol.scenarios.SCENARIOS``. Ten are
registered: ``A``, the wild-type ethanol strain on a fed-batch baseline;
``B``, an engineered ethanol + isobutanol co-producer run as a batch; and
eight ``opt_*`` scenarios, each freezing the best design of one
strain-design campaign (:doc:`kinetic_optimization`), among them the
flagship co-producer ``opt_PI_TRY_informed``.

.. dropdown:: All ten scenarios: feeding and pinned MPSPs

   .. list-table:: Registered scenarios and their pinned MPSPs (purity-adjusted, $/kg, at a 15 % IRR)
      :header-rows: 1
      :widths: 30 34 16 10 10
      :class: ibo-scenarios

      * - Scenario
        - What it is
        - Feeding: threshold / target (g/L), spike cap
        - Ethanol MPSP
        - Isobutanol MPSP
      * - ``A``
        - Wild-type ethanol strain (isobutanol pathway off), **fed-batch
          baseline**; enzyme burden on but inert.
        - 217.125 / 221.25, cap 16 (7 spikes fire)
        - 0.86021
        - NaN
      * - ``B``
        - Engineered ethanol + isobutanol co-producer, **batch**; enzyme
          burden off (its pathway constants exceed the proteome cap).
        - batch: cap 0, target 140 (threshold 34.25 inert)
        - 0.59530
        - 1.2800
      * - ``opt_PI_uninformed``
        - Best trial #103 of the profitability (PI) campaign started from
          scratch (seed 350); stayed ethanol-only. PI 0.0761.
        - 230.77 / 235.77, cap 43
        - 0.82091
        - NaN
      * - ``opt_PI_TRY_informed``
        - Best trial #1912 of the PI campaign preloaded with the six TRY
          campaigns' designs — the **flagship co-producer** (isobutanol 39.6
          + ethanol 40.6 g/L). PI 0.8559.
        - 189.19 / 194.19, cap 1
        - 0.46066
        - 1.1123
      * - ``opt_IBO_titer``
        - Best trial #84 of the isobutanol-titer campaign: 46.04 g/L.
        - 295 / 300, cap 46
        - 207.46\ :sup:`*`
        - 2.4625
      * - ``opt_IBO_yield``
        - Best trial #841 of the isobutanol-yield campaign: 0.3814 g/g.
        - 35.06 / 40.06, cap 21
        - 5,760.9\ :sup:`*`
        - 1.8533
      * - ``opt_IBO_productivity``
        - Best trial #1472 of the isobutanol-productivity campaign:
          2.063 g/L/h.
        - 295 / 300, cap 21
        - 491.16\ :sup:`*`
        - 2.2694
      * - ``opt_EtOH_titer``
        - Best trial #1820 of the ethanol-titer campaign: 168.45 g/L.
        - 295 / 300, cap 35
        - 0.85672
        - NaN
      * - ``opt_EtOH_yield``
        - Best trial #1923 of the ethanol-yield campaign: 0.4847 g/g.
        - 85.18 / 107.54, cap 8
        - 0.86838
        - NaN
      * - ``opt_EtOH_productivity``
        - Best trial #365 of the ethanol-productivity campaign: 6.897 g/L/h.
        - 295 / 300, cap 15
        - 0.90813
        - NaN

   :sup:`*` A product made in near-zero amounts gets an absurd purity-adjusted
   MPSP; these are artifacts of the metric, not bugs. NaN means the product is
   not made at all.

   .. admonition:: Reading the pins
      :class: note

      - **A and B** are pinned on build 6d4776f. Scenario A's ethanol MPSP,
        0.86021 $/kg, is 3.835 $/GGE (quoted as 3.83); its IRR at market
        prices is 12.79 %.
      - **The eight** ``opt_*`` **scenarios** freeze the best trial (by each
        campaign's own objective) of the eight strain-design campaigns
        described in :doc:`kinetic_optimization`. Their pins are live-model
        reproductions on build 95e9362: every fermentation metric matches the
        recorded trial to < 1e-5, but TEA metrics differ from the
        campaign-recorded values (build 2b5b27d) by about 1–5 % — for example
        ``opt_PI_TRY_informed`` reads PI 0.8559 live vs 0.8083 recorded.
      - **Trial numbers** are the 0-based ``trial_number`` of the campaign
        CSVs. Some figure captions on :doc:`kinetic_optimization` count
        simulations from 1 instead, so #841 appears there as trial 842.

.. note::

   **Scenario B is a reference build, not a feasible strain.** Its
   pathway constants exceed the proteome cap (so it runs without the
   enzyme burden), and its workbook still sets the Adh6 rate constant
   ``k_17 = 44``, while the corrected model value is 0.1077, so Adh6
   (r17) runs about 400× too fast in B. Treat B's numbers as provisional
   until that row is corrected and re-pinned, and do not read its IRR as
   evidence that co-production pays; the feasible co-production result is
   ``opt_PI_TRY_informed`` (:doc:`kinetic_optimization`).


Run it
------

The case study needs the first NSKinetics release that ships
``nskinetics.biorefineries``, which is not on PyPI yet, plus optional
dependencies that plain NSKinetics does not (the Bioindustrial-Park
``biorefineries`` package, ``hensmith`` and a few analysis libraries). Until
that release is out, install it from source with the ``isobutanol`` extra,
then check the installation against scenario A's pins.
:ref:`isobutanol-installation` has the full instructions and the validated
dependency versions.

.. code-block:: bash

   $ git clone https://github.com/sarangbhagwat/nskinetics
   $ cd nskinetics
   $ git checkout dev
   $ pip install -e ".[isobutanol]"
   $ python -m nskinetics.biorefineries.isobutanol.smoke 1

Then build the flowsheet, load a scenario and solve its TEA:

.. code-block:: python

   from nskinetics.biorefineries import isobutanol

   isobutanol.load()                 # build + baseline-simulate (once per process)
   from nskinetics.biorefineries.isobutanol import scenarios

   bundle = scenarios.load_scenario('A')   # workbook kinetics, feeding strategy,
                                           # burden policy; simulates 2-3 times
   results = bundle['solve_TEA']()         # TEA on the current flowsheet state
   print(results)
   # {'IRR': ..., 'MPSPs': {'ethanol': ..., 'isobutanol': nan}}
   # ethanol MPSP: purity-adjusted $/kg at a 15 % IRR (pinned at 0.86021 for 'A')

Every change goes through the same simulate-then-solve loop. For example,
to turn scenario A into a pure batch at a 230 g/L glucose charge:

.. code-block:: python

   bundle['fbs_spec'].max_n_spikes = 0                # no glucose spikes
   bundle['model_specification'](target_conc=230.0)   # re-simulate the flowsheet
   print(bundle['solve_TEA']())

   V406 = bundle['V406']                              # the NSKinetics fermentor
   print(V406.tau, V406.nsk_results_specific_tau_dict['[s_EtOH]'])

Things to know before you script with it:

- :func:`isobutanol.load() <nskinetics.biorefineries.isobutanol.load>`
  builds the flowsheet and runs a full baseline simulation, so expect
  tens of seconds rather than seconds (a pinned smoke case — load plus three
  scenario simulations — takes about 2–4 minutes). It is idempotent: a second
  call in the same process is a no-op, and the separation configuration
  (``separation_processes``, a non-empty subset of
  ``('IBO_EtOH', 'ethanol')``) cannot be changed without a new process.
- ``load_scenario`` applies the workbook through one baseline ``Model``
  evaluation and then runs one baseline ``model_specification``; the first
  call with the enzyme burden on also evaluates scenario A once, to
  calibrate the burden. It mutates the shared model state, so to reproduce
  a pin exactly, load **one scenario per fresh process**, as the smoke
  cases do.
- ``bundle`` also carries ``model`` (the biosteam uncertainty
  ``Model``), ``tea``, ``fbs_spec`` (the fed-batch strategy specification),
  ``namespace_dict``, ``feeding_kwargs``, ``expected`` (the pins), the
  active ``burden_model`` and ``burden_on``.
- Run simulating scripts **one at a time**: concurrent processes share an
  on-disk numba cache and can corrupt it.
- The analysis scripts write their outputs under one results root
  (``analyses/results`` in the package, or the folder named by the
  ``NSKINETICS_IBO_RESULTS_DIR`` environment variable). No analysis results
  are distributed with NSKinetics: regenerate them, or point the variable at
  a folder that already holds them. :doc:`running_analyses` lists every
  script, what it writes and how long it takes.


.. _isobutanol-conventions:

Conventions used on these pages
-------------------------------

.. dropdown:: Units, metrics, build tags and trial numbers
   :open:

   - **MPSP / MESP.** MPSP is the minimum product selling price at a 15 %
     IRR, solved for one product while the other is sold at its market price,
     and *purity-adjusted* ($ per kg of the product's main chemical). MESP is
     the ethanol MPSP in $ per gallon of gasoline equivalent (GGE):

     .. math::

        \mathrm{MESP}\ [\$/\mathrm{GGE}]
        &= \mathrm{MPSP}\ [\$/\mathrm{kg}] \times \frac{0.789 \times 3.785411784}{0.67} \\
        &= 4.45775 \times \mathrm{MPSP}\ [\$/\mathrm{kg}]

     (ethanol density 0.789 kg/L, 3.785 L/gal, 0.67 GGE per gallon of
     ethanol). Because each MPSP is a break-even *given the other product's
     revenue*, MESP is not a stable headline for co-producing designs; those
     pages lead with IRR and break-even prices instead.
   - **IRR** is the internal rate of return at market prices (ethanol
     0.835 $/kg, isobutanol 1.500 $/kg, 2023$). A design that loses money
     outright reports IRR = −inf.
   - **PI** (profitability index) = NPV at the 15 % hurdle ÷ total capital
     investment (TCI). It is dimensionless, breaks even at 0 (PI > 0 ⇔
     IRR > 15 %), and ranks designs exactly as IRR does wherever an IRR
     exists — but not as MESP once isobutanol is co-produced.
   - **Basis.** Titer, cell density and productivity are per litre of
     *water*, not broth; yield is g product per g glucose *added*;
     productivity = titer / τ.
   - **Model versions ("builds").** The analyses ran over several weeks on
     three versions of the model that differ only in the
     heat-exchanger-network synthesis package,
     `hensmith <https://github.com/BioSTEAMDevelopmentGroup/hensmith>`__,
     used for ``HXN1001``. Fermentation results are identical across them;
     TEA results shift by up to a few percent. Each number on these pages
     names its version by a short tag, and **the tag is the hensmith git
     commit** it ran on; the validated stack's hensmith clone is at
     ``95e9362``. To re-run an analysis on its own build, see
     :ref:`isobutanol-switch-build`.

     .. list-table::
        :header-rows: 1
        :widths: 14 40 46

        * - Version
          - Scenario A on this version
          - Used by
        * - 2b5b27d
          - 0.86233 $/kg = 3.844 $/GGE, IRR 12.60 %, TCI 139.3 MM$
          - every optimization campaign, the TEA breakdowns, the 2-D kinetic
            sweeps, the Sobol' runs, the objective-landscape 1-D slices; the
            feeding grid ran on code equivalent to it (2b5b27d-equivalent)
        * - 6d4776f
          - 0.86021 $/kg = 3.835 $/GGE, IRR 12.79 %, TCI 139.09 MM$
          - the 6,000-sample uncertainty run, the baseline trajectory figure,
            the A and B smoke-test pins
        * - 95e9362
          - —
          - the live ``opt_*`` scenario pins (TEA metrics 1–5 % off the
            2b5b27d campaign record)

     Campaign IRRs (build 2b5b27d) and the uncertainty run's MESPs (build
     6d4776f) therefore come from different builds and are not directly
     comparable, which is why the pages never set one beside the other.
   - **Trial numbers.** ``#842`` is a 0-based campaign-CSV
     ``trial_number``, the form used in prose, in ``scenarios.py`` and in
     the objective-landscape figures. The two campaign figures on
     :doc:`kinetic_optimization` count simulations from 1 instead
     ("trial 843", "simulation 939"), so the same design appears there one
     number higher (for the relay, simulation = ``trial_number`` − 999).
   - **Reproducibility margin.** The objective is deterministic, but a
     design re-simulated in a fresh process can land slightly differently
     because the flowsheet converges from a different starting state. Each
     page that compares close PI values states which margin it uses
     (0.0105 PI cross-process, or 0.0075 PI between campaigns).


.. _isobutanol-references:

References
----------

The base kinetic model, the isobutanol benchmarks behind the scope of the
co-production result, and the sources of the enzyme-burden constants, as
cited on the case-study pages.

.. _ibo-ref-jach2022:

Jach, M. E., Serefko, A., Ziaja, M. & Kieliszek, M. (2022). Yeast protein as
an easily accessible food source. *Metabolites* **12**, 63.
https://doi.org/10.3390/metabo12010063

.. _ibo-ref-kondo2012:

Kondo, T., Tezuka, H., Ishii, J., Matsuda, F., Ogino, C. & Kondo, A.
(2012). Genetic engineering to enhance the Ehrlich pathway and alter carbon
flux for increased isobutanol production from glucose by *Saccharomyces
cerevisiae*. *J. Biotechnol.* **159**, 32–37.
https://doi.org/10.1016/j.jbiotec.2012.01.022

.. _ibo-ref-lane2020:

Lane, S., Zhang, Y., Yun, E. J., Ziolkowski, L., Zhang, G., Jin, Y.-S. &
Avalos, J. L. (2020). Xylose assimilation enhances the production of
isobutanol in engineered *Saccharomyces cerevisiae*. *Biotechnol. Bioeng.*
**117**, 372–381. https://doi.org/10.1002/bit.27202

.. _ibo-ref-lei2001:

Lei, F., Rotbøll, M. & Jørgensen, S. B. (2001). A biochemically structured
model for *Saccharomyces cerevisiae*. *J. Biotechnol.* **88**, 205–221.
https://doi.org/10.1016/S0168-1656(01)00269-3

.. _ibo-ref-metzlraz2017:

Metzl-Raz, E., Kafri, M., Yaakov, G., Soifer, I., Gurvich, Y. & Barkai, N.
(2017). Principles of cellular resource allocation revealed by
condition-dependent proteome profiling. *eLife* **6**, e28034.
https://doi.org/10.7554/eLife.28034

.. _ibo-ref-sanchez2017:

Sánchez, B. J., Zhang, C., Nilsson, A., Lahtvee, P.-J., Kerkhoven, E. J. &
Nielsen, J. (2017). Improving the phenotype predictions of a yeast
genome-scale metabolic model by incorporating enzymatic constraints.
*Mol. Syst. Biol.* **13**, 935. https://doi.org/10.15252/msb.20167411

.. _ibo-ref-scott2010:

Scott, M., Gunderson, C. W., Mateescu, E. M., Zhang, Z. & Hwa, T. (2010).
Interdependence of cell growth and gene expression: origins and
consequences. *Science* **330**, 1099–1102.
https://doi.org/10.1126/science.1192588

.. _ibo-ref-su2021:

Su, Y., Shao, W., Zhang, A. & Zhang, W. (2021). Improving isobutanol
tolerance and titers through EMS mutagenesis in *Saccharomyces cerevisiae*.
*FEMS Yeast Res.* **21**, foab012. https://doi.org/10.1093/femsyr/foab012

.. _ibo-ref-wess2019:

Wess, J., Brinek, M. & Boles, E. (2019). Improving isobutanol production
with the yeast *Saccharomyces cerevisiae* by successively blocking competing
metabolic pathways as well as ethanol and glycerol formation.
*Biotechnol. Biofuels* **12**, 173.
https://doi.org/10.1186/s13068-019-1486-8
