Kinetic optimization
====================

Strain design can be posed as an optimization over kinetic parameters. Each
candidate strain is a set of enzyme capacities and product-tolerance
strengths. Paired with a fed-batch feeding policy, it runs through the fermentor
``V406`` and then the whole biorefinery, and the result is scored. This page
asks two questions. First, if a strain program optimizes **titer, rate or yield
(TRY)** directly, does it reach the design that **earns the most**? Second, can
those TRY campaigns help a profitability search? To answer them,
``nskinetics.biorefineries.isobutanol`` ran Bayesian optimization campaigns over
the same 12 strain and feeding decisions under different objectives, and froze
each campaign's best design as a reproducible ``opt_*`` scenario.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 11 of 12
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      TRY optima (six objectives × two seeds) failed the 15 % IRR hurdle.
      None exceeded 15.5 %, against the starting strain's 12.6 %.

   .. grid-item-card:: 1.7 → 22.9 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      IRR across two **consecutive** trials (#841 → #842) of the
      isobutanol-yield campaign, for 7 % less yield. The campaign returned
      the 1.7 % design.

   .. grid-item-card:: 15.8 vs 27.3 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Final IRR of two profitability campaigns that differed **only in random
      seed**. One froze on an ethanol-only design for 95 % of its budget.

   .. grid-item-card:: 52 and 67
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Simulations that two TRY-informed relay campaigns needed to reach
      PI ≥ 0.75. The successful unseeded campaign needed 939; the other never
      got there.

.. admonition:: In one sentence
   :class: tip

   TRY optima miss the profit optimum, a profitability search alone is
   seed-dependent, and seeding it with the TRY campaigns' trials makes it
   fast and repeatable. The best design found co-produces about 40 g/L each of
   ethanol and isobutanol with Pdc knocked down, not deleted. It returns
   27.3 % IRR against 12.6 % at essentially the same capital (137.7 vs
   139.3 MM$). That design is a **prospective, bounds-limited target**; see
   :ref:`the scope of the co-production result <ibo-opt-scope>`.

.. _ibo-opt-conventions:

.. note:: **Conventions on this page**

   - **Model build.** Every campaign number (IRR, PI, TCI, titers, Φ\ :sub:`M`)
     is the value the campaign recorded under the hensmith heat-exchanger-network
     build **2b5b27d**. Under that build the starting strain (scenario A) has
     IRR 12.6 %; it is 12.8 % under build 6d4776f, which the current scenario-A
     pin uses. The ``opt_*`` pins in
     :ref:`ibo-opt-scenarios` are live-model reproductions (build **95e9362**).
     Their fermentation metrics match the campaign records to < 10\ :sup:`−5`,
     but their TEA metrics differ by 1–5 %. The two sets are never mixed in one
     statement.
   - **Trial numbers.** ``#841`` is the 0-based ``trial_number`` of the
     trajectory CSVs, as used in ``scenarios.py``. The two
     figures, and phrases like "simulation 52" or "trial 18", use the 1-based
     **simulated index**: trial = ``trial_number`` + 1, or ``trial_number`` − 999
     for the TRY-informed relay. So #841/#842 are trials 842/843 in the
     figures, and #1912 is trial 913.
   - **Reproducibility margin.** The objective is deterministic, but it drifts
     slightly with the converged state each simulation starts from
     (load-path drift). This page treats PI differences below
     **0.0105** as a tie. That value is the cross-process load-path drift
     recorded by the first relay campaign.
   - **Basis.** Titers and productivities are per litre of **water**, not broth.
     Yields are g per g of glucose added. IRR is at default prices (ethanol
     0.835 $/kg, isobutanol 1.500 $/kg, 2023$).

At a glance
-----------

.. list-table::
   :widths: 26 74
   :stub-columns: 1

   * - Question
     - Do TRY optima and the profit optimum agree, and can TRY campaigns speed
       up the profitability search?
   * - Decision space
     - 12 variables: four enzyme capacities, two capacity groups, three
       inhibition-strength multipliers and three feeding variables
       (preset ``ethanol_isobutanol`` × ``metabolic_split_12d``, scenario-A
       anchored)
   * - Constraints
     - Enzyme burden (Φ\ :sub:`M` < F\ :sub:`flex` = 0.245 g/gDCW; see
       :doc:`enzyme_burden`) and a fed-batch volume-ratio cap of 20
   * - Objectives
     - PI (log-tail) for profitability; ethanol and isobutanol titer, yield
       and productivity for the six TRY "scouts"
   * - Engine
     - Optuna Gaussian process (feasibility-aware, Matérn-5/2 ARD, log-EI),
       2,000 trials per campaign (relays: 1,000 preloaded + 1,000 simulated)
   * - Outputs
     - One trajectory CSV per campaign (not shipped); eight ``opt_*``
       scenarios; the two figures below
   * - Model build
     - 2b5b27d for every campaign number; 95e9362 for the ``opt_*`` pins

Method
------

The optimizer
~~~~~~~~~~~~~

The optimizer is the
:mod:`~nskinetics.biorefineries.isobutanol.kinetic_optimization` module
(``ko``). It is driven by ``analyses/kinetic_optimization/optimize_kinetics_BO.py``
and launched in production through the crash- and stall-resilient supervisor
``optimize_kinetics_BO_supervised.py`` in the same folder. Three
engines, ``OPTIMIZATION_METHODS = ('tpe', 'gp', 'dual_annealing')``, share one
evaluation site, ``evaluate_decision_point``. Every trial goes through the
same steps: burden prune → volume prune → simulate → a ``COMPLETE``, ``NAN``
or ``FAIL`` row in the trajectory CSV. Infeasible proposals are logged as
``INFEASIBLE`` and pruned **before** any simulation.

.. tab-set::

   .. tab-item:: GP (the publication campaigns)

      Optuna's Gaussian-process sampler, wrapped as the feasibility-aware
      ``FeasibleGPSampler``. It uses a Matérn-5/2 kernel with automatic
      relevance determination and a log expected-improvement acquisition.
      Start-up points come from a feasibility-filtered Latin-hypercube design
      by default. An infeasible GP proposal is swapped for the best feasible
      point among quasi-random candidates. The GP accepts at most
      ``GP_MAX_DIMENSIONS = 15`` decision variables, so it runs only on the
      compact presets (``metabolic_minimal_subset``, ``metabolic_14d`` and the
      split presets). The publication campaigns used a 50-point LHS start-up,
      a deterministic objective and **no learned constraints**.

   .. tab-item:: TPE (the default)

      Optuna's Tree-structured Parzen Estimator in a resumable SQLite study.
      A feasibility-aware variant (``FeasibleTPESampler``) checks every
      start-up draw and TPE candidate against the constraints before
      proposing it, so no sampled trial is logged ``INFEASIBLE``.

   .. tab-item:: Dual annealing

      SciPy's ``dual_annealing`` over the unit cube of the search space. The
      energy is the objective divided by a per-objective ``energy_scale``,
      with a finite penalty for failed or infeasible points. A restart resumes
      from the best ``COMPLETE`` row of the CSV. There is no SQLite store.

The decision space
~~~~~~~~~~~~~~~~~~

Search spaces come from named **presets** on two axes:

- ``study_target_products`` ∈ {``ethanol_only``, ``ethanol_isobutanol``};
- ``study_type`` ∈ {``metabolic``, ``metabolic_protein`` (default),
  ``metabolic_minimal``, ``metabolic_minimal_subset``, ``metabolic_14d``,
  ``metabolic_split_14d``, ``metabolic_split_12d``}.

Every preset starts at the scenario-A (wild-type) kinetics. The publication
campaigns use ``ethanol_isobutanol`` × ``metabolic_split_12d``, anchored on
scenario A (study-name tag ``_aA``), with the enzyme-burden constraint on. Its
12 variables:

.. dropdown:: The 12 decision variables of ``metabolic_split_12d`` and their campaign bands
   :open:

   .. list-table::
      :header-rows: 1
      :widths: 5 27 40 28

      * - No.
        - Variable
        - Meaning
        - Band
      * - 1
        - ``k_3``
        - Pdc capacity (pyruvate → acetaldehyde)
        - log, 10\ :sup:`−3`\ ×–4× A (A = 5.81 g/L/h)
      * - 2
        - ``k_6``
        - Adh1 capacity (acetaldehyde → ethanol)
        - log, 10\ :sup:`−3`\ ×–4× A (A = 2.82)
      * - 3
        - ``k_13``
        - ALS (Ilv2/Ilv6), the branch point into isobutanol
        - **absolute**, log 10\ :sup:`−3`–4.0 g/L/h (A = 0)
      * - 4
        - ``k_17``
        - Adh6, the isobutyraldehyde → isobutanol step
        - log, 10\ :sup:`−3`\ ×–20× A (A = 0.1077; cap 2.154)
      * - 5
        - ``glycolysis``
        - one multiplier on ``k_1l``, ``k_1h`` and ``k_1e``
        - log, 0.2×–4×
      * - 6
        - ``ehrlich_downstream``
        - ``k_14`` (Ilv5), absolute, with ``k_15`` (Ilv3) = 1.015× and
          ``k_16`` (Aro10) = 0.879× of it (stoichiometric ratios)
        - log, 10\ :sup:`−3`–4.0 g/L/h
      * - 7–9
        - ``inhib_ethanol``, ``inhib_isobutanol``, ``inhib_acetate``
        - one multiplier per effector family on every inhibition coefficient
          of that product (lower = more tolerant)
        - log, 0.75×–1.5×
      * - 10–12
        - ``threshold_conc``, ``target_delta``, ``max_n_spikes``
        - fed-batch feeding: spike-trigger threshold, target − threshold,
          spike **cap**
        - threshold 0–295 g/L; target − threshold 5–500 g/L (target
          ≤ 300 g/L); cap 0–50

   The spike feed is pinned at 600 g/L. ALS and the downstream Ehrlich group
   are absolute capacities because the wild-type strain has neither enzyme;
   trial 0 clips scenario A's zero Ehrlich rates up to the band floor. The
   ethanol-only twin of the preset has 8 variables (no isobutanol pathway).

Two **feasibility constraints** are checked before each simulation:

- **Enzyme burden.** Every capacity costs proteome. A design whose modelled
  metabolic pool reaches the flexible protein sector (Φ\ :sub:`M` ≥
  F\ :sub:`flex` = 0.245 g/gDCW) is infeasible. Below that, growth is derated
  once Φ\ :sub:`M` passes the penalty-free budget of 0.1347 g/gDCW. See
  :doc:`enzyme_burden`.
- **Fed-batch volume ratio.** A closed-form upper bound on the volume growth
  from spiking,
  ((600 − threshold)/(600 − target))\ :sup:`max_n_spikes` ≤ 20. It never
  admits a blow-up, but it may prune a high-cap point that the runtime guard
  would have accepted.

Objectives
~~~~~~~~~~

``OBJECTIVE_REGISTRY`` lists every objective a campaign can maximize or
minimize: IRR, PI, PI (log-tail), ethanol and isobutanol MPSP, titers, yields,
productivities, TCI, cell density, isobutanol yield × titer, and broth-basis
twins of the titers and cell density. The profitability campaigns optimize the
**profitability index**:

.. math::

   \mathrm{PI} = \frac{\mathrm{NPV}_{15\,\%}}{\mathrm{TCI}},
   \qquad
   \mathrm{PI}_{\text{log-tail}} =
   \begin{cases}
     \mathrm{PI} & \mathrm{PI} \ge 0 \\
     -\log(1 - \mathrm{PI}) & \mathrm{PI} < 0
   \end{cases}

where NPV\ :sub:`15 %` is the net present value at the 15 % hurdle rate and
default prices. Why not IRR?

- IRR is −∞ over a dead zone that covers about **94 %** of the domain, and an
  MPSP is NaN whenever its product has zero flow. Neither gives the optimizer a
  gradient to follow.
- PI is defined at every simulated point. It is dimensionless, breaks even at
  0, and **PI > 0 exactly when IRR > 15 %**.
- Wherever a genuine IRR exists, PI ranks designs exactly as IRR does: Spearman
  1.000 over 13,494 rows, with 0 discordant pairs.
- PI is **not** rank-equivalent to ethanol MESP once isobutanol is co-produced.
  Along the ``k_3`` slice, 11 % of pairs are discordant.
- The log-tail form is C\ :sup:`1`-continuous and strictly increasing, so its
  argmax is PI's. It only stops one very negative PI from swamping the GP's
  target standardization.

The campaigns
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Set
     - Label
     - What it is
   * - ``_rs350`` (2026-09-23)
     - uninformed + six TRY scouts
     - Seven GP campaigns, seed 350, 2,000 trials each. Objectives: PI
       (log-tail), and isobutanol and ethanol yield, titer and productivity.
       The 50-point LHS start-up is shared by the seven campaigns of one seed.
   * - Relay ``_rl15c111dc`` (2026-09-24)
     - TRY-informed
     - A fresh GP study maximizing PI (log-tail), seed 350. It was preloaded
       with 1,000 of the six scouts' 12,000 trials, which were not
       re-simulated, and then ran 1,000 simulated trials. The preload exceeded
       the start-up count, so there was no LHS.
   * - Seed 33960 (2026-09-16/17)
     - reference uninformed
     - The same seven objectives plus a price-weighted-yield campaign, on a
       different seed. Its PI campaign is the second uninformed campaign
       (best #1602).
   * - First relay ``_rlba1b2315`` (2026-09-23)
     - TRY-informed (2)
     - Seed 20260923, preloaded from the seven seed-33960 process-level
       campaigns.

.. dropdown:: How a relay is preloaded

   Every donor trial already recorded ``PI`` and ``PI (log-tail)`` as tracked
   metrics. ``ko.select_relay_rows`` therefore picks rows, and
   ``ko.relay_frozen_trials`` inserts them into the optuna store as
   ``COMPLETE`` trials that carry the donor's recorded value, with no
   re-simulation. For ``_rl15c111dc`` the 12,000 scout rows were filtered in
   four steps:

   1. Drop 32 ``FAIL`` and 214 unconverged rows, leaving 11,754.
   2. De-duplicate 274, leaving 11,480 unique rows.
   3. Keep all **327** rows with PI ≥ PI\ :sub:`A`, the starting strain's PI.
   4. Fill to 1,000 with **673** maximin (space-filling) rows.

   The keep step selects on PI, so the seed set is **not** blind to
   profitability. Relays run only under the GP engine. Optuna's TPE would
   class preloaded trials without a constraint record as infeasible, so they
   would become its "bad" density. The relay study name carries ``_rl`` plus
   a hash of the donor set and relay settings. ``n_trials`` counts simulated
   trials only.

.. dropdown:: The supervisor: crash and stall resilience

   ``optimize_kinetics_BO_supervised.py`` is stdlib-only. It relaunches the
   resumable study after the known CVODE segfault. It kills a child that logs
   no trial for ``--stall-timeout-min`` minutes (default 3), because a hung
   native integrator call cannot be interrupted in-process. The killed
   in-flight trial is logged at once as a ``LOST`` row, recovered from the
   ``<study>_inflight.json`` sidecar the child writes before each simulation,
   so there is one row per budgeted trial and no gaps. Under the
   **empty-attempt rule**, an attempt that logs no new trial is resumed if its
   child had started a simulation. It aborts at once if the child never
   reached one, and after ``--max-empty-attempts`` (default 5) consecutive
   empty attempts in any case. Separately, the evaluation site restores the
   flowsheet to the last converged snapshot after every ``FAIL`` trial, so a
   failure cannot contaminate the next trial.

Results
-------

The two figures below summarize the eight seed-350 campaigns. The text adds
the second seed set and the first relay where a claim needs them.

.. figure:: /_static/images/isobutanol/optimization/kinbo_main_scouts_campaigns.png
   :target: ../../_static/images/isobutanol/optimization/kinbo_main_scouts_campaigns.png
   :class: pub-figure
   :alt: Six-panel figure. (a) Best IRR so far against simulated trials on a
         log axis: the TRY-informed campaign jumps to 19.8 % on trial 2 and
         climbs to 27.3 %, while the uninformed campaign plateaus at 15.8 %
         from trial 104. (b) IRR against the isobutanol share of alcohol titer:
         every design above 15.8 % makes isobutanol. (c) Per-campaign strips
         of IRR, with open best-visit and filled returned-design markers; the
         isobutanol scouts visit 19–23 % but return low or no IRR. (d) Titers
         and (e) metabolic proteome of the returned designs, with the
         penalty-free budget line. (f) Share of trials making isobutanol
         against IRR.
   :width: 100%

   **Isobutanol TRY scouts visit designs above the 15.8 % ethanol-only plateau
   of the uninformed profitability campaign but return low-IRR designs.
   Seeded with their trials, a profitability search co-produces above that
   plateau on trial 2, beats every seed by trial 13 and reaches 27.3 %.**
   Eight GP campaigns (seed 350) searched the same 12-d space. Two maximized
   profitability (PI, plotted as IRR): **uninformed** (cyan, 2,000 trials) and
   **TRY-informed** (teal, 1,000 simulated trials after 1,000 preloaded
   seeds). Six TRY scouts each maximized one titer, rate or yield of ethanol
   (amber) or isobutanol (violet). Their 12,000 simulations supplied the seeds:
   327 scout trials at least as profitable as the starting strain, plus 673
   space-filling ones. Trials are 1-based simulated indices. Dashed cyan: the
   uninformed plateau (tinted above). Dotted grey: the starting strain,
   12.6 %. Dotted violet: the best seed, 22.9 %. "Loss": IRR < 0 or no IRR.
   "Making isobutanol" means ≥ 5 g/L; "co-production" means ≥ 5 g/L of each
   alcohol.
   **(a)** Best IRR so far; the left strip shows the seeds. In the TRY-informed
   campaign every improvement from trial 7 on co-produced. It reached 25.3 % on
   trial 18 and refined 26.2 % (trial 25) to 27.3 % (trial 913 = #1912). The
   uninformed campaign stayed at 15.8 % (ethanol only) from trial 104 (#103)
   to the end of its budget.
   **(b)** IRR against the isobutanol share of alcohol titer, for every design
   making ≥ 1 g/L alcohol. No design making < 5 g/L isobutanol exceeds 15.8 %.
   All 296 designs above 22.9 % are TRY-informed co-producers.
   **(c)** Each campaign's trials (strip), its highest-IRR visit (open marker)
   and the design it returned (filled; the argmax of its own objective). The
   isobutanol scouts visited 110, 1 and 1 designs above 15.8 % but returned
   1.7 % and two with no IRR. The ethanol scouts never exceeded 15.4 % and
   returned 8.1–12.7 %. **Look at the long horizontal gaps:** each one is a
   design that was visited but not returned.
   **(d, e)** Titers and metabolic proteome of the returned designs. Dotted:
   the penalty-free budget, 0.1347 g/(g DCW). ×: growth factor beyond it.
   Hatched: Adh1, Adh6.
   **(f)** Share of trials making isobutanol against best-visit (open) and
   returned (filled) IRR. The uninformed campaign made isobutanol in 4.9 % of
   its trials; the TRY-informed one in 65.9 %. The hexagon marks the seeds.
   IRRs are at default prices under build 2b5b27d.

TRY optima miss the profit optimum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first runs were the objectives a strain program would choose. **Eleven of
twelve** process-level optima (six objectives × two seeds) failed the 15 %
hurdle. None exceeded 15.5 % IRR; the maximum was 15.45 %, from the seed-33960
isobutanol-yield optimum. Each optimum fails through a different channel. The
seed-350 returned designs:

.. list-table:: Returned designs of the seed-350 campaigns (build 2b5b27d)
   :header-rows: 1
   :widths: 24 10 15 12 39

   * - Campaign (objective)
     - Trial
     - Objective
     - IRR
     - How it fails, or wins
   * - Starting strain (scenario A)
     - —
     - —
     - 12.6 %
     - TCI 139.3 MM$; ethanol 114.5 g/L, yield 0.460 g/g
   * - Isobutanol yield
     - #841
     - 0.3814 g/g
     - 1.7 %
     - 93 % of theoretical yield (0.411 g/g), but **2.3× the capital**
       (316 vs 139 MM$). Growth is derated to ×0.30 and cells fall to
       2.2 g/L.
   * - Isobutanol titer
     - #84
     - 46.0 g/L
     - none
     - Yield 0.15 g/g (37 % of theoretical). Annual sales fall below
       operating cost, so IRR = −∞.
   * - Isobutanol productivity
     - #1472
     - 2.06 g/L/h
     - none
     - Same channel: yield 0.15 g/g, sales below operating cost.
   * - Ethanol yield
     - #1923
     - 0.485 g/g
     - 12.1 %
     - Proteome overload (growth ×0.63, 4.2 g/L cells); MESP +0.9 %
   * - Ethanol titer
     - #1820
     - 168.5 g/L
     - 12.7 %
     - MESP −0.1 % despite 168 g/L
   * - Ethanol productivity
     - #365
     - 6.90 g/L/h
     - 8.1 %
     - MESP +5.9 % despite 6.9 g/L/h
   * - **Profitability, uninformed**
     - #103
     - PI 0.049
     - 15.8 %
     - Ethanol only; Pdc at its 4× bound, isobutanol pathway off
   * - **Profitability, TRY-informed**
     - #1912
     - PI 0.808
     - **27.3 %**
     - Co-produces 39.6 g/L isobutanol + 40.6 g/L ethanol at TCI 137.7 MM$

- **The yield optima overload the proteome to throttle growth.** They run at
  157 % (isobutanol) and 130 % (ethanol) of the penalty-free budget, and cell
  density falls to 2.2 and 4.2 g/L against 13.0. The dilute, slow broth that
  results is what costs the capital (:doc:`enzyme_burden`).
- **The 2.3× capital belongs to the campaign shown; the direction is
  systematic.** Across seeds, the isobutanol-yield optimum reproduced its
  yield but not its capital. The seed-33960 optimum reached 15.5 % IRR at a
  TCI of 245 MM$.
- **A revenue-weighted objective does not fix it.** One seed-33960 campaign
  maximized 1.5002·Y\ :sub:`IBO` + 0.835·Y\ :sub:`EtOH`. It still returned PI
  −0.93, IRR −6.0 % and TCI 347 MM$. This is a single seed.

Campaigns pass through better designs than they return
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**112 of 11,968** process-level trials beat the best design of the trapped
uninformed campaign, and 110 of them came from the isobutanol-yield campaign.
The sharpest case is a pair of consecutive trials:

.. list-table:: Isobutanol-yield campaign, consecutive trials
   :header-rows: 1
   :widths: 40 30 30

   * - Outcome
     - #841 (returned)
     - #842 (next trial)
   * - Isobutanol yield (g/g)
     - **0.381**
     - 0.354 (−7 %)
   * - IRR
     - 1.7 %
     - **22.9 %**
   * - TCI (MM$)
     - 316
     - 183
   * - Fermentation time τ (h)
     - 102
     - 31
   * - Metabolic proteome Φ\ :sub:`M` (g/gDCW)
     - 0.212
     - 0.173
   * - Growth factor
     - 0.30
     - 0.65

The yield gain at #841 was bought with proteome, and the derated culture is the
expensive one. The campaign scored only yield, so it ranked #841 first. A
campaign that reports only its best-objective strain discards its most
valuable observation. The lesson is to **score the economics of every trial
while it runs**.

The scouts do not, however, find the profit optimum themselves. **None** of the
23,944 process-level trials across both seeds reached PI 0.809. The best
process-level trial anywhere had PI 0.658 (IRR 25.1 %). None of the 11,968
seed-350 scout trials lies within 0.5 unit-cube units of the TRY-informed
optimum #1912.

Profitability search is seed-dependent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two uninformed PI campaigns with identical settings ended at **15.8 %** and
**27.3 %** IRR:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * -
     - Seed 350
     - Seed 33960
   * - Final best
     - #103: PI 0.049, ethanol only
     - #1602: PI 0.809, co-production
   * - Trajectory
     - Frozen from #103 on (95 % of the budget)
     - First passed PI 0.75 at #938
   * - Best point of the 50-point start-up
     - PI −0.975, ethanol only (144 g/L ethanol)
     - PI −0.299, co-production (28.0 g/L isobutanol + 62.1 g/L ethanol)
   * - Trials making isobutanol after start-up
     - 4.9 %
     - 39.2 %
   * - Trials near both the ALS and Adh6 caps
     - 2
     - 636

The divergence traces to **each seed's own** start-up sample. The six TRY
objectives, by contrast, reproduced across the two seeds to within 0.01–6 %.
The mechanism is conjunctive. Isobutanol pays only when ALS, the downstream
Ehrlich enzymes and Adh6 rise **together**, so the profitable region is narrow
and high-order inside a broad ethanol basin (see :doc:`objective_landscape`).
A black-box optimizer on an economic objective can therefore report
"ethanol-only is optimal" for stochastic reasons alone, at least under the
sampler settings used here (deterministic objective, no learned constraints).
Two seeds are an **existence proof, not a failure rate**.

Scouts make the search fast and repeatable, not higher
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two relay campaigns started from the TRY campaigns' trials, so each began
already holding its best donor row: PI 0.504 (IRR 22.9 %, the TRY-informed
relay's #842) and PI 0.658 (IRR 25.1 %, the first relay).

.. list-table:: First simulated trial at which best-so-far PI reaches a threshold
   :header-rows: 1
   :widths: 16 21 21 21 21

   * - PI
     - TRY-informed relay
     - First relay
     - Uninformed, seed 33960
     - Uninformed, seed 350
   * - 0.65
     - 18
     - 19
     - 104
     - never
   * - 0.75
     - **52**
     - **67**
     - 939
     - never
   * - 0.80
     - 913
     - never (0.7999)
     - 1603
     - never

- IRR ≥ 25 % came at simulation 18 and 19 for the two relays.
- The relays' finals, PI **0.808** and **0.800**, tie the reference uninformed
  campaign's **0.809** within the 0.0105 reproducibility margin. **The benefit
  is speed and reliability, not a higher ceiling.**
- The seeded route used 13,000 simulations (12,000 scout + 1,000). The
  uninformed campaign used 2,000.
- **No warm-start, random-search or single-donor control was run.** These
  comparisons therefore bound the value of reusing donor campaigns that were
  run anyway, not the minimum information needed. Of the 128 preloaded rows
  with PI > 0, 124 came from the isobutanol-yield donor. Only one preloaded
  row exceeded PI 0.5 (#842).
- Two 4,091-design Sobol' samples of the same space bracket the result. Under
  the campaigns' own sampling measure no design exceeded the hurdle (best
  13.7 %). Under a linear engineering prior 19 did, the best at 21.5 %, which
  is below the best seed (:doc:`sensitivity`).

The ceiling was set by the bounds, not the search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All **48** designs with PI ≥ 0.75 sit at four search bounds: the ALS cap
(4.0 g/L/h), the Adh6 cap (2.154 = 20× native), the isobutanol-inhibition floor
(0.75) and the target − threshold floor (5 g/L). #1912 sits at five bounds and
#1602 at six, so higher IRR may exist outside the searched ranges. Wide-band
``metabolic_14d`` campaigns reached PI 0.99, but in a different decision space
that is not comparable.

The flagship design inverts the textbook isobutanol strain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The TRY-informed optimum #1912 (scenario ``opt_PI_TRY_informed``) does not
delete Pdc and does not drop ethanol. It **knocks Pdc down** and keeps ethanol
as a co-product.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: +19 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Sales; revenue per kg of sugar +22 % (0.384 → 0.469 $/kg)

   .. grid-item-card:: 39.6 + 40.6 g/L
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Isobutanol + ethanol. Isobutanol supplies 55 % of sales.

   .. grid-item-card:: 0.33×
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Pdc capacity relative to wild type: knocked down, not deleted

   .. grid-item-card:: 47 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Share of the metabolic proteome in the isobutanol pathway
      (0.0607 of Φ\ :sub:`M` 0.130 g/gDCW)

Total alcohol yield **falls** to 0.403 g/g. Revenue still rises, because
isobutanol sells at 1.80× ethanol per kg; the economic optimum is not the
carbon-yield optimum. Capital is essentially unchanged (TCI 137.7 vs
139.3 MM$).

.. dropdown:: The flagship's full decision vector (#1912)

   .. list-table::
      :header-rows: 1
      :widths: 34 26 40

      * - Variable
        - Value
        - Note
      * - Pdc ``k_3``
        - 1.91 g/L/h
        - 0.33× wild type
      * - Adh1 ``k_6``
        - 7.44 g/L/h
        - 2.6× wild type
      * - ALS ``k_13``
        - 4.0 g/L/h
        - at the cap
      * - Adh6 ``k_17``
        - 2.154 g/L/h
        - at the cap (20× native)
      * - ``glycolysis``
        - 0.98×
        - baseline-like
      * - ``ehrlich_downstream``
        - 0.575 g/L/h
        - Ilv5 anchor; Ilv3 and Aro10 at fixed ratios
      * - ``inhib_ethanol``, ``inhib_isobutanol``
        - 0.75, 0.75
        - both at the floor (most tolerant allowed)
      * - ``inhib_acetate``
        - 1.02
        -
      * - Feeding
        - threshold 189.19 g/L, target 194.19 g/L
        - target − threshold at the 5 g/L floor; 1 spike (near-batch)
      * - Outcome
        - τ 40.8 h, cells 14.2 g/L
        - growth not derated (Φ\ :sub:`M` 97 % of the penalty-free budget)

- **Pdc is knocked down, not deleted.** Across the 48 designs with PI ≥ 0.75,
  which come from three campaigns, Pdc falls in **two clusters: 0.08–0.40×
  and 0.51–0.93×** wild type. Among 31,905 pooled designs, the best with Pdc
  at the 10\ :sup:`−3`\ × band floor reached only **PI 0.52** (IRR 23.1 %,
  TCI 169 MM$). On the ``k_3`` slice through #1912, wild-type Pdc still gives
  PI 0.653 (IRR 25.1 %). The knock-down is therefore worth about +0.16 PI
  (:doc:`objective_landscape`).
- **A ridge, not a point.** Ten designs from three campaigns have
  PI ≥ 0.79. ALS and Adh6 sit at their caps in all ten, with both tolerance
  multipliers at maximum, near-batch feeding and a Pdc knock-down. The
  isobutanol:ethanol titer ratio (0.88–3.98) and Adh1 (0.39–4.0×) are **free**.
  The product split is a design choice within the margin.
- **MESP is not a stable headline here.** Two equally profitable designs,
  #1912 and #1602 (both 2b5b27d), have ethanol MESPs of 2.12 and 1.60 $/GGE. The
  robust statement is a **break-even**. The design holds the 15 % hurdle
  while isobutanol sells above **1.13 $/kg** (25 % below the model price) with
  ethanol at 0.835 $/kg. Alternatively, it holds while ethanol sells above
  **0.475 $/kg** (below its 2021–2025 low of 0.518 $/kg) with isobutanol at
  1.50 $/kg.

.. figure:: /_static/images/isobutanol/optimization/kinbo_s2_design_fingerprints.png
   :target: ../../_static/images/isobutanol/optimization/kinbo_s2_design_fingerprints.png
   :class: pub-figure
   :alt: Heatmap with one row per design (starting strain, the returned design
         of each of the eight campaigns, and the isobutanol scouts' best
         visits) and twelve columns: native capacities as fold change, the
         isobutanol-pathway capacities in absolute units, three inhibition
         multipliers and three feeding variables. Outlined cells sit at a
         search bound. Right-hand columns give IRR, each campaign's own
         objective, the metabolic proteome and the growth factor.
   :width: 100%

   **Design fingerprints of the returned designs and of the isobutanol
   scouts' best visits.** Rows: the starting strain (scenario A); the design
   each campaign returned (filled markers); and, in italics with open markers,
   each isobutanol scout's highest-IRR trial. Trial numbers are 1-based
   simulated indices. **Native capacities** are fold change against the
   starting strain (Pdc ``k_3`` 5.81, Adh1 ``k_6`` 2.82 and Adh6 ``k_17``
   0.1077 g/L/h; glycolysis multiplier 1). Blue is lower and red higher, with
   the colour clipped at 1/16× and 16×. **Isobutanol-pathway capacities** (ALS
   ``k_13`` and the Ilv5/Ilv3/Aro10 group) are absolute, in g/L/h (purple).
   **Inhibition multipliers** run 0.75–1.5, where lower means a more tolerant
   strain (green). **Feeding** is shaded grey within its search range; spikes
   are shown as actual/cap. Cells with a black outer and white inner outline
   sit at a search bound. **Right:** IRR (bold = above the 15.8 % plateau),
   the campaign's own objective, Φ\ :sub:`M` (dark = above the 0.1347 budget)
   and the growth factor.
   **What to look at:** the TRY-informed row throttles the entry to the
   ethanol branch (Pdc 0.33×, Adh1 2.6×) at baseline-like glycolysis (0.98×).
   It puts ALS and Adh6 at their caps and stays within the proteome budget
   (Φ\ :sub:`M` 0.130, growth ×1.00). The uninformed design instead puts Pdc
   at its 4.0× bound and switches the isobutanol pathway off (Φ\ :sub:`M`
   0.153, growth ×0.84). The isobutanol-yield scout's returned design
   (trial 842 = #841) and best visit (trial 843 = #842) are consecutive
   trials. The 7.6 % higher yield it was rewarded for cost 21 IRR points
   (1.68 % vs 22.90 %): TCI 316 vs 183 MM$, batch 102 vs 31 h, growth ×0.30
   vs ×0.65. The inhibition multipliers carry no proteome cost, so the searches
   use them freely: 12 of the 24 in returned designs sit at the 0.75 floor. A
   multiplier on a product the design does not make is inert, so this block
   says less about strain design than the capacity columns do. IRRs are under
   build 2b5b27d (starting strain 12.6 %).

.. _ibo-opt-scope:

.. warning:: **Scope of the co-production result: read before quoting 27.3 %**

   - **Prospective titers.** 39.6 g/L isobutanol is about **15–20×** the best
     reported *S. cerevisiae* titers: 2.09 g/L from glucose (Wess et al. 2019)
     and 2.6 g/L from xylose (Lane et al. 2020). It is also **2.5×** the best
     engineered tolerance, 16 g/L (Su et al. 2021). The optimum ran at the
     most tolerant isobutanol setting allowed, so this result is
     **prospective target-setting**, not a near-term strain. A retained,
     partially knocked-down Pdc nevertheless inverts the academic standard of
     *PDC1*/*ADH1* deletion (Kondo et al. 2012; Wess et al. 2019).
   - **Bounds-limited.** All 48 designs with PI ≥ 0.75 sit at the ALS and
     Adh6 caps and at the isobutanol-tolerance and target − threshold floors.
     The ceiling reflects the chosen bands.
   - **One isobutanol price.** The model uses a single price, 1.500 $/kg
     (2023$), indexed from a 2025 product-listing price of 1.49 $/kg. No price
     sensitivity of the optimum has been run. Quote the break-even (1.13 $/kg
     isobutanol) alongside the IRR.
   - **Two seeds.** The seed dependence of the uninformed search and the speed
     of the relays come from single runs per setting. They are an existence
     proof, not a rate.

The economically optimal strain was neither the most productive, the most
tolerant nor the highest-yielding. It was the one that spent its proteome
where the facility earned the most per unit of sugar.

.. _ibo-opt-scenarios:

The eight ``opt_*`` scenarios
-----------------------------

Each campaign's returned design, the best trial by its **own** objective, is
frozen as a scenario in ``scenarios.SCENARIOS``. Each has its own parameter
workbook (``parameter-distributions_corn_IBO_EtOH_<scenario>.xlsx``), feeding
strategy and pinned MPSPs. They make the campaign optima reproducible
baselines for further analysis; for example, ``evaluate_EtOH_k13_k7ii.py``
anchors its sweep on ``opt_PI_TRY_informed`` by default. All eight run with
the enzyme burden **on**
(``ScenarioSpec.burden_default``) and are checked by smoke cases 9–16
(``python -m nskinetics.biorefineries.isobutanol.smoke 9`` … ``16``; see
:doc:`running_analyses`).

.. list-table:: ``opt_*`` scenario pins (live model, build 95e9362)
   :header-rows: 1
   :widths: 28 19 15 17 21

   * - Scenario
     - Source trial
     - Pinned objective
     - Feeding (threshold / target g/L; spike cap)
     - EtOH / IBO MPSP ($/kg)
   * - ``opt_PI_uninformed``
     - #103, ``_rs350`` PI campaign
     - PI 0.0761
     - 230.77 / 235.77; 43
     - 0.82091 / nan
   * - ``opt_PI_TRY_informed``
     - #1912, relay ``_rl15c111dc``
     - PI 0.8559
     - 189.19 / 194.19; 1
     - 0.46066 / 1.11234
   * - ``opt_IBO_titer``
     - #84
     - 46.04 g/L
     - 295 / 300; 46
     - 207.46\ :sup:`*` / 2.4625
   * - ``opt_IBO_yield``
     - #841
     - 0.3814 g/g
     - 35.06 / 40.06; 21
     - 5760.9\ :sup:`*` / 1.8533
   * - ``opt_IBO_productivity``
     - #1472
     - 2.063 g/L/h
     - 295 / 300; 21
     - 491.16\ :sup:`*` / 2.2694
   * - ``opt_EtOH_titer``
     - #1820
     - 168.45 g/L
     - 295 / 300; 35
     - 0.85672 / nan
   * - ``opt_EtOH_yield``
     - #1923
     - 0.4847 g/g
     - 85.18 / 107.54; 8
     - 0.86838 / nan
   * - ``opt_EtOH_productivity``
     - #365
     - 6.897 g/L/h
     - 295 / 300; 15
     - 0.90813 / nan

:sup:`*` The ethanol MPSPs of the isobutanol scenarios are purity-adjusted
MPSP artifacts of a near-zero co-product, not bugs. ``nan`` marks a product
that is not made.

.. important::

   The pins are **live-model** values (build 95e9362), not the campaign
   records. Fermentation metrics reproduce the recorded trials to
   < 10\ :sup:`−5`, but TEA metrics move by 1–5 %. For example, PI moves
   from 0.8083 (campaign) to **0.8559** (pin) for ``opt_PI_TRY_informed``, and
   from 0.0493 to **0.0761** for ``opt_PI_uninformed``. Quote campaign numbers
   from the campaign, scenario numbers from the scenario, and name the one you
   use.

How the scenarios are built and verified:

.. tab-set::

   .. tab-item:: Build: ``build_opt_split12d_workbooks.py``

      One fresh process per scenario
      (``build_opt_split12d_workbooks.py <scenario>``; it simulates). The
      script:

      1. Builds the biorefinery (the campaigns' build: both trains, splitter
         at 1.0) and snapshots scenario A's full live kinetics and isobutanol
         price.
      2. Reproduces the campaign trial with ``ko.reproduce_split12d_trial``
         (mode ``'both'``, scenario-A-referenced burden on). It aborts on an
         error, a cross-check mismatch, or any **fermentation-level** metric
         off its recorded value. TEA-level warnings are printed and tolerated,
         because they move with the hensmith build.
      3. Merges the trial's applied kinetics over the scenario-A snapshot.
         ``k_7``/``k_8`` keep their *intended* values, because the burden
         derates them at simulation time exactly as during the campaign.
      4. Writes the workbook: a copy of the scenario-B workbook with every
         kinetic row's baseline replaced (triangular ±20 %). It also records a
         JSON reproduction under the results folder (``opt_scenario_reproductions``).
      5. Reloads the workbook, re-simulates, and aborts unless the MPSPs and
         **every** tracked metric match the reproduction within 10\ :sup:`−3`.
      6. Prints the ``ScenarioSpec`` block to paste into ``scenarios.py``.

   .. tab-item:: Verify: ``verify_opt_scenarios.py``

      ``verify_opt_scenarios.py <scenario>`` is read-only and simulates. In a
      fresh kernel it runs the scenario's baseline exactly as its smoke test
      does. It then resets the kinetics to scenario A's reference, because the
      A workbook cannot undo the opt_* ``k_13``–``k_17`` and
      isobutanol-inhibition values. Finally it reproduces the campaign trial
      and compares both MPSPs and every tracked metric at a relative tolerance
      of 10\ :sup:`−3`, with identical feeding. Exit 0 and ``MATCH`` mean the
      scenario's baseline **is** the reproduced trial. All eight scenarios
      returned ``MATCH`` on 2026-10-08 (max relative difference 3.9e-4).

To load one in a session:

.. code-block:: python

   from nskinetics.biorefineries import isobutanol

   isobutanol.load()                     # build + baseline-simulate once per process
   bundle = isobutanol.load_scenario('opt_PI_TRY_informed')   # burden on by default
   isobutanol.solve_TEA()

Caveats and limitations
-----------------------

- **Model builds.** The campaigns, the figures, the TEA breakdowns and every
  campaign number on this page are under build 2b5b27d.
  The ``opt_*`` pins are under 95e9362. The current scenario-A/B smoke pins
  are under 6d4776f. The build shift is larger at co-production designs:
  #1602's IRR moves from 27.30 % to 28.02 % between 2b5b27d and 6d4776f.
  Never place a campaign IRR beside an uncertainty-analysis MESP
  (:doc:`uncertainty`).
- **Single runs.** Every campaign is one run at one seed. The seed-dependence
  and relay-speed claims rest on two seeds and two relays, under the sampler
  settings used (deterministic objective, no learned constraints).
- **No compute-matched control.** No warm start from the single best donor,
  no random-seed preload and no single-donor ablation were run.
- **Reproducibility margin.** This page uses 0.0105 PI. The caption notes
  shipped with the figure scripts (``plots/optimization_figure_ideas``)
  compute a stricter, design-pair-based margin (0.0075 PI, from 49 pairs of
  identical designs simulated after different trial histories). Statements
  here do not use it.
- **Superseded items.** The ``opt_IRR`` scenario (#1602 of the seed-33960 PI
  campaign) was removed on 2026-10-08; use ``opt_PI_TRY_informed``. The
  ``metabolic_minimal_subset`` optima of 2026-09-07/08 (with ``k_17`` = 44
  workbooks) are removed too. Figures of the parameter-set plotter
  (``plots/plot_kin_opt_parameter_sets.py``) are not shown on this page.
- **Basis.** Titers and productivities are per litre of water, not broth.

Reproduce it
------------

.. important::

   No analysis results ship with NSKinetics: the campaign trajectory CSVs,
   optuna stores, relay manifests and reproduction JSONs must be
   regenerated, or the ``NSKINETICS_IBO_RESULTS_DIR`` environment variable
   pointed at a results folder that already holds them. The campaigns are long-running, simulating
   scripts and are **opt-in**. Exact campaign numbers reproduce only under
   the hensmith build the campaigns used (2b5b27d); under the current build,
   TEA metrics shift by 1–5 %.

.. list-table::
   :header-rows: 1
   :widths: 30 44 26

   * - Script
     - What it does / writes
     - Cost
   * - ``analyses/kinetic_optimization/optimize_kinetics_BO.py``
     - Driver: ``run(objective=..., study_type=..., n_trials=2000, ...)``.
       Writes ``<study>_trajectory.csv`` (one row per trial, with every
       tracked metric and the burden columns) and the optuna store. Rerunning
       a study name resumes it. Importing the file builds and
       baseline-simulates the biorefinery.
     - Long: 2,000 simulated trials per campaign (~8–17 h)
   * - ``analyses/kinetic_optimization/optimize_kinetics_BO_supervised.py``
     - The production launcher. Flags include ``--method gp``,
       ``--study-target-products``, ``--study-type``, ``--objective``,
       ``--seed``, ``--n-startup-trials`` and, for relays, ``--relay-from``.
       ``--help`` lists them all.
     - As above; survives crashes and stalls
   * - ``analyses/kinetic_optimization/reproduce_split12d_trial.py``
     - Re-simulates one recorded trial and compares it, metric by metric,
       with its trajectory row. It writes nothing.
     - One load plus one trial
   * - ``analyses/opt_scenarios/build_opt_split12d_workbooks.py <scenario>``
     - Writes one ``opt_*`` workbook and its reproduction JSON, and prints
       the registry pin
     - One fresh process per scenario
   * - ``analyses/opt_scenarios/verify_opt_scenarios.py <scenario>``
     - Read-only ``MATCH`` check of a scenario against its campaign trial
     - One fresh process per scenario
   * - ``plots/optimization_figure_ideas/``: ``check_facts.py``,
       ``fig_main_scouts.py``, ``fig_s2_fingerprints.py``
     - Sim-safe (pandas only). They check every figure number against the
       campaign CSVs, then draw the two figures above as
       ``kinBO_main_scouts_<stamp>`` / ``kinBO_S2_fingerprints_<stamp>``
       (``.png`` + ``.pdf``, plus ``*_latest`` copies) and write the caption
       files.
     - Seconds, but they need the eight campaign CSVs

Related pages: :doc:`objective_landscape` explains *why* profitability search
is seed-dependent; :doc:`enzyme_burden` gives the proteome constraint behind
the growth factors; :doc:`sensitivity` ranks the same 12 decisions by variance
share.
