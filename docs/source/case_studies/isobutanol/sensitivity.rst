Global sensitivity analysis
===========================

Twelve strain and feeding decisions shape the co-production biorefinery: six
enzyme capacities, three product-inhibition strengths and three fed-batch
feeding settings. This page asks which of them carry the variance in
**profitability**, and whether those are the same decisions that carry the
variance in **titer, rate and yield (TRY)**. The answer comes from a two-stage
variance-based global sensitivity analysis (GSA): Sobol' closed indices and
Shapley effects, computed on surrogates fitted to 4,096 simulated designs per
sampling measure. It runs over the same decision space that the
:doc:`kinetic_optimization` campaigns search. Read it as a *pre-optimization screen*.
It ranks levers by how much they move the typical design, which here is a
failing one; it cannot say which levers the rare profitable designs need.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 3.5 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Share of profitability variance carried by ALS (``k_13``) under the
      screening prior, against 36 % of ethanol-titer and 23 % of
      isobutanol-titer variance. ALS carries ≤ 3.5 % under every measure.

      +++
      Shapley effects, screening prior (SCR), build 2b5b27d

   .. grid-item-card:: 10–43 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Share of profitability variance carried by Adh1 (``k_6``) across the
      three measures. Adh1 carries ≤ 1.6 % of isobutanol titer or yield.

      +++
      Shapley effects, all three measures, build 2b5b27d

   .. grid-item-card:: 0.46 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Share of the 4,091 screening designs that clear the 15 % hurdle
      (PI > 0). Every one of them co-produces ethanol and isobutanol.

      +++
      Screening prior (SCR) sample, build 2b5b27d

   .. grid-item-card:: ≥ 95 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Share of Var(PI) carried by designs that do worse than the baseline,
      under every measure. The indices therefore rank **failure modes**.

      +++
      All three measures, build 2b5b27d

.. admonition:: The robust finding: titer's driver is not profit's driver
   :class: tip

   The branch-point enzyme ALS (``k_13``) is a top-three driver of isobutanol
   titer and yield under every sampling measure (15–29 % of their variance).
   It carries only 1.3–3.5 % of the variance in the profitability index, and
   its closed first-order index for PI is ≤ 0.9 %. Adh1 (``k_6``) is the mirror
   image: 10–43 % of Var(PI) but ≤ 1.6 % of isobutanol titer or yield. This
   divergence holds under **all three** sampling measures. The
   :ref:`proteome-sink ranking <ibo-gsa-measure-dependent>` below does
   **not**: it changes with the prior.

At a glance
-----------

.. list-table::
   :widths: 28 72
   :stub-columns: 1

   * - Question
     - Which of the 12 split-12d decisions carry the variance in profitability
       (PI), and are they the ones that carry the variance in TRY?
   * - Decision space
     - The 12 variables of the ``ethanol_isobutanol`` × ``metabolic_split_12d``
       campaigns, starting from scenario A, with the enzyme-burden and
       fed-batch volume constraints enforced
   * - Sampling
     - Scrambled Sobol' stream, rejection-sampled to the feasible set. Three
       input measures (campaign, screening, custom), 4,096 designs each
   * - Headline output
     - Raw PI = NPV at the 15 % hurdle ÷ TCI. IRR is never analysed.
   * - Estimators
     - Closed subset indices and exact Shapley effects on a surrogate (GP or
       gradient boosting, chosen by 5-fold CV Q²)
   * - Model build
     - The 2b5b27d-era model, the same as the optimization campaigns
       (scenario-A PI = −0.1383). See :ref:`ibo-gsa-caveats`.
   * - Scripts
     - ``analyses/sensitivity/evaluate_sobol_split12d.py`` (stage 1,
       simulations) → ``analyses/sensitivity/analyze_sobol_split12d.py``
       (stage 2, sim-safe analysis)

Method
------

The decision space
~~~~~~~~~~~~~~~~~~

The GSA samples the 12 decision variables of the split-12d kinetic
optimization. The spike feed is fixed at 600 g/L. ``stage_1_max_x`` is pinned
at 5 g/L. Every design starts from scenario A and uses the default
both-trains separation build. Four of the six enzyme axes scale a wild-type
(scenario-A) value. The other two, ALS (``k_13``) and the downstream Ehrlich
lever, are **absolute capacities**, because the wild-type strain has neither
enzyme.

.. dropdown:: The 12 decision variables and their bands under each measure

   .. list-table::
      :header-rows: 1
      :widths: 25 31 15 15 14

      * - Variable
        - What it sets
        - Campaign (CAM)
        - Screening (SCR)
        - Custom (CUS)
      * - ``k_3``
        - Pdc capacity (r3, pyruvate → acetaldehyde); A = 5.81
        - log 10⁻³×–4× A
        - log 0.1×–4× A
        - linear 0–4× A
      * - ``k_6``
        - Adh1 capacity (r6, acetaldehyde → ethanol); A = 2.82
        - log 10⁻³×–4× A
        - log 0.1×–4× A
        - linear 0–4× A
      * - ``k_13``
        - ALS (Ilv2/Ilv6) capacity (r13, pyruvate → acetolactate). A = 0, so
          the band is absolute
        - log 10⁻³–4.0
        - linear 0–4.0
        - linear 0–4.0
      * - ``k_17``
        - Adh6 capacity (r17, isobutyraldehyde → isobutanol); A = 0.1077
        - log 10⁻³×–20× A
        - linear 0–2.154 (0–20× A)
        - linear 0–2.154
      * - ``glycolysis``
        - One multiplier on the three r1 capacities ``k_1l``, ``k_1h`` and
          ``k_1e``
        - log 0.2×–4×
        - log 0.1×–4×
        - linear 0–4×
      * - ``ehrlich_downstream``
        - One lever that sets ``k_14`` (Ilv5) absolutely. ``k_15`` (Ilv3) =
          1.015× it and ``k_16`` (Aro10) = 0.8790× it, fixed stoichiometric
          weights. A = 0
        - log 10⁻³–4.0
        - linear 0–4.0
        - linear 0–4.0
      * - ``inhib_ethanol``
        - Multiplier on ``k_1ie``, ``k_4ie``, ``k_7ie``, ``k_10ie``,
          ``k_17ie``. **Higher = less tolerant**
        - log 0.75×–1.5×
        - log 0.75×–1.5×
        - linear 0.75×–1.5×
      * - ``inhib_isobutanol``
        - Multiplier on ``k_1ii``, ``k_4ii``, ``k_6ii``, ``k_7ii``, ``k_10ii``
        - log 0.75×–1.5×
        - log 0.75×–1.5×
        - linear 0.75×–1.5×
      * - ``inhib_acetate``
        - Multiplier on ``k_1ia``, ``k_4ia``, ``k_6ia``, ``k_7ia``,
          ``k_10ia``, ``k_17ia``. Acetate is a by-product, not a product
        - log 0.75×–1.5×
        - log 0.75×–1.5×
        - linear 0.75×–1.5×
      * - ``threshold_conc``
        - Sugar concentration that triggers a glucose spike
        - linear 0–295 g/L
        - linear 0–295 g/L
        - linear 0–295 g/L
      * - ``target_delta``
        - Target − threshold. The target is min(300, threshold + Δ) g/L
        - linear 5–500 g/L
        - linear 5–500 g/L
        - linear 5–500 g/L
      * - ``max_n_spikes``
        - **Cap** on the number of spikes, not the realized count
        - 0–50
        - 0–50
        - 0–50

**Feasibility.** Each design must satisfy two constraints, both shared with
the campaigns.

- **Enzyme burden** (see :doc:`enzyme_burden`). The modelled enzyme load Φ_M
  must stay below the flexible proteome F_flex = 0.245 g/gDCW. Growth is
  derated once Φ_M exceeds 0.135 g/gDCW. The wild type sits at 0.070.
- **Fed-batch volume.** The spikes may grow the working volume by at most a
  factor of 20:
  ((600 − threshold)/(600 − target))^max_n_spikes ≤ 20.

About two-thirds of the campaign box is infeasible. A classical Saltelli design
needs every point of each base sample, so it cannot be used here. Stage 1
instead **rejection-samples a scrambled Sobol' stream** down to the feasible
set, which leaves the design uniform over that set in the chosen measure.

Three sampling measures
~~~~~~~~~~~~~~~~~~~~~~~

Variance shares depend on how widely and on what scale each input is sampled.
The analysis was therefore run under three measures over the same 12
variables:

.. list-table::
   :header-rows: 1
   :widths: 14 38 10 16 22

   * - Tag
     - Measure
     - Seed
     - Rows complete / failed (feasible fraction of candidates)
     - PI surrogate (CV Q²); PI > 0 rows
   * - **CAM**
     - The campaign box itself, log-uniform. 10⁻³× native knock-downs sit
       inside it
     - 20260920
     - 4,091 / 5 (32.2 %)
     - gradient boosting 0.926; 0
   * - **SCR**
     - "Screening" engineering prior. ``k_13``, ``k_17`` and the downstream
       Ehrlich lever are linear on [0, ceiling]. ``k_3``, ``k_6`` and
       glycolysis are log-uniform on 0.1×–4×
     - 20260920
     - 4,091 / 5 (27.4 %)
     - Gaussian process 0.827; 19 (0.46 %)
   * - **CUS**
     - Custom, all-linear: every enzyme axis on 0–4× (pathway axes on
       [0, ceiling]). Bands from ``linear_0-4x_2026-09-22.json``, which
       ships in the package's ``data/sobol_bands``
     - 20260922
     - 4,095 / 1 (17.9 %)
     - gradient boosting 0.857; 15

The campaign measure answers a narrower question: what separates a dead strain
from a merely bad one? Its 10⁻³× native knock-downs *are* the failure tail.
Of its feasible designs, 89 % have ``k_3`` or ``k_6`` below 0.25× A, and 59 %
make less than 5 g/L ethanol. The co-production corner (``k_13`` ≥ 2,
``k_17`` ≥ 1, downstream Ehrlich ≥ 0.4, no native knock-down below 0.1×) holds
only 0.07 % of CAM rows, against 23.1 % (SCR) and 18.7 % (CUS). The two prior
measures (SCR and CUS) were built to give that corner real mass.

.. note::

   The GSA and the campaigns share the 12 variables but **not** the measure.
   The SCR bands differ from the campaign bands in three ways: the glycolysis
   floor is 0.1× instead of 0.2×; the pathway axes run linearly from zero; and
   there are no native knock-downs below 0.1×.

Surrogates and estimators
~~~~~~~~~~~~~~~~~~~~~~~~~

Scoring every subset of inputs needs a conditional expectation per subset.
That is affordable only on a surrogate fitted to the simulated rows. Stage 2
(``analyze_sobol_split12d.py``) works in five steps:

1. **Fit a surrogate per metric.** A Gaussian process (constant × Matérn-5/2
   with automatic relevance determination) competes against a
   ``HistGradientBoostingRegressor``, and 5-fold cross-validated Q² picks the
   winner. A metric is flagged unreliable below Q² = 0.8; its indices are
   still computed. By default only the two PI metrics are offered a GP,
   because GP prediction is expensive. Every index is an index of the
   surrogate, so Q² bounds how much of the true variance it carries.
2. **Closed subset indices.** For every subset *u* of the 12 inputs:

   .. math::

      S_u = \frac{\operatorname{Var}\!\left(\mathbb{E}[\,Y \mid X_u\,]\right)}{\operatorname{Var}(Y)}

   Each one is estimated by pick-freeze, with the partner points drawn by
   **conditional rejection sampling** so that they stay feasible.
3. **Exact Shapley effects.** These come from all 2¹² closed indices with the
   standard weights. They sum to 100 %, and they are the attribution to quote:

   .. math::

      \mathrm{Sh}_i = \sum_{u \subseteq D \setminus \{i\}}
      \frac{|u|!\,(d-|u|-1)!}{d!}\,\bigl(S_{u \cup \{i\}} - S_u\bigr),
      \qquad d = 12

4. **First-order and total indices, for reference only.** The feasibility
   constraint makes the inputs *dependent*, so S1 > ST is possible and a
   small S1 can come out negative (estimation error).
5. **Replicates.** Three replicate estimates give the means and standard
   deviations. Convergence is checked at 1,024, 2,048 and the full row
   count.

**Output metric.** The headline is raw **PI** (profitability index): NPV at
the 15 % hurdle divided by TCI. PI is dimensionless and breaks even at 0, and
PI > 0 exactly when IRR > 15 %. **IRR is never analysed.** It has no finite
value for 78.8–93.1 % of the sampled designs (depending on the measure), so
its variance is undefined. Wherever IRR is finite, PI ranks designs exactly as
IRR does: Spearman 1.0000 on the finite-IRR rows, and sign(PI) =
sign(IRR − 15 %) in all 12,277 completed rows. Ethanol MESP is not a tracked
column, so it is not analysed. Every other finite, non-constant tracked metric
(titers, yields, productivities, cell density, TCI) also gets indices.

Stage 2 also computes **target indices** on the indicator "PI > PI_A", which
asks what decides whether a design beats scenario A (PI_A = −0.1383). Where
enough rows lie above that threshold, it also computes **conditional
indices** on PI within that zone.

Results
-------

Titer's driver is not profit's driver (robust)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cleanest divergence is at the pyruvate branch point. The table gives
Shapley effects (% of each metric's variance) for the two divergent levers
under each measure:

.. list-table:: Shapley effects (%) of the two divergent levers, by measure
   :header-rows: 1
   :stub-columns: 1
   :widths: 40 20 20 20

   * - Lever → metric
     - CAM
     - SCR
     - CUS
   * - ALS ``k_13`` → **PI**
     - **2.6**
     - **3.5**
     - **1.3**
   * - ALS ``k_13`` → isobutanol titer / yield
     - 27.7 / 28.7
     - 23.1 / 22.6
     - 15.3 / 16.6
   * - ALS ``k_13`` → ethanol titer / yield
     - 1.1 / 1.3
     - 35.8 / 33.3
     - 13.8 / 11.3
   * - Adh1 ``k_6`` → **PI**
     - **43.1**
     - **10.2**
     - **12.2**
   * - Adh1 ``k_6`` → isobutanol titer / yield
     - 0.2 / 0.2
     - 0.8 / 0.2
     - 1.6 / 0.5

ALS's closed first-order index for PI is ≤ 0.9 % under every measure (SCR
0.009 ± 0.018). By closed first-order index, the variance a lever explains
on its own, ALS is the single best lever under SCR for ethanol titer and
yield and for isobutanol titer (0.321, 0.296 and 0.199). By Shapley effect,
which also credits interactions, it ranks first for ethanol titer and yield
but third for isobutanol titer (see the full profile below). The mechanism
is transparent in the raw samples. Under SCR,
from the lowest to the highest ``k_13`` quintile:

.. list-table:: Median outcomes by ALS (``k_13``) quintile, screening measure
   :header-rows: 1
   :stub-columns: 1
   :widths: 40 30 30
   :class: ibo-num ibo-compact

   * - Quantity
     - Lowest quintile (0–0.77)
     - Highest quintile (3.17–4.0)
   * - Ethanol titer (g/L)
     - 52.1
     - 6.6
   * - Isobutanol titer (g/L)
     - 5.3
     - 21.0
   * - Total alcohol (g/L)
     - 61.1
     - 29.6
   * - PI
     - −1.87
     - −2.30

ALS **reallocates carbon** between two products coupled through one enzyme
set. Each single-product metric sees one side of the trade, while PI sees only
the small net. The branch-point enzyme that dominates the product titers does
not decide profitability variance.

The same split appears in the closed-subset view under SCR. All six TRY
metrics share their best triple, {``k_13``, glycolysis, downstream Ehrlich}.
PI's best triple swaps ALS for Adh1: {``k_6``, glycolysis, downstream
Ehrlich}.

.. dropdown:: Full Shapley profile under the screening measure: PI vs TRY metrics

   Shapley effects (%), SCR, from gradient-boosting surrogates (a GP for PI).
   Only metrics with a reliable surrogate (Q² ≥ 0.8) are shown. The
   productivity and TCI surrogates fell below 0.8 and are omitted. The
   inhibition multipliers are not tabulated for the TRY metrics in the source,
   so the columns do not sum to 100.

   .. list-table::
      :header-rows: 1
      :stub-columns: 1
      :widths: 22 13 13 13 13 13 13

      * - Lever
        - PI
        - IBO titer
        - IBO yield
        - EtOH titer
        - EtOH yield
        - Cell density
      * - glycolysis
        - **40.0**
        - 27.1
        - 26.7
        - 27.1
        - 28.7
        - **74.7**
      * - downstream Ehrlich
        - **29.9**
        - **30.3**
        - 20.2
        - 22.6
        - 16.9
        - 19.0
      * - ``k_6`` (Adh1)
        - **10.2**
        - 0.8
        - 0.2
        - 3.1
        - 4.8
        - 0.3
      * - ``k_13`` (ALS)
        - 3.5
        - 23.1
        - 22.6
        - **35.8**
        - **33.3**
        - 1.1
      * - ``k_3`` (Pdc)
        - 3.3
        - 9.8
        - 8.7
        - 6.7
        - 7.6
        - 0.7
      * - ``k_17`` (Adh6)
        - 2.4
        - 4.8
        - 4.2
        - 2.3
        - 1.1
        - 1.9
      * - threshold / Δ / cap
        - 3.6 / 3.5 / 1.5
        - 0.8 / 0.3 / 0.0
        - 7.7 / 6.9 / 1.5
        - 0.5 / 0.2 / 0.1
        - 2.8 / 3.0 / 1.1
        - —
      * - surrogate Q²
        - 0.827
        - 0.898
        - 0.852
        - 0.880
        - 0.845
        - 0.865

   The raw data tell the same story without a surrogate. Under SCR, the
   design with the highest isobutanol titer (41.7 g/L) has PI −1.25. The one
   with the highest ethanol titer (142.9 g/L) has PI −0.81. The design with
   the highest cell density (309 g/L) has the lowest PI of all 4,091
   (−5.432). The profitable designs are high-yield on **both** products, not
   TRY-maximal on any one.

.. _ibo-gsa-measure-dependent:

Which enzyme leads depends on the prior (measure-dependent)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ranking of the **profitability** levers changes with the sampling
measure. This page shows all three side by side instead of picking one.

.. grid:: 1 1 3 3
   :gutter: 3

   .. grid-item-card:: 73 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      {glycolysis, downstream Ehrlich, Adh1}, closed index 0.728 ± 0.011
      (Shapley sum 80 %). Stable in sample size: 0.715 / 0.733 / 0.728.

      +++
      Screening prior (SCR)

   .. grid-item-card:: 83 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Same triple, closed index 0.827 ± 0.010. **Not converged**: still
      rising, 0.776 / 0.804 / 0.827.

      +++
      Custom, all-linear (CUS)

   .. grid-item-card:: 87 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      A different triple, {Pdc, Adh1, glycolysis}, closed index
      0.873 ± 0.006. Stable: 0.874 / 0.871 / 0.873.

      +++
      Campaign box (CAM)

The convergence values are at 1,024 rows, 2,048 rows and the full sample.

.. list-table:: Shapley effects (%) on PI for all 12 levers, by measure (sum = 100; replicate s.d. ≤ 1.0 point)
   :header-rows: 1
   :stub-columns: 1
   :widths: 34 22 22 22

   * - Lever
     - CAM
     - SCR
     - CUS
   * - glycolysis
     - 8.3
     - **40.0**
     - **43.8**
   * - downstream Ehrlich (``k_14``–``k_16``)
     - 3.3
     - **29.9**
     - **31.1**
   * - ``k_6`` (Adh1)
     - **43.1**
     - **10.2**
     - **12.2**
   * - ``k_3`` (Pdc)
     - **38.9**
     - 3.3
     - 3.0
   * - ``k_13`` (ALS)
     - 2.6
     - 3.5
     - 1.3
   * - ``k_17`` (Adh6)
     - 3.4
     - 2.4
     - 0.4
   * - ``threshold_conc``
     - 0.2
     - 3.6
     - 2.8
   * - ``target_delta``
     - 0.2
     - 3.5
     - 3.3
   * - ``max_n_spikes``
     - 0.1
     - 1.5
     - 1.1
   * - ``inhib_isobutanol``
     - −0.1
     - 1.2
     - 0.2
   * - ``inhib_ethanol``
     - 0.1
     - 0.8
     - 0.8
   * - ``inhib_acetate``
     - −0.1
     - 0.1
     - −0.0

Small negative entries are estimation error.

.. list-table:: Closed-subset indices for PI, :math:`S_u` (± replicate s.d.)
   :header-rows: 1
   :stub-columns: 1
   :widths: 34 22 22 22

   * - Subset *u*
     - CAM
     - SCR
     - CUS
   * - {glycolysis, Ehrlich, ``k_6``}
     - 0.367 ± 0.031
     - **0.728 ± 0.011**
     - **0.827 ± 0.010**
   * - {``k_3``, ``k_6``, glycolysis}
     - **0.873 ± 0.006**
     - 0.388 ± 0.020
     - 0.422 ± 0.013
   * - {glycolysis, Ehrlich}
     - 0.043 ± 0.019
     - 0.575 ± 0.006
     - 0.698 ± 0.009
   * - {``k_3``, ``k_6``}
     - 0.750 ± 0.009
     - 0.051 ± 0.010
     - 0.105 ± 0.008
   * - ``k_13`` alone
     - 0.000 ± 0.017
     - 0.009 ± 0.018
     - 0.002 ± 0.013
   * - feeding {threshold, Δ, cap}
     - 0.006 ± 0.022
     - 0.073 ± 0.016
     - 0.060 ± 0.002
   * - inhibition {EtOH, IBO, acetate}
     - 0.005 ± 0.002
     - 0.030 ± 0.012
     - −0.000 ± 0.014

**Why the two prior measures put the proteome sinks first.** Glycolysis and
the downstream Ehrlich enzymes are the model's two largest proteome sinks.
Together they account for 95.7 % (SCR), 95.4 % (CUS) and 97.0 % (CAM) of the
variance of the sampled enzyme pools. Since designs worse than the baseline
carry ≥ 95 % of Var(PI), the variance mostly measures how badly a design
fails, and in the prior measures that is set largely by how much proteome
it spends:

- The enzyme load Φ_M on its own explains η² = **0.39** of Var(PI) under SCR
  and **0.47** under CUS (0.012 under CAM).
- PI is hump-shaped in glycolysis. Median PI by SCR glycolysis quintile is
  −2.40 (0.10–0.19×), −1.92 (0.19–0.34×), −1.67 (0.34–0.64×), −1.73
  (0.64–1.25×) and −3.01 (1.25–4.0×). CUS shows the same shape.
- PI falls at high downstream-Ehrlich expression (SCR quintiles −1.88, −1.75,
  −1.82, −2.15, −2.80).
- **Every** design with PI > 0 (19 under SCR plus 15 under CUS) keeps
  Φ_M ≤ 0.17 g/gDCW, which is under 70 % of the flexible proteome. Their
  median growth-derating factor is 1.0 (no derating).

**Why the campaign measure puts Pdc and Adh1 first.** CAM's variance is
dominated by "native ethanol pathway knocked out or not", because its bands
reach down to 10⁻³× native. Under CAM, ``k_6`` alone gives η² = 0.25 of PI.

.. important::

   A variance share in this section is a property of the sampling measure
   it was computed under, so the same lever can carry 3 % under one measure
   and 39 % under another; every share on this page names its measure.
   *Measure-dependent:* the size of ``k_3`` (39 → 3 %), glycolysis (8 → 40–44
   %), downstream Ehrlich (3 → 30–31 %), which triple is best, and the closed
   index of any triple. *Robust across all three:* ``k_6`` is in PI's top
   three; ``k_13`` and ``k_17`` each carry ≤ 3.5 % of Var(PI); the three
   inhibition groups carry ≤ 1.2 % each and ≤ 3 % jointly; and designs worse
   than baseline carry ≥ 95 % of Var(PI). This page therefore reports the
   proteome-sink triple under every measure (73 % SCR, 83 % CUS, 87 % CAM
   for a different triple) rather than a single headline number.

Profitable designs are rare, and all of them co-produce
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Outcome distribution of the sampled designs (completed rows; baseline PI_A = −0.1383)
   :header-rows: 1
   :stub-columns: 1
   :widths: 40 20 20 20

   * - Quantity
     - CAM
     - SCR
     - CUS
   * - PI median / max
     - −4.870 / −0.077
     - −2.101 / **0.411**
     - −2.081 / 0.308
   * - Rows beating PI_A
     - 2 (0.05 %)
     - 33 (0.81 %)
     - 46 (1.12 %)
   * - Rows with PI > 0 (IRR > 15 %)
     - 0
     - 19 (0.46 %)
     - 15 (0.37 %)
   * - Share of Σ(PI − mean)² from rows ≤ PI_A
     - 99.5 %
     - 95.7 %
     - 95.5 %
   * - Best row
     - PI −0.077, IRR 13.7 %, ethanol-only
     - PI 0.411, IRR 21.5 %, IBO 32.3 / EtOH 42.6 g/L
     - PI 0.308, IRR 19.9 %, IBO 26.5 / EtOH 45.8 g/L

All 34 sampled designs with PI > 0 co-produce, at isobutanol 10.2–34.1 g/L and
ethanol 16.1–125.1 g/L (IRR 15.1–21.5 %). The best ethanol-only rows reach IRR
13.7 % (CAM), 12.6 % (SCR) and 14.3 % (CUS). This holds for
**space-filling samples only**: the ethanol-only profitability campaign
reached IRR 15.8 % (:doc:`kinetic_optimization`).

What a variance-based screen cannot say
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Entry into the profitable zone is unlearnable from these samples.** The
  surrogate for the indicator "PI > PI_A" reaches a cross-validated Q² of
  only **0.024** (SCR) and 0.064 (CUS). Conditional indices on "how much
  better" could not be estimated. The zone surrogates have negative Q², and
  at the observed zone frequency about 25,000 SCR points would be needed to
  supply the 200 zone rows the estimator requires.
- **Feeding looks unimportant globally but gates the good region.** The
  three feeding variables explain only 0.6–7.3 % of Var(PI). Yet 62 % of the
  sampled designs sit at the 300 g/L target ceiling, against 1 of 33 (SCR)
  and 0 of 46 (CUS) of the designs that beat the baseline (binomial
  p ≤ 1e-12). See :doc:`feeding_strategies` for the feeding policy itself.
- **Search, not sampling, finds the good designs.** The best of the 4,091
  SCR designs reaches PI 0.41 (IRR 21.5 %). The Bayesian-optimization
  campaigns on the same model build reached PI ≈ 0.80–0.81
  (:doc:`kinetic_optimization`). A variance-based screen ranks what moves
  the typical, failing design, not what to pursue. Finding the good designs
  took search.
- **A small share does not make a lever dispensable.** Every design with
  PI ≥ 0.75 found by search sits at the ALS cap and at the most tolerant
  isobutanol setting (:doc:`kinetic_optimization`), and along the ALS slice
  through the flagship design #1912, PI is positive only above ``k_13`` ≈
  0.39 (:doc:`objective_landscape`). ALS and isobutanol tolerance pay only
  together with the downstream Ehrlich enzymes and Adh6, inside a
  profitable corner (0.46 % of the SCR designs clear the hurdle) too small
  to register in a variance that failing designs dominate (≥ 95 %).

Smaller results
~~~~~~~~~~~~~~~

- **Tolerance multipliers carry little variance within their 0.75–1.5×
  band** (≤ 1.2 % each, ≤ 3 % jointly), yet the best designs sit at its
  most tolerant end: all 48 with PI ≥ 0.75 at the isobutanol floor, and
  the ten with PI ≥ 0.79 at both the ethanol and isobutanol floors
  (:doc:`kinetic_optimization`). The band is only a
  2-fold range against the 40-fold expression bands. :doc:`kinetic_sweeps`
  explores a much wider inhibition range.
- **Adh6 saturates at about 2× native.** In the SCR subset with
  ``k_13`` ≥ 1, 0.4 ≤ Ehrlich ≤ 2 and derating factor d > 0.9 (n = 879),
  median isobutanol titer rises from 3.5 g/L (``k_17`` ≤ 0.05) to 21.9 g/L
  (0.1–0.2). It then plateaus at 24.6–26.2 g/L for ``k_17`` > 0.2. This is
  a property of the model, not a measured saturation: Adh6's native
  capacity (0.1077 g/L/h) is not fitted to data but derived from the fitted
  Adh1 capacity through the two enzymes' measured abundance ratio, so
  "2× native" is relative to that estimate.
- **The response is interaction-dominated.** Under SCR, S1 / ST are 23.5 /
  57.1 % for glycolysis, 15.6 / 44.4 % for downstream Ehrlich and 4.5 / 17.4 %
  for ``k_6``.

.. _ibo-gsa-caveats:

Caveats and limitations
-----------------------

.. warning::

   **Model build.** The three Sobol' runs (2026-09-20 to 09-23) used the
   2b5b27d-era heat-exchanger-network build, the same build as every
   optimization campaign. Under that build the scenario-A baseline is PI
   −0.1383 (MESP 3.844 $/GGE, IRR 12.60 %). The :doc:`uncertainty` analysis
   used the later 6d4776f build, so a number on this page and an
   uncertainty-page MESP differ by the build shift as well as by the
   analysis. The GSA has not been re-run
   under the later build. Its expected effect on the indices is small but
   unverified.

- **Measure dependence.** Variance shares are properties of the chosen
  measure and bands. Only the statements listed as robust above survive all
  three.
- **Dependent inputs.** First-order and total indices do not bracket the
  Shapley effect, and S1 > ST occurs for the feeding variables. Quote Shapley
  effects or closed indices.
- **Surrogate-based.** Every index is an index of the fitted surrogate. The
  PI surrogate's Q² is 0.827 (SCR), 0.857 (CUS) and 0.926 (CAM). The
  log-tail PI surrogate under SCR (Q² 0.797, just below the 0.8 flag) gives
  the same ranking.
- **Burden-model contingency.** The proteome-sink ranking depends on the
  burden constants (σ_eff, F_flex, kcat/MW). No sensitivity of the indices to
  those constants has been run; it is an open follow-up. Treat the
  proteome story as a model-contingent prediction.
- **No MESP.** Ethanol MESP is not tracked, so the GSA says nothing about
  MESP ordering. PI and MESP are not rank-equivalent once isobutanol is
  co-produced.
- **Figures.** The GSA bar charts and Shapley heatmaps are not reproduced
  here; this page tabulates the numbers instead. Stage 2 regenerates the
  figures from the stage-1 outputs.

Reproduce it
------------

.. tab-set::

   .. tab-item:: Stage 1: simulate (opt-in, long)

      ``analyses/sensitivity/evaluate_sobol_split12d.py`` simulates 4,096
      feasible designs by default, each through the campaigns' own
      evaluation path. It is a **long, opt-in run** of about 20 hours or
      more. Run one simulation process at a time. From the ``analyses/``
      folder:

      .. code-block:: bash

         # campaign measure (default)
         python sensitivity/evaluate_sobol_split12d.py
         # screening prior
         IBO_SOBOL_MEASURE=screening python sensitivity/evaluate_sobol_split12d.py
         # custom all-linear bands (the band file ships in data/sobol_bands)
         IBO_SOBOL_MEASURE=custom IBO_SOBOL_BANDS=linear_0-4x_2026-09-22.json \
             IBO_SOBOL_SEED=20260922 python sensitivity/evaluate_sobol_split12d.py

      .. rst-class:: ibo-long-literals

      - Environment overrides: ``IBO_SOBOL_N_POINTS`` (default 4096) and
        ``IBO_SOBOL_SEED`` (default 20260920).
      - The run is **resumable**: a logged row is never re-simulated. It is
        also **nested**: raising ``IBO_SOBOL_N_POINTS`` extends the same
        design. It prints ``STUDY_NAME=`` first, so a long run can be
        supervised: ``python supervise_sweep.py
        sensitivity/evaluate_sobol_split12d.py --stem <STUDY_NAME>
        --checkpoint-suffix _trajectory.csv``.
      - It writes ``<study>_trajectory.csv`` (the campaign column set),
        ``<study>_design.json`` (the single source of truth for the space,
        bands, burden and ``meta.baseline_PI``) and
        ``<study>_candidates.csv``.
      - Study stems of the three runs:
        ``kin_sobol_ethanol_isobutanol_metabolic_split_12d_rb0.001-4_ib0.75-1.5_aA_burden_seed20260920``
        (CAM),
        ``..._split_12d_screening_rb0.1-4_lin0_ib0.75-1.5_aA_burden_seed20260920``
        (SCR) and ``..._split_12d_custom9f1df84b_aA_burden_seed20260922``
        (CUS).

   .. tab-item:: Stage 2: analyse (sim-safe)

      ``analyses/sensitivity/analyze_sobol_split12d.py`` reads the stage-1
      CSV and design JSON only. It loads its helpers by file path and never
      imports the biorefinery, so it is safe to run next to a simulation.

      .. code-block:: bash

         python sensitivity/analyze_sobol_split12d.py --self-test          # synthetic, no data
         python sensitivity/analyze_sobol_split12d.py --study-name <STUDY_NAME>
         python sensitivity/analyze_sobol_split12d.py --study-name <STUDY_NAME> --replot

      - It writes the index, subset, target, conditional and convergence
        tables (``*_sobol_indices.csv``, ``*_sobol_PI_subsets.csv``, ...),
        ``*_sobol_summary.txt``, the PI bars, the Shapley heatmap and
        ``*_sobol_bars_stack.{png,pdf}``.
      - Useful flags: ``--n-base`` (4096), ``--n-replicates`` (3),
        ``--min-rows`` (200), ``--gp-metrics`` (default: the two PI metrics),
        ``--pi-threshold`` (default: the record's ``meta.baseline_PI``) and
        ``--min-zone-rows`` (200).
      - Cost: one GP prediction of the 4,096-point design takes ~0.6 s, and
        each GP metric needs 4,094 subsets × the replicates, about 43 min per
        GP metric per replicate. A GP-backed stage 2 therefore runs for hours.
        Gradient-boosting surrogates predict in milliseconds, and
        ``--replot`` redraws the figures from the saved tables in seconds.

.. note::

   **Where results go.** No analysis results ship with NSKinetics. Both
   stages read and write under the results root, ``analyses/results``
   inside the package by default. To reuse an earlier run's trajectories,
   point ``NSKINETICS_IBO_RESULTS_DIR`` at a folder that holds them.
   :doc:`running_analyses` covers the installation, the results root and
   the one-simulation-at-a-time rule.

See also
--------

- :doc:`kinetic_optimization`: the Bayesian-optimization campaigns over the same 12
  decisions, which find the profitable designs this screen cannot locate.
- :doc:`enzyme_burden`: the proteome-allocation model behind the feasibility
  constraint and the growth derating.
- :doc:`uncertainty`: Monte Carlo uncertainty of the *baseline* ethanol plant,
  where contextual inputs, not kinetics, dominate cost risk.

Next: :doc:`running_analyses` shows how to install the case study and
reproduce every result on these pages.
