Uncertainty analysis
====================

How risky is the baseline corn-ethanol plant (scenario A, the wild-type
strain) once prices, plant context and the fitted kinetics are all uncertain?
And how much of that cost risk comes from the strain? This page answers by
propagating 57 uncertain inputs (40 of them kinetic parameters of the
NSKinetics fermentation model) through 6,000 Monte Carlo simulations of the
full flowsheet and TEA. The answer: **context, not kinetics, sets cost risk.**
Corn and DDGS prices, plant uptime and scale dominate the spread in the
minimum ethanol selling price (MESP). The kinetics act only through titer,
rate and yield (TRY), and among those only yield has real per-unit leverage.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 3.83 $/GGE
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      MESP at baseline. The 5th–95th percentile range is 3.64–4.07 $/GGE
      and the median is 3.85. Capital is 139 MM$ [125–151].

      +++
      Scenario A and the 6,000-sample run, build 6d4776f

   .. grid-item-card:: ~92 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Share of the variance in MESP from contextual inputs: corn price 41 %,
      DDGS price 24 %, operating days 10.5 %, capacity 7.9 %.

      +++
      6,000-sample uncertainty run, build 6d4776f

   .. grid-item-card:: ~5 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Share from all 40 kinetic parameters together. Their effect is fully
      mediated by TRY: context + TRY explains MESP with linear R² 0.993.

      +++
      6,000-sample uncertainty run, build 6d4776f

   .. grid-item-card:: −0.54
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      MESP elasticity to ethanol yield, close to corn price's +0.65. For
      titer it is −0.08, and for productivity only −0.014.

      +++
      6,000-sample uncertainty run, build 6d4776f

.. admonition:: Key insights
   :class: tip

   - **A small variance share is not small leverage.** Per percent, ethanol
     yield is worth almost as much as corn price. Titer and productivity
     barely register: fermentation is ~10 % of installed equipment cost,
     while corn is ~86 % of material cost.
   - **A strain program cannot outrun the commodity spread.** Across their
     sampled ranges, ethanol production and DDGS revenue move the
     binned-median MESP by 0.21 and 0.23 $/GGE. The three TRY metrics move it
     by only 0.07–0.09 $/GGE.
   - **One kinetic lever shows through.** Once contextual variance is
     regressed out, ethanol inhibition of glycolysis (``k_1ie``) correlates
     ρ = 0.70 with the residual MESP. :doc:`kinetic_sweeps` maps that lever.

At a glance
-----------

.. list-table::
   :widths: 28 72
   :stub-columns: 1

   * - Scenario
     - A: wild-type ethanol strain, isobutanol pathway off. Baseline
       fed-batch policy: threshold 217.125 g/L, target 221.25 g/L, 16-spike
       cap.
   * - Run
     - ``_IBO_2026.9.29-23.45``: 6,000 samples, **0 failed**
   * - Inputs
     - 57 uncertain inputs (40 kinetic + 17 contextual), plus one dummy
       "Blank parameter" as a null reference
   * - Outputs
     - MESP (purity-adjusted ethanol MPSP at a 15 % IRR, in $/GGE), TCI,
       ethanol production, titer, yield, productivity, IRR at default prices
   * - Model build
     - 6d4776f (scenario-A baseline 0.86021 $/kg = 3.835 $/GGE, IRR
       12.79 %). See :ref:`ibo-unc-caveats`.
   * - Scripts
     - ``analyses/uncertainty/uncertainties_IBO_EtOH.py`` (simulation) →
       ``plots/plot_uncertainty_MPSP_vs_TCI.py`` and
       ``plots/plot_uncertainty_spearman_tornado.py`` (sim-safe figures)

Method
------

Sampling
~~~~~~~~

``uncertainties_IBO_EtOH.py`` loads the biorefinery and fixes scenario A's
feeding policy. It then loads the scenario's parameter-distribution workbook,
``parameter-distributions_corn_IBO_EtOH_A.xlsx``, and draws 6,000
Latin-hypercube samples (``model.sample(N=6000, rule='L')``, NumPy seed 3221).
The seed makes the samples identical across launches. Each sample is a full
flowsheet simulation followed by a TEA solve. The script sets scenario A's
workbook and feeding policy directly rather than through ``load_scenario``,
so it runs **without** the enzyme-burden constraint (:doc:`enzyme_burden`).
At scenario A ± 20 % the constraint would be inert anyway: the modelled
enzyme load stays at or below about 0.084 g/gDCW, well under the point at
which growth derating would start (0.1347 g/gDCW at wild-type growth
capacity).

**Kinetic inputs (40).** Each fitted kinetic parameter of the NSKinetics
*S. cerevisiae* model gets an independent **triangular ±20 %** distribution
(0.8 / 1.0 / 1.2 × the fitted value). Physically constrained parameters are
excluded by *type*, not by value: five in scenario A, which takes it from 45
to 40 kinetic inputs, plus ``k_16r`` and ``P_10i`` in B:

.. list-table::
   :header-rows: 1
   :widths: 26 30 44

   * - Excluded
     - Category
     - Why a flat ±20 % is wrong for it
   * - ``k_6r`` (and ``k_16r`` in B)
     - Thermodynamic reverse ratio
     - It encodes the apparent equilibrium constant. Perturbing it alone
       violates the Haldane relationship and can drive Adh backwards.
   * - ``K_2``, ``K_9``
     - Structural saturation constants
     - These are zero-order devices, not measured Michaelis constants.
   * - ``P_10e``, ``P_10a`` (and ``P_10i`` in B)
     - Piecewise death thresholds
     - These are switches, not smooth gains.

**Contextual inputs (17).** Twelve TEA prices and finance terms, corn starch
content, plant capacity, the aeration-rate safety factor and two storage
times. Examples:

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - Input
     - Distribution
     - Values
   * - Corn price
     - triangular
     - 0.1591 / 0.1741 / 0.1879 $/wet-kg (−8.6 / +7.9 %)
   * - DDGS price
     - uniform
     - 0.1144–0.1717 $/kg (±20 %)
   * - Annual operating days
     - triangular
     - 297 / 330 / 363 d
   * - Plant capacity (corn feed)
     - triangular
     - 45,578 / 56,972 / 68,366 kg/h (±20 %)
   * - Corn starch content
     - uniform
     - 0.684–0.743
   * - Federal corporate tax rate
     - uniform
     - 0.15–0.28

The ethanol price is **not** sampled. MESP is the quantity solved for, at the
fixed 15 % IRR. IRR is removed from the sampled set and reported as a metric,
solved at default product prices (ethanol 0.835 $/kg, 2023$). The
"Blank parameter" is a dummy uniform 0.9–1.1 input that changes nothing. Its
ρ(MESP) = −0.0004 serves as a built-in null.

Metrics and statistics
~~~~~~~~~~~~~~~~~~~~~~

**MESP** is the purity-adjusted ethanol MPSP at a 15 % IRR (product purity
0.95046), converted from $/kg to $ per gallon of gasoline equivalent:

.. math::

   \text{MESP}\;[\$/\text{GGE}] &= \text{MPSP}\;[\$/\text{kg}]
   \times \frac{0.789 \times 3.785411784}{0.67} \\
   &= \text{MPSP}\;[\$/\text{kg}] \times 4.45775

Titer and productivity are per litre of **water**, not broth. Yield is g per g
of sugar **added**, and productivity = titer/τ. The analyses on this page are:

- **Spearman rank correlations** (with p-values) of every output against every
  input. These are the workbook's own ``Spearman`` and ``Spearman p-values``
  sheets. With 6,000 samples, the null 95 % band is :math:`|\rho| < 0.0253`.
- **Variance shares.** MESP is regressed on the standardized inputs, and each
  input's share is β²·Var(z)/Var(MESP).
- **Elasticities.** log MESP is regressed on log of the 17 contextual inputs
  and log yield, titer and productivity (R² 0.993).
- **Context residual.** MESP is regressed on the 17 contextual inputs
  (R² 0.934, residual s.d. 0.034 $/GGE), and the residual is correlated with
  each kinetic input.

Results
-------

Outcome distributions
~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Baseline and Monte Carlo percentiles (6,000 samples)
   :header-rows: 1
   :stub-columns: 1
   :widths: 34 16 16 16 18

   * - Metric
     - Baseline
     - 5th
     - Median
     - 95th
   * - MESP ($/GGE)
     - **3.835**
     - 3.637
     - 3.852
     - 4.075
   * - TCI (MM$)
     - 139.1
     - 125.1
     - 137.8
     - 151.0
   * - Ethanol production (MM gal/y, contained)
     - 47.02
     - 39.33
     - 46.38
     - 54.08
   * - Titer (g/L water)
     - 114.5
     - 104.7
     - 114.0
     - 122.6
   * - Yield (g/g sugar added)
     - 0.4602
     - 0.4510
     - 0.4596
     - 0.4659
   * - Productivity (g/L/h)
     - 1.859
     - 1.503
     - 2.057
     - 2.674
   * - Fermentation time τ (h)
     - 61.6
     - 44.9
     - 55.3
     - 72.2
   * - Glucose spikes
     - 7
     - 2
     - 7
     - 12
   * - IRR at default prices
     - 12.79 %
     - 8.00 %
     - 12.48 %
     - 16.66 %

The spread is narrow. The coefficients of variation are 3.4 % for MESP,
5.7 % for TCI, 1.0 % for yield, 4.6 % for titer and 17.7 % for productivity.
Every sample lies inside the 2021–2025 ethanol market price range
(2.31–5.15 $/GGE) and the reported capital range (110.3–195.6 MM$, 2023$).
92.5 % of samples lie inside the 40–60 MM gal/y production range. These
range checks are sanity checks only: the market range is far wider than the
MESP spread. At the model's ethanol price, the midpoint of the 2021–2025
range (3.722 $/GGE = 0.835 $/kg), **17.0 %** of the samples clear the 15 %
hurdle. The baseline MESP sits 0.112 $/GGE (2.9 %) above that price.

.. figure:: /_static/images/isobutanol/uncertainty/uncertainty_bivariate_mesp_scenario_a.png
   :target: ../../_static/images/isobutanol/uncertainty/uncertainty_bivariate_mesp_scenario_a.png
   :alt: Six-panel figure of kernel-density clouds from 6,000 Monte Carlo samples. MESP sits on the shared y-axis, against total capital investment, ethanol production, DDGS revenue, ethanol yield, ethanol titer and ethanol productivity. Each cloud is compact around a white baseline diamond near 3.8 dollars per GGE. Panel A carries a red Pareto staircase. Panels B to F carry black binned-median lines that slope only gently downward.
   :width: 100%
   :class: pub-figure

   **MESP against capital, scale, coproduct revenue and TRY** (run
   ``_IBO_2026.9.29-23.45``, build 6d4776f). Follow the black binned-median
   lines: across their sampled ranges, ethanol production (B) and DDGS
   revenue (C) move MESP by 0.21 and 0.23 $/GGE, while yield, titer and
   productivity (D–F) move it by only 0.07–0.09 $/GGE. The red staircase in
   A is **not** a design trade-off; see below. Click the figure to open it
   at full resolution.

.. dropdown:: How to read this figure

   - **Clouds and markers.** Filled contours are kernel-density estimates of
     the 6,000 samples. The white diamond is the baseline. The black lines
     are binned medians over ten equal-count bins.
   - **Margins.** Marginal box plots show the 25th–75th percentiles (box),
     the 5th–95th (whiskers) and the minimum and maximum samples (dots).
     Boxed values are Spearman ρ with MESP.
   - **Grey reference boxes.** Light-grey boxes span typical corn-ethanol
     ranges of each abscissa (TCI, production, DDGS revenue and high-gravity
     TRY ranges) × the 2021–2025 ethanol market price range. The darker grey
     band is the 2021–2025 retail gasoline range, 2.16–4.84 $/GGE.
   - **Panel C.** The model's DDGS price (~130 $/short ton) lies below its
     2021–2025 low, so the samples hug the low edge of the DDGS box.

Where the variance comes from
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Shares of Var(MESP) (standardized linear regression)
   :header-rows: 1
   :widths: 50 25 25

   * - Group / input
     - Share of Var(MESP)
     - Type
   * - **All 17 contextual inputs**
     - **91.8 %**
     - context
   * - corn price
     - 41.1 %
     - context
   * - DDGS price
     - 23.8 %
     - context
   * - annual operating days
     - 10.5 %
     - context
   * - plant capacity
     - 7.9 %
     - context
   * - corn starch content
     - 4.8 %
     - context
   * - corporate tax rate
     - 2.2 %
     - context
   * - **All 40 kinetic parameters**
     - **5.4 %**
     - kinetics
   * - ``k_1ie`` (ethanol inhibition of glycolysis)
     - 3.0 %
     - kinetics
   * - ``k_10`` (active-biomass decay)
     - 0.7 %
     - kinetics
   * - nonlinearity / interaction residual
     - 1.1 %
     - —

The listed group shares sum to 98.3 %. The linear fit on all 57 inputs
explains 98.9 % of Var(MESP) (the remaining 1.1 % is the residual row),
while the individual β²·Var(z) shares add up to only 97.2 %; the
difference, about 1.7 %, is carried by the covariance (cross) terms of the
regression, from the weak chance correlations between independently sampled
inputs in a finite sample.

The kinetic share acts **entirely through TRY**. TRY alone explains little
of MESP (linear R² 0.068), but adding TRY to the contextual inputs raises the
linear R² from 0.934 to **0.993**. In the other direction, TRY are almost
purely kinetic outcomes: the 40 kinetic inputs give linear R² of 0.854
(yield), 0.984 (titer) and 0.814 (productivity), against ≤ 0.006 for the
contextual inputs.

**Capital and MESP correlate only through plant size.** ρ(TCI, MESP) = −0.26,
but the partial rank correlation given capacity is −0.01. The red Pareto
staircase in panel A (15 points) is a "luck" frontier traced by favourable
price draws. Its points have a median corn-price percentile of 7, and 14 of
15 sit above the median DDGS price. It shows a best-case envelope, not a
TCI-versus-MESP design trade-off.

Small variance share ≠ small leverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variance shares depend on the assumed input ranges. Elasticities do not:
they measure leverage per percent of change in each input.

.. list-table:: MESP elasticities, d ln MESP / d ln x (log-log regression, R² 0.993)
   :header-rows: 1
   :widths: 40 20 40

   * - Driver
     - Elasticity
     - Reading
   * - Corn price
     - +0.65
     - +1 ¢/kg corn → +0.145 $/GGE
   * - **Ethanol yield**
     - **−0.54**
     - +0.01 g/g → −0.047 $/GGE
   * - Corn starch content
     - −0.32
     - more starch also means less DDGS to sell
   * - Annual operating days
     - −0.27
     - output without added capital
   * - DDGS price
     - −0.14
     - coproduct credit
   * - Plant capacity
     - −0.12
     - economy of scale
   * - **Ethanol titer**
     - **−0.08**
     - +10 g/L → −0.027 $/GGE
   * - **Ethanol productivity**
     - **−0.014**
     - +1 g/L/h → −0.026 $/GGE

Yield already sits at **90 %** of the 0.511 g/g stoichiometric maximum. A
linear extrapolation beyond the sampled range says that pushing yield, titer
and productivity all to their ceilings would lower MESP by roughly
**0.4 $/GGE (~10 %)**. That is about the size of the contextual
5th–95th-percentile span. Improving ethanol TRY within the incumbent process
cannot deliver a step change in economics. :doc:`kinetic_sweeps`, next,
checks this far beyond the ±20 % neighbourhood, and
:doc:`kinetic_optimization` then turns to a higher-value co-product.

.. dropdown:: TRY ceilings, item by item (indicative linear extrapolation)

   .. list-table::
      :header-rows: 1
      :widths: 60 40
      :class: ibo-num ibo-compact

      * - Change from baseline
        - ΔMESP ($/GGE)
      * - Yield 0.460 → 0.511 g/g (stoichiometric maximum)
        - −0.238
      * - Titer 114.5 → 142 g/L
        - −0.073
      * - Productivity 1.86 → 4.4 g/L/h
        - −0.066
      * - All three combined
        - ≈ −0.38 (−9.8 %)

   The three TRY metrics are collinear outcomes of the same kinetics, so the
   partial effects are model-dependent. Treat the ceiling as indicative.
   Raising the baseline productivity from 1.86 to the 2.0 g/L/h bottom of the
   industrial range would change MESP by only −0.004 $/GGE.

Spearman drivers of cost and of TRY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: /_static/images/isobutanol/uncertainty/spearman_tornado_scenario_a.png
   :target: ../../_static/images/isobutanol/uncertainty/spearman_tornado_scenario_a.png
   :alt: Five tornado charts of Spearman rank correlations. MESP is driven by corn price (+0.64), DDGS price (−0.48), operating days, capacity, starch, the kinetic parameter k_1ie (+0.18) and tax. TCI is driven mainly by plant capacity (+0.93). Ethanol titer, productivity and yield are driven entirely by kinetic parameters such as k_7ie, k_1ie, k_7 and k_10.
   :width: 100%
   :class: pub-figure

   **Spearman tornados for MESP, TCI and the three ethanol TRY metrics**
   (run ``_IBO_2026.9.29-23.45``). Each panel shows only the inputs with
   :math:`|\rho| \geq 0.10` and p < 0.05. Bars are sorted by :math:`|\rho|`
   and coloured by category: prices and tax, plant scale and uptime, corn
   composition, fermentation kinetics. **What to look at:** the colour
   change between columns. The economic outcomes on the left (A, B) are
   driven mostly by context. Kinetics enter MESP through a single bar
   (``k_1ie``) and TCI through two weak ones (``k_7ie`` and ``k_7``, both at
   :math:`|\rho|` = 0.15). The fermentation outcomes on the right (C–E) are
   driven entirely by kinetics. The kinetic parameters that dominate TRY
   (``k_7ie``, ``k_7``) barely register in MESP: ρ(MESP) is +0.04 for
   ``k_7ie`` and −0.07 for ``k_7``.

The kinetic drivers of TRY differ from the kinetic driver of cost:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Outcome
     - Strongest inputs (Spearman ρ)
   * - MESP
     - corn price +0.64, DDGS price −0.48, operating days −0.31, capacity
       −0.28, starch −0.22, ``k_1ie`` +0.18, tax +0.16
   * - TCI
     - capacity +0.93, starch +0.23, ``k_7ie`` +0.15, ``k_7`` −0.15
   * - Titer
     - ``k_7ie`` −0.57, ``k_1ie`` −0.50, ``k_7`` +0.42, ``k_10`` −0.41
   * - Productivity
     - ``k_7`` +0.55, ``k_1ie`` −0.48, ``k_7ie`` −0.46
   * - Yield
     - ``k_1ie`` −0.52, ``k_7ie`` +0.45, ``k_1e`` +0.35, ``K_1i`` −0.32,
       ``k_7`` −0.30

Ten of the 40 kinetic parameters lie outside the null band for MESP, and 8
pass a Bonferroni correction (p < 0.05/57). Among the contextual inputs, 9 of
17 pass. Isobutanol storage time is a structural null in scenario A, which
makes no isobutanol.

The hidden kinetic signal: ``k_1ie``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the raw ranking, ethanol inhibition of glycolysis (``k_1ie``) is only 6th
of 57 (ρ +0.18). Once the contextual variance is regressed out, it
correlates **ρ = 0.70** with the residual MESP. No other kinetic input comes
close:

.. list-table:: Correlation with the context-residual MESP vs with raw MESP
   :header-rows: 1
   :stub-columns: 1
   :widths: 40 30 30
   :class: ibo-num ibo-compact

   * - Input
     - ρ (residual MESP)
     - ρ (raw MESP)
   * - ``k_1ie``
     - **+0.704**
     - +0.182
   * - ``k_10``
     - +0.332
     - +0.078
   * - ``k_1e``
     - −0.246
     - −0.053
   * - ``K_1i``
     - +0.231
     - +0.062
   * - ``k_4ie``
     - −0.213
     - −0.061
   * - ``k_7``
     - −0.166
     - −0.065
   * - Blank parameter (null)
     - −0.012
     - −0.0004

``k_1ie`` lowers titer, yield and productivity together, with no trade-off
among them. Across its deciles, the MESP median rises 3.821 → 3.893 $/GGE
while productivity falls 2.39 → 1.81 g/L/h. :doc:`kinetic_sweeps` maps this
lever on a much wider range, jointly with glycolytic capacity. Ethanol
inhibition of **growth** (``k_7ie``) is MESP-neutral near scenario A
(elasticity +0.008): its yield gain almost exactly cancels its titer and
productivity losses.

.. _ibo-unc-caveats:

Caveats and limitations
-----------------------

.. warning::

   **Model build.** This run used the 6d4776f heat-exchanger-network build
   (baseline 3.835 $/GGE, IRR 12.79 %). The Sobol' GSA, the optimization
   campaigns and the kinetic sweeps used the earlier 2b5b27d build (baseline
   3.844 $/GGE, IRR 12.60 %), a 0.25 % offset in the baseline. A campaign
   IRR and an MESP from this page therefore come from different builds, and
   the offset grows to a few percent at co-production designs, so the two
   are not directly comparable.

.. caution::

   **Earlier runs are superseded.** Two earlier Monte Carlo runs of this
   analysis (``_IBO_2026.9.20-22.24`` and ``_IBO_2026.9.29-18.23``) are
   superseded by the run on this page: the first used an older, wider
   starch band and a starch setter that was later fixed, and the second ran
   before that fix. Cite only run ``_IBO_2026.9.29-23.45``.

- **Shares are conditional on the assumed distributions.** Kinetic
  parameters span ±20 %, while corn price spans only −8.6 / +7.9 %. Prices are
  sampled independently, with no corn–DDGS co-movement. The elasticities are
  the range-independent robustness check.
- **Parametric, not strain-to-strain, variation.** The TRY spread reflects
  uncertainty in one fitted strain model. It covers only part of the
  industrial ranges (yield 25 %, titer 42 %, productivity 28 % of their
  0.42–0.48 g/g, 100–142 g/L and 2.0–4.4 g/L/h ranges).
- **Regression-based leverage.** The elasticities come from regressions on
  Monte Carlo outcomes. TRY are collinear, so the partial effects are
  model-dependent, and the TRY ceilings are linear extrapolations.
- **Productivity sits low.** The baseline productivity (1.86 g/L/h) and
  43.8 % of the samples lie below the 2.0–4.4 g/L/h high-gravity range. Given
  the −0.014 elasticity, this conservatism does not bias MESP.
- **Local results.** The ``k_7ie`` neutrality and the ``k_1ie`` signal hold in
  the ±20 % neighbourhood of the fitted model. :doc:`kinetic_sweeps` and the
  :doc:`sensitivity` analysis explore far wider ranges.
- **Significance workbook.** The companion
  ``significant_parameters_from_sensitivity_analysis_*.xlsx`` file the run
  writes came out **empty** for this run; use the ``Spearman`` and
  ``Spearman p-values`` sheets of the results workbook instead.

Reproduce it
------------

.. tab-set::

   .. tab-item:: Monte Carlo run (opt-in, long)

      ``analyses/uncertainty/uncertainties_IBO_EtOH.py`` runs the full
      analysis for scenario A: it loads the biorefinery, draws 6,000 samples
      and simulates each one. It is a **long, opt-in run** of about 30
      hours. Run it alone: never alongside another simulation process.

      - **Inputs:** the parameter-distribution workbook
        ``parameter-distributions_corn_IBO_EtOH_A.xlsx``, shipped in the
        package's ``data/parameter_distributions`` folder.
      - **Outputs:** written to the results root and stamped with the launch
        time, ``_IBO_<stamp>_['A']_6000sims_A_0_baseline.xlsx`` and
        ``..._A_1_full_evaluation.xlsx``. The second has the sheets
        ``Parameters``, ``TEA results``, ``TEA percentiles``, ``Spearman``,
        ``Spearman p-values`` and ``Raw data``.
      - **Crash recovery:** the fermentation-coupled flowsheet can
        occasionally die inside native integrator code. Progress is
        therefore checkpointed every 10 samples to
        ``uncertainties_A_autosave.pickle``. A relaunch resumes from it,
        since the fixed seed regenerates the same samples. The checkpoint is
        deleted once the workbooks are written.

   .. tab-item:: Figures (sim-safe)

      Both plotters read workbooks only and never import the biorefinery, so
      they are safe to run at any time. With no argument, each picks the
      newest scenario-A ``*_1_full_evaluation.xlsx`` in the results root.

      .. code-block:: bash

         # from the package folder
         python plots/plot_uncertainty_MPSP_vs_TCI.py [<results.xlsx>] [--baseline-file F] [--out-dir DIR] [--stem STEM] [--dpi N]
         python plots/plot_uncertainty_spearman_tornado.py [<results.xlsx>] [--out-dir DIR] [--stem STEM] [--dpi N]

      They write ``uncertainty_bivariate_A_<tag>.{png,pdf}`` and
      ``spearman_tornado_A_<tag>.{png,pdf,csv}`` to
      ``publication/Uncertainty`` under the results root. The tornado CSV
      lists every bar shown.

.. note::

   **What must be regenerated.** No analysis results ship with
   NSKinetics, so the 6,000-sample workbooks behind this page must be
   regenerated, or ``NSKINETICS_IBO_RESULTS_DIR`` pointed at a folder that
   holds an earlier run's workbooks, before re-plotting. The tornado
   plotter also writes a CSV of every bar it draws (26 rows for this run),
   a ready-made table of the significant correlations shown above.
   :doc:`running_analyses` covers the installation, the results root and
   the one-simulation-at-a-time rule.

See also
--------

- :doc:`sensitivity`: a variance-based GSA over the *strain-design*
  decisions, where the question is profitability rather than cost risk.
- :doc:`kinetic_sweeps`: the glycolytic capacity × ethanol inhibition plane
  that follows up on the ``k_1ie`` signal.
- :doc:`feeding_strategies`: the fed-batch policy held fixed in this run.

Next: :doc:`kinetic_sweeps` follows the ``k_1ie`` signal far outside the
± 20 % band and asks how much economic headroom the classic ethanol-strain
levers have.
