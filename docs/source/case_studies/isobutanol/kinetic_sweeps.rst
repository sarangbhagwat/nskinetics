Kinetic sweeps
==============

.. currentmodule:: nskinetics

Strain engineers usually reach for two classic ethanol levers: more
glycolytic capacity and more ethanol tolerance. NSKinetics runs the kinetic
model inside the flowsheet, so any combination of kinetic parameters can be
priced directly. What does the plant earn if the strain's glycolysis is twice
as fast, or if ethanol stops inhibiting it? This page asks how much economic
headroom these levers have at the scenario-A process (wild-type ethanol
strain, fixed fed-batch policy), together with their siblings: ``k_1e`` alone,
pyruvate decarboxylase (Pdc) and alcohol dehydrogenase 1 (Adh1). It also asks
whether more tolerance is always better. To answer, it sweeps 40 × 40 grids of
kinetic-parameter multipliers. Each grid point is a full fed-batch
fermentation, flowsheet simulation and techno-economic analysis.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 1.46 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      How much the best design in the glycolysis × ethanol-inhibition plane
      improves on the baseline strain (3.788 vs 3.844 $/GGE); only 8.4 % of
      the plane beats the baseline at all

      +++
      Glycolysis × ethanol-inhibition plane, build 2b5b27d

   .. grid-item-card:: 13.7 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Highest IRR anywhere in the plane: no design clears the 15 % hurdle

      +++
      Glycolysis × ethanol-inhibition plane, build 2b5b27d

   .. grid-item-card:: +21 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      MESP of a strain with ethanol inhibition fully relieved, despite 3.7×
      the productivity; IRR falls from 12.6 % to −9.4 %

      +++
      Glycolysis × ethanol-inhibition plane, build 2b5b27d

   .. grid-item-card:: 19 of 22
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      glycolysis columns without growth derating whose cost-optimal
      inhibition is exactly where the 16-spike cap starts to bind

      +++
      Glycolysis × ethanol-inhibition plane, build 2b5b27d

.. admonition:: Key insights
   :class: tip

   - **Ethanol-only strain levers have little headroom.** At the baseline
     process, the best design in the plane is only 1.46 % cheaper than the
     baseline strain. Four distinct planes (glycolysis, ``k_1e``, Pdc and
     Adh1), plus two earlier runs of the glycolysis plane, six runs in all,
     agree: best improvements of 0.80–1.53 % and highest IRRs of
     13.2–13.7 %. A step change requires changing what the plant sells
     (:doc:`kinetic_optimization`).
   - **More ethanol tolerance can cost more.** The MESP optimum sits at
     about 0.77× the fitted inhibition coefficients (0.77–0.79× in all six
     runs). In every run whose inhibition axis reaches zero, relieving
     inhibition entirely raises MESP 18–23 % above the baseline at the
     baseline enzyme levels, because the extra tolerance feeds growth, not
     ethanol (raising Adh1 4× shrinks the penalty to 5.6 %; see
     `Sibling planes: k_1e, Pdc and Adh1`_).
   - **The process, not the cell, sets the tolerance optimum.** Lower
     inhibition pays only until the fed-batch spike cap binds. A tolerance
     target quoted without a feeding policy is meaningless. The *location*
     of the optimum depends on the cap, which was not varied.
   - **Single-metric optima sit on cliff edges or outside the valley.** The
     titer optimum matches the baseline MESP but needs 4.1× the cost
     optimum's aeration, and it sits one grid step from a +19 % cliff. The
     yield and productivity optima cost 5.7 % and 13.1 % more. The
     cell-density and TCI "optima" make no ethanol.
   - **The landscape is a flat, curved valley bounded by three steep
     cliffs**, each with a different mechanism.

Method
------

Why this plane
~~~~~~~~~~~~~~

In the uncertainty analysis, the one kinetic lever that showed through the
commodity-price spread was ethanol inhibition of glycolysis, with
ρ = 0.70 against context-adjusted MESP (:doc:`uncertainty`; build
``6d4776f``). The sweeps follow that lever far outside the ±20 % band, but
bundled: the y-axis scales ``k_1ie`` together with the other three
ethanol-inhibition coefficients (AcDH, growth and death; see **What the axes
mean** below), against glycolytic capacity on the x-axis. ``k_1ie`` has not
been scanned on its own.

What one sweep point is
~~~~~~~~~~~~~~~~~~~~~~~

Each point scales one or two parameter *families* on the fermentor's kinetic
model (``V406.nsk_kinetic_model``). Every member of a family is set to its
baseline value times the same multiplier, so the ratios within the family
are preserved (this is ``kinetic_optimization.expand_grouped_values`` for a
single group). The biorefinery's model specification then runs the fed-batch
fermentation inside the converged flowsheet, the TEA is solved, and the
metrics are recorded. The following are held fixed at every point:

- **Feeding.** Scenario A's policy: threshold 217.125 / target 221.25 g/L,
  600 g/L spike feed, at most **16 spikes**, ``tau_max`` 120 h. Feeding is
  **not** re-optimized per point (``perform_feeding_strategy_opt =
  False``). :doc:`feeding_strategies` explains these settings.
- **Spike count and harvest time.** Both emerge from the model as described
  on :doc:`feeding_strategies`. A point that ends at **16 spikes** is a
  culture that could have taken more sugar. A point that ends at **1 spike**
  exhausted the spike-reduction retry
  (:class:`~nskinetics.units.SpikeReduceRetry`) with sugar left over.
- **Enzyme burden.** The constraint is on, referenced to scenario A. It
  derates growth (``k_7`` and ``k_8``) by
  :math:`d = \mathrm{clip}\big((F_\mathrm{flex} - \Phi_M)/\phi_T,\ 0,\ 1\big)`,
  where :math:`\Phi_M` is the metabolic proteome the enzyme capacities
  require (:doc:`enzyme_burden`). Higher glycolytic capacity raises
  :math:`\Phi_M`. Derating starts at 2.35× glycolysis. :math:`\Phi_M` at 4×
  is 0.214 g/gDCW, below :math:`F_\mathrm{flex}` = 0.245, so no point in the
  glycolysis plane is infeasible. Points that fail are left empty.

The planes
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 18 23 20 26 13

   * - Plane
     - x-axis (40 values)
     - y-axis (40 values)
     - Anchor, spike cap
     - Status
   * - Glycolysis × ethanol inhibition
     - Glycolysis multiplier 0.001–4× on ``k_1l``, ``k_1h``, ``k_1e``
       (scenario A: 1.43 / 0.584 / 47.1 g/L/h)
     - Ethanol-inhibition multiplier 0.001–2×
     - A, 16
     - Current (build ``2b5b27d``)
   * - ``k_1e`` × ethanol inhibition
     - ``k_1e`` alone, 0.2–4× (9.42–188.4 g/L/h)
     - 0.75–1.5× (the campaigns' band)
     - A, 16
     - Current (``2b5b27d``)
   * - ``k_3`` (Pdc) × ethanol inhibition
     - 0.001–4×
     - 0.001–1.5×
     - A, 16
     - Current (``2b5b27d``), no figure
   * - ``k_6`` (Adh1) × ethanol inhibition
     - 0.001–4×
     - 0.001–1.5×
     - A, 16
     - Current (``2b5b27d``), no figure
   * - ``k_13`` (ALS) × isobutanol inhibition
     - 0.001–4× of 4.0 g/L/h (0.004–16 g/L/h)
     - Isobutanol-inhibition multiplier 0.75–1.5× on ``k_1ii``, ``k_4ii``,
       ``k_6ii``, ``k_7ii``, ``k_10ii``
     - ``opt_PI_TRY_informed``, 1
     - **No current result** (see `The isobutanol analogue`_)

**What the axes mean.** The x-axis of the glycolysis plane is the capacity
of the lumped glycolytic reaction r1. Its three members are the low-affinity,
high-affinity and acetaldehyde-activated glucose-uptake capacities. The
ethanol-inhibition multiplier scales four active coefficients, so it bundles
four mechanisms (the fifth group member, ``k_16ie``, is zero and has no
effect):

- r1 glycolysis × exp(−\ ``k_1ie``\ ·[EtOH]), with ``k_1ie`` = 0.02 L/g
- r4 acetaldehyde → acetate (AcDH) × exp(−\ ``k_4ie``\ ·[EtOH]), with 0.04 L/g
- r7 growth on glucose × exp(−\ ``k_7ie``\ ·[EtOH]), with 0.04 L/g
- r10 active-cell death × exp(``k_10ie``·([EtOH] − 100 g/L)) above
  100 g/L ethanol, with 0.04 L/g

**A larger multiplier means stronger inhibition, that is, less ethanol
tolerance.** At 115 g/L ethanol the multiplier changes the rates as follows:

.. list-table::
   :header-rows: 1
   :stub-columns: 1

   * - Multiplier
     - Glycolysis factor
     - Growth factor
     - Death-rate factor
   * - 2×
     - 0.010
     - 1e-4
     - ×3.32
   * - 1× (fitted)
     - 0.100
     - 0.0100
     - ×1.82
   * - 0.77×
     - 0.170
     - 0.0290
     - ×1.59
   * - 0.001×
     - ≈ 1
     - ≈ 1
     - ≈ 1

The y-axis is therefore **not** "growth inhibition". It moves glycolysis,
AcDH, growth and lethality together. Note too that growth (r7) draws glucose
directly, in parallel with glycolysis, not downstream of it. That structural
fact drives both the left-hand cliff and the "tolerance feeds growth" result.

.. note::

   **Conventions.** MESP is the purity-adjusted ethanol MPSP at a 15 % IRR
   × 4.45775 ($/kg → $/GGE). IRR is at default prices (ethanol 0.835 $/kg),
   and −∞ means the plant loses money outright. Titers are per litre of
   water, and yield is g per g of glucose added. The point (1, 1) is **not**
   a grid node, so the baseline values are interpolated. One grid point
   failed (glycolysis 1.847×, inhibition 0.206×). The figures fill it with
   the mean of its four neighbours. The 0.001× glycolysis column makes no
   ethanol, so its MESP is undefined. All values come from
   heat-exchanger-network build ``2b5b27d``, on which the scenario-A
   baseline MESP is 3.844 $/GGE.

Results
-------

The glycolysis × ethanol-inhibition plane
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tab-set::

   .. tab-item:: IRR

      .. figure:: /_static/images/isobutanol/kinetic_sweeps/glycolysis_x_ethanol_inhibition_irr.png
         :target: ../../_static/images/isobutanol/kinetic_sweeps/glycolysis_x_ethanol_inhibition_irr.png
         :alt: Heat map of internal rate of return over glycolysis capacity
               (x, 0 to 4 times baseline) and ethanol inhibition (y, 0 to 2
               times baseline). An orange profitable region, outlined by a
               white IRR = 0 contour, stretches diagonally from low glycolysis
               and low inhibition up to high glycolysis and high inhibition.
               Dark grey fills the lower-right triangle (high glycolysis, low
               inhibition), the far left edge and the top-left corner, where
               IRR is below -15 percent or undefined. A cyan star labelled
               profitability sits near 1.3 times glycolysis and 0.77 times
               inhibition, just below and right of the white baseline diamond
               at (1, 1). The white titer triangle sits on the lower-right
               edge of the profitable region near (2.6, 1.0), the yield
               pentagon near (3.2, 1.6), the productivity square near
               (1.0, 0.15), and the cell-density circle and cyan TCI pentagon
               on the left axis at 0 times glycolysis.
         :class: pub-figure
         :width: 100%
         :figclass: ibo-fig-measure

         **Profitability over glycolytic capacity × ethanol inhibition.**
         The colour is IRR at default prices, from −15 % to 15 %. The white
         contour is IRR = 0. Dark grey marks IRR below −15 % or undefined
         (−∞: the design loses money outright). The cyan star is the
         highest IRR (13.7 %) and the cyan pentagon the lowest TCI. The white
         markers are the baseline strain (diamond, interpolated at (1, 1)),
         titer (triangle), yield (pentagon), productivity (square) and cell
         density (circle). The TCI and cell-density markers sit on the
         0.001× glycolysis column, where essentially no ethanol is made, so
         treat them as degenerate. Look at how close the star sits to the
         diamond, and how the titer and productivity markers sit on the
         plateau's lower edge, right at the cliff.

   .. tab-item:: MESP

      .. figure:: /_static/images/isobutanol/kinetic_sweeps/glycolysis_x_ethanol_inhibition_mesp.png
         :target: ../../_static/images/isobutanol/kinetic_sweeps/glycolysis_x_ethanol_inhibition_mesp.png
         :alt: The same glycolysis capacity by ethanol inhibition grid
               coloured by MESP from 3.5 to 8 dollars per gallon of gasoline
               equivalent. A broad orange-yellow valley with diagonal white
               hatching covers most of the plane between a thin dark band at
               the far left and a dark grey triangle in the lower right; the
               hatched region, bounded by a white contour, is where MESP lies
               inside the ethanol market price range. A cyan star labelled
               MESP sits near (1.33, 0.77), close to the white baseline
               diamond at (1, 1). White markers for titer, yield and
               productivity, and a cyan TCI pentagon and white cell-density
               circle on the left axis, are placed as in the IRR version.
         :class: pub-figure
         :width: 100%
         :figclass: ibo-fig-measure

         **The same plane as MESP** (3.5–8 $/GGE; values above 8 are dark
         grey). The hatched region inside the white contour is MESP within
         the 2021–2025 ethanol market range (2.31–5.15 $/GGE). No node falls
         below 2.31, so the hatch is simply MESP ≤ 5.15 and covers 58.5 % of
         the nodes. The 0.001× glycolysis column is blank (MESP undefined).
         The star is the MESP optimum, 3.788 $/GGE at (1.334×, 0.770×).

The baseline sits near the bottom of a broad valley. Only 134 of 1,600 nodes
(8.4 %) beat it, and IRR is positive at only 689 nodes (43.1 %). The
marked points compare as follows:

.. list-table:: Marked points of the glycolysis × ethanol-inhibition plane (build ``2b5b27d``)
   :header-rows: 1
   :widths: auto
   :class: ibo-num

   * - Point
     - (glycolysis ×, inhibition ×)
     - MESP ($/GGE)
     - vs optimum
     - IRR
     - Titer (g/L) / yield (g/g) / productivity (g/L/h)
     - Spikes
     - Aeration (kmol O₂/h)
   * - Baseline (interpolated)
     - (1, 1)
     - 3.844
     - +1.5 %
     - 12.6 %
     - 114.6 / 0.460 / 1.86
     - 7
     - 13.6
   * - **MESP optimum**
     - (1.334, 0.770)
     - **3.788**
     - 0
     - **13.7 %**
     - 130.7 / 0.464 / 3.96
     - 16 (cap)
     - 21.0
   * - Titer maximum
     - (2.564, 1.026)
     - 3.845
     - +1.5 %
     - 12.6 %
     - **133.8** / 0.475 / 4.33
     - 16
     - 86.3
   * - Yield maximum
     - (3.180, 1.590)
     - 4.003
     - +5.7 %
     - 10.0 %
     - 108.5 / **0.481** / 1.44
     - 1
     - 168.7
   * - Productivity maximum
     - (1.026, 0.155)
     - 4.285
     - +13.1 %
     - 3.1 %
     - 112.5 / 0.400 / **7.15**
     - 16
     - 12.6
   * - Cell density max, TCI min
     - (0.001, 0.052) and (0.001, 0.872)
     - undefined
     - n/a
     - −∞
     - no ethanol
     - 16
     - 3.2

- **The baseline strain is already near the cost optimum.** The MESP
  optimum improves on the baseline by 1.46 %. Its worst neighbour one grid
  step away is only +0.4 %, so it is robust.
- **Single-metric optima are fragile or capital-heavy.** The titer optimum
  matches the baseline MESP, but it needs 4.1× the cost optimum's aeration
  and 11.6 MM$ more capital (TCI 148.3 vs 136.7 MM$). One grid step further,
  MESP is +19 %. The yield optimum needs 27 MM$ more capital than the
  baseline and about 12× its aeration. Its high yield comes from growth
  suppressed by burden derating (d = 0.64 there), not from better
  conversion. Capital
  tracks aeration across the plane (ρ = +0.84), because high glycolytic flux
  raises the oxygen demand of NADH re-oxidation during the aerobic stage.

More ethanol tolerance can cost more
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The response along the inhibition axis at baseline glycolytic capacity
(the x = 1.026 column) is U-shaped:

.. list-table:: MESP versus ethanol inhibition at 1.026× glycolysis (build ``2b5b27d``)
   :header-rows: 1
   :stub-columns: 1
   :class: ibo-num

   * - Inhibition ×
     - MESP ($/GGE)
     - IRR
     - Yield (g/g)
     - Titer (g/L)
     - Productivity (g/L/h)
     - Cells (g/L)
     - τ (h)
     - Spikes
     - TCI (MM$)
   * - 0.001
     - **4.667**
     - **−9.4 %**
     - 0.365
     - 102.9
     - 6.96
     - 29.5
     - 14.8
     - 16
     - 132.3
   * - 0.309
     - 4.071
     - 7.8 %
     - 0.423
     - 119.0
     - 6.93
     - 22.6
     - 17.2
     - 16
     - 132.1
   * - 0.514
     - 3.901
     - 11.4 %
     - 0.443
     - 124.8
     - 6.01
     - 18.8
     - 20.8
     - 16
     - 132.6
   * - **0.770**
     - **3.809**
     - **13.3 %**
     - 0.458
     - 128.9
     - 3.56
     - 15.2
     - 36.2
     - 16
     - 134.7
   * - 1.026
     - 3.843
     - 12.6 %
     - 0.461
     - 113.0
     - 1.74
     - 12.4
     - 65.0
     - 6
     - 139.9
   * - 1.539
     - 4.482
     - 0.9 %
     - 0.360
     - 81.2
     - 0.74
     - 9.0
     - 109.2
     - 1
     - 140.2
   * - 2.000
     - 5.190
     - −59.7 %
     - 0.282
     - 63.4
     - 0.83
     - 7.4
     - 76.9
     - 1
     - 130.5

Full de-inhibition, the classic tolerance-engineering target, makes the
strain 3.7× more productive (1.86 → 6.96 g/L/h) and cuts capital by
7.6 MM$. Yet it raises MESP by 21 % (3.84 → 4.67 $/GGE) and turns the plant
unprofitable (IRR 12.6 % → −9.4 %). The carbon accounting shows where the
gain goes. Cell density rises from 13 to 29.5 g/L and yield falls from 0.460
to 0.365 g/g: tolerance buys cell-hours rather than per-cell production.
Relief helps only down to about 0.77× the fitted coefficients. Between 1.0×
and 0.77×, lower inhibition *helps* (−0.9 %). Within 1 % of the column
minimum the inhibition spans about 0.62–1.03×, so the bottom of the U is
flat. The same U appears in every run of the four planes: the optimum sits
at 0.77–0.79× in all six, and in the four runs whose inhibition axis
reaches zero, full relief costs 18–23 % over the baseline.

The process sets the optimum: the 16-spike cap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each glycolysis column, compare the MESP-minimizing inhibition with the
inhibition at which the realized spike count first reaches its cap of 16.
In **19 of the 22** columns without growth derating (glycolysis < 2.35×),
the two are identical, and in the other three they are one grid step
apart. Across all 39 columns that make ethanol, the two lie within one grid
step in 24; every larger mismatch lies at 2.46× or above, where derating
begins. The mechanism can be read off the table above:

- **Before the cap binds** (1.026× → 0.770×), lower inhibition lets the
  culture accept more sugar. The spike count goes from 6 to 16, titer from
  113 to 129 g/L and τ from 65 to 36 h, while yield holds steady (0.461 →
  0.458).
- **Once the cap binds**, no more sugar can be delivered. Extra tolerance
  then mainly relaxes growth and AcDH, so cells rise and yield falls.

At the baseline process, then, tolerance pays only while the process can
deliver more sugar. Tolerance targets and fed-batch capacity have to be
designed together. The *location* of the optimum depends on the cap chosen
(823 of the 1,600 nodes sit at the cap, the optimum among them). The plane
has not been re-run at other caps (for example 8 and 32), so whether the
same cap-onset mechanism sets the optimum at other caps is untested; only
its location at cap 16 is shown here.

The MESP-optimal inhibition **rises** with glycolytic capacity: gently
while growth is not derated (0.77× at 0.9–1.2×, about 0.97× by 2.36×; see
the valley below), then steeply past the derating onset:

.. list-table::
   :header-rows: 1
   :stub-columns: 1

   * - Glycolysis ×
     - 0.9–1.2
     - 2.36
     - 2.5
     - 2.8
     - 3.1
     - 3.4
   * - MESP-optimal inhibition ×
     - 0.77
     - ≈ 0.97
     - 1.08
     - 1.23
     - 1.49
     - 1.69

**An iso-MESP valley.** From 1.13× to 2.36× glycolysis, raising the
inhibition from 0.77× to 0.97× keeps MESP within 0.54 % of the optimum
(3.788–3.808 $/GGE). Along the same valley, TCI rises from 135.4 to
145.2 MM$ and aeration from 14.5 to 67.6 kmol O₂/h (4.7×). Glycolytic
capacity and tolerance are thus interchangeable at equal MESP but not at
equal capital. The valley's flatness is reproducible (replicate runs agree
to about 0.05 %), but the *ordering* of designs along it is not: the
differences are tenths of a percent, comparable to the 0.25 % offset between
model builds. Do not single out one best glycolytic capacity. Beyond the
derating onset the valley floor climbs:
+1.2 % at 2.46×, +2.9 % at 2.77×, +5.7 % at 3.18×, +25 % at 3.49× and
+153 % at 4.0×.

Three cliffs, three mechanisms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 18 36 32

   * - Cliff
     - Where
     - What the sweep outputs show
     - Mechanism (evidence)
   * - Left edge
     - Glycolysis ≤ about 0.2×
     - At 0.104×: MESP 6.53–8.79 $/GGE, yield 0.18–0.24 g/g, cells
       66–125 g/L. At 0.001×: cells 205.6 g/L and essentially no ethanol.
     - Growth (r7), which draws glucose directly, out-competes a crippled
       glycolysis (read from the outputs plus the rate laws). Below 0.2× this
       is a **model-structure artifact**: a yeast without glycolysis could not
       grow on glucose, but the lumped r7 bypasses r1. Treat this region as
       outside the model's valid range.
   * - Diagonal
     - Low inhibition with high glycolysis. MESP first exceeds 5.15 $/GGE
       at 1.23× (inhibition 0.001×) and at 2.46× (inhibition 0.770×),
       along a roughly linear boundary.
     - Fast fermentations (τ 14–18 h) at the 16-spike cap whose yield
       collapses (row 0.309×: 0.42 g/g at 1.0× → 0.25 at 2.05×) without a
       biomass sink (cells about 10–14 g/L). Aeration climbs to
       199 kmol O₂/h and TCI to 170 MM$. Along the boundary itself
       (1.23–2.46×) burden derating is inactive or nearly so (d = 1 below
       2.35×, 0.95 at 2.46×); beyond 2.35× the burden also starts to derate
       growth, so the far right of each row mixes both mechanisms.
     - "Kinetic overflow": once ethanol stops throttling glycolysis, flux
       outruns the ethanol branch (**inferred**). *Where* the carbon goes is
       **not established**. One reading, a leak to acetate through the
       de-inhibited AcDH, rests on indirect evidence. The counterfactual
       probe
       ``analyses/kinetic_sweeps/probe_glycolysis_cliff_mechanism.py`` is
       built to test it, but its results have not been analysed. Treat the
       diagonal's mechanism as an open question.
   * - Right edge
     - Glycolysis ≈ 3.3–3.6× for inhibition ≳ 1.2×
     - Growth derated (d 0.46–0.60), cells 1–3 g/L, τ 73–119 h (up to
       ``tau_max``), 1–5 spikes (retry exhausted, sugar left over), yield
       0.17–0.35 g/g.
     - Proteome burden. Overexpressing glycolysis about 3.5× leaves too
       little of the flexible proteome sector for the translation
       machinery, which roughly halves growth capacity, so the
       fermentation cannot finish within 120 h. The values of d, τ and the
       spike counts are read directly from the outputs; the causal chain is
       inferred.

Sibling planes: k_1e, Pdc and Adh1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Four distinct planes, six runs in all, from the same model state tell the
same story:

.. list-table:: Sibling planes (all scenario A, 16-spike cap, build ``2b5b27d``)
   :header-rows: 1
   :widths: 30 14 26 14 16

   * - Plane (x-band × inhibition band)
     - Baseline MESP ($/GGE)
     - Plane optimum ($/GGE at x, y)
     - Best vs baseline
     - Max IRR
   * - Glycolysis × inhibition, 0.001–4 × 0.001–2 (shown above)
     - 3.844
     - 3.788 at (1.334×, 0.770×)
     - −1.46 %
     - 13.7 %
   * - Glycolysis × inhibition, 0.2–4 × 0.75–1.5 (earlier run)
     - 3.8441
     - 3.7858 at (1.272×, 0.788×)
     - −1.52 %
     - 13.7 %
   * - Glycolysis × inhibition, 0.2–4 × 0.001–2 (earlier run)
     - 3.8443
     - 3.7865 at (1.272×, 0.770×)
     - −1.50 %
     - 13.7 %
   * - ``k_1e`` × inhibition, 0.2–4 × 0.75–1.5
     - 3.8436
     - 3.7858 at (1.467×, 0.788×)
     - −1.50 %
     - 13.7 %
   * - ``k_3`` (Pdc) × inhibition, 0.001–4 × 0.001–1.5
     - 3.8437
     - 3.8128 at (0.821×, 0.770×)
     - −0.80 %
     - 13.2 %
   * - ``k_6`` (Adh1) × inhibition, 0.001–4 × 0.001–1.5
     - 3.8439
     - 3.7852 at (3.897×, 0.770×)
     - −1.53 %
     - 13.7 %

- **Pdc (**\ ``k_3``\ **) is not a lever for ethanol economics above the baseline.**
  Across 1–4× ``k_3``, the per-row MESP spread has a median of 0.008 % and
  a maximum of 0.58 %. Knocking Pdc down below about 0.3–0.9× (depending
  on inhibition) is a cliff.
- **Adh1 (**\ ``k_6``\ **) matters only for de-inhibited strains.** At 0.001×
  inhibition, raising Adh1 from 1× to 4× lowers MESP from 4.694 to
  3.998 $/GGE and raises yield from 0.363 to 0.423 g/g. It cuts the full
  de-inhibition penalty from 23.1 % to 5.6 %. At 0.77× inhibition the same
  change moves MESP only from 3.811 to 3.786 $/GGE.

.. figure:: /_static/images/isobutanol/kinetic_sweeps/k1e_x_ethanol_inhibition_mesp.png
   :target: ../../_static/images/isobutanol/kinetic_sweeps/k1e_x_ethanol_inhibition_mesp.png
   :alt: Heat map of MESP over the glycolytic capacity k_1e (x, about 10 to
         190 grams per litre per hour) and ethanol inhibition (y, 0.75 to 1.5
         times baseline). Almost the whole plane is orange-yellow, between
         about 3.8 and 4.3 dollars per gallon of gasoline equivalent. A dark
         blue band and a dark grey region at k_1e above about 150 to 175 mark
         a steep cliff on the right. A cyan star labelled MESP sits near
         k_1e 69 and inhibition 0.79. White markers for titer near (138, 0.88),
         yield near (157, 1.23) and productivity near (110, 0.75), and a
         white cell-density circle and cyan TCI pentagon in the bottom-left
         corner.
   :class: pub-figure
   :width: 100%
   :figclass: ibo-fig-measure

   **MESP over** ``k_1e`` **× ethanol inhibition.** Only ``k_1e`` is scaled;
   ``k_1l`` and ``k_1h`` stay at baseline. The x-axis is absolute (0.2–4× of
   47.1 g/L/h) and the inhibition axis covers the campaigns' 0.75–1.5× band.
   The colour runs from 3.75 to 5.25 $/GGE, with values above that in dark
   grey. Markers are as in the glycolysis plane. The plane optimum (star) is 3.7858 $/GGE
   at 1.467× ``k_1e`` and 0.788× inhibition, 1.50 % below this plane's
   baseline. The baseline itself (47.1 g/L/h, 1×) is not marked. Note how
   flat the orange plateau is before the right-hand cliff.

The isobutanol analogue
~~~~~~~~~~~~~~~~~~~~~~~

The same question for the second product sweeps the Ehrlich-entry capacity
``k_13`` (acetolactate synthase, ALS) against the grouped isobutanol
inhibition. Currently there is **no result to report**. The script
``evaluate_EtOH_k13_inhib_isobutanol.py`` is now anchored on the
co-production design ``opt_PI_TRY_informed`` (trial #1912; trial numbers on
these pages are the 0-based CSV ``trial_number`` used by ``scenarios.py``),
and it has not been run on that anchor yet. Earlier runs around a
since-removed anchor design are superseded and are not reported.
:doc:`kinetic_optimization` covers the co-production designs themselves.

Caveats and limitations
-----------------------

- **Model build.** Every number on this page comes from heat-exchanger-network
  build ``2b5b27d`` (baseline 3.844 $/GGE). The current scenario-A pin is
  from build ``6d4776f`` (0.86021 $/kg, or 3.835 $/GGE), an offset of about
  0.25 %. The offset is smaller than every claim here except the ordering
  along the 0.54 % valley. Do not compare these values with the uncertainty
  analysis's MESP as if they were one basis.
- **Fixed feeding policy.** Every point uses scenario A's policy with a
  16-spike cap, and feeding was not re-optimized. "Little headroom" is
  therefore conditional on the baseline feeding policy. The location of the
  tolerance optimum depends on the cap, which was not varied.
- **A grouped inhibition axis.** The y-axis moves glycolytic, AcDH, growth
  and lethality coefficients together and cannot separate them. Scanning the
  four coefficients one at a time has not been done.
- **Model structure.** Growth (r7) consumes glucose directly, so the left
  edge below 0.2× is not biologically meaningful.
- **Diagonal-cliff mechanism.** The sweep CSVs contain no pyruvate,
  acetaldehyde or acetate columns. The overflow reading is inferred, and the
  counterfactual probe that could settle it has not been analysed.
- **Discreteness.** Integer spike counts and discrete τ rows cause jitter of
  about 0.1–1 % in the stalled, 1-spike region. Do not read meaning into
  those wiggles. "A flat, curved valley bounded by three steep cliffs"
  describes the surface better than "non-monotonic".
- **Reproducibility margin.** Interpolating the glycolysis plane onto the two earlier
  glycolysis runs gives a median relative difference of −0.05 % on plateau
  points (larger only at cliff edges, where it is an interpolation effect).
  **Do not interpret MESP differences below about 0.1 % between runs.**
- **Degenerate markers.** The cell-density and TCI "optima" of the glycolysis plane sit on
  the no-ethanol column. Which grid row each lands on is a numerical
  tie-break.

.. dropdown:: What this plane does and does not show

   - The x-axis is **glycolytic capacity**, not Pdc expression. The Pdc
     plane is flat above the baseline.
   - The MESP optimum lies at *less* inhibition than fitted (0.77×), and
     full relief is the costly case; stronger ethanol inhibition does not
     improve the economics.
   - Growth-specific inhibition has no current sweep result. The
     ``k_1e`` × ``k_7ie`` plane (``evaluate_EtOH_k1e_k7ie.py``) has only
     been run on an earlier version of the kinetic model, and that run is
     qualitative only.
   - ``evaluate_EtOH_k13_k7ii.py``, ``evaluate_EtOH_k1e_k1ie.py`` and
     ``evaluate_EtOH_k1e_k7ie.py`` are exploratory planes; no results from
     them are presented here.

Reproduce it
------------

No analysis results ship with NSKinetics: every sweep CSV must be
regenerated, or ``NSKINETICS_IBO_RESULTS_DIR`` must point to a results folder
that already holds it. The sweeps are long, opt-in simulations. Each of the
1,600 points is a converged flowsheet simulation plus a TEA solve, so a
plane takes roughly 7–11 hours (plan for overnight). Run one simulation
process at a time (:doc:`running_analyses`). The sweep scripts are in
``analyses/kinetic_sweeps/`` and the plotters in ``plots/`` of the package
folder.

.. list-table::
   :header-rows: 1
   :widths: 44 38 18

   * - Script
     - What it does and writes
     - Checkpointed
   * - ``evaluate_EtOH_glycolysis_inhib_ethanol.py``
     - Glycolysis × ethanol-inhibition sweep. One CSV per metric in the results root, e.g.
       ``ibo_(40, 40, 1)_glyco_inhib_Spike_rb0.001-4_ib0.001-2_opt=False_max_n=16__<metric>.csv``.
     - yes
   * - ``evaluate_EtOH_k1e_inhib_ethanol.py``
     - ``k_1e`` plane.
     - **no**
   * - ``evaluate_EtOH_k3_inhib_ethanol.py``, ``evaluate_EtOH_k6_inhib_ethanol.py``
     - Pdc and Adh1 planes.
     - yes
   * - ``evaluate_EtOH_k13_inhib_isobutanol.py``
     - ``k_13`` × isobutanol inhibition on ``opt_PI_TRY_informed`` (writes
       ``max_n=1`` CSVs). Needed before its figure can be drawn.
     - yes
   * - ``plot_glycolysis_inhib_ethanol_A.py``, ``plot_k1e_inhib_ethanol_A.py``
     - The figures above (both glycolysis-plane twins come from the first). They read
       the CSVs only and never import the package, so they are safe to run
       beside a simulation. Drawing code: ``_kinetic_sweep_figure.py``.
       Output: ``publication/Kinetic-sweeps/*.{png,pdf}``.
     - n/a (seconds)
   * - ``plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py``
     - Figure for the rerun ``k_13`` sweep. It fails until that sweep has run.
     - n/a

Two environment variables control a sweep script.
``IBO_SWEEP_REPLOT_FROM_CSV=1`` skips every simulation and redraws from
saved CSVs. ``IBO_SWEEP_FRESH=1`` discards a leftover checkpoint and starts
over.

Checkpointed sweeps and the supervisor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A native integrator crash (a CVODE segfault: exit code 5, no traceback)
kills a sweep process outright. One such crash killed the first ``k_13`` run
at about point 764 of 1,600. The checkpointed sweeps survive this:

- **Checkpoint CSV** ``<script stem>_checkpoint.csv``. The sweep writes one
  row per evaluated point, holding the grid indices, the axis values, a state
  (``OK``, ``ERROR`` or ``LOST``) and every metric. On relaunch, the script
  checks the checkpoint's axis values against its own grid (a checkpoint from
  other bounds raises an error) and resumes. The checkpoint is deleted when
  the script finishes, plots included.
- **In-flight sidecar** ``<script stem>_inflight.json``. The sweep writes it
  just before each point's simulation and clears it afterwards. If a relaunch
  finds it, the previous process died in that point. The relaunch logs the
  point as ``LOST`` (all NaN) and steps past it, so one bad point cannot
  crash the sweep forever.
- **Supervisor** ``supervise_sweep.py``. This stdlib-only script never
  imports the model. It relaunches the sweep in a fresh process each time the
  child exits non-zero, strictly one process at a time. It kills a child
  whose in-flight sidecar is older than ``--stall-timeout-min`` (default 10;
  a hung integrator call). Loads and plotting write no sidecar, so they are
  never killed. It gives up after ``--max-empty-attempts`` (default 3)
  consecutive attempts that add no checkpoint row, or after
  ``--max-attempts`` (default 25) attempts in all. The Sobol' sampler uses
  the same supervisor with ``--stem`` and ``--checkpoint-suffix
  _trajectory.csv`` (:doc:`sensitivity`).

.. tab-set::

   .. tab-item:: Bash

      .. code-block:: bash

         # from <package folder>/analyses: run under the supervisor (one sweep at a time)
         python supervise_sweep.py kinetic_sweeps/evaluate_EtOH_glycolysis_inhib_ethanol.py

         # redraw the figures from existing CSVs (no simulation)
         python ../plots/plot_glycolysis_inhib_ethanol_A.py
         python ../plots/plot_k1e_inhib_ethanol_A.py

   .. tab-item:: PowerShell

      .. code-block:: powershell

         # from <package folder>\analyses
         python supervise_sweep.py kinetic_sweeps\evaluate_EtOH_glycolysis_inhib_ethanol.py

         python ..\plots\plot_glycolysis_inhib_ethanol_A.py
         python ..\plots\plot_k1e_inhib_ethanol_A.py

The sweep prints a message each time it logs a point as ``LOST``. A lost
point becomes a NaN cell in the CSVs and figures, just like a point that
fails the enzyme-burden or fed-batch-volume guards.

Next: :doc:`kinetic_optimization` opens the design space to isobutanol
co-production and searches it with Bayesian optimization.
