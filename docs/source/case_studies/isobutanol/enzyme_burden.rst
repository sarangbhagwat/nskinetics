Enzyme burden
=============

A kinetic model's enzyme capacities look like free dials, but each one is a
specific activity per unit biomass, that is, enzyme amount × turnover. Raising
a capacity means expressing more enzyme out of a finite proteome. Left
unconstrained, an optimizer raises every capacity at once and proposes a strain
no cell can make. This page asks **which strain designs are physiologically
expressible, and how proteome spending shapes the economics.** The isobutanol
biorefinery answers with a closed-form **proteome-allocation constraint**. The
enzyme pools a design needs are summed and closed against the cell's flexible
protein sector, and growth is derated as that sector fills. The constraint is
enforced at a single choke point, so every analysis inherits it. These
analyses include the feeding and kinetic sweeps, the uncertainty runs, the
Sobol' analysis and the optimization campaigns.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 0.1347
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      g/gDCW: the **penalty-free metabolic budget**. The wild-type strain
      uses 0.0701. Above the budget, growth slows linearly.

   .. grid-item-card:: 157 % / 130 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Budget used by the isobutanol- and ethanol-**yield** optima. Growth is
      derated to ×0.30 and ×0.63, and cells fall to 2.2 and 4.2 g/L.

   .. grid-item-card:: 97 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Budget used by the most profitable design found (#1912, IRR 27.3 %),
      with **no** growth derating

   .. grid-item-card:: η² = 0.39
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Share of the profitability variance that Φ\ :sub:`M` explains on its
      own (screening prior). Every PI > 0 design in those samples has
      Φ\ :sub:`M` ≤ 0.17.

.. admonition:: In one sentence
   :class: tip

   Yield optima win by overloading the proteome to throttle growth. The dilute,
   slow broth that results is what costs the capital. The most profitable
   designs spend their proteome close to the penalty-free budget, and every
   profitable design in the space-filling samples stays at Φ\ :sub:`M`
   ≤ 0.17 g/gDCW.

At a glance
-----------

.. list-table::
   :widths: 26 74
   :stub-columns: 1

   * - Question
     - Which designs can a cell express, and what does proteome spending do
       to growth, titer and profit?
   * - Module
     - :mod:`~nskinetics.biorefineries.isobutanol.enzyme_burden`
       (stdlib-only: imports neither numpy, biosteam nor nskinetics, so it
       can also be loaded by file path)
   * - Inputs
     - The kinetic capacities of a design (``k_1*``, ``k_2``, ``k_3``,
       ``k_5*``, ``k_6``, ``k_13``–``k_17``) and the growth capacities
       ``k_7``/``k_8``
   * - Output
     - The metabolic pool Φ\ :sub:`M`, the translation sector φ\ :sub:`T`, the
       growth-derating factor *d*, and feasibility (Φ\ :sub:`M` <
       F\ :sub:`flex`)
   * - Where it acts
     - Inside ``system.load_simulate``, once per simulation; in the optimizer,
       also as a pre-simulation feasibility check
   * - Calibration
     - Always against scenario A (wild type), where it is exactly inert
       (Φ\ :sub:`M` 0.0701, *d* = 1)

Method
------

How the proteome is divided
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cell's protein (0.49 g per g dry cell weight) is split in half. One half
is a growth-invariant **housekeeping** sector. The other is a **flexible**
sector shared by the translation machinery and the modelled metabolic
enzymes. The bar below is drawn roughly to scale:

.. grid:: 1
   :gutter: 1

   .. grid-item-card:: Housekeeping
      :columns: 12 12 6 6
      :text-align: center

      **0.245** g/gDCW

      Growth-invariant half of cell protein. It cannot be reallocated.

   .. grid-item-card:: Translation φ\ :sub:`T`
      :columns: 12 6 3 3
      :text-align: center

      **0.110** g/gDCW

      Active ribosomes at wild-type growth

   .. grid-item-card:: Metabolic budget
      :columns: 12 6 3 3
      :text-align: center

      **0.1347** g/gDCW

      Penalty-free room for the modelled enzymes, Φ\ :sub:`M`

Translation plus the metabolic budget make up the flexible sector,
F\ :sub:`flex` = 0.49 × (1 − 0.50) = **0.245** g/gDCW. A design may spend
beyond the budget, but only by taking room from translation, which slows
growth. At Φ\ :sub:`M` = F\ :sub:`flex` growth stops and the design is
**infeasible**.

From capacities to grams of enzyme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each modelled reaction step is converted to a pool in g enzyme/gDCW by one of
two routes:

.. tab-set::

   .. tab-item:: Native steps (ratio route)

      A native step is charged by ratio to its measured wild-type pool, so the
      constraint is exactly inert at the reference capacities. A
      co-regulated block is sized by its most-raised term:

      .. math::

         \mathrm{pool}_i = \mathrm{pool}_{\mathrm{wt},i}\;
         \max_j \frac{k_{ij}}{k_{ij,\mathrm{ref}}}

      .. list-table:: Wild-type pools at protein content 0.49
         :header-rows: 1
         :widths: 12 38 22 28

         * - Step
           - Enzyme(s)
           - Pool (g/gDCW)
           - Capacities charged
         * - r1
           - Glycolysis lump (17 genes)
           - 0.0479
           - ``k_1h``, ``k_1l``, ``k_1e``
         * - r2
           - PDH complex + TCA
           - 0.0035
           - ``k_2``
         * - r3
           - Pdc1 (+ Pdc5/6)
           - 0.0093
           - ``k_3``
         * - r4
           - Ald6 (fixed pool)
           - 0.0035
           - none
         * - r5
           - Acs2 (+ Acs1)
           - 0.00087
           - ``k_5``, ``k_5e``
         * - r6
           - Adh1 (+ Adh2–5)
           - 0.0044
           - ``k_6``
         * - r17
           - Adh6 (native, constitutive)
           - 0.00066
           - ``k_17``

      The pools come from quantitative proteomics (per-gene molecules per
      cell, cross-checked against PaxDB). They were tabulated at a protein
      content of 0.45 and are rescaled to 0.49 at import. The r4 pool is fixed
      because the model already carries the acetaldehyde-dehydrogenase pool
      as its own state variable. Adh6 is charged by abundance, not flux,
      because its expression is constitutive.

   .. tab-item:: Heterologous Ehrlich steps (kcat route)

      The engineered Ehrlich steps have no wild-type capacity to scale from.
      They are converted with turnover numbers and molecular masses at an
      average in-vivo saturation σ\ :sub:`eff` = 0.5:

      .. math::

         \mathrm{pool}_i = \frac{k_i\, M_{\mathrm{enz}}}
         {k_{\mathrm{cat}}\,\sigma_{\mathrm{eff}}\cdot 3600\cdot M_{\mathrm{sub}}}

      .. list-table:: Enzyme cost per unit capacity (g enzyme/gDCW per g/gDCW/h)
         :header-rows: 1
         :widths: 12 48 40

         * - Step
           - Enzyme
           - Cost per unit ``k``
         * - r13
           - ALS, Ilv2 + Ilv6
           - 0.0056
         * - r14
           - Ilv5 (*N. crassa* kcat as surrogate)
           - 0.0136
         * - r15
           - Ilv3 (*A. fumigatus* kcat as surrogate)
           - 0.0138
         * - r16
           - Aro10
           - 0.0180

      The campaigns move the three downstream steps together, as
      ``ehrlich_downstream``. Within its band, that group alone can take up to
      0.174 g/gDCW (see Results). σ\ :sub:`eff` scales every Ehrlich pool
      inversely, which makes it the largest single lever on the cost of a
      high-isobutanol strain.

Closure and growth derating
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pools sum to Φ\ :sub:`M`. The translation sector scales with the larger of
the two demanded growth capacities, and growth is derated linearly once
Φ\ :sub:`M` crowds into it:

.. math::

   \Phi_M = \sum_i \mathrm{pool}_i, \qquad
   \varphi_T = 0.110 \cdot \max\!\left(\frac{k_7}{k_{7,\mathrm{ref}}},
   \frac{k_8}{k_{8,\mathrm{ref}}}\right), \qquad
   d = \mathrm{clip}\!\left(\frac{F_\mathrm{flex} - \Phi_M}{\varphi_T},\,0,\,1\right)

.. math::

   k_{7,\mathrm{eff}} = d\,k_7, \qquad k_{8,\mathrm{eff}} = d\,k_8, \qquad
   \text{feasible} \iff \Phi_M < F_\mathrm{flex}

In the optimization campaigns ``k_7`` and ``k_8`` are not decision variables,
so φ\ :sub:`T` stays at 0.11025. Derating then starts at Φ\ :sub:`M` > 0.1347,
and *d* reaches 0 at Φ\ :sub:`M` = 0.245. Only capacities are charged.
Inhibition and affinity parameters carry no proteome cost.

.. list-table:: The four constants
   :header-rows: 1
   :widths: 30 12 58

   * - Constant
     - Value
     - Basis
   * - ``PROTEIN_CONTENT``
     - 0.49
     - Crude protein of brewery *S. cerevisiae* biomass (Jach et al. 2022)
   * - ``HOUSEKEEPING_FRACTION``
     - 0.50
     - Growth-invariant half of the proteome (Scott et al. 2010)
   * - ``TRANSLATION_FRACTION_WT``
     - 0.225
     - Active ribosomal sector of fast-growing yeast: ≈ 30 % ribosomal
       protein, of which ≈ 25 % is a non-translating reserve (Metzl-Raz et al.
       2017). This is the **looser** of the two options considered.
   * - ``SIGMA_EFF``
     - 0.50
     - In-vivo saturation fitted in the GECKO enzyme-constrained yeast model
       (Sánchez et al. 2017)

Where it is enforced
~~~~~~~~~~~~~~~~~~~~

There is exactly **one** derating site: the simulate choke point
``system.load_simulate``.

- The ``_apply_enzyme_burden()`` context manager wraps the whole convergence
  loop. It derates ``k_7``/``k_8`` once on entry and restores them in a
  ``finally``. It is re-entrant: the outermost application owns the
  derate/restore.
- A guard on the fermentor ``V406`` routes every kinetic run through the same
  context. That covers a bare ``V406.simulate()``, a system re-simulation and
  the fed-batch specification's own runs, so each is derated exactly once.
- An over-cap point raises ``BurdenInfeasibleError`` (a ``ValueError``
  subclass, re-exported as ``system.EnzymeBurdenInfeasibleError``). Sweeps
  turn it into a blank (NaN) cell. The optimizer checks the same predicate
  **before** simulating and logs the trial as ``INFEASIBLE``.
- Burden studies record the pools, Φ\ :sub:`M`, φ\ :sub:`T`,
  F\ :sub:`flex`, the derating factor and the effective ``k_7``/``k_8`` in
  every trajectory row (``BURDEN_COLUMNS``).

**Policy per scenario** lives in ``ScenarioSpec.burden_default``, and
``load_scenario(..., burden=None)`` applies it. The ``BurdenModel`` is
**always** built from scenario A's live kinetics:

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - Scenario
     - Burden
     - Why
   * - A (wild type)
     - on
     - Exactly inert: Φ\ :sub:`M` 0.0701, *d* = 1
   * - B (engineered co-producer)
     - **off**
     - Its Ehrlich constants are infeasible (Φ\ :sub:`M` 0.286 > 0.245).
       The cell could carry at most 0.81× of them.
   * - ``opt_*``
     - on
     - As during the campaigns that produced them

Results
-------

Yield optima overload the proteome to throttle growth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The returned designs of the seven seed-350 campaigns and the TRY-informed
relay (see :doc:`kinetic_optimization`) span the whole range, from well inside
the budget to far beyond it:

.. list-table:: Proteome use of the campaign optima
   :header-rows: 1
   :widths: 28 14 14 14 14 16

   * - Design
     - Φ\ :sub:`M` (g/gDCW)
     - Φ\ :sub:`M` / budget
     - Growth factor *d*
     - Cells (g/L)
     - IRR
   * - Starting strain (scenario A)
     - 0.0701
     - —
     - 1.00
     - 12.97
     - 12.6 %
   * - **TRY-informed #1912**
     - 0.1301
     - **0.97**
     - **1.00**
     - 14.19
     - **27.3 %**
   * - Uninformed #103
     - 0.1526
     - 1.13
     - 0.84
     - 7.63
     - 15.8 %
   * - **Isobutanol yield #841**
     - 0.2121
     - **1.57**
     - **0.30**
     - 2.21
     - 1.7 %
   * - Isobutanol titer #84
     - 0.1156
     - 0.86
     - 1.00
     - 14.45
     - none
   * - Isobutanol productivity #1472
     - 0.1351
     - 1.00
     - 0.997
     - 10.91
     - none
   * - **Ethanol yield #1923**
     - 0.1751
     - **1.30**
     - **0.63**
     - 4.19
     - 12.1 %
   * - Ethanol titer #1820
     - 0.1233
     - 0.91
     - 1.00
     - 13.56
     - 12.7 %
   * - Ethanol productivity #365
     - 0.1126
     - 0.84
     - 1.00
     - 13.37
     - 8.1 %

Trial numbers are 0-based CSV indices. IRRs are as recorded under the
campaigns' build, 2b5b27d. Φ\ :sub:`M`, *d* and cell density are
fermentation-level and do not depend on that build.

- **Yield optima trade growth for yield.** Less carbon goes to biomass and more
  to product, which is how the yield optima reach 93–95 % of theoretical
  yield. The cost is a dilute, slow broth: τ is 102 h and 51 h, and
  isobutanol-yield #841 needs **2.3×** the starting strain's capital (316 vs
  139 MM$).
- **The next trial shows the price.** In the very next trial, #842, the
  isobutanol-yield campaign spent less proteome (Φ\ :sub:`M` 0.212 → 0.173,
  growth ×0.30 → ×0.65). It gave up 7 % of the yield (0.381 → 0.354 g/g),
  ran in 31 h instead of 102 h, and returned 22.9 % IRR instead of 1.7 %. The
  campaign scored only yield, so it returned #841.
- **Profit optima sit at the budget.** The two equally profitable co-production
  designs from different campaigns both sit at the budget: #1912 at 0.97
  (no derating) and the seed-33960 optimum #1602 at 1.006 (growth ×0.993).
  The uninformed ethanol optimum overshoots mildly: 113 % of the budget, 16 %
  growth derating and 7.6 g/L cells, and it still beats the starting strain.
- **The flagship spends its proteome on the right pathway.** At #1912,
  0.0607 g/gDCW (**47 %** of Φ\ :sub:`M` 0.130) is in the isobutanol pathway,
  of which Adh6 is 0.0133. Glycolysis takes 0.047 and the ethanol branch
  (Pdc + Adh1) 0.0145. On the Pdc slice through #1912, raising Pdc pushes the
  design over the budget from ``k_3`` 4.85 g/L/h (0.835× wild type), and the
  derating factor falls to 0.734 at the top of the band
  (:doc:`objective_landscape`).

.. figure:: /_static/images/isobutanol/burden/returned_designs_irr_titer_proteome.png
   :target: ../../_static/images/isobutanol/burden/returned_designs_irr_titer_proteome.png
   :class: pub-figure
   :alt: Three aligned panels, one row per campaign. Left: IRR strips of each
         campaign's trials, with open best-visit and filled returned-design
         markers. Middle: stacked ethanol and isobutanol titers of the
         returned designs. Right: stacked metabolic proteome of the returned
         designs by sector (glycolysis, TCA/acetate, Pdc, Adh1, ALS to Aro10,
         Adh6), with a dotted line at the 0.1347 budget and growth factors
         x0.84, x0.63, x0.30 and x0.997 printed beside the bars that cross it.
   :width: 100%

   **What each campaign returned, and what it cost in proteome** (panels c–e
   of the campaign figure on :doc:`kinetic_optimization`; seed 350 and the
   TRY-informed relay). **(c)** IRR of every trial (strip), the highest-IRR
   visit (open) and the returned design (filled). **(d)** Ethanol (amber) and
   isobutanol (violet) titers of the returned designs. **(e)** Their metabolic
   proteome Φ\ :sub:`M` by sector: glycolysis and TCA/acetate (grey), Pdc and
   Adh1 (amber; Adh1 hatched), ALS→Aro10 and Adh6 (violet; Adh6 hatched).
   Dotted: the penalty-free budget, 0.1347 g/(g DCW). The translation sector
   (φ\ :sub:`T` = 0.110) and the housekeeping sector are omitted. The four
   designs beyond the budget grow at ×0.84, ×0.63, ×0.30 and ×0.997 of their
   penalty-free rate. **What to look at:** the two longest proteome bars
   belong to the yield optima, the two most growth-derated returned designs.
   The isobutanol-yield optimum also returns only 1.7 % IRR in (c). The
   TRY-informed design (★) splits its proteome between the two branches and
   stays inside the budget.
   Trial numbers in the source figure are 1-based; IRRs are under build
   2b5b27d.

Proteome spending separates failing designs from merely poor ones
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Sobol' samples of the :doc:`sensitivity` analysis draw designs without an
optimizer. On them, the **single scalar** Φ\ :sub:`M` explains much of the
profitability variance under the two engineering-prior measures:

.. list-table:: Φ\ :sub:`M` and profitability on space-filling samples
   :header-rows: 1
   :widths: 40 20 20 20

   * - Statistic
     - Screening (SCR)
     - Custom (CUS)
     - Campaign (CAM)
   * - η²(PI | Φ\ :sub:`M`), 20 quantile bins
     - **0.39**
     - 0.47
     - 0.012
   * - Designs growth-derated (Φ\ :sub:`M` > 0.135)
     - 63.9 %
     - 81.1 %
     - 25.7 %

- **Every PI > 0 design is frugal.** All 34 PI > 0 designs (SCR and CUS) have
  Φ\ :sub:`M` ≤ 0.170, under 70 % of the flexible proteome, and their median
  *d* is 1.0. Their glycolysis lies between 0.45× and 1.87× wild type, with a
  median at or below wild type.
- **Two levers dominate the pool.** Glycolysis and the downstream Ehrlich
  group carry 95–97 % of the variance of the sampled pools. Each can use
  about 70–78 % of F\ :sub:`flex` on its own; their maximum contributions over
  the bands are 0.191 and 0.174 g/gDCW. PI is hump-shaped in glycolysis and
  falls at high downstream-Ehrlich expression. The classic target of "more
  glycolytic capacity" has its largest economic effect as a *penalty*.
- **The ranking depends on the prior.** Under the campaign measure (CAM),
  variance is dominated by whether the native ethanol pathway is knocked out,
  and Φ\ :sub:`M` explains almost nothing (η² 0.012). Always name the measure
  next to a variance share.

The same mechanism shapes the ethanol-only plane of
:doc:`kinetic_sweeps`. Along glycolytic capacity at scenario A's fixed feeding,
burden derating starts at **2.35×** and causes the right-edge cliff near
3.5× (growth factor ≈ 0.46–0.60, cells 1–3 g/L, τ → 120 h). Even at 4×,
Φ\ :sub:`M` is 0.214 < 0.245, so no point on that plane is infeasible.

Caveats and limitations
-----------------------

- **Model-contingent.** The mechanism rests on a linear derating law and four
  constants. σ\ :sub:`eff` is the largest lever on the Ehrlich cost, and
  ``TRANSLATION_FRACTION_WT`` = 0.225 was chosen as the **looser** option. A
  conservative alternative (housekeeping 0.58, translation 0.22) would lower
  the hard cap from 0.245 to 0.206, but it was not adopted. **No
  sensitivity** of the results to σ\ :sub:`eff`, F\ :sub:`flex` or the
  housekeeping fraction has been run; it is an open follow-up.
- **One-directionally permissive.** Φ\ :sub:`M` covers only the modelled
  steps. Amino-acid, nucleotide and lipid synthesis, the full respiratory
  chain and the transporters carry no burden term. The constraint therefore
  lets an engineered pathway expand into proteome that a real cell has
  occupied. Isobutanol-specific couplings are not charged either: NADPH
  supply for Ilv5/Adh6, the competing valine/isoleucine use of Ilv2/Ilv5/Ilv3,
  and thiamine-diphosphate demand.
- **Surrogate kinetics.** The Ilv5 and Ilv3 costs use fungal ortholog
  turnover numbers (*N. crassa*, *A. fumigatus*); no *S. cerevisiae* values
  exist.
- **Capacity units.** The model file and the parameter-set figures write
  capacities in g/L/h, while the burden treats them as biomass-specific. A
  harmonization of the two is still open.
- **Scenario B runs burden-free** because its literature Ehrlich point is
  infeasible. B's numbers are therefore not comparable with burden-on
  results without saying so.
- **Other figures.** The proteome-bar and Voronoi-treemap plotters
  (``plots/plot_kin_opt_parameter_sets.py``, ``plots/plot_proteome_voronoi.py``)
  have only been run on older campaigns, so their figures are not shown. The
  figure above reflects the current campaigns.

Reproduce it
------------

The burden module itself needs no simulation and no results files. Importing
it through the package needs the optional ``isobutanol`` dependencies
(``pip install "nskinetics[isobutanol]"``; see :doc:`running_analyses`), but
builds nothing:

.. code-block:: python

   from nskinetics.biorefineries.isobutanol import enzyme_burden as eb

   print(eb.F_FLEX, eb.PHI_T_WT)               # 0.245 0.11025
   print(round(eb.F_FLEX - eb.PHI_T_WT, 4))    # 0.1347, the penalty-free budget
   print(round(eb.ehrlich_unit_cost('r13'), 4))   # 0.0056 g enzyme/gDCW per unit k_13

:meth:`BurdenModel.from_reference(...)
<nskinetics.biorefineries.isobutanol.enzyme_burden.BurdenModel.from_reference>`
builds the model from scenario-A kinetics; ``evaluate`` and
``describe_point`` then report a design's pools, Φ\ :sub:`M`, φ\ :sub:`T`
and *d*. In a full biorefinery session,
:func:`isobutanol.load_scenario(name)
<nskinetics.biorefineries.isobutanol.scenarios.load_scenario>` installs the
scenario's burden policy automatically (see :doc:`/API/biorefineries`).

The tables above come from the per-trial burden columns of the campaign and
Sobol' trajectory CSVs. No analysis results ship with NSKinetics: regenerate
them (:doc:`kinetic_optimization`, :doc:`sensitivity`), or point
``NSKINETICS_IBO_RESULTS_DIR`` at a results folder that holds them. The
figure above is a crop of the campaign figure drawn by
``plots/optimization_figure_ideas/fig_main_scouts.py``.
