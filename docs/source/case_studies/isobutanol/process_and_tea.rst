.. _isobutanol-process-and-tea:

Process and techno-economics
============================

What does a corn dry-grind biorefinery with a *kinetic* fed-batch
fermentation cost and earn — for the wild-type ethanol strain (scenario A)
and for an engineered strain that co-produces isobutanol (scenario B and the
``opt_*`` designs) — and where do its capital and operating costs go? This
page describes the flowsheet and the TEA basis every other analysis in this
case study inherits, then reads the baseline numbers and the cost and revenue
breakdowns of the strain-design optima.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 19 unit groups
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      reconciled so every priced feed and utility cash flow is counted
      exactly once; together they sum to the annual operating cost.

      +++
      Cost and utility bookkeeping

   .. grid-item-card:: ~10 %
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      of installed equipment cost is fermentation, while corn is ~86 % of
      material cost — which is why productivity barely moves MESP.

      +++
      Scenario A, build 6d4776f

   .. grid-item-card:: +19 % sales
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      for the flagship co-producer (#1912) vs scenario A, at −1.2 % TCI and
      +0.2 % operating cost: the gain is all revenue.

      +++
      Campaign optima, build 2b5b27d

   .. grid-item-card:: 2.27× capital
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      for the isobutanol-yield optimum (TCI 315.6 vs 139.3 MM$), despite
      the highest sales of any design.

      +++
      Campaign optima, build 2b5b27d


The flowsheet
-------------

The biorefinery is assembled on top of the BioSTEAM ``biorefineries.corn``
dry-grind ethanol model. Corn handling, milling, liquefaction and
saccharification come from that model unchanged; from the saccharified slurry
onward the process is replaced by the NSKinetics fed-batch section, a
solvent-free ethanol + isobutanol separation train, and a high-rate
wastewater treatment, boiler-turbogenerator and facilities layer.

.. figure:: /_static/images/isobutanol/overview/block_flow_diagram.png
   :target: ../../_static/images/isobutanol/overview/block_flow_diagram.png
   :class: pub-figure
   :width: 100%
   :alt: Block flow diagram. Panel A, fermentation process: saccharified corn slurry is split into two multi-effect evaporation and dilution trains giving an initial feed and a spike feed to the fed-batch fermentor, which receives compressed air and yeast; vented gases are scrubbed and the recovered alcohols join the beer sent to distillation. Panel B, separation process: beer column, ethanol rectifier, molecular sieves and denaturant mixing give ethanol; stripper, decanter and drying column give isobutanol above 99.9 wt%; stillage goes through centrifugation, evaporation, oil separation and drying to DDGS and corn oil.

   Block flow diagram at the baseline separation gate (all broth to the
   ethanol + isobutanol train). **A, fermentation:** follow the slurry
   through the splitter into the two conditioning trains — the upper one
   makes the initial charge, the lower one the concentrated spike feed —
   into the fermentor; the vent scrubber recovers stripped alcohols into the
   beer sent to distillation. **B, separation:** the top row makes fuel ethanol, the dashed
   middle row (it runs only when isobutanol is made) makes isobutanol, and
   the bottom rows turn the stillage into DDGS and corn oil. Not drawn: the
   idle ethanol-primary train, heat exchangers, pumps and tanks, wastewater
   treatment, the boiler-turbogenerator and utilities. Drawn by
   ``plots/plot_block_flow_diagram.py``.

.. list-table:: Process sections
   :header-rows: 1
   :widths: 24 34 42

   * - Section
     - Main units
     - Role
   * - Feedstock handling and saccharification (corn model)
     - grain handling ``MH101``, storage ``V102``, hammer mill, slurry
       mixing, jet cooker, liquefaction, saccharification
     - Receives 56,972 kg/h of wet corn (72 % starch, dry basis) and turns
       the starch into a glucose-rich slurry.
   * - Sugar-solution preparation
     - splitter ``S301``; multi-effect evaporators ``F301`` (initial feed)
       and ``F302`` (spike feed), dilution-water mixers, coolers
     - Splits the slurry between the two trains and concentrates or
       dilutes each to its setpoint: the initial charge at the target
       concentration, the spike feed at 600 g/L glucose.
   * - Fermentation
     - fermentor ``V406``; air compressor ``K330`` and let-down valve
       ``V330``; seed-yeast tank; vent scrubber ``V409``
     - Runs the kinetic fed-batch model (below); the scrubber captures the
       alcohols stripped by the off-gas (in scenario B about 1.2 % of the
       ethanol and 1.9 % of the isobutanol produced) and recycles them into
       the beer sent to distillation.
   * - Alcohol recovery and purification
     - gating splitter ``S201``; beer column ``D101``, rectifier ``D102``,
       molecular sieve ``MS201``, stripper ``D103``, decanter, drying column
       ``D104``; ethanol day tank and denaturant chain, isobutanol storage
       ``V514``
     - Recovers fuel ethanol (denatured) and isobutanol; see
       :ref:`isobutanol-separation`.
   * - DDGS and corn-oil recovery
     - centrifugation, multi-effect evaporator ``Ev607``, oil separation,
       drum dryer ``D610``, thermal oxidizer ``X611``
     - Turns whole stillage into DDGS and corn oil; broth acetic acid leaves
       with the dryer exhaust and is burned rather than sold as DDGS.
   * - Wastewater treatment (area 700)
     - aqueous-waste mixer ``M501`` and a high-rate treatment train
     - Treats evaporator condensates, the thin-stillage purge and other
       aqueous wastes.
   * - Boiler-turbogenerator (area 800) and facilities (area 900)
     - ``BT801``; cooling tower ``CT901``, chilled-water package ``CWP901``,
       process-water center ``PWC901``, clean-in-place, air distribution,
       fire-water tank
     - Supply steam, power, cooling and process water; the
       process-water center reuses stripper bottoms.
   * - Heat integration
     - heat-exchanger network ``HXN1001``
     - Recovers process heat; the version of the package that synthesizes
       it (the *build tag*) shifts TEA results by up to a few percent.

Separation units are named by their separation-factory IDs (``D101`` …
``D104``, ``MS201``, ``S201``); the flowsheet renames them, so reach them
through ``isobutanol.sep_udct`` after :func:`~nskinetics.biorefineries.isobutanol.load`.


Fermentation section (``V406``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sugar-preparation trains, the fermentor and its aeration loop are built
by one NSKinetics factory,
:obj:`~nskinetics.processes.create_sugar_prep_and_fermentation_system`, which
the biorefinery's ``system.py`` calls. ``V406`` is a
:class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol` unit
wired to the shipped *S. cerevisiae* model (``te_r`` from
``nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo``). On every run:

- **Fed-batch feeding.** A :class:`~nskinetics.FeedSpike` fires whenever
  glucose falls to the *threshold* and tops it back up to the *target* with
  600 g/L feed, at most *spike cap* times. A cap of 0 makes the run a batch
  (scenario B); the threshold is then inert.
- **Two-stage aeration.** The culture is aerobic until cell density reaches
  5.0 g/L or a time cap, then anaerobic; the compressor is sized to twice the
  stage-1 oxygen demand and re-simulated after every fermentation run.
- **Fermentation time τ** is the first time at or after the last spike at
  which glucose falls below 0.5 g/L, capped at τ\ :sub:`max` = 120 h.
- **Closures.** The fermentor compares the spike-feed sugar its kinetics
  added with what the spike train delivers and re-derives the splitter
  setting until they agree (absolute relative residual ≤ 1e-3; exactly 0
  for a batch). The enzyme-burden constraint (:doc:`enzyme_burden`) derates
  growth for strains whose enzyme load exceeds the proteome budget; it is on
  but inert in scenario A and off in scenario B.

.. figure:: /_static/images/isobutanol/overview/sugar_prep_feeding_aeration_concept.png
   :target: ../../_static/images/isobutanol/overview/sugar_prep_feeding_aeration_concept.png
   :class: pub-figure
   :width: 100%
   :alt: Three-panel concept diagram. Panel a: splitter feeding two evaporator, dilution-mixer and cooler trains into the fed-batch fermentor, with a fed-batch strategy specification setting the actuators and a compressor and valve supplying air during stage 1. Panel b: schematic glucose, ethanol, isobutanol and cell trajectories in one batch, with threshold-triggered spikes, a stage-1 aerobic window and the batch end at the first glucose minimum. Panel c: flowchart of one simulation, from specifying the feeding strategy through simulating the kinetics, closing feed and air balances and simulating the full flowsheet until the state stops drifting, with a failure-recovery cascade.

   How the fermentation section is controlled and simulated. **a** — the
   fed-batch strategy specification (pink) root-solves each train's
   evaporator and dilution mixer so the initial charge lands on the target
   and the spike feed on its set concentration; aeration (orange) is sized
   by demand pull during stage 1 only. **b** — the in-batch logic,
   schematic only: spikes restore the target whenever glucose hits the
   threshold, up to the cap, and the batch ends when glucose is exhausted
   after the last spike. **c** — one simulation is a fixed-point sweep
   (kinetics → feed and air balances → full flowsheet) repeated until the
   flowsheet state stops drifting, with a recovery cascade (orange) for
   failed runs. Drawn by ``plots/conceptual_feeding_aeration_sugar_prep.py``.

.. figure:: /_static/images/isobutanol/overview/baseline_A_fermentation_trajectory.png
   :target: ../../_static/images/isobutanol/overview/baseline_A_fermentation_trajectory.png
   :class: pub-figure
   :width: 100%
   :figclass: ibo-fig-measure
   :alt: Concentration versus time for the scenario A baseline fermentation over 80 hours: glucose starts near 221 g/L with small ripples from seven feed spikes between 9 and 11.5 hours, then falls to near zero by about 60 hours; ethanol rises to about 115 g/L; cell density plateaus near 13 g/L; acetate and isobutanol stay near zero. Dashed lines mark the end of aeration near 12 hours and the end of fermentation near 62 hours.

   The scenario-A baseline fermentation inside ``V406`` (threshold 217.125 /
   target 221.25 g/L, cap 16). The ripples on the glucose curve are the
   **7 spikes**, all between 9.0 and 11.5 h, inside the aerobic stage, which
   ends at 12.13 h. Fermentation ends at τ = **61.62 h**, the first time
   glucose drops below 0.5 g/L after the last spike, with ethanol at
   **114.54 g/L** (it would still creep to 114.95 g/L by 120 h), cells
   12.97 g/L and acetate 1.25 g/L; the isobutanol pathway is off in A.
   Concentrations per litre of water. Panel B of the published feeding
   figure ``feeding_strategy_baseline_A`` (see :doc:`feeding_strategies`);
   trajectory from build 6d4776f.


.. _isobutanol-separation:

Alcohol separation: solvent-free heteroazeotropic distillation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Isobutanol and water form a minimum-boiling heteroazeotrope (89.8 °C,
66.6 wt% isobutanol) that splits into two liquid phases on cooling — at
40 °C an organic phase of 87.5 wt% isobutanol and an aqueous phase of
10.7 wt%. The train exploits this instead of an extraction solvent:

1. The **beer column** ``D101`` takes both alcohols overhead (dilute
   isobutanol is far more volatile than water) and sends the stillage to
   DDGS recovery.
2. The **rectifier** ``D102`` makes near-azeotropic ethanol overhead, which
   the **molecular sieve** ``MS201`` polishes (~99.2 wt%) before
   denaturing; isobutanol leaves in the rectifier bottoms.
3. The **stripper** ``D103`` concentrates the isobutanol overhead (its
   nearly isobutanol-free bottoms go to process-water reuse); the
   **decanter** crosses the azeotrope, returning the aqueous phase to the
   rectifier feed.
4. The **drying column** ``D104`` distils the organic phase; on the
   isobutanol-rich side water is the light key, so dry isobutanol
   (> 99.9 wt%) is the bottoms product and the azeotropic overhead returns
   to the decanter.

Every column design specification is re-derived from the current feed on
each simulation pass, so the same train handles roughly 0–200 g/L of either
alcohol — from an ethanol-only broth (stripper, decanter and drying column
idle and uncosted) to an isobutanol-rich one.

**Ethanol-primary branch.** The gating splitter ``S201`` can instead route
broth to a conventional ethanol train adapted from ``biorefineries.ethanol``
(beer column, rectifier, molecular sieve). At the baseline split of 1.0 it
idles at exactly zero flow and is not designed or costed. It cannot recover
isobutanol — what reaches it ends up in wastewater treatment — so it is the
right comparison only for ethanol-only strains. Choose which trains are
built with ``isobutanol.load(separation_processes=...)`` (``('IBO_EtOH',)``,
``('ethanol',)`` or both, the default); intermediate splits
(``0 < split < 1``) have not been exercised in the full flowsheet.


TEA basis
---------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Item
     - Setting
   * - TEA model
     - The corn package's ``ConventionalEthanolTEA`` with its default
       financial parameters (MACRS-7 depreciation, 40/60 construction
       schedule, working capital 5 % of fixed capital), except
       ``lang_factor=None``.
   * - Hurdle rate
     - 15 % IRR (the rate at which MPSPs and PI are computed).
   * - Operating time
     - 330 days per year (7,920 h/y).
   * - Income tax
     - 21 % (scenario workbook baseline).
   * - Price year
     - **2023 dollars.** Every stream and utility price is converted once
       from its source year with the BLS chemicals producer price index
       (WPU06).
   * - Cost index
     - CEPCI 797.9 (2023 annual average).
   * - Capital costing
     - **Bare-module, unit by unit, with no Lang factor.** Each unit's
       installed cost is its purchase cost × (F\ :sub:`BM` + F\ :sub:`D`
       F\ :sub:`P` F\ :sub:`M` − 1); fixed capital is the sum, with no
       separate indirect, field, contingency or site costs. A unit without
       bare-module data counts at purchase cost.
   * - Utilities
     - Steam, cooling water and chilled water are regenerated on site at
       zero price; the only utility cash flows are purchased natural gas
       and sold excess electricity.
   * - Feedstock
     - Corn, 56,972 kg/h (wet), 0.17414 $/wet-kg.
   * - Products
     - Ethanol 0.835 $/kg (midpoint of its 2021–2025 range); isobutanol
       1.500 $/kg (a 1.49 $/kg 2025 listing, indexed to 2023$); DDGS
       0.14304 $/kg; corn oil.

**Selling prices and returns.**
:func:`~nskinetics.biorefineries.isobutanol.system.solve_TEA` (a wrapper of
:func:`~nskinetics.biorefineries.isobutanol.system.solve_TEA_at_IRR`) works
on the current flowsheet state without
re-simulating and returns:

- **MPSP** for each product: the minimum selling price at a 15 % IRR while
  the other product sells at its market price, **purity-adjusted** — the
  solved stream price divided by the main chemical's mass fraction, so it is
  in $ per kg of ethanol or isobutanol. (Scenario A's denatured ethanol is
  95.046 % ethanol; its unadjusted MPSP would read 3.64 $/GGE instead of
  3.83.) A product with no flow gets NaN. MESP is the ethanol MPSP in $/GGE
  (× 4.45775; see :ref:`isobutanol-conventions`).
- **IRR** with every product at its market price. When no genuine root
  exists on IRR > −100 % and the plant loses money outright, IRR is
  reported as **−inf**; NaN is reserved for "not solved".

The optimization and sensitivity analyses rank designs by the
**profitability index**, PI = NPV at 15 % ÷ TCI, because IRR is −inf over
most of the design space and MPSP is NaN wherever a product is absent (see
:doc:`kinetic_optimization`).


Unit groups
~~~~~~~~~~~

For cost and utility breakdowns the flowsheet is partitioned into 19 unit
groups, available as ``isobutanol.unit_groups_dict`` after
:func:`~nskinetics.biorefineries.isobutanol.load`. The bookkeeping rule is
that every priced feed and every utility cash flow is counted in exactly one
group, so the groups' operating cost × operating hours equals the TEA's
annual operating cost (AOC), and their installed costs sum to the installed
equipment cost. Units shared by both products are not allocated; they get
their own group.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Group
     - Contents
   * - feedstock acquisition
     - grain handling and corn storage (the corn purchase is booked here,
       once)
   * - feedstock saccharification
     - milling through saccharification and the slurry cooler
   * - sugar solution preparation
     - splitter ``S301`` and both feed-conditioning trains
   * - fermentation
     - ``V406``, aeration ``K330``/``V330``, yeast tank and pump
   * - alcohol recovery
     - the shared units: broth pump, vent-scrubber recycle, gating splitter
       ``S201``, beer column, rectifier (which performs the
       ethanol/isobutanol split), stillage merge
   * - ethanol purification
     - ethanol superheater, molecular sieve and condenser; the ethanol-primary
       train when built; the ethanol merge and day-tank pump
   * - isobutanol purification
     - stripper, decanter loop and drying column
   * - storage and handling
     - storage tanks, DDGS handling, denaturant
   * - DDGS recovery
     - stillage processing through the dryer, including the thermal
       oxidizer's fuel
   * - wastewater treatment
     - waste mixer and high-rate treatment train
   * - heat exchanger network
     - ``HXN1001`` (its savings appear as credits)
   * - boiler; turbogenerator
     - ``BT801`` split into two groups: the boiler (with its solids mixer)
       and the turbogenerator cost item
   * - cooling utility facilities
     - cooling tower and chilled-water package
   * - other facilities
     - process-water center, clean-in-place, air distribution, fire water
   * - four cost-only groups
     - natural gas for steam generation, natural gas for product drying,
       fixed operating cost, excess electricity (a credit)

Each group reports installed equipment cost, cooling duty, heating duty,
electricity consumption and operating cost (USD/h):

.. code-block:: python

   import biosteam as bst
   from nskinetics.biorefineries import isobutanol
   # ... after isobutanol.load() and scenarios.load_scenario(...)

   groups = list(isobutanol.unit_groups_dict.values())
   print(bst.UnitGroup.df_from_groups(groups))


Results
-------

Baseline scenarios A and B
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Scenario baselines (build 6d4776f; per litre of water)
   :header-rows: 1
   :stub-columns: 1
   :widths: 30 35 35

   * -
     - A — wild-type, fed-batch
     - B — engineered co-producer, batch
   * - Feeding
     - threshold 217.125 / target 221.25 g/L, cap 16 → 7 spikes
     - batch, target 140 g/L
   * - Enzyme burden
     - on, inert (Φ\ :sub:`M` 0.0701 g/gDCW, no derating)
     - off (Φ\ :sub:`M` 0.286 would exceed the 0.245 cap)
   * - Fermentation time τ
     - 61.62 h
     - 66.67 h
   * - Ethanol titer
     - 114.54 g/L
     - 26.30 g/L
   * - Isobutanol titer
     - —
     - 29.06 g/L
   * - Ethanol yield / productivity
     - 0.4602 g/g / 1.859 g/L/h
     - —
   * - Ethanol MPSP
     - 0.86021 $/kg (3.835 $/GGE)
     - 0.59530 $/kg (≈ 2.654 $/GGE)
   * - Isobutanol MPSP
     - NaN (not produced)
     - 1.2800 $/kg
   * - IRR at market prices
     - 12.79 %
     - 21.76 %
   * - TCI / AOC
     - 139.09 MM$ / 115.67 MM$/y
     - —

A dash means the source reports no value. **Scenario B is provisional:** its
workbook still uses ``k_17 = 44`` for Adh6 instead of the corrected 0.1077,
so r17 runs about 400× too fast (see `Caveats and limitations`_). B is a
reference build, not a feasible strain: its pathway constants also exceed
the proteome cap. Do not read its 21.76 % IRR as evidence that
co-production pays; the feasible co-production result is
``opt_PI_TRY_informed`` (:doc:`kinetic_optimization`).

**Reading the baseline.** Scenario A's MESP, 3.83 $/GGE, is about
0.11 $/GGE (2.9 %) above the model's ethanol price (3.72 $/GGE =
0.835 $/kg, the midpoint of the 2021–2025 range), so at that price the
plant earns 12.79 % rather than the 15 % hurdle. Its titer and yield sit
inside the US high-gravity ranges (100–142 g/L, 0.42–0.48 g/g); its
productivity, 1.86 g/L/h, is below the 2.0–4.4 g/L/h range — but that
conservatism hardly matters: raising productivity to 2.0 g/L/h would lower
MESP by only 0.004 $/GGE.

Where the money goes
~~~~~~~~~~~~~~~~~~~~

In scenario A the fermentation section is **13.74 of 132.46 MM$ of installed
equipment cost (10.4 %)**, while corn is **~86 % of material cost**. A faster
fermentation shrinks only the fermentor, so productivity moves only a small
slice of cost; yield, which sets how much product each kilogram of corn
becomes, moves MESP far more (the per-unit leverage is quantified on
:doc:`uncertainty`). Build 6d4776f.

TEA breakdowns at the campaign optima
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The table below re-simulates scenario A and the optimum of each of the eight
strain-design campaigns (:doc:`kinetic_optimization`) and reads the TEA
(every trial reproduced within 0.17 %). It is on **build 2b5b27d**, the build
the campaigns ran under, so its scenario-A row (0.8623 $/kg, IRR 12.60 %)
differs slightly from the 6d4776f pin above — compare rows within this
table, not across tables. MM$/y unless noted; corn-oil revenue is 4.0 MM$/y
in every row.

.. list-table:: Cost and revenue at the campaign optima (build 2b5b27d)
   :header-rows: 1
   :stub-columns: 1
   :widths: auto
   :class: ibo-num

   * - Design
     - IRR
     - EtOH MPSP $/kg ($/GGE)
     - IBO MPSP $/kg
     - TCI MM$
     - AOC
     - Sales
     - EtOH rev.
     - IBO rev.
     - DDGS rev.
   * - Baseline A
     - 12.60 %
     - 0.8623 (3.84)
     - —
     - 139.3
     - 115.9
     - 138.3
     - 117.3
     - 0
     - 17.0
   * - Uninformed PI #103
     - 15.83 %
     - 0.8258 (3.68)
     - —
     - 137.0
     - 114.8
     - 142.0
     - 122.1
     - 0
     - 15.9
   * - **TRY-informed PI #1912**
     - **27.26 %**
     - 0.4749 (2.12)
     - 1.127
     - **137.7**
     - 116.2
     - **164.6**
     - 52.1
     - 90.4
     - 18.0
   * - IBO yield #841
     - 1.68 %
     - (5,769)\ :sup:`*`
     - 1.854
     - 315.6
     - 173.9
     - 192.8
     - 0
     - 171.6
     - 17.2
   * - IBO titer #84
     - −inf
     - (210)\ :sup:`*`
     - 2.474
     - 148.2
     - 123.0
     - 106.4
     - 0.2
     - 68.6
     - 33.6
   * - IBO productivity #1472
     - −inf
     - (499)\ :sup:`*`
     - 2.282
     - 139.3
     - 122.3
     - 113.3
     - 0.1
     - 67.6
     - 41.7
   * - EtOH yield #1923
     - 12.06 %
     - 0.8703 (3.88)
     - —
     - 156.0
     - 118.8
     - 142.9
     - 123.5
     - 0
     - 15.4
   * - EtOH titer #1820
     - 12.68 %
     - 0.8614 (3.84)
     - —
     - 140.6
     - 116.3
     - 139.0
     - 118.7
     - 0
     - 16.3
   * - EtOH productivity #365
     - 8.07 %
     - 0.9130 (4.07)
     - —
     - 141.1
     - 116.2
     - 132.3
     - 112.0
     - 0
     - 16.3

:sup:`*` Near-zero ethanol co-product: a purity-adjusted MPSP artifact.

.. admonition:: Flagship #1912 vs baseline A — the gain is all on the revenue side
   :class: tip

   - **Capital and cost barely move:** TCI −1.2 % (137.7 vs 139.3 MM$),
     AOC +0.2 %. **Sales rise 26.2 MM$/y (+19 %)**, and isobutanol supplies
     55 % of them.
   - **Why:** revenue per kg of sugar rises from 0.384 to 0.469 $/kg
     (+22 %) — baseline 0.460 g/g × 0.835 $/kg vs 0.199 g/g isobutanol ×
     1.500 + 0.204 g/g ethanol × 0.835. Total alcohol yield *falls* to
     0.403 g/g, but that is stoichiometry, not lost carbon (88 % vs 90 % of
     the theoretical maximum): isobutanol sells at 1.80× ethanol per kg
     (1.45× per kg of glucose at theoretical yield). This revenue-per-sugar
     change accounts for essentially all of the 131 MM$ gain in NPV.
   - **Equipment is nearly unchanged** (installed, MM$): ethanol
     purification 6.52 → 3.98, isobutanol purification 0 → 1.76,
     fermentation 13.7 → 12.6, boiler 46.0 → 46.3, DDGS 14.8 → 16.8.
   - **Utilities:** the heat-exchanger network recovers 26 % of gross
     heating demand (vs 21 %), so net heating rises only 4.5 %; cooling duty
     and electricity each fall 14 %.
   - **Why both alcohols:** in this model isobutanol inhibits glycolysis
     about 3.75× more strongly per gram than ethanol, and its cell-death
     term switches on above 25 g/L (against 100 g/L for ethanol) and then
     grows with concentration. Isobutanol-only optima therefore stall near
     46 g/L, while the flagship splits its product and accumulates about
     80 g/L of total alcohol.

**Single-metric optima fail through different channels:**

- **Isobutanol yield (#841)** reaches 93 % of theoretical yield and the
  highest sales of any design (192.8 MM$/y), yet returns IRR 1.7 %. Its
  slow (102 h), dilute broth (20.9 g/L isobutanol) carries 4.5× the
  baseline's water per kg of glucose, so capital is 2.27× (boiler 46.0 →
  108.3, fermentation 13.7 → 56.1, DDGS 14.8 → 47.3, wastewater 12.0 →
  30.5 MM$ installed) and natural gas for steam costs 52.2 vs
  11.7 MM$/y.
- **Isobutanol titer and productivity (#84, #1472)** lose money outright:
  sales (106.4 / 113.3 MM$/y) fall below operating cost (123.0 /
  122.3 MM$/y), because yield is only 0.15 g/g (37 % of theoretical) —
  consistent with the culture stalling under isobutanol toxicity with much
  of the loaded glucose unconverted (residual sugar is not a saved metric).
- **Ethanol titer, yield and productivity (#1820, #1923, #365)** change MESP
  by −0.1 %, +0.9 % and +5.9 % despite reaching 168 g/L, 0.485 g/g and
  6.9 g/L/h. Concentrating sugar for their high-gravity feeds raises
  sugar-preparation heating to 58–69 GJ/h (vs 16).

**Read co-production economics as break-evens.** MESP is not a stable
headline for co-producing designs: two equally profitable designs from
different campaigns read 2.12 and 1.60 $/GGE (#1912, and #1602, the best
design of a second, independent profitability campaign run with random seed
33960; see :doc:`kinetic_optimization`; both build 2b5b27d), because each
MPSP depends on the other product's revenue. The
robust statement for #1912 is that it clears the 15 % hurdle while
isobutanol sells above **1.13 $/kg** (25 % below the modeled price, ethanol
at 0.835 $/kg), or while ethanol sells above **0.475 $/kg** (below its
2021–2025 low of 0.518 $/kg, isobutanol at 1.500 $/kg).

.. note::

   The live ``opt_PI_TRY_informed`` scenario reproduces #1912's fermentation
   exactly but, on build 95e9362, reads ethanol MPSP 0.46066 $/kg,
   isobutanol MPSP 1.1123 $/kg and PI 0.8559 — vs 0.4749, 1.127 and 0.8083
   in the campaign record above. Use one build per comparison.


Caveats and limitations
-----------------------

.. admonition:: Build provenance
   :class: warning

   Three builds of the heat-exchanger-network package appear on this page:
   **6d4776f** (scenario A/B pins, the baseline trajectory, the cost-structure
   numbers), **2b5b27d** (the campaign-optima breakdowns, break-evens and
   flagship comparison) and **95e9362** (the live ``opt_*`` pins). They move
   TEA results by up to a few percent — more at co-production designs, where
   the IRR of #1602 (the seed-33960 profitability optimum) moved from
   27.30 % to 28.02 % between builds — while fermentation metrics reproduce
   across them to < 1e-5. A difference between numbers from two builds
   therefore mixes the build shift into the result.

- **Scenario B's Adh6 constant.** The B workbook still sets ``k_17 = 44``
  while the corrected model value is 0.1077, so Adh6 (r17) runs about 400×
  too fast in B. Correcting it is a kinetic change that re-pins B; until
  then, quote B's numbers as provisional. (The ``opt_*`` workbooks carry
  each trial's own ``k_17`` and are not affected.)
- **The baselines are not the feeding optima.** Scenario A's 217/221 g/L
  policy was the optimum of an older, coarser sweep on a superseded model;
  on the current 40 × 40 grid the optimum is 226/230 g/L with 4 spikes, and
  A sits 0.54 % above it (build 2b5b27d; :doc:`feeding_strategies`). B's
  batch policy is likewise the IRR optimum of an older 25 × 25 sweep.
- **Capital basis.** Fixed capital is the bare-module installed-cost sum
  with no separate indirect, field, contingency or site costs; compare
  designs within this basis rather than against full-estimate TEAs.
- **Prices.** The isobutanol price rests on a single 2025 product listing,
  with no limit on how much isobutanol the market absorbs; the co-production
  results are therefore stated with break-even prices. The modeled DDGS
  price (about 130 $ per short ton) lies below its 2021–2025 low.
- **Basis.** Titers and productivities are per litre of water, which
  overstates them relative to the per-litre-of-beer figures industry
  usually quotes; the broth/water conversion (≈ 0.9014 × (1 − vent
  fraction)) has been verified for scenario B only.
- **Near-zero co-products** produce absurd purity-adjusted MPSPs (the
  starred entries above); they are artifacts, not results.
- **Separation configuration.** The ethanol-primary train cannot recover
  isobutanol, and intermediate ``S201`` splits have not been exercised in the
  full flowsheet.
- **Figure.** The 3 × 3 TEA-breakdown figure is not reproduced on this
  page; regenerate it as below.


Reproduce it
------------

Scripts that simulate are opt-in: run each explicitly, in a fresh process,
and never two at once (they share an on-disk numba cache). No analysis
results ship with NSKinetics; every input or output lives under the results
root, which you can point at an existing results folder with the
``NSKINETICS_IBO_RESULTS_DIR`` environment variable. Script paths are
relative to the package folder; :doc:`running_analyses` lists every script.

.. list-table::
   :header-rows: 1
   :widths: 30 44 26

   * - Script or command
     - What it does
     - Writes / runtime
   * - ``python -m nskinetics.biorefineries.isobutanol.smoke 1`` (and ``2``)
     - Scenarios A and B: load, then three simulate + ``solve_TEA`` cycles,
       asserting the pins (within 1 %), stability (within 0.5 % of the first
       simulation) and the spike-feed residual. In Python,
       :func:`~nskinetics.biorefineries.isobutanol.smoke.load_simulate_baseline`
       returns the three result dicts.
     - nothing; ~2–4 min per case
   * - ``... smoke 3`` to ``... smoke 8``
     - The same checks in other separation configurations: 3/4 build only
       the ethanol + isobutanol train (reproducing 1/2), 5/6 only the
       ethanol-primary train, 7/8 both trains with all broth gated to the
       ethanol-primary one.
     - nothing; ~2–4 min per case
   * - ``... smoke 9`` to ``... smoke 16``
     - The eight ``opt_*`` scenarios, additionally re-deriving each
       campaign's objective.
     - nothing; ~2–4 min per case
   * - ``analyses/opt_scenarios/verify_opt_scenarios.py <scenario>``
     - Checks that an ``opt_*`` scenario's baseline equals a fresh
       reproduction of its campaign trial. Needs the campaign trajectory
       CSVs.
     - nothing; minutes
   * - ``analyses/tea_breakdowns/collect_tea_breakdowns_split12d.py``
     - Re-simulates scenario A and the eight campaign optima, checks the
       unit-group closures, and records the TEA breakdowns. Needs the
       campaign trajectory CSVs. To reproduce the table above, switch
       hensmith to the campaigns' build, commit ``2b5b27d``
       (:ref:`isobutanol-switch-build`); on another build the script still
       runs, but the TEA numbers shift.
     - ``tea_breakdowns_split12d_<stamp>.json``; ~4 min
   * - ``plots/plot_tea_breakdowns_split12d.py``
     - Draws the 3 × 3 breakdown figure from the newest JSON (never
       simulates).
     - PNG/PDF + one CSV per scenario; seconds
   * - ``analyses/opt_scenarios/build_opt_split12d_workbooks.py <scenario>``
     - Regenerates an ``opt_*`` workbook and its registry pin from a recorded
       trial. Needs the campaign trajectory CSVs.
     - a parameter workbook; minutes
   * - ``plots/plot_block_flow_diagram.py``,
       ``plots/conceptual_feeding_aeration_sugar_prep.py``
     - The two schematics on this page (matplotlib only; never simulate).
     - ``publication/Schematics`` (PNG/PDF/SVG); seconds
   * - ``plots/plot_feeding_strategy_baseline_A.py``
     - The feeding figure whose panel B is shown above. Panel A needs the
       ``analyses/feeding_strategies/evaluate_feeding_strategies.py`` sweep
       CSVs; panel B simulates scenario A once if its trajectory cache is
       missing.
     - ``publication/Feed-strat``; seconds with cached data

Next: :doc:`feeding_strategies` asks which feeding policy minimizes MESP
for the baseline strain.
