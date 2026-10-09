Objective landscape
===================

The :doc:`kinetic_optimization` campaigns raise two puzzles. TRY campaigns
pass through profitable designs but do not return them, and an uninformed
profitability search (one not preloaded with TRY designs) succeeds or fails
depending on its random seed. This
page looks at the *shape* of the objectives over the 12-d strain-design and
feeding space to explain both, and compares the profitability index (PI)
with the titer, rate and yield (TRY) objectives. The answer is about
**topology, not roughness**. Locally, PI is not intrinsically rougher than the
TRY objectives. What sets it apart is that its profitable region is tiny and
split into **separated islands**. One-dimensional slices through the best
design show what each island is made of.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: 0 of 4,091
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Designs with PI > 0 in a uniform sample of the campaigns' own space.
      The profitable region is tiny.

      +++
      Sobol' campaign-box sample, build 2b5b27d

   .. grid-item-card:: −0.07 vs −0.39
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      PI's fitness–distance correlation in the optimizer-sampled pool vs a
      uniform Sobol' sample. The "needle" look comes from the sampling.

      +++
      Optimizer pool vs Sobol' sample, build 2b5b27d

   .. grid-item-card:: +0.50 → −0.68
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Range of PI along the isobutanol-yield plateau (0.37–0.38 g/g). High
      yield does not imply profit.

      +++
      Optimizer pool, build 2b5b27d

   .. grid-item-card:: 0.32–0.33×
      :text-align: center
      :class-card: ibo-stat
      :shadow: sm

      Pdc capacity, relative to wild type, at the PI peak of the Pdc slice
      through the best design (PI 0.8086). Wild-type Pdc gives 0.653.

      +++
      Pdc slice through #1912, build 2b5b27d

.. admonition:: What this page does *not* claim
   :class: important

   The figures below were first drawn to contrast a "needle-like" PI with
   "smooth" TRY objectives, and the plotting scripts' docstrings and original
   file names (``…_needle_plateau``) still say so. **That contrast does not
   survive unbiased sampling.** It holds only in
   the pool of trials the optimizers chose to sample. This page therefore reads
   every landscape figure as a picture of the *optimizer-sampled pool's
   topology*. Do not cite them as evidence that PI is intrinsically rougher
   than the TRY objectives. Likewise, the sharp trench in the Pdc
   slice is a separation-design regime, **not** evidence of objective
   roughness.

.. note:: **Conventions on this page.** Trial numbers (#841, #842, #1882,
   #1912) are the 0-based CSV ``trial_number``, as in the figures. The
   campaign figures on :doc:`kinetic_optimization` count from 1, so the same
   design appears there one number higher (#842 is "trial 843"). Every
   dataset here ran under the hensmith build **2b5b27d**, the campaigns'
   build. PI = NPV at the 15 % hurdle ÷ TCI, so PI > 0 exactly when
   IRR > 15 % (see :doc:`kinetic_optimization`).

At a glance
-----------

.. list-table::
   :widths: 26 74
   :stub-columns: 1

   * - Question
     - Why do TRY campaigns traverse profitable designs without returning
       them, and why is profitability search seed-dependent?
   * - Data
     - The 13,622 unique ``COMPLETE`` trials of the seven seed-350
       ``metabolic_split_12d`` GP campaigns (``_rs350``). As an unbiased
       control, two uniform-feasible Sobol' samples of the same space.
       1-D sweeps through the TRY-informed optimum #1912.
   * - Distance
     - Euclidean in the campaigns' own 12-d unit cube (log axes
       log-normalized, linear axes min–max), the distance the GP itself sees.
       The cube diagonal is √12 = 3.464.
   * - Statistics
     - Fitness–distance correlation, exact rank variogram, nearest-neighbour
       Spearman, distance profiles, a t-SNE map, and 1-D slices
   * - Model build
     - 2b5b27d throughout

Method
------

Landscape statistics on the optimizer pool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The seven ``_rs350`` campaigns (PI (log-tail) and the six TRY objectives,
2,000 trials each) are pooled and de-duplicated by the shared data layer
``_objective_landscape_data.py``. That gives 13,965 ``COMPLETE`` rows with a
finite PI, of which 13,622 are unique. Each trial's decision vector is mapped
into the unit cube. For PI and each TRY objective *f*, the analysis
computes:

- **Fitness–distance correlation** (FDC; Jones & Forrest 1995): the Spearman
  correlation of *f* with the distance from *f*'s own best trial *x*\ :sup:`*`.
  Strongly negative means the objective falls off steadily away from its
  optimum, a searchable hill. Near 0 means the optimum gives no long-range
  signal.

  .. math::

     \mathrm{FDC} = \rho_S\!\left(f(x_i),\; \lVert x_i - x^{*} \rVert\right)

- **Rank variogram**: the normalized semivariance of *percentile ranks* r
  against pair distance h, computed **exactly** over all ~92.8 million trial
  pairs. γ = 1 means a pair at that distance is no more alike than two random
  trials. Ranks make the seven objectives commensurate and neutralize PI's
  heavy loss tail.

  .. math::

     \gamma(h) = \frac{\tfrac12\,\mathrm{E}\!\left[(r_i - r_j)^2 \,\middle|\, d_{ij} \approx h\right]}{\mathrm{Var}(r)}

- **Nearest-neighbour Spearman**: the rank correlation of each trial's value
  with that of its nearest 12-d neighbour.
- **Distance profiles** around the most profitable pooled trial (#842), and
  a **t-SNE map** of the whole pool.

The unbiased control
~~~~~~~~~~~~~~~~~~~~

An optimizer pool is not a uniform sample. Each GP campaign exploits around
its own optimum, so short-distance pairs crowd onto whatever plateau that
campaign climbed. The same statistics were therefore recomputed on the two
**uniform-feasible Sobol' samples** from the :doc:`sensitivity` analysis
(no optimizer involved). One is drawn under the campaign measure and one
under the screening (engineering-prior) measure. Each has 4,091 completed
designs; the screening-measure statistics use the 4,086 on which all seven
objectives are finite (five have no finite ethanol productivity).

1-D slices through the best design
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``analyses/objective_landscape/evaluate_axis_try_informed_optimum.py`` holds every decision at the
TRY-informed optimum #1912 (PI 0.808, IRR 27.3 %, one glucose spike) and
sweeps one variable. It runs in two modes:

- a **coarse screen** of the 11 decision variables other than ``k_13`` (357
  rows: 354 OK, 3 infeasible);
- a **full-resolution** sweep of one variable. For Pdc (``k_3``) this is 219
  points: 218 log-spaced over the campaign band 0.00581–23.24 g/L/h
  (10\ :sup:`−3`\ × to 4× wild type) plus the anchor value.

``evaluate_k13_try_informed_optimum.py`` (same folder) adds a full ALS (``k_13``) sweep of
161 points. Every point is restored to the anchor's converged flowsheet
snapshot before it is simulated, so the sweep order can neither create nor
hide a spike. The enzyme burden is on, as in the campaigns. Each sweep first
re-simulates #1912 and checks it against the recorded row (the ``k_3`` anchor
reproduced PI 0.80855 against the recorded 0.80832).

Results
-------

In the optimizer pool, PI looks needle-like…
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pooled PI runs from −5.56 to 0.504 (median −1.47). The most profitable
trial is #842 of the isobutanol-yield campaign (PI 0.504, IRR 22.9 %): a
Pdc-null, pure-isobutanol design at 0.354 g/g. The highest-yield trial, #841
(0.381 g/g), is about 0.40 away and has PI −0.683. Of the 166 trials with
PI > 0, 124 come from the isobutanol-yield campaign. Forty are ethanol-only,
all with PI ≤ 0.049, and every trial with PI > 0.05 makes isobutanol (yield
≥ 0.114 g/g).

.. tab-set::

   .. tab-item:: Optimizer pool (13,622 trials)

      .. list-table:: Fitness–distance correlation in the optimizer pool
         :header-rows: 1
         :widths: 60 40

         * - Objective
           - FDC
         * - **Profitability index (PI)**
           - **−0.069**
         * - Isobutanol yield
           - −0.568
         * - Isobutanol titer
           - −0.693
         * - Isobutanol productivity
           - −0.711
         * - Ethanol yield
           - −0.806
         * - Ethanol titer
           - −0.798
         * - Ethanol productivity
           - −0.695

      PI is also the roughest objective at **every lag** of the rank
      variogram up to h ≈ 1.5:

      .. list-table:: Normalized rank semivariance γ(h), optimizer pool
         :header-rows: 1
         :widths: 30 35 35

         * - Pair distance h
           - PI
           - Six TRY objectives
         * - 0.07
           - 0.097
           - 0.003–0.066
         * - 0.36
           - 0.150
           - 0.035–0.059
         * - 0.72
           - 0.349
           - 0.162–0.197
         * - 1.11
           - 0.607
           - 0.350–0.375

      Nearest-neighbour Spearman: PI 0.753 vs 0.857–0.889 for the TRY
      objectives.

   .. tab-item:: Unbiased Sobol' samples (the control)

      .. list-table:: Fitness–distance correlation on uniform-feasible samples
         :header-rows: 1
         :widths: 40 30 30

         * - Sampling measure
           - PI
           - Six TRY objectives
         * - Campaign measure (4,091 designs)
           - **−0.39**
           - −0.35 to −0.48
         * - Screening measure (4,086 designs)
           - **−0.21**
           - −0.02 to −0.35

      Under the campaign measure, PI's rank semivariance at h ≈ 0.73 is
      0.467, against 0.414–0.481 for the TRY objectives. On samples that no
      optimizer chose, **PI is as smooth as the TRY objectives.**

…but not under unbiased sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pool contrast is a sampling effect. Six of the seven pooled campaigns
exploited around a TRY objective's optimum, so short-distance pairs are
over-represented on TRY plateaus. The isobutanol-yield campaign also
contributed most of the trials near #841 and #842. Unbiased samples remove
that bias, and with it PI's apparent roughness. The figures below should be
read with this in mind.

.. figure:: /_static/images/isobutanol/landscape/fitness_distance_optimizer_pool.png
   :target: ../../_static/images/isobutanol/landscape/fitness_distance_optimizer_pool.png
   :class: pub-figure
   :alt: Three panels. (a) Profitability index of 13,622 trials against
         unit-cube distance from trial 842, with the 166 profitable trials in
         red at the top and FDC rho = -0.07. (b) Isobutanol yield against
         distance from trial 841, with FDC rho = -0.57 and the 40 ethanol-only
         profitable designs at zero yield labelled. (c) Rank semivariance
         against pair distance for PI and six TRY objectives; the PI curve
         lies above the others.
   :width: 100%

   **Fitness–distance structure of the optimizer-sampled pool** (13,622
   unique trials of the seven seed-350 campaigns). **(a)** PI of every trial
   against its unit-cube distance from the most profitable pooled trial (#842,
   star). The axis is linear above 0 and log-tail compressed below 0, with
   ticks in raw PI; no trial is dropped. **(b)** Isobutanol yield against the
   distance from the highest-yield trial (#841, diamond). Trials with PI > 0
   are red; the 40 ethanol-only ones sit at zero yield. **(c)** Normalized
   rank semivariance against pair distance (h ≤ 1.6) for PI (red) and the six
   TRY objectives, computed over all pairs. The FDC values above (a) and (b)
   are each objective's Spearman correlation with distance from its own best.
   **How to read it:** in *this pool* PI has a near-zero FDC and the highest
   semivariance at every distance. That is a property of where the
   optimizers sampled, not of PI. On uniform Sobol' samples, PI's FDC
   (−0.39) matches the TRY objectives' (−0.35 to −0.48). What the figure shows
   robustly is that the profitable (red) trials sit in separate groups, not on
   one hill.

.. figure:: /_static/images/isobutanol/landscape/distance_profiles_around_842.png
   :target: ../../_static/images/isobutanol/landscape/distance_profiles_around_842.png
   :class: pub-figure
   :alt: Two panels sharing an x axis of distance from trial 842. (a) PI of
         each trial with binned median and interquartile band; the median
         drops below zero within about 0.2 of trial 842. (b) Isobutanol yield
         of the same trials; the median stays near 0.37 out to about 0.8 and
         then falls. The neighbour trial 1882 is ringed in both panels.
   :width: 100%

   **Distance profiles around the most profitable pooled trial (#842).**
   Same anchor and same trials in both panels. **(a)** PI and **(b)**
   isobutanol yield against unit-cube distance from #842, with the binned
   median (0.1-wide bins, n ≥ 10) and interquartile range. Solid line:
   break-even, PI = 0. Dashed: the scenario-A strain (PI −0.14). Ringed: the
   neighbour #1882 at d = 0.135, which has a *higher* yield than #842 (0.379
   vs 0.354 g/g) but PI −0.43 instead of +0.50. Trials with PI < −3.5
   (n = 295) are drawn below the dotted floor; 11,010 trials beyond d = 1.5
   (44 with PI > 0) are not shown. **How to read it:** around this anchor,
   PI's median collapses by d ≈ 0.2 while yield stays on its plateau to
   d ≈ 0.8. Along the isobutanol-yield plateau (0.37–0.38 g/g), PI runs from
   +0.50 to −0.68. This explains why a yield campaign can step off a
   profitable design without noticing. It does **not** show that PI is
   intrinsically rougher than yield (see the Sobol' control above).

.. figure:: /_static/images/isobutanol/landscape/tsne_embedding_optimizer_pool.png
   :target: ../../_static/images/isobutanol/landscape/tsne_embedding_optimizer_pool.png
   :class: pub-figure
   :alt: Two t-SNE maps of the same 13,622 trials with identical coordinates.
         (a) Coloured by profitability index, with high values in small
         patches at the lower left. (b) Coloured by isobutanol yield, with a
         high-yield arm at the left; profitable trials are ringed in red and
         cluster on that arm, plus a group on the zero-yield region at the
         right.
   :width: 100%

   **t-SNE map of the optimizer pool,** coloured by PI **(a)** and by
   isobutanol yield **(b)**, with identical coordinates (perplexity 40;
   trustworthiness 0.99 at k = 10). t-SNE preserves *neighbourhoods*, not
   distances: map distances and cluster sizes are not metric, and the axes
   have no units. The ρ above each panel is the 12-d nearest-neighbour rank
   correlation (PI 0.75, isobutanol yield 0.89). Two caveats temper that
   contrast. The yield ρ is lifted by the split between the two-thirds of
   trials that make no isobutanol (a flat zero floor) and the co-production
   third. Within co-production trials alone, the 12-d ρ is 0.59 for PI and
   0.67 for yield. Also, 40 of the 166 PI > 0 trials (rings) are ethanol-only
   designs with PI ≤ 0.05, which sit on the zero-yield region. PI colours are
   clipped at −1 (NPV = −TCI); 60 % of trials lie below. **How to read it:**
   as topology. Profitable designs occupy small, separate patches, most of
   them on the high-yield arm.

What is robust: islands
~~~~~~~~~~~~~~~~~~~~~~~

The profitable region is tiny and fragmented. With distances measured from
the TRY-informed optimum #1912:

.. list-table:: The profitable islands of the 12-d space
   :header-rows: 1
   :widths: 30 20 50

   * - Island
     - PI
     - Where it is
   * - Uniform background
     - PI > 0 at **0 of 4,091**
     - Uniform campaign-measure sample (best PI −0.077)
   * - Ethanol-only shelf
     - ≤ 0.05
     - 2.0 unit-cube units from #1912
   * - Pdc-null, pure-isobutanol island
     - ≤ 0.50
     - 1.15 units from #1912 (contains #842)
   * - Co-production island
     - ≈ 0.81
     - Contains #1912. It is a ridge rather than a point: #1602 of the second
       seed lies 1.08 units away at equal PI.

No seed-350 scout trial came within 0.5 units of #1912; the nearest of their
11,968 trials is 0.505 away. Isobutanol pays only when ALS, the downstream Ehrlich
enzymes and Adh6 rise together, so the co-production island is narrow and
high-order inside a broad ethanol basin. An optimizer whose 50-point start-up
lands only in the ethanol basin can stay there. This is the seed dependence
seen in :doc:`kinetic_optimization`.

Yield and profit co-locate, except along Pdc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Locally, the picture differs from the global pool, where high yield does
not imply profit (PI runs from +0.50 to −0.68 along the yield plateau). The
1-D slices below run through #1912, where ALS, Adh6 and both tolerances
already sit at their bounds. Along them the isobutanol-yield maximum and
the PI maximum coincide, to within one grid step, along **10 of the 11
screened directions plus** ``k_13``; ``inhib_acetate`` is the one-step case
(0.976 vs 1.023). Two caveats apply. First, 6 of the
12 coordinates sit at a band bound, where coincidence is partly trivial. The
interior coincidences are threshold 189.2 g/L, glycolysis 0.98×, ``k_6`` 7.44
and ``ehrlich_downstream`` 0.575, on a 31-point grid. Second, ``k_3`` (Pdc) is
**the only direction where they separate**. The yield maximum there is the
Pdc-null end (0.228 g/g at ``k_3`` 0.0077, PI −0.007), while PI peaks at
``k_3`` 1.906. Along ``k_13``, PI and isobutanol yield are both monotone within
the band and peak together at the 4.0 g/L/h ceiling (PI 0.8084). Just beyond
the band, PI is 0.8113 at 4.063, and PI is positive from ``k_13`` 0.389
upward.

The Pdc slice
~~~~~~~~~~~~~

.. figure:: /_static/images/isobutanol/landscape/pdc_slice_through_1912.png
   :target: ../../_static/images/isobutanol/landscape/pdc_slice_through_1912.png
   :class: pub-figure
   :alt: Line chart against Pdc capacity k3 on a log axis from about 0.006 to
         23 g/L/h. Profitability index (red) starts near zero, drops to about
         -0.8 in a flat-bottomed band between 0.03 and 0.1, recovers, dips to
         about -0.95 near 0.45, then rises steeply to a peak of about 0.81
         near 1.9 before declining to about 0.42 at the top of the band. IRR
         (purple) tracks PI on an outer axis. Ethanol yield (blue) rises
         monotonically from zero to about 0.32 g/g. A dotted vertical line
         near 4.85 marks the onset of burden derating.
   :width: 85%
   :align: center

   **PI (red, left axis), ethanol yield (blue, right axis) and IRR (purple,
   outer axis) along Pdc capacity** ``k_3`` **through #1912.** The sweep has
   219 points over the campaign band, 0.00581–23.24 g/L/h (10\ :sup:`−3`\ × to
   4× wild type, where wild type is 5.81). The other eleven decisions stay at
   #1912, with the burden on. Dotted: growth derated by enzyme burden from
   ``k_3`` 4.85 (0.835× wild type). **What to look at:** the peak, PI 0.8086
   (IRR 27.3 %) at ``k_3`` 1.865–1.906 (0.32–0.33× wild type). PI ≥ 0.75 only
   at 0.31–0.37× wild type, and wild-type Pdc still gives PI 0.653
   (IRR 25.1 %). The **flat-bottomed drop between** ``k_3`` **≈ 0.031 and
   0.098 g/L/h** is a *trace-ethanol recovery regime*. There, broth ethanol is
   0.003–2.3 g/L, enough to switch on the ethanol recovery train (boiler
   +25.7 MM$, TCI +33 MM$). The model always recovers trace ethanol and never
   weighs recovery against purging, so this is a separation-design regime,
   not biology, and not evidence of objective roughness. The small sawtooth
   follows the discrete fermentation-time selection (PI vs τ Spearman 0.93 on
   the Pdc-null plateau). Ethanol yield, a single-product metric, is trivially
   monotone here.

.. list-table:: Regions of the Pdc slice through #1912
   :header-rows: 1
   :widths: 30 20 22 28

   * - Region
     - ``k_3`` (× wild type)
     - PI
     - Notes
   * - Pdc-null plateau
     - ≤ 0.005×
     - −0.095 to +0.03
     - IRR 13.4–15.5 %; ethanol yield ≤ 2 × 10\ :sup:`−6` g/g
   * - Trace-ethanol recovery regime
     - 0.005–0.017×
     - −0.71 to −0.84
     - Broth ethanol 0.003–2.3 g/L; boiler +25.7 MM$, TCI +33 MM$
   * - Smooth trough
     - 0.078×
     - −0.954
     - Minimum of the slice
   * - **Peak**
     - **0.32–0.33×**
     - **0.8086**
     - IRR 27.3 %; total alcohol yield 0.403 g/g
   * - Near wild type
     - 1.01×
     - 0.653
     - IRR 25.1 %
   * - Band top
     - 4×
     - 0.420
     - Burden derating from 0.835× wild type

The trace-ethanol regime is not specific to this slice. Across the pooled
seed-350 campaigns and the TRY-informed relay, **0 of 1,538**
isobutanol-producing designs with broth ethanol between 0.002 and 0.5 g/L are
profitable.

.. dropdown:: Profitability of isobutanol-producing designs by broth ethanol (seed-350 pool + relay)

   Designs with isobutanol yield > 0.1 g/g in the pooled seed-350 campaigns
   plus the TRY-informed relay, grouped by broth ethanol titer:

   .. list-table::
      :header-rows: 1
      :widths: 30 20 25 25

      * - Broth ethanol (g/L)
        - n
        - Profitable (PI > 0)
        - Max PI
      * - ≤ 0.002
        - 1,228
        - 11.6 %
        - 0.504
      * - 0.002–0.5
        - 1,538
        - **0 %**
        - −0.30
      * - 0.5–2.5
        - 439
        - 0.2 %
        - 0.21
      * - 2.5–10
        - 232
        - 18 %
        - 0.68
      * - 10–30
        - 231
        - 60 %
        - 0.77
      * - > 30
        - 391
        - 94 %
        - 0.81

The qualitative design rule is defensible: **make either essentially no
ethanol or a real co-product.** Its magnitude is specific to this
separation design.

Caveats and limitations
-----------------------

- **The pool is biased by construction.** The FDC, variogram, distance-profile
  and t-SNE results describe where seven GP campaigns chose to sample. The
  needle-vs-smooth contrast built into the three landscape plotting scripts'
  docstrings and file names fails on unbiased samples. Use the figures for
  topology (islands, plateaus), not roughness.
- **The slices are anchor-specific.** A 1-D slice fixes the other eleven
  decisions at #1912. At the same Pdc capacity, pooled designs span PI −5.1 to
  +0.79, so the trough and peak positions belong to this anchor. The
  Pdc-null benchmark is therefore the pooled best, PI 0.52 among 31,905
  designs (:doc:`kinetic_optimization`), not this slice's Pdc-null plateau.
- **Bounds.** Six of #1912's twelve coordinates sit at a band bound, so
  yield–PI coincidence along those axes is partly trivial.
- **Unexplained right edge.** The left edge of the trace-ethanol regime
  matches the switch-on of trace-ethanol recovery. The right edge (around
  ``k_3`` 0.1 g/L/h, where the boiler drops back) is **not explained**. One
  candidate is a heat-exchanger-network topology change, which would be
  build-sensitive. Heating duties were not recorded.
- **Build.** All data are 2b5b27d. Under the current build, the TEA metrics of
  co-production designs shift by a few percent (see :doc:`kinetic_optimization`).
- **Figures.** The ``k_13`` sweep figure is not shown; its results are
  given here in text.

Reproduce it
------------

.. important::

   No analysis results ship with NSKinetics: the seven ``_rs350``
   trajectory CSVs and the sweep CSVs must be regenerated (the campaigns are
   long, opt-in runs; see :doc:`kinetic_optimization`), or
   ``NSKINETICS_IBO_RESULTS_DIR`` set to a results folder that holds them.
   The slices reproduce #1912 only under the campaigns' model version
   (2b5b27d). Script paths are relative to the package folder.

.. list-table::
   :header-rows: 1
   :widths: 34 44 22

   * - Script
     - What it does / writes
     - Cost
   * - ``plots/plot_objective_landscape_fitness_distance.py``
     - FDC figure (``objective_landscape_fitness_distance.png``).
       ``--pi-scale linear`` draws a plain linear PI axis.
     - Sim-safe (CSV reads)
   * - ``plots/plot_objective_landscape_needle_plateau.py``
     - Distance-profile figure (``objective_landscape_needle_plateau.png``);
       ``--out-dir``, ``--stem``, ``--dpi``
     - Sim-safe
   * - ``plots/plot_objective_landscape_embedding.py``
     - t-SNE figure (``objective_landscape_embedding.png``). It reuses the
       cached fit, keyed by a hash of the coordinates, settings and sklearn
       version; ``--recompute`` forces a fresh fit.
     - Sim-safe
   * - ``analyses/objective_landscape/evaluate_axis_try_informed_optimum.py``
     - 1-D sweeps through #1912. Set the mode with the ``IBO_AXIS_SWEEP``
       environment variable: ``screen`` (default; 11 variables, coarse) or a
       variable name such as ``k_3``. ``IBO_AXIS_N_POINTS`` (default 160)
       sets the resolution; the ``k_3`` figure used 218 plus the anchor.
       Writes ``evaluate_axis_try_informed_optimum_<name>.csv`` plus an
       anchor JSON. Run it in a fresh process under ``supervise_sweep.py``.
     - Simulates; opt-in
   * - ``analyses/objective_landscape/evaluate_k13_try_informed_optimum.py``
     - The 161-point ``k_13`` sweep through #1912
     - Simulates; about 15–25 min
   * - ``plots/plot_axis_try_informed_sweep.py``
     - ``--param k_3`` draws the Pdc-slice figure above. ``--screen`` prints
       the slice-roughness ranking and draws a small-multiples overview.
     - Sim-safe
   * - ``plots/plot_k13_try_informed_sweep.py``
     - PI and isobutanol yield along ``k_13``
     - Sim-safe

All plotting scripts read CSVs through the data layer, which they load by file
path. They never import ``biosteam`` or build the biorefinery, so they can run
alongside a simulation.
