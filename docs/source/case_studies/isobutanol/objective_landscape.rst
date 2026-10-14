Objective landscape
===================

The shape of the objectives over the 12-d design space, from the pooled
trials of the seed-350 campaigns (:doc:`kinetic_optimization`) and from 1-D
slices through the flagship design #1912. Commands run from the package
folder (:doc:`running_analyses`).

Landscape of the optimizer pool
-------------------------------

.. code-block:: bash

   # sim-free; they read the seven _rs350 campaign CSVs
   python plots/plot_objective_landscape_fitness_distance.py
   python plots/plot_objective_landscape_needle_plateau.py
   python plots/plot_objective_landscape_embedding.py   # the t-SNE fit is cached after the first run

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

   Fitness–distance structure of the 13,622 pooled trials: profitable (red)
   trials sit in separate groups, not on one hill.

.. figure:: /_static/images/isobutanol/landscape/distance_profiles_around_842.png
   :target: ../../_static/images/isobutanol/landscape/distance_profiles_around_842.png
   :class: pub-figure
   :alt: Two panels sharing an x axis of distance from trial 842. (a) PI of
         each trial with binned median and interquartile band; the median
         drops below zero within about 0.2 of trial 842. (b) Isobutanol yield
         of the same trials; the median stays near 0.37 out to about 0.8 and
         then falls. The neighbour trial 1882 is ringed in both panels.
   :width: 100%

   Around the most profitable pooled trial (#842), PI collapses within a short
   distance while isobutanol yield stays on its plateau.

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

   t-SNE map of the pool coloured by PI **(a)** and isobutanol yield **(b)**.

Profitable designs form small, separated islands (an ethanol-only shelf, a
Pdc-null isobutanol island and the co-production island). PI looks
needle-like only in this optimizer-sampled pool; on unbiased Sobol' samples
it is as smooth as the TRY objectives.

The Pdc slice through #1912
---------------------------

.. code-block:: bash

   # simulates; requires hensmith 2b5b27d (otherwise the anchor check aborts)
   python analyses/objective_landscape/evaluate_axis_try_informed_optimum.py \
       --sweep k_3 --n-points 218
   python plots/plot_axis_try_informed_sweep.py --param k_3

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
   :width: 100%
   :figclass: ibo-fig-measure

   PI, ethanol yield and IRR along Pdc capacity ``k_3`` through #1912.

PI peaks at 0.32–0.33× wild-type Pdc, so the best design knocks Pdc down
rather than deleting it; the trench at ``k_3`` ≈ 0.03–0.1 is a trace-ethanol
regime of the separation design, not objective roughness.
