Feeding strategies
==================

A 40 × 40 sweep of the fed-batch policy (glucose threshold × target, with
the spike cap optimized for MESP at every point) for scenario A. Commands run
from the package folder (:doc:`running_analyses`).

.. code-block:: bash

   # the sweep: thousands of simulations, more than a day; not checkpointed
   IBO_SWEEP_SCENARIO=A IBO_SWEEP_STEPS=40 \
       python analyses/feeding_strategies/evaluate_feeding_strategies.py
   # the figure (simulates scenario A once for panel B, then caches it)
   python plots/plot_feeding_strategy_baseline_A.py

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

   **(A)** MESP over the valid threshold × target policies; the star is the
   optimum (226/230 g/L, 4 spikes, 3.82 $/GGE). **(B)** The scenario-A
   baseline fermentation.

Sugar loading, not feed mode, sets cost: a target 10 g/L above the slurry's
own concentration switches on feed evaporation (+4.6 % MESP), while a plain
batch at 230 g/L is only 0.41 % worse than the best fed-batch policy.
