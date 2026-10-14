Uncertainty analysis
====================

A 6,000-sample Monte Carlo of scenario A over prices, plant context, corn
composition and all 40 kinetic parameters, with Spearman analysis. Commands
run from the package folder (:doc:`running_analyses`).

.. code-block:: bash

   # the Monte Carlo run: ~30 h, resumable (autosaves every 10 samples)
   python analyses/uncertainty/uncertainties_IBO_EtOH.py
   # the figures (read the newest *_1_full_evaluation.xlsx in the results root)
   python plots/plot_uncertainty_MPSP_vs_TCI.py
   python plots/plot_uncertainty_spearman_tornado.py

.. figure:: /_static/images/isobutanol/uncertainty/uncertainty_bivariate_mesp_scenario_a.png
   :target: ../../_static/images/isobutanol/uncertainty/uncertainty_bivariate_mesp_scenario_a.png
   :alt: Six-panel figure of kernel-density clouds from 6,000 Monte Carlo samples. MESP sits on the shared y-axis, against total capital investment, ethanol production, DDGS revenue, ethanol yield, ethanol titer and ethanol productivity. Each cloud is compact around a white baseline diamond near 3.8 dollars per GGE. Panel A carries a red Pareto staircase. Panels B to F carry black binned-median lines that slope only gently downward.
   :width: 100%
   :class: pub-figure

   MESP against capital, scale, coproduct revenue and titer, rate and yield;
   black lines are binned medians.

.. figure:: /_static/images/isobutanol/uncertainty/spearman_tornado_scenario_a.png
   :target: ../../_static/images/isobutanol/uncertainty/spearman_tornado_scenario_a.png
   :alt: Five tornado charts of Spearman rank correlations. MESP is driven by corn price (+0.64), DDGS price (−0.48), operating days, capacity, starch, the kinetic parameter k_1ie (+0.18) and tax. TCI is driven mainly by plant capacity (+0.93). Ethanol titer, productivity and yield are driven entirely by kinetic parameters such as k_7ie, k_1ie, k_7 and k_10.
   :width: 100%
   :class: pub-figure

   Spearman correlations (:math:`|\rho| \geq 0.10`, p < 0.05) for MESP, TCI
   and the three ethanol TRY metrics, coloured by input category.

Context, not kinetics, sets cost risk: prices, uptime, scale and starch
carry ~92 % of the MESP variance and all kinetic parameters ~5 %. The kinetic
parameters that drive TRY barely register in MESP.
