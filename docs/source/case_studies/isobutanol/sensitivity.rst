Global sensitivity analysis
===========================

Sobol' closed indices and Shapley effects of the 12 ``metabolic_split_12d``
decisions, computed on surrogates fitted to 4,096 feasible simulated designs
per sampling measure (campaign, screening and custom priors). Commands run
from the package folder (:doc:`running_analyses`).

.. code-block:: bash

   # stage 1: simulate the Sobol' design (~20 h per measure, resumable)
   python analyses/sensitivity/evaluate_sobol_split12d.py                                   # campaign
   IBO_SOBOL_MEASURE=screening python analyses/sensitivity/evaluate_sobol_split12d.py       # screening
   IBO_SOBOL_MEASURE=custom IBO_SOBOL_BANDS=linear_0-4x_2026-09-22.json \
       IBO_SOBOL_SEED=20260922 python analyses/sensitivity/evaluate_sobol_split12d.py      # custom

   # stage 2: surrogates, indices, tables and figures (sim-free)
   python analyses/sensitivity/analyze_sobol_split12d.py --study-name <STUDY_NAME>

.. list-table:: Shapley effects (% of variance) of the two divergent levers
   :header-rows: 1
   :stub-columns: 1
   :widths: 40 20 20 20

   * - Lever → metric
     - Campaign
     - Screening
     - Custom
   * - ALS ``k_13`` → **PI**
     - 2.6
     - 3.5
     - 1.3
   * - ALS ``k_13`` → isobutanol titer
     - 27.7
     - 23.1
     - 15.3
   * - Adh1 ``k_6`` → **PI**
     - 43.1
     - 10.2
     - 12.2
   * - Adh1 ``k_6`` → isobutanol titer
     - 0.2
     - 0.8
     - 1.6

Titer's driver is not profit's driver: ALS carries up to 29 % of isobutanol
titer and yield variance but ≤ 3.5 % of profitability variance, and Adh1 the
reverse, under every measure. Since ≥ 95 % of Var(PI) comes from designs worse
than the baseline, the indices rank failure modes, not the rare profitable
designs.
