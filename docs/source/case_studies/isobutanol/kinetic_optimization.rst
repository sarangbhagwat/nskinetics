Kinetic optimization
====================

Bayesian optimization (GP, seed 350) of 12 strain and feeding decisions
(``metabolic_split_12d``: six enzyme capacities, three product inhibitions,
three feeding settings) for six titer, rate and yield (TRY) objectives and
for profitability (PI), with the enzyme burden on. Commands run from the
package folder (:doc:`running_analyses`); each campaign takes ~8–17 h, and
the published numbers are on hensmith build ``2b5b27d``
(:ref:`isobutanol-switch-build`).

.. code-block:: bash

   KO=analyses/kinetic_optimization/optimize_kinetics_BO_supervised.py
   P=kin_opt_ethanol_isobutanol_metabolic_split_12d
   S=gp_rb0.001-4_ib0.75-1.5_aA_rs350_burden
   COMMON="--study-target-products ethanol_isobutanol --study-type metabolic_split_12d
           --method gp --gp-deterministic --gp-no-learned-constraints --seed 350"

   # the uninformed profitability campaign and the six TRY scouts (2,000 trials each)
   python $KO --objective "PI (log-tail)" $COMMON --n-startup-trials 50 --n-trials 2000 \
       --study-name ${P}_pi_log-tail_${S}
   for obj in "IBO yield" "IBO titer" "IBO productivity" \
              "EtOH yield" "EtOH titer" "EtOH productivity"; do
     slug=$(echo "$obj" | tr 'A-Z ' 'a-z_')
     python $KO --objective "$obj" $COMMON --n-startup-trials 50 --n-trials 2000 \
         --study-name ${P}_${slug}_${S}
   done

   # the TRY-informed relay: preload the scouts' trials, then 1,000 more
   python $KO --objective "PI (log-tail)" $COMMON --n-trials 1000 \
       --relay-from ${P}_ibo_yield_${S} ${P}_ibo_titer_${S} ${P}_ibo_productivity_${S} \
                    ${P}_etoh_yield_${S} ${P}_etoh_titer_${S} ${P}_etoh_productivity_${S} \
       --relay-max-rows 1000 --relay-keep-above -0.12953

   # the figures (sim-free; they read the campaign CSVs)
   python plots/optimization_figure_ideas/fig_main_scouts.py
   python plots/optimization_figure_ideas/fig_s2_fingerprints.py

.. figure:: /_static/images/isobutanol/optimization/kinbo_main_scouts_campaigns.png
   :target: ../../_static/images/isobutanol/optimization/kinbo_main_scouts_campaigns.png
   :class: pub-figure
   :alt: Six-panel figure. (a) Best IRR so far against simulated trials on a
         log axis: the TRY-informed campaign jumps to 19.8 % on trial 2 and
         climbs to 27.3 %, while the uninformed campaign plateaus at 15.8 %
         from trial 104. (b) IRR against the isobutanol share of alcohol titer:
         every design above 15.8 % makes isobutanol. (c) Per-campaign strips
         of IRR, with open best-visit and filled returned-design markers; the
         isobutanol scouts visit 19–23 % but return low or no IRR. (d) Titers
         and (e) metabolic proteome of the returned designs, with the
         penalty-free budget line. (f) Share of trials making isobutanol
         against IRR.
   :width: 100%

   Best IRR so far in each campaign and the designs they return. In the
   figure, "seeds" are the 1,000 preloaded donor trials, not random seeds.

.. figure:: /_static/images/isobutanol/optimization/kinbo_s2_design_fingerprints.png
   :target: ../../_static/images/isobutanol/optimization/kinbo_s2_design_fingerprints.png
   :class: pub-figure
   :alt: Heatmap with one row per design (starting strain, the returned design
         of each of the eight campaigns, and the isobutanol scouts' best
         visits) and twelve columns: native capacities as fold change, the
         isobutanol-pathway capacities in absolute units, three inhibition
         multipliers and three feeding variables. Outlined cells sit at a
         search bound. Right-hand columns give IRR, each campaign's own
         objective, the metabolic proteome and the growth factor.
   :width: 100%

   Decision vectors of the returned designs; outlined cells sit at a search
   bound.

Eleven of twelve TRY optima miss the 15 % IRR hurdle. The uninformed
profitability search stalled on an ethanol-only design (15.8 % IRR), while
the relay returned a co-producer at 27.3 % IRR (#1912: 39.6 g/L isobutanol +
40.6 g/L ethanol, Pdc knocked down to 0.33×) at the starting strain's
capital.

.. _ibo-opt-scope:

.. warning::

   #1912 is a **prospective, bounds-limited target**: its isobutanol titer is
   15–20× the best reported in *S. cerevisiae*
   (:ref:`Wess et al. 2019 <ibo-ref-wess2019>`,
   :ref:`Lane et al. 2020 <ibo-ref-lane2020>`;
   :ref:`Su et al. 2021 <ibo-ref-su2021>` for tolerance), it sits at the
   search bounds, and the result rests on one isobutanol price and two random
   seeds.

.. _ibo-opt-scenarios:

The ``opt_*`` scenarios
-----------------------

Each campaign's best trial is frozen as an ``opt_*`` scenario (pins on
:doc:`index`) and checked by smoke cases 9–16:

.. code-block:: bash

   python analyses/opt_scenarios/build_opt_split12d_workbooks.py opt_PI_TRY_informed   # rewrites the workbook
   python analyses/opt_scenarios/verify_opt_scenarios.py opt_PI_TRY_informed          # read-only MATCH check

.. code-block:: python

   from nskinetics.biorefineries import isobutanol

   isobutanol.load()
   bundle = isobutanol.load_scenario('opt_PI_TRY_informed')
   print(bundle['solve_TEA']())
