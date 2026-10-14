.. _isobutanol-process-and-tea:

Process and techno-economics
============================

The flowsheet and TEA basis that every other analysis inherits. Commands run
from the package folder (:doc:`running_analyses`).

Flowsheet
---------

.. code-block:: bash

   python plots/plot_block_flow_diagram.py
   python plots/conceptual_feeding_aeration_sugar_prep.py

.. figure:: /_static/images/isobutanol/overview/block_flow_diagram.png
   :target: ../../_static/images/isobutanol/overview/block_flow_diagram.png
   :class: pub-figure
   :width: 100%
   :alt: Block flow diagram. Panel A, fermentation process: saccharified corn slurry is split into two multi-effect evaporation and dilution trains giving an initial feed and a spike feed to the fed-batch fermentor, which receives compressed air and yeast; vented gases are scrubbed and the recovered alcohols join the beer sent to distillation. Panel B, separation process: beer column, ethanol rectifier, molecular sieves and denaturant mixing give ethanol; stripper, decanter and drying column give isobutanol above 99.9 wt%; stillage goes through centrifugation, evaporation, oil separation and drying to DDGS and corn oil.

   Block flow of the biorefinery: **A**, the slurry is split into an
   initial-feed and a spike-feed train and fermented in ``V406``; **B**,
   solvent-free distillation recovers ethanol and isobutanol (dashed units run
   only when isobutanol is made) and the stillage becomes DDGS and corn oil.

.. figure:: /_static/images/isobutanol/overview/sugar_prep_feeding_aeration_concept.png
   :target: ../../_static/images/isobutanol/overview/sugar_prep_feeding_aeration_concept.png
   :class: pub-figure
   :width: 100%
   :alt: Three-panel concept diagram. Panel a: splitter feeding two evaporator, dilution-mixer and cooler trains into the fed-batch fermentor, with a fed-batch strategy specification setting the actuators and a compressor and valve supplying air during stage 1. Panel b: schematic glucose, ethanol, isobutanol and cell trajectories in one batch, with threshold-triggered spikes, a stage-1 aerobic window and the batch end at the first glucose minimum. Panel c: flowchart of one simulation, from specifying the feeding strategy through simulating the kinetics, closing feed and air balances and simulating the full flowsheet until the state stops drifting, with a failure-recovery cascade.

   How the fermentation section is controlled: **a**, the feeding
   specification root-solves each train's evaporator and dilution mixer;
   **b**, glucose spikes restore the target until the cap; **c**, one
   simulation is a fixed-point sweep with a recovery cascade.

The corn dry-grind base comes from ``biorefineries.corn``; the NSKinetics
fermentor ``V406`` and its two feed trains are built by
:obj:`~nskinetics.processes.create_sugar_prep_and_fermentation_system`.

Baseline fermentation (scenario A)
----------------------------------

.. code-block:: bash

   python plots/plot_feeding_strategy_baseline_A.py   # simulates A once, then caches the trajectory

.. figure:: /_static/images/isobutanol/overview/baseline_A_fermentation_trajectory.png
   :target: ../../_static/images/isobutanol/overview/baseline_A_fermentation_trajectory.png
   :class: pub-figure
   :width: 100%
   :figclass: ibo-fig-measure
   :alt: Concentration versus time for the scenario A baseline fermentation over 80 hours: glucose starts near 221 g/L with small ripples from seven feed spikes between 9 and 11.5 hours, then falls to near zero by about 60 hours; ethanol rises to about 115 g/L; cell density plateaus near 13 g/L; acetate and isobutanol stay near zero. Dashed lines mark the end of aeration near 12 hours and the end of fermentation near 62 hours.

   Scenario A in ``V406``: 7 glucose spikes during the aerobic stage, then
   batch to τ = 61.6 h with 114.5 g/L ethanol (per litre of water).

TEA basis
---------

Capital is costed per unit from bare-module factors (no Lang factor), at
CEPCI 797.9 and 2023 prices (ethanol 0.835 $/kg, isobutanol 1.500 $/kg), with
a 15 % IRR hurdle. Scenario A's MESP is 3.83 $/GGE (IRR 12.8 %); fermentation
is only ~10 % of installed equipment cost while corn is ~86 % of material
cost.

TEA breakdowns at the campaign optima
-------------------------------------

.. code-block:: bash

   python analyses/tea_breakdowns/collect_tea_breakdowns_split12d.py   # ~4 min; needs the campaign CSVs
   python plots/plot_tea_breakdowns_split12d.py                         # 3 × 3 breakdown figure

The flagship co-producer (#1912, :doc:`kinetic_optimization`) lifts sales by
19 % at 1.2 % lower capital than scenario A, so its gain is all revenue; the
isobutanol-yield optimum needs 2.27× the capital because its broth is slow and
dilute.
