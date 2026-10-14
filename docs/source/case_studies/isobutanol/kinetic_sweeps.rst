Kinetic sweeps
==============

40 × 40 sweeps of kinetic-parameter multipliers at the scenario-A process;
every grid point is a full fermentation, flowsheet simulation and TEA.
Commands run from the package folder (:doc:`running_analyses`).

Glycolysis × ethanol inhibition
-------------------------------

.. code-block:: bash

   # ~7–11 h, checkpointed; the supervisor relaunches it after a crash
   python analyses/supervise_sweep.py \
       analyses/kinetic_sweeps/evaluate_EtOH_glycolysis_inhib_ethanol.py
   python plots/plot_glycolysis_inhib_ethanol_A.py

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

         IRR over glycolytic capacity × ethanol inhibition; the star is the
         best design (13.7 %), the diamond the baseline strain.

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

         The same plane as MESP; the optimum is 3.788 $/GGE at (1.334×,
         0.770×).

The baseline strain sits within 1.5 % of the best design on this plane, and
removing ethanol inhibition entirely *raises* MESP by 21 %, because the extra
tolerance feeds growth once the 16-spike feeding cap binds.

``k_1e`` × ethanol inhibition
-----------------------------

.. code-block:: bash

   # ~7–11 h, not checkpointed
   python analyses/kinetic_sweeps/evaluate_EtOH_k1e_inhib_ethanol.py
   python plots/plot_k1e_inhib_ethanol_A.py

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

   MESP over ``k_1e`` alone × ethanol inhibition.

Scaling ``k_1e`` alone tells the same story: the plane optimum is only 1.50 %
below the baseline. The Pdc and Adh1 planes
(``evaluate_EtOH_k3_inhib_ethanol.py``, ``evaluate_EtOH_k6_inhib_ethanol.py``)
agree: ethanol-only strain levers have little headroom.
