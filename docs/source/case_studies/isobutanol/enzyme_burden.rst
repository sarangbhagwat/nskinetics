Enzyme burden
=============

Each enzyme capacity costs proteome. A closed-form proteome-allocation
constraint sums the enzyme pools a design needs, closes them against the
cell's flexible protein sector and derates growth as that sector fills
(constants from :ref:`Metzl-Raz et al. 2017 <ibo-ref-metzlraz2017>`,
:ref:`Scott et al. 2010 <ibo-ref-scott2010>`,
:ref:`Sánchez et al. 2017 <ibo-ref-sanchez2017>` and
:ref:`Jach et al. 2022 <ibo-ref-jach2022>`). It is on for scenario A (where it
is inert), the ``opt_*`` scenarios, the kinetic sweeps, the optimization
campaigns and the Sobol' analysis.

.. code-block:: python

   from nskinetics.biorefineries.isobutanol import enzyme_burden as eb

   print(round(eb.F_FLEX - eb.PHI_T_WT, 4))       # 0.1347 g/gDCW, the penalty-free budget
   print(round(eb.ehrlich_unit_cost('r13'), 4))   # 0.0056 g enzyme/gDCW per unit k_13

.. code-block:: bash

   # the figure: panels c–e of the campaign figure (sim-free; reads the campaign CSVs)
   python plots/optimization_figure_ideas/fig_main_scouts.py

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

   What each campaign returned and what it cost in proteome; the
   TRY-informed design (★) stays inside the budget.

The isobutanol- and ethanol-yield optima run at 157 % and 130 % of the
penalty-free budget, throttling growth into a slow, dilute and
capital-intensive broth; the flagship profit optimum uses 97 % with no
derating.
