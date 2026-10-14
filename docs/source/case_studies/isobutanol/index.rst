.. _isobutanol-case-study:

Isobutanol biorefinery
======================

.. toctree::
   :maxdepth: 1
   :hidden:

   process_and_tea
   feeding_strategies
   uncertainty
   kinetic_sweeps
   kinetic_optimization
   objective_landscape
   enzyme_burden
   sensitivity
   running_analyses

``nskinetics.biorefineries.isobutanol`` is a corn dry-grind biorefinery,
built in `BioSTEAM <https://biosteam.readthedocs.io/>`__, that makes fuel
ethanol, isobutanol, DDGS and corn oil. Its fermentor, ``V406``, integrates
an NSKinetics fed-batch model of *Saccharomyces cerevisiae* (the
:ref:`Lei et al. 2001 <ibo-ref-lei2001>` core plus an engineered Ehrlich
isobutanol pathway), so titer, yield, rate and fermentation time are
*outputs* of the strain kinetics and the feeding policy, and they set the
size and cost of everything downstream. Every analysis below ships with the
package and can be re-run (:doc:`running_analyses`).

Analyses
--------

- :doc:`process_and_tea`: the flowsheet, the TEA basis and the cost
  breakdowns of the baseline and the optimized strains.
- :doc:`feeding_strategies`: which fed-batch policy minimizes the ethanol
  selling price.
- :doc:`uncertainty`: how much of the baseline plant's cost risk comes from
  the strain rather than from prices and plant context.
- :doc:`kinetic_sweeps`: how much the classic ethanol-strain levers
  (glycolysis, ethanol tolerance) can pay.
- :doc:`kinetic_optimization`: Bayesian optimization of 12 strain and
  feeding decisions for titer, rate, yield or profitability.
- :doc:`objective_landscape`: why profitable designs are hard to find.
- :doc:`enzyme_burden`: the proteome-allocation constraint on strain designs.
- :doc:`sensitivity`: a Sobol'/Shapley ranking of the same 12 decisions.

The kinetic model and the process
---------------------------------

.. tab-set::

   .. tab-item:: Kinetic model inside V406

      .. figure:: /_static/images/isobutanol/overview/kinetic_model_conceptual_diagram.png
         :target: ../../_static/images/isobutanol/overview/kinetic_model_conceptual_diagram.png
         :class: pub-figure
         :width: 100%
         :alt: Reaction network of the S. cerevisiae fed-batch model: glucose to pyruvate (r1), pyruvate to acetaldehyde (r3) and ethanol (r6), acetate and respiration branches (r2, r4, r5, r8), growth and biomass states (r7, r9 to r11), and the engineered isobutanol pathway from pyruvate through alpha-acetolactate, 2,3-dihydroxyisovalerate, alpha-ketoisovalerate and isobutyraldehyde to isobutanol (r13 to r17), with product-inhibition, glucose-repression and oxygen-dependence markers, plus boxes for fed-batch glucose feeding and two-stage aeration control.

         The *S. cerevisiae* fed-batch model in ``V406``: glycolysis (r1),
         the ethanol branch (Pdc1, Adh1; r3, r6), respiration and the
         engineered isobutanol pathway (r13–r17), with product inhibition
         and two-stage aeration. Drawn by
         ``nskinetics/models/s_cerevisiae_ferm_fb_inhib_mod_ibo/conceptual_diagram.py``.

   .. tab-item:: Process around it

      .. figure:: /_static/images/isobutanol/overview/block_flow_diagram.png
         :target: ../../_static/images/isobutanol/overview/block_flow_diagram.png
         :class: pub-figure
         :width: 100%
         :alt: Block flow diagram. Panel A, fermentation process: saccharified corn slurry is split into two multi-effect evaporation and dilution trains giving an initial feed and a spike feed to the fed-batch fermentor, which receives compressed air and yeast; vented gases are scrubbed and the recovered alcohols join the beer sent to distillation. Panel B, separation process: beer column, ethanol rectifier, molecular sieves and denaturant mixing give ethanol; stripper, decanter and drying column give isobutanol above 99.9 wt%; stillage goes through centrifugation, evaporation, oil separation and drying to DDGS and corn oil.

         Block flow of the biorefinery: **A**, sugar preparation and the
         fed-batch fermentor; **B**, ethanol and isobutanol recovery and the
         DDGS train. Drawn by ``plots/plot_block_flow_diagram.py``.

Scenarios
---------

A *scenario* bundles a parameter workbook (baseline kinetics plus
uncertainty distributions), a feeding strategy and pinned selling prices,
registered in ``nskinetics.biorefineries.isobutanol.scenarios.SCENARIOS``:
``A`` (wild-type ethanol strain, fed-batch), ``B`` (an engineered
ethanol + isobutanol co-producer, batch) and eight ``opt_*`` scenarios, each
freezing the best design of one optimization campaign
(:ref:`ibo-opt-scenarios`).

.. dropdown:: Pinned MPSPs of the ten scenarios ($/kg, purity-adjusted, 15 % IRR)

   .. list-table::
      :header-rows: 1
      :widths: 34 46 10 10

      * - Scenario
        - What it is
        - Ethanol
        - Isobutanol
      * - ``A``
        - Wild-type ethanol strain, fed-batch baseline
        - 0.86021
        - NaN
      * - ``B``
        - Engineered co-producer, batch (reference build; its Adh6 constant
          ``k_17`` is still too high, so treat it as provisional)
        - 0.59530
        - 1.2800
      * - ``opt_PI_uninformed``
        - Profitability campaign started from scratch
        - 0.82091
        - NaN
      * - ``opt_PI_TRY_informed``
        - Profitability campaign preloaded with the TRY campaigns' trials
          (the flagship co-producer)
        - 0.46066
        - 1.1123
      * - ``opt_IBO_titer`` / ``_yield`` / ``_productivity``
        - Isobutanol titer, yield and productivity optima
        - —\ :sup:`*`
        - 2.4625 / 1.8533 / 2.2694
      * - ``opt_EtOH_titer`` / ``_yield`` / ``_productivity``
        - Ethanol titer, yield and productivity optima
        - 0.85672 / 0.86838 / 0.90813
        - NaN

   :sup:`*` A near-zero co-product gets an absurd purity-adjusted MPSP.
   NaN means the product is not made.

Quick start
-----------

Install from source with the ``isobutanol`` extra
(:ref:`isobutanol-installation`), then build the flowsheet, load a scenario
and solve its TEA:

.. code-block:: python

   from nskinetics.biorefineries import isobutanol

   isobutanol.load()                         # build + baseline-simulate (once per process)
   bundle = isobutanol.load_scenario('A')    # workbook kinetics + feeding strategy
   print(bundle['solve_TEA']())              # {'IRR': ..., 'MPSPs': {'ethanol': ..., 'isobutanol': nan}}

   # re-simulate as a plain batch at 230 g/L, then re-solve
   bundle['model_specification'](threshold_conc=226.1, target_conc=230.,
                                 max_n_spikes=0)
   print(bundle['solve_TEA']())

Run simulating scripts one at a time, and one ``load()`` per process.

.. _isobutanol-conventions:

Conventions
-----------

- **MPSP** is the purity-adjusted minimum product selling price ($/kg) at a
  15 % IRR. **MESP** is the ethanol MPSP in $/GGE (= 4.45775 × $/kg).
- **IRR** is at market prices (ethanol 0.835 $/kg, isobutanol 1.500 $/kg,
  2023$). **PI** = NPV at the 15 % hurdle ÷ TCI, so PI > 0 ⇔ IRR > 15 %.
- Titers and productivities are per litre of **water**; yields are g per g
  of glucose added.
- **Build tags** (``2b5b27d``, ``6d4776f``, ``95e9362``) are commits of the
  heat-exchanger-network package hensmith. They shift TEA results by up to a
  few percent, so numbers from different builds are not compared
  (:ref:`isobutanol-switch-build`).
- **Trial numbers** such as ``#1912`` are the 0-based ``trial_number`` of the
  campaign CSVs.

.. _isobutanol-references:

References
----------

.. _ibo-ref-jach2022:

Jach, M. E., Serefko, A., Ziaja, M. & Kieliszek, M. (2022). Yeast protein as
an easily accessible food source. *Metabolites* **12**, 63.
https://doi.org/10.3390/metabo12010063

.. _ibo-ref-kondo2012:

Kondo, T., Tezuka, H., Ishii, J., Matsuda, F., Ogino, C. & Kondo, A.
(2012). Genetic engineering to enhance the Ehrlich pathway and alter carbon
flux for increased isobutanol production from glucose by *Saccharomyces
cerevisiae*. *J. Biotechnol.* **159**, 32–37.
https://doi.org/10.1016/j.jbiotec.2012.01.022

.. _ibo-ref-lane2020:

Lane, S., Zhang, Y., Yun, E. J., Ziolkowski, L., Zhang, G., Jin, Y.-S. &
Avalos, J. L. (2020). Xylose assimilation enhances the production of
isobutanol in engineered *Saccharomyces cerevisiae*. *Biotechnol. Bioeng.*
**117**, 372–381. https://doi.org/10.1002/bit.27202

.. _ibo-ref-lei2001:

Lei, F., Rotbøll, M. & Jørgensen, S. B. (2001). A biochemically structured
model for *Saccharomyces cerevisiae*. *J. Biotechnol.* **88**, 205–221.
https://doi.org/10.1016/S0168-1656(01)00269-3

.. _ibo-ref-metzlraz2017:

Metzl-Raz, E., Kafri, M., Yaakov, G., Soifer, I., Gurvich, Y. & Barkai, N.
(2017). Principles of cellular resource allocation revealed by
condition-dependent proteome profiling. *eLife* **6**, e28034.
https://doi.org/10.7554/eLife.28034

.. _ibo-ref-sanchez2017:

Sánchez, B. J., Zhang, C., Nilsson, A., Lahtvee, P.-J., Kerkhoven, E. J. &
Nielsen, J. (2017). Improving the phenotype predictions of a yeast
genome-scale metabolic model by incorporating enzymatic constraints.
*Mol. Syst. Biol.* **13**, 935. https://doi.org/10.15252/msb.20167411

.. _ibo-ref-scott2010:

Scott, M., Gunderson, C. W., Mateescu, E. M., Zhang, Z. & Hwa, T. (2010).
Interdependence of cell growth and gene expression: origins and
consequences. *Science* **330**, 1099–1102.
https://doi.org/10.1126/science.1192588

.. _ibo-ref-su2021:

Su, Y., Shao, W., Zhang, A. & Zhang, W. (2021). Improving isobutanol
tolerance and titers through EMS mutagenesis in *Saccharomyces cerevisiae*.
*FEMS Yeast Res.* **21**, foab012. https://doi.org/10.1093/femsyr/foab012

.. _ibo-ref-wess2019:

Wess, J., Brinek, M. & Boles, E. (2019). Improving isobutanol production
with the yeast *Saccharomyces cerevisiae* by successively blocking competing
metabolic pathways as well as ethanol and glycerol formation.
*Biotechnol. Biofuels* **12**, 173.
https://doi.org/10.1186/s13068-019-1486-8
