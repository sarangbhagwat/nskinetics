Demonstrations
==============

The tutorials build kinetic models and run them inside a single process unit.
The demonstrations go one step further: they plug an NSKinetics model into a
complete process model with a techno-economic analysis (TEA), so a change to
an enzyme capacity, a product tolerance or a feeding setpoint is priced all
the way to the plant's selling price and rate of return. Each demonstration
ships with NSKinetics as a subpackage of ``nskinetics.biorefineries``, with
its flowsheet, scenarios and the scripts behind every analysis shown.

A second kind of demonstration reproduces a published kinetic study. Its
models ship as a subpackage of ``nskinetics.models``, with the reference
data and the functions that redraw the study's figures, and its page states
what agrees with the publication and what does not.

.. grid:: 1
   :gutter: 3

   .. grid-item-card:: Isobutanol biorefinery
      :link: isobutanol/index
      :link-type: doc
      :img-top: /_static/images/isobutanol/overview/kinetic_model_conceptual_diagram.png
      :img-alt: Reaction network of the S. cerevisiae fed-batch model with the engineered isobutanol pathway, fed-batch glucose feeding and two-stage aeration.
      :class-card: ibo-case-card
      :shadow: md

      A corn dry-grind biorefinery whose fermentor runs the shipped
      *S. cerevisiae* fed-batch model with an engineered isobutanol
      pathway. **Strains that maximize titer, rate or yield miss the
      profitable design; the best strain found co-produces isobutanol and
      ethanol and more than doubles the plant's internal rate of return
      (12.6 % → 27.3 %) at the same capital.** That strain is a
      prospective, bounds-limited in-silico target, not a near-term one: its
      39.6 g/L isobutanol is 15–20× the best reported *S. cerevisiae* titers,
      and the result rests on one isobutanol price and two random seeds
      (:ref:`scope <ibo-opt-scope>`).

      +++
      ``nskinetics.biorefineries.isobutanol`` · flowsheet, TEA, 10 scenarios,
      8 analyses · IRRs: campaign record, build 2b5b27d

   .. grid-item-card:: Pd/C oxidation of 3-HP to malonic acid
      :link: pd_c_3hp_oxidation
      :link-type: doc
      :img-top: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11.png
      :img-alt: Heat maps of malonic acid selectivity and yield over reaction temperature and time, reproduced with the shipped models, with contour lines at 0.1, 0.3 and 0.5 and a ridge of maxima running from short times at high temperature to long times at low temperature.
      :class-card: ibo-case-card
      :shadow: md

      Four kinetic models of the oxidation of 3-hydroxypropionic acid to
      malonic acid over Pd/C with O2 or H2O2 (Kim et al., 2026) that redraw
      the model curves of the paper's Figs. 9, 10 and 11, with a record of
      where the paper, its data deposit and the models disagree. **The
      published 50 °C model, extended to other temperatures as the paper
      describes, reproduces neither the published heat maps nor the
      experiments at 25, 75 and 100 °C; the paper's Figs. 10 and 11 were
      made with two other models that are not in the paper.**

      +++
      ``nskinetics.models.pd_c_3hp_oxidation`` · 4 models, deposited and
      digitized reference data, 7 figures

.. toctree::
   :maxdepth: 2
   :hidden:

   isobutanol/index
   pd_c_3hp_oxidation
