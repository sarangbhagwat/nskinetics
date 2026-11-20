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
      :img-top: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig9.png
      :img-alt: Normalized concentrations of 3-HP, malonic acid and by-products against time with H2O2 and with O2 as the oxidant, reproduced with the shipped models: solid model lines on the black dotted deposited curves, circles for the experimental points, and the reaction list next to each panel.
      :class-card: ibo-case-card
      :shadow: md

      Two kinetic models of the oxidation of 3-hydroxypropionic acid to
      malonic acid over Pd/C with O2 or H2O2 (Kim et al., 2026) that redraw
      the model curves of the paper's Fig. 9, with a record of where the
      paper, its data deposit and the models disagree. **The published
      50 °C model, extended to other temperatures as the paper describes,
      does not reproduce the experiments at 25, 75 and 100 °C; the model
      curves of the paper's Fig. 10a–c and the heat maps of its Fig. 11 are
      not produced by these models and are not reproduced here.** The
      points and lines of its Arrhenius plot (Fig. 10d) are compared with
      lines of the O2 model's activation energies.

      +++
      ``nskinetics.models.pd_c_3hp_oxidation`` · 2 models, deposited and
      digitized reference data, 4 figures

.. toctree::
   :maxdepth: 2
   :hidden:

   isobutanol/index
   pd_c_3hp_oxidation
