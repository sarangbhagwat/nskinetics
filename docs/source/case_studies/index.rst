Case Studies
============

The tutorials build kinetic models and run them inside a single process unit.
The case studies go one step further: they plug an NSKinetics model into a
complete process model with a techno-economic analysis (TEA), so a change to
an enzyme capacity, a product tolerance or a feeding setpoint is priced all
the way to the plant's selling price and rate of return. Each case study
ships with NSKinetics as a subpackage of ``nskinetics.biorefineries``, with
its flowsheet, scenarios and the scripts behind every analysis shown.

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

.. toctree::
   :maxdepth: 2
   :hidden:

   isobutanol/index
