(Non-)Steady State Kinetics Simulation
======================================

.. currentmodule:: nskinetics

.. toctree::
   :maxdepth: 2
   :hidden:

   tutorial/index
   concepts
   API/api
   contributing/index

.. grid:: 1 1 2 2

    .. grid-item::
    
        .. image:: _static/images/demo/loop_dark.png
           :class: only-dark
           :align: center
           :height: 200

        .. image:: _static/images/demo/loop_light.png
           :class: only-light
           :align: center
           :height: 200

    .. grid-item::

        NSKinetics is a fast, flexible, and convenient package in Python for simulating steady- and non-steady-state
        reaction kinetics — including microbial fermentation and enzyme kinetics — and connecting them to
        techno-economic analysis (TEA) and life-cycle assessment (LCA) under uncertainty. Kinetic models are
        declared as SBML — most easily authored as `Antimony <https://tellurium.readthedocs.io/en/latest/antimony.html>`__
        text, or imported from an existing SBML file — and wrapped in a :class:`KineticModel`, which
        adds unit-aware value access and a Python event API (:class:`Event`, and the higher-level
        :class:`FeedSpike` for fed-batch feeding) on top of a Tellurium RoadRunner engine that performs the actual
        ODE integration. The same kinetic model can then drive a biosteam process unit through the
        :class:`~nskinetics.units.NSKBatchReactor` bridge, coupling kinetics directly to TEA —
        and ready-made flowsheet sections built around these units ship as system factories in
        ``nskinetics.processes`` — see the :doc:`quickstart tutorial <tutorial/01_quickstart>` and the
        :obj:`~nskinetics.processes.create_sugar_prep_and_fermentation_system` factory.

Quickstart
----------

One factory call builds a complete, industrially configured process around a
real kinetic model:
:obj:`~nskinetics.processes.create_sugar_prep_and_fermentation_system`
assembles the sugar-solution preparation and fed-batch fermentation section of
an actual biorefinery model — a splitter feeding parallel initial-feed and
spike-feed conditioning trains (multi-effect evaporator, pumps, dilution-water
mixer, heat exchanger), a
:class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol` fermentor
driven by the shipped *S. cerevisiae* kinetic model, and a compressed-air
aeration loop. The interactive demo below runs it end to end: build, simulate,
change the fed-batch strategy, and inspect the reactor.

.. raw:: html

   <div style="margin:1.5rem auto;max-width:1120px;">
     <iframe src="_static/quickstart_demo.html"
             title="NSKinetics quickstart — interactive demo"
             loading="lazy"
             style="width:100%;height:840px;border:1px solid rgba(128,128,128,0.25);border-radius:14px;display:block;"></iframe>
     <p style="text-align:center;font-size:0.85rem;margin-top:0.4rem;opacity:0.8;">
       Interactive quickstart demo &mdash;
       <a href="_static/quickstart_demo.html" target="_blank" rel="noopener">open in a new tab</a>.
     </p>
   </div>

See the full :doc:`tutorial <tutorial/index>` for the rest of the workflow —
writing and simulating a kinetic model from scratch, the
:class:`Event`/:class:`FeedSpike` API used inside the fermentor here, a tour
of the shipped *S. cerevisiae* kinetic model, and the kinetics-to-biosteam
bridge behind ``V406``.


.. grid:: 1 2 3 4
    :class-row: sd-align-major-center


    .. grid-item-card:: Getting Started
       :text-align: center
       :link: tutorial/index
       :link-type: doc
       :padding: 1

       .. image:: _static/images/icons/getting-started_dark.png
          :height: 100
          :class: only-dark
          :align: center

       .. image:: _static/images/icons/getting-started_light.png
          :height: 100
          :class: only-light
          :align: center

       Tutorials on NSKinetics


    .. grid-item-card:: Key Concepts
       :text-align: center
       :link: concepts
       :link-type: doc
       :padding: 1

       .. image:: _static/images/icons/concepts_dark.png
          :height: 100
          :class: only-dark
          :align: center

       .. image:: _static/images/icons/concepts_light.png
          :height: 100
          :class: only-light
          :align: center

       The ideas behind NSKinetics


    .. grid-item-card:: API Reference
       :text-align: center
       :link: API/api
       :link-type: doc
       :padding: 1

       .. image:: _static/images/icons/api_dark.png
          :height: 100
          :class: only-dark
          :align: center

       .. image:: _static/images/icons/api_light.png
          :height: 100
          :class: only-light
          :align: center

       Detailed documentation


Installation
------------

Get the latest version of NSKinetics from `PyPI <https://pypi.org/project/nskinetics/>`__. If you have an installation of Python with pip, simply install it with:

.. code-block:: bash
    
    $ pip install nskinetics

To get the git version, use:

.. code-block:: bash

    $ git clone git://github.com/sarangbhagwat/nskinetics

Or download directly from the `GitHub page <https://github.com/sarangbhagwat/nskinetics>`__.

Common Issues
-------------

* **Cannot install/update NSKinetics:**

  If you are having trouble installing or updating NSKinetics, it may be due to dependency issues. You can bypass these using:
  
  .. code-block:: bash

     $ pip install --user --ignore-installed nskinetics

  You can make sure you install the right version by including the version number:

  .. code-block:: bash

     $ pip install nskinetics==<version>

  E.g., for version 0.1.4:

  .. code-block:: bash

     $ pip install nskinetics==0.1.4

