(Non-)Steady State Kinetics Simulation
======================================

.. currentmodule:: nskinetics

.. toctree::
   :maxdepth: 2
   :hidden:

   tutorial/index
   concepts
   API/api

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
        :class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol` bridge, coupling kinetics directly to TEA —
        and ready-made flowsheet sections built around these units ship as system factories in
        ``nskinetics.processes`` (see the :doc:`quickstart <tutorial/01_quickstart>`).

Quickstart
----------

Write a tiny kinetic model in Antimony, wrap it in a :class:`KineticModel`, and simulate it:

.. code-block:: python

   import tellurium as te
   import nskinetics as nsk

   model = """
   model demo()
     compartment env; species S in env, P in env;
     S = 10; P = 0; env = 1; k = 0.3;
     J: S => P; k*S*env;
   end
   """
   r = te.loadAntimonyModel(model)
   km = nsk.KineticModel(r, units={'time': 'h', 'conc': 'g/L'})
   km.validate_units()
   km.reset()

   result = km.simulate(0, 10, 101, ['time', '[S]', '[P]'])
   print('t=10:', result[-1])
   km.plot_simulation_results(labels=['S', 'P'])

This prints ``t=10: [10.     0.498  9.502]`` — species ``S`` has decayed from an initial concentration of 10 g/L
to about 0.498 g/L by ``t=10`` h under first-order decay (``k=0.3``), while ``P`` has risen to about 9.502 g/L.

.. figure:: /_static/images/examples/index_landing.png
   :width: 400

   ``km.simulate`` records the trajectory; ``plot_simulation_results`` plots it.

See :doc:`tutorial/index` for the full walkthrough, including events, fed-batch feeding, and the biosteam/TEA
bridge.


.. grid:: 1 2 3 4


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

