NSKinetics: Non-Steady State Kinetics Simulation
================================================

.. toctree::
   :maxdepth: 2
   :hidden:
   
   API/api
   user_guide

**NSKinetics** is a fast, flexible, and convenient package to simulate non-steady state reaction kinetics, especially for systems involving enzymatic conversions. Models for multiple inhibitory phenomena (competitive, non-competitive, uncompetitive, and "mechanism-based") are also included. NSKinetics enables the construction, simulation, and analysis of reaction systems governed by mass action kinetics or other user-defined rate laws. It supports features such as species concentration spikes, event triggers, inverse modeling (parameter fitting to experimental data), and design of experiments.

.. note::

   Autosummary-generated `API <https://nskinetics.readthedocs.io/en/latest/API/api.html>`__ is available. Full documentation is expected soon.

Installation
------------

Get the latest version of NSKinetics from `PyPI <https://pypi.org/project/nskinetics/>`__. If you have an installation of Python with pip, simply install it with:

.. code-block:: bash
    
    $ pip install nskinetics

To get the git version, run:

.. code-block:: bash

    $ git clone git://github.com/sarangbhagwat/nskinetics


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

.. toctree::
   :maxdepth: 2

