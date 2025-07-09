Non-Steady State Kinetics Simulation
====================================

.. toctree::
   :maxdepth: 2
   :hidden:
   
   API/api
   user_guide

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

        **NSKinetics** is a fast, flexible, `open-source <https://github.com/sarangbhagwat/nskinetics>`__ platform in Python to simulate and evaluate
        non-steady state reaction kinetics, especially for systems involving enzymatic conversions 
        with inhibitory phenomena. NSKinetics enables the construction, simulation, and analysis 
        of reaction systems governed by mass action kinetics or other user-defined rate laws, 
        as well as the optimal design of experiments for parameter identifiability.


.. grid:: 1 2 3 4
   
   
    .. grid-item-card:: Getting Started
       :text-align: center
       :link: https://nskinetics.readthedocs.io/en/latest/user_guide/Getting_started.html
       :link-type: url
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


    .. grid-item-card:: Inverse Modeling
       :text-align: center
       :link: https://nskinetics.readthedocs.io/en/latest/API/inverse_modeling.html
       :link-type: url
       :padding: 1
       
       .. image:: _static/images/icons/inverse-modeling_dark.png
          :height: 100
          :class: only-dark
          :align: center
          
       .. image:: _static/images/icons/inverse-modeling_light.png
          :height: 100
          :class: only-light
          :align: center
       
       Fit parameters to data

       
    .. grid-item-card:: Design of Experiments
       :text-align: center
       :link: https://nskinetics.readthedocs.io/en/latest/API/doe.html
       :link-type: url
       :padding: 1
       
       .. image:: _static/images/icons/doe_dark.png
          :height: 100
          :class: only-dark
          :align: center
          
       .. image:: _static/images/icons/doe_light.png
          :height: 100
          :class: only-light
          :align: center
    
       Computationally driven


    .. grid-item-card:: API Reference
       :text-align: center
       :link: https://nskinetics.readthedocs.io/en/latest/API/api.html
       :link-type: url
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

Or download directly from the `GitHub page <https://github.com/sarangbhagwat/nskinetics>`.

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

