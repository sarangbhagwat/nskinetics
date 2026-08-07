Quickstart
==========

The fastest way to see NSKinetics work is to load a real, shipped model,
simulate it, and look at two species. NSKinetics ships an SBML model of a
fed-batch *S. cerevisiae* fermentation alongside the package; this snippet
loads it, runs it for 15 simulated hours, and prints glucose and ethanol
concentrations at the start and end:

.. code-block:: python

   import os, nskinetics as nsk
   sbml = os.path.join(os.path.dirname(nsk.__file__),
                       'examples', 's_cerevisiae_ferm_fb_inhib_mod_ibo_sbml.xml')
   trs = nsk.TelluriumReactionSystem.from_sbml(sbml)
   trs._units = {'time': 'h', 'conc': 'g/L'}
   trs.reset()
   res = trs.simulate(0, 15, 16, ['time', '[s_glu]', '[s_EtOH]'])
   print(res[0], res[-1])
   trs.plot_simulation_results(labels=['glucose', 'ethanol'])

Running it in the ``HP_2024`` environment prints:

.. code-block:: text

   [ 0. 15.  0.] [ 1.500e+01 -8.661e-15  4.545e+00]

The first row is ``t=0``: glucose (``s_glu``) starts at 15 g/L, ethanol
(``s_EtOH``) at 0. The second row is ``t=15`` h: glucose has been consumed
down to essentially zero (``-8.661e-15`` is solver-level numerical noise,
not a real negative concentration), and ethanol has risen to about
4.545 g/L.

.. figure:: /_static/images/examples/tutorial_01_quickstart.png
   :width: 400

   Glucose is consumed to ~0 over 15 h while ethanol accumulates to ~4.5 g/L.

That is the whole loop — load a model, simulate, read an array. The next
page, :doc:`02_build_and_simulate`, builds a model from scratch instead of
loading one, and explains each step (units, value access, ``reset()``)
along the way; :doc:`../concepts` covers the ideas in more depth.
