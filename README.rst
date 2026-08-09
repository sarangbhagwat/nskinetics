.. image:: docs/source/_static/images/logo/logo_nskinetics_light_white-circle.png
  :width: 250

===============================================================
The (Non-)Steady state Kinetics simulation package (NSKinetics)
===============================================================

.. image:: http://img.shields.io/pypi/v/nskinetics.svg?style=flat
   :target: https://pypi.python.org/pypi/nskinetics
   :alt: Version_status
.. image:: http://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: https://nskinetics.readthedocs.io/en/latest/
   :alt: Documentation
.. image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
   :alt: license
.. image:: https://img.shields.io/pypi/pyversions/nskinetics.svg
   :target: https://pypi.python.org/pypi/nskinetics
   :alt: Supported_versions
.. image:: https://coveralls.io/repos/github/sarangbhagwat/nskinetics/badge.svg?cachebuster=202507072
   :target: https://coveralls.io/github/sarangbhagwat/nskinetics?branch=main

**Contents**

.. contents::
   :local:

What is NSKinetics?
-------------------

NSKinetics is a fast, flexible, and convenient package in Python for simulating steady- and non-steady-state reaction kinetics — including microbial fermentation and enzyme kinetics — and connecting them to techno-economic analysis (TEA) and life-cycle assessment (LCA) under uncertainty. Kinetic models are declared as SBML — most easily authored as `Antimony <https://tellurium.readthedocs.io/en/latest/antimony.html>`__ text, or imported from an existing SBML file — and wrapped in a ``KineticModel``, which adds unit-aware value access and a Python event API on top of a Tellurium RoadRunner engine that performs the actual ODE integration. ``Event`` covers a single trigger/assignment pair (a parameter switch, a control action); the higher-level ``FeedSpike`` builds on it for fed-batch feeding, topping a species back up to a target concentration whenever it drops below a threshold. The same kinetic model can then drive a `BioSTEAM <https://biosteam.readthedocs.io/en/latest/>`__ process unit through the ``NSKBatchReactor`` bridge, coupling kinetics directly to TEA.

Installation
------------

Get the latest version of NSKinetics from `PyPI <https://pypi.org/project/nskinetics/>`__. If you have an installation of Python with pip, simply install it with:

.. code-block:: bash

    $ pip install nskinetics

To get the git version, run:

.. code-block:: bash

    $ git clone git://github.com/sarangbhagwat/nskinetics

For help on common installation issues, please visit the `documentation <https://nskinetics.readthedocs.io/en/latest/>`__.

Documentation
-------------

NSKinetic's `full documentation <https://nskinetics.readthedocs.io/en/latest/>`__ includes a staged tutorial, starting from a minimal model and building up to a full process/TEA-coupled fed-batch fermentation. The quickstart example below runs that full picture end to end.

One factory call builds a complete, industrially configured process around a real kinetic model: ``nskinetics.processes.create_sugar_prep_and_fermentation_system`` assembles the sugar-solution preparation and fed-batch fermentation section of an actual biorefinery model — a splitter feeding parallel initial-feed and spike-feed conditioning trains (multi-effect evaporator, pumps, dilution-water mixer, heat exchanger), a ``FermentationSaccharomycesEthanolIsobutanol`` fermentor driven by the shipped *S. cerevisiae* kinetic model, and a compressed-air aeration loop. This example runs it with default settings, end to end: build, impose the fed-batch strategy, simulate, inspect the reactor, and finish with a kinetic parameter sweep that carries a strain property all the way through to process economics. Total simulation time is a few seconds once the imports are loaded; every number and figure below is the output of the code shown.

**Step 1: Build the process**

The factory needs a flowsheet, a chemical set (``set_thermo=True`` activates
its own, shipped set), and material for its two inlets: a saccharified sugar
slurry and a seed culture. The seed's ``Yeast`` mass matters — the fermentor
maps it onto the kinetic model's initial cell density, so a yeast-free seed
would mean nothing grows. biosteam's cost correlations warn when a transient
solver state falls outside their validated ranges; those warnings are benign
here and silenced up front, as is a known-benign ``RuntimeWarning`` from
re-running biosteam's design logic:

.. code-block:: python

   import warnings
   import biosteam as bst
   import nskinetics as nsk
   from biosteam.exceptions import CostWarning, DesignWarning

   warnings.filterwarnings('ignore', category=CostWarning)
   warnings.filterwarnings('ignore', category=DesignWarning)
   warnings.filterwarnings('ignore', message=r'.*method added unit results.*',
                           category=RuntimeWarning)
   warnings.filterwarnings('ignore', message=r'.*has been replaced in registry.*',
                           category=RuntimeWarning)

   bst.main_flowsheet.set_flowsheet('quickstart')
   sugar_ferm_sys = nsk.processes.create_sugar_prep_and_fermentation_system(
       set_thermo=True)
   f = bst.main_flowsheet
   S301, F301, F302, V406, K330, V330 = (f.S301, f.F301, f.F302,
                                         f.V406, f.K330, f.V330)

   slurry = S301.ins[0]              # saccharified_slurry
   seed = V406.ins[1]
   slurry.imass['Water'] = 100000.   # kg/hr
   slurry.imass['Glucose'] = 16000.  # kg/hr — 160 g glucose per L water
   seed.imass['Water'] = 500.        # kg/hr
   seed.imass['Yeast'] = 40.         # kg/hr — sets the initial cell density
   sugar_ferm_sys.operating_hours = 330 * 24

   sugar_ferm_sys.diagram()

.. figure:: docs/source/_static/images/examples/tutorial_01_quickstart_flowsheet.png
   :width: 700

   The factory's flowsheet: ``S301`` splits the slurry between the
   initial-feed train (``F301`` → ``M301`` → ``H301``) and the spike-feed
   train (``F302`` → ``M302`` → ``H302``), both feeding the ``V406``
   fermentor; ``K330``/``V330`` supply compressed air for aeration.

**Step 2: Impose the fed-batch strategy**

The factory also builds a ``FedBatchStrategySpecification`` wired to its own
units and attaches it to the fermentor as ``V406.fbs_spec`` (alias
``V406.fed_batch_strategy_specification``) — built but *not* imposed. The
spec holds three concentrations and a maximum batch time: the initial feed
is conditioned to ``target_conc``, a fed-batch spike fires whenever the
reactor's glucose falls to ``threshold_conc``, and each spike draws from a
concentrated ``spike_conc`` feed. Calling ``load_specifications()`` imposes
the strategy: it writes the concentrations into the kinetic model, then tunes
the physical flowsheet to match — solving each train's evaporator vapor
fraction (``V``) or dilution-water ratio until the streams actually reach the
specified concentrations, and setting the splitter so the two trains supply
glucose in the ratio the simulated fermentation consumed it:

.. code-block:: python

   fbs = V406.fbs_spec
   print(fbs)
   print(fbs.current_specifications)

   fbs.load_specifications()
   print('feed conc [g/L]:', fbs.get_feed_conc())
   print('spike conc [g/L]:', fbs.get_spike_conc())
   print('S301 split to initial feed:', S301.split[0])
   print('F301.V:', F301.V, ' F302.V:', F302.V)

.. code-block:: text

   FedBatchStrategySpecification(target_conc=220.0, threshold_conc=210.0, spike_conc=600.0)
   {'target_conc': 220.0, 'threshold_conc': 210.0, 'spike_conc': 600.0, 'tau_max': 72}
   feed conc [g/L]: 217.37708020633303
   spike conc [g/L]: 592.8465856369048
   S301 split to initial feed: 0.5819806503003419
   F301.V: 0.2763304965688685  F302.V: 0.7346545168658999

The 160 g/L slurry cannot meet either specification as-is, so the solve put
real numbers on both trains: the initial-feed evaporator boils off 28% of its
water to reach ~217 g/L (within ~1% of the 220 g/L target) and the spike
train boils off 73% to reach ~593 g/L, with 58% of the slurry routed to the
initial feed. Change a specification (say,
``fbs.load_specifications(target_conc=240.)``) and the same call re-tunes
evaporators, dilution water, and splitter to deliver it.

**Step 3: Simulate and inspect the reactor**

``sugar_ferm_sys.simulate()`` runs the whole flowsheet; inside it, ``V406``
hands its mixed feed to the kinetic model, integrates the fed-batch
fermentation over the full ``tau_max`` window (spikes, aeration switching,
and a spike-count retry included), picks the harvest time ``tau`` by its
``tau_update_policy`` (default: the time ethanol peaks), and converts that
one time point back into its effluent stream. The aeration demand propagates
*backward* — ``V406`` writes it onto the valve's outlet — so the compressor
train is re-simulated afterwards to refresh it:

.. code-block:: python

   sugar_ferm_sys.simulate()
   V330.simulate()   # refresh the aeration loop at the new air demand
   K330.simulate()

   d = V406.nsk_results_specific_tau_dict
   print('tau [h]:', V406.tau)
   print('ethanol titer [g/L]:', d['[s_EtOH]'])
   print('residual glucose [g/L]:', d['[s_glu]'])
   print('cell density [g/L]:', d['[x]'])
   print('fed-batch spikes:', d['curr_n_glu_spikes'])
   print('yield [g EtOH / g glucose]:', d['y_EtOH_glu_added'])

   V406.plot_simulation_results(variables=['[s_glu]', '[s_EtOH]', '[x]'],
                                labels=['glucose', 'ethanol', 'cells'])

.. code-block:: text

   tau [h]: 51.6036036036036
   ethanol titer [g/L]: 138.5123245288004
   residual glucose [g/L]: 0.25470496269632203
   cell density [g/L]: 14.53238208170289
   fed-batch spikes: 9.0
   yield [g EtOH / g glucose]: 0.4629172623016083

.. figure:: docs/source/_static/images/examples/tutorial_01_quickstart_kinetics.png
   :width: 500

   The full kinetic trajectory behind the process result. Early on, 9 spikes
   hold glucose near the 210–220 g/L band while aerobic growth builds cells
   (aeration ends at ~5.3 h, when the cell density hits its stage-1 cutoff);
   glucose is then drawn down to ~0 as ethanol climbs to ~139 g/L, and the
   dashed ``fermentation end`` line marks the selected harvest time
   ``tau`` ≈ 51.6 h.

The same simulation drives a real process unit. The effluent stream carries
the kinetic result mapped onto chemicals, and ``V406.results()`` is
biosteam's standard design-and-cost table — sized and costed from the kinetic
``tau``:

.. code-block:: python

   V406.show(N=100)
   print(V406.results())

.. code-block:: text

   FermentationSaccharomycesEthanolIsobutanol: V406
   ins...
   [0] glucose_initial_feed  from  HXutility-H301
       phase: 'l', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water    2.34e+03
                       Glucose  51.7
   [1] seed
       phase: 'l', T: 298.15 K, P: 101325 Pa
       flow (kmol/hr): Water  27.8
                       Yeast  1.77
   [2] glucose_spike_feed  from  HXutility-H302
       phase: 'l', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water    616
                       Glucose  37.1
   [3] s5  from  IsenthalpicValve-V330
       phase: 'g', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): O2  3.64
                       N2  13.7
   outs...
   [0] fermentation_vent
       phase: 'g', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water       7.25
                       AceticAcid  0.00386
                       CO2         172
                       O2          1.82
                       N2          13.7
   [1] fermentation_effluent
       phase: 'l', T: 305.15 K, P: 101325 Pa
       flow (kmol/hr): Water       2.98e+03
                       Ethanol     163
                       AceticAcid  0.886
                       Glucose     0.0765
                       Yeast       34.8
   FermentationSaccharomycesEthanolIsobutanol         Units      V406
   Electricity         Power                             kW       201
                       Cost                          USD/hr      15.7
   Chilled water       Duty                           kJ/hr -2.75e+07
                       Flow                         kmol/hr  1.85e+04
                       Cost                          USD/hr       138
   Design              Reactor volume                    m3  3.89e+03
                       Batch time                        hr       109
                       Loading time                      hr      54.6
                       Number of reactors                           2
                       Recirculation flow rate        m3/hr      32.1
                       Reactor duty                   kJ/hr -9.18e+06
                       Cleaning and unloading time       hr         3
                       Working volume fraction                    0.9
   Purchase cost       Heat exchangers (x2)             USD  2.92e+04
                       Reactors (x2)                    USD  1.86e+06
                       Agitators (x2)                   USD  1.16e+05
                       Cleaning in place                USD  4.66e+05
                       Recirculation pumps (x2)         USD  5.08e+04
   Total purchase cost                                  USD  2.52e+06
   Utility cost                                      USD/hr       153

**Step 4: From kinetics to process economics — a strain-tolerance sweep**

Ethanol inhibits its own production — the kinetic model applies exponential
ethanol-inhibition terms (parameters ``k_1ie``, ``k_4ie``, ``k_7ie``) to
glucose uptake, acetate formation, and growth. Scaling all three together is
a one-knob model of strain ethanol tolerance: a multiplier below 1 is a more
tolerant (engineered) strain, above 1 a more sensitive one. Because the whole
process responds through the same objects used above, asking "what is
ethanol tolerance worth?" is just a loop — at each multiplier, re-impose the
fed-batch strategy and re-simulate, then read outcomes at every level:
kinetic (titer), strategy (spike count), process (production rate), and
economic (installed cost per annual tonne of ethanol capacity):

.. code-block:: python

   from matplotlib import pyplot as plt

   km = V406.nsk_kinetic_model
   inhibition_params = ('k_1ie', 'k_4ie', 'k_7ie')
   baseline = {p: km.get_value(p) for p in inhibition_params}

   multipliers = [0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 1.75, 2.0]
   records = []
   for m in multipliers:
       for p, v in baseline.items():
           km.set_value(p, m * v)
       V406.fbs_spec.load_specifications()   # re-impose the strategy
       sugar_ferm_sys.simulate()
       V330.simulate()
       K330.simulate()
       d = V406.nsk_results_specific_tau_dict
       ethanol_flow = V406.outs[1].imass['Ethanol']  # kg/hr
       annual_ethanol = ethanol_flow * sugar_ferm_sys.operating_hours / 1000.
       records.append((d['[s_EtOH]'], d['curr_n_glu_spikes'], ethanol_flow,
                       sugar_ferm_sys.installed_equipment_cost / annual_ethanol))

   for p, v in baseline.items():   # restore the shipped model's kinetics
       km.set_value(p, v)

   plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
   plt.rcParams['font.size'] = '11'
   panels = ('Ethanol titer [g/L]', 'Fed-batch spikes [count]',
             'Ethanol production [kg/hr]', 'Installed cost [USD per t/yr]')
   fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.6), sharex=True)
   for i, (ax, label) in enumerate(zip(axes.flat, panels)):
       ax.plot(multipliers, [r[i] for r in records],
               '-o', color='#1f77b4', linewidth=1.5, markersize=5)
       ax.axvline(1.0, linestyle='--', linewidth=1.0, color='gray')
       ax.set_ylabel(label)
       ax.tick_params(axis='both', which='both', direction='inout')
   for ax in axes[1]:
       ax.set_xlabel('Ethanol inhibition strength [x baseline]')
   fig.tight_layout()

.. figure:: docs/source/_static/images/examples/tutorial_01_quickstart_sweep.png
   :width: 650

   One kinetic property, read out at four levels of the model. More tolerant
   strains (left of the dashed baseline) sustain more fed-batch spikes and
   higher titers; more sensitive strains collapse to a single spike and low
   titer. Production and capital intensity follow — mostly through yield and
   batch time — with steps where discrete design decisions (spike counts,
   reactor count, the selected ``tau``) shift.

The whole nine-point sweep simulates in a few seconds. At 0.7× inhibition the
strain sustains 25 spikes and 187.5 g/L ethanol at ~89 USD of installed
equipment per annual tonne of capacity; at 2× it manages one spike, 72.5 g/L,
and ~155 USD per annual tonne — a strain property, priced in process terms.
That is NSKinetics' core loop: kinetics you can edit at one end, equipment
and economics that respond at the other.

See the `full tutorial <https://nskinetics.readthedocs.io/en/latest/tutorial/index.html>`__ for the rest of the workflow — writing and simulating a kinetic model from scratch, the ``Event``/``FeedSpike`` API used inside the fermentor here, a tour of the shipped *S. cerevisiae* kinetic model, and the kinetics-to-biosteam bridge behind ``V406``.


Bug reports
-----------

To report bugs, please use NSKinetics's Bug Tracker at:

    https://github.com/sarangbhagwat/nskinetics

Contributing
------------
For guidelines on how to contribute, visit:

    [link to be added]


License information
-------------------

See ``LICENSE.txt`` for information on the terms & conditions for usage
of this software, and a DISCLAIMER OF ALL WARRANTIES.

Although not required by the NSKinetics license, if it is convenient for you,
please cite NSKinetics if used in your work. Please also consider contributing
any changes you make back, and benefit the community.


About the authors
-----------------

NSKinetics was created and developed by `Sarang S. Bhagwat <https://github.com/sarangbhagwat>`__ as part of the `Scown Group <https://cscown.com/>`__ and the `Energy & Biosciences Institute <https://energybiosciencesinstitute.org/>`__ at the `University of California, Berkeley (UC Berkeley) <https://www.berkeley.edu/>`__. 


References
----------
.. [1] ` To be added <link to be added>`__.


