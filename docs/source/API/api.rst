API Reference
=============

The current public API of NSKinetics centers on the Tellurium/SBML workflow,
reproducible runs and batch scans (:class:`~nskinetics.Condition`,
:meth:`KineticModel.run <nskinetics.KineticModel.run>`,
:func:`~nskinetics.run_many`, :func:`~nskinetics.scan`), the
:mod:`nskinetics.analysis` subpackage for sensitivity, identifiability,
screening and calibration, the biosteam process bridge, and ready-made
process section factories. The
last page documents the isobutanol biorefinery that ships as a case study
(``nskinetics.biorefineries.isobutanol``; optional dependencies).

.. toctree::
   :maxdepth: 2

   kinetic_model
   events
   units
   processes
   flux_analysis
   analysis
   models
   exceptions
   biorefineries
