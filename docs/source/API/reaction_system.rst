TelluriumReactionSystem
=======================

.. currentmodule:: nskinetics

:class:`TelluriumReactionSystem` is a thin wrapper around a Tellurium
extended RoadRunner model. It adds unit-aware value access and a Python
event API. It has **no** ``solve``/``simulate`` method of its own — you run
the model through the underlying RoadRunner (``trs._te.simulate(...)``) or
through the :doc:`process bridge <units>`.

.. autoclass:: TelluriumReactionSystem
   :members:
   :show-inheritance:
