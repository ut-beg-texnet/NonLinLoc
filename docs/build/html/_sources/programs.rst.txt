NonLinLoc Programs
==================

NonLinLoc contains core and utilitiy programs. Core contains the
main programs for setting up velocity models, calculating travel-times and locating.
Utilitiess contains a range of programs for format conversion, grid manipulation, etc.

All NonLinLoc core programs and many utiltiy programs use a control file to specify program run parameters.


Core
----

Core routines of the NonLinLoc project.

.. toctree::
   :maxdepth: 1

   programs/core.Vel2Grid
   programs/core.Vel2Grid3D
   programs/core.Grid2Time
   programs/core.Time2EQ
   programs/core.NLLoc
   programs/core.LocSum
   programs/core.Grid2GMT
   programs/core.Loc2ssst


Utils
-----

Various utility functions.

.. toctree::
   :maxdepth: 1

   programs/utils.oct2grid
   

Control
-------

Core routines of the NonLinLoc project.

.. toctree::
   :maxdepth: 1

   programs/control.ControlFile
   