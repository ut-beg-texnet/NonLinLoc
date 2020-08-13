| 

Time2EQ - travel-time grid to synthetic observations
====================================================

Given a hypocenter location and a set of travel-time grids, Time2EQ
calculates predicted travel-times.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

`Overview <#_overview_>`__ - `Running the program-Input <#_running_>`__
- `Output <#_output_>`__ - `Processing and Display of
results <#_processing_>`__ - `[NonLinLoc Home] <./index.html>`__

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Overview
--------

The Time2EQ program calculates predicted travel-times between one or
more synthetic events and one or more stations. Predicted take-off
angles at the source are also calculated if an event mechanism is given
and the corresponidng take-off angles grids are available.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Running the program - Input
---------------------------

Synopsis: ``Time2EQ InputControlFile``

The Time2EQ program takes a single argument ``InputControlFile`` which
specifies the complete path and filename for an `Input Control
File <./control.html>`__ with certain required and optional statements
specifying program parameters and input/output file names and locations.
See the `Time2EQ Statements section <./control.html#_Time2EQ_>`__ of the
Input Control File for more details. Note that to run Time2EQ the
`Generic Statements section <./control.html#_generic_>`__ of the Input
Control File must contain the ``CONTROL`` and ``TRANS`` (Geographic
Transformation) statements.

In addition, the Time2EQ program requires:

#. Files containing a 2D or a 3D **Travel-time grids** (and optionally
   **Angles grids**) created by the program
   `Grid2Time <./Grid2Time.html>`__ for each phase type at each station.

The names, locations and other information for these files is specified
in the `Time2EQ Statements section <./control.html#_Time2EQ_>`__ of the
Input Control File.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Output
------

The predicted travel-times are written to an observation file in
`NonLinLoc phase file format <formats.html#_phase_nlloc_>`__.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Processing and Display of results
---------------------------------

The predicted travel-time files can be used as input phase/observation
files for location with the program `NLLoc <./NLLoc.html>`__.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Back to `the NonLinLoc site Home page <./index.html>`__.

*Anthony Lomax*
