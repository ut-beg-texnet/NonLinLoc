| 

Vel2Grid3D - velocity model description to 3D model grid
========================================================

Given an existing 3D velocity model defined by interpolation between
control point nodes and created by velocity inversion programs such as
SimulPS, Simul2000 and FDTomo to a 3D model grid. (outputs a 3D Grid) ,
**Vel2Grid3D** generates a 3D model Grid header and buffer files
containing velocities, slownesses or other model specification.

 

\ `Overview <#_overview_>`__ - `Running the
program-Input <#_running_>`__ - `Output <#_output_>`__ - `Processing and
Display of results <#_processing_>`__ - `[NonLinLoc
Home] <index.html>`__

 

Overview
--------

The Vel2Grid3D program converts an existing 3D velocity model defined by
interpolation between control point nodes and created by velocity
inversion programs such as SimulPS, Simul2000 and FDTomo into a `3D
Grid <formats.html#_grid_>`__ file containing velocity or slowness
values.

The Vel2Grid3D program uses a "flat earth", rectangular, left-handed,
*x,y,z,t* coordinate system (positive X = East, positive Y = North,
positive Z = down). Distance units are kilometers.

**Important note:**

| Because of there are many variations of the SimulPS, Simul2000 and
  FDTomo type model formats, it is necessary to edit the Vel2Grid3D.c
  source code to select a specific sub-format using C ``#define``
  statements.
| The current options are:

1. default, early SimulPs format with %5f (fixed width) format for
reading velocity values and right-handed, East-negative coordinate
system

2. Zhang %7f (fixed width) format for reading velocity values and
left-handed, East-positive coordinate system

3. ETH format (\*f variable width, space delimited) for reading velocity
values

4. A left-handed, East-positive coordinate system

 

Running the program - Input
---------------------------

Synopsis: \ ``Vel2Grid3D3D``\ ``InputControlFile``\ 

The Vel2Grid3D program takes a single argument \ ``InputControlFile``\ 
which specifies the complete path and filename for an `Input Control
File <control.html>`__ with certain required and optional statements
specifying program parameters and input/output file names and locations.
See the `Vel2Grid3D Statements section <control.html#_Vel2Grid3D_>`__ of
the Input Control File for more details. Note that to run Vel2Grid3D the
`Generic Statements section <control.html#_generic_>`__ of the Input
Control File must contain the \ ``CONTROL``\  and \ ``TRANS``\ 
(Geographic Transformation) statements.

In addition, the Vel2Grid3D program requires a set of `Vel2Grid
Statements <control.html#_Vel2Grid_>`__ and `Vel2Grid3D
Statements <control.html#_Vel2Grid3D_>`__ in the Input Control File:

#. 
#. 

 

Output
------

The velocity or slowness values throughout the requested grid are
written to a new `3D Grid File <formats.html#_grid_>`__. For a
description of the naming convention for these grid files, see the
\ ``VGOUT``\  statement in the Vel2Grid Statements section of the Input
Control File.

 

Processing and Display of results
---------------------------------

The 3D model grids can be post-processed with the program
`Grid2GMT <Grid2GMT.html>`__ to produce a GMT command script for
plotting with the `GMT plotting
package <http://gmt.soest.hawaii.edu/>`__.

 

Back to `the NonLinLoc site Home page <index.html>`__.

*Anthony Lomax*
