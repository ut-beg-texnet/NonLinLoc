| 

Vel2Grid - velocity model description to 3D model grid
======================================================

Given a velocity model description, **Vel2Grid** generates a 3D model
Grid header and buffer files containing velocities, slownesses or other
model specification.


Overview
--------

The Vel2Grid program converts analytic or other velocity model
specifications into a `3D Grid <./formats.html#_grid_>`__ file
containing velocity or slowness values.

The Vel2Grid program uses a "flat earth", rectangular, left-handed,
*x,y,z,t* coordinate system (positive X = East, positive Y = North,
positive Z = down). Distance units are kilometers.


Running the program - Input
---------------------------

Synopsis: ``Vel2Grid InputControlFile``

The Vel2Grid program takes a single argument ``InputControlFile`` which
specifies the complete path and filename for an `Input Control
File <./control.html>`__ with certain required and optional statements
specifying program parameters and input/output file names and locations.
See the `Vel2Grid Statements section <./control.html#_Vel2Grid_>`__ of
the Input Control File for more details. Note that to run Vel2Grid the
`Generic Statements section <./control.html#_generic_>`__ of the Input
Control File must contain the ``CONTROL`` and ``TRANS`` (Geographic
Transformation) statements.

In addition, the Vel2Grid program requires a set of `Vel2Grid
Statements <./control.html#_Vel2Grid_>`__ in the Input Control File that
specify a layered model or a 3D velocity model. The velocity model can
be specified in the control file by:

#. A set of ```LAYER`` <./control.html#_Vel2Grid_layer_>`__ statements
   defining a horizontally layered model with constant or
   constant-gradient velocity and density in each layer.
#. A set of ```VERTEX`` <./control.html#_Vel2Grid_vertex_>`__,
   ```EDGE`` <./control.html#_Vel2Grid_edge_>`__, and
   ```POLYGON2`` <./control.html#_Vel2Grid_polygon2_>`__ statements
   defining a 2D polygon model and a
   ```2DTO3DTRANS`` <./control.html#_Vel2Grid_2d3dtrans_>`__ statement
   to convert this 2D model into a 3D model. Optionally, there may be a
   set of ```LAYER`` <./control.html#_Vel2Grid_layer_>`__ statements
   defining a horizontally layered background model. This background
   model must be defined if the transformed 2D polygon model does not
   completeley fill the requested 3D grid.


Output
------

The velocity or slowness values throughout the requested grid are
written to a new `3D Grid File <formats.html#_grid_>`__. For a
descrition of the naming convention for these grid files, see the
```VGOUT`` <./control.html#_Vel2Grid_vgout_>`__ statement in the
Vel2Grid Statements section of the Input Control File.


Processing and Display of results
---------------------------------

The 3D model grids can be post-processed with the program
`Grid2GMT <./Grid2GMT.html>`__ to produce a GMT command script for
plotting with the `GMT plotting
package <http://gmt.soest.hawaii.edu>`__.

