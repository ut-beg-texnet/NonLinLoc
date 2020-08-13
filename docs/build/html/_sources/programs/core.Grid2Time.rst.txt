| 

Grid2Time - 3D model grid to travel-time and angles grids
=========================================================

Given a velocity model grid, **Grid2Time** calculates the travel-time
from a source point in a 3D grid to all other points in the grid.
Optionally, the program also estimates the take-off angles for rays from
each point in the grid to the source.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

`Overview <#_overview_>`__ - `Podvin and Lecomte
Algorithm <#_pl_algorithm_>`__ - `Take-Off Angles
Algorithm <#_angles_algorithm_>`__ - `Running the
program-Input <#_running_>`__ - `Output <#_output_>`__ - `Processing and
Display of results <#_processing_>`__ - `[NonLinLoc
Home] <./index.html>`__

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Overview
--------

The Grid2Time program calculates the travel-times between a station and
all nodes of an *x,y,z* spatial grid using the Eikonal finite-difference
scheme of `Podvin and Lecomte (1991) <./references.html>`__. The results
are stored on disk in a **travel-time** `3D
Grid <./formats.html#_grid_>`__ files.

Optionally, the Grid2Time program also estimates the take-off angles for
rays from each point in the grid to the source by examining the
gradients of the travel-time field. These results are stored in an
**angles** grid.

The 3D travel-time computation and the size of the output time-grid
files grow rapidly with grid dimension. However, for location in
horizontally layered models the travel-times can be stored on compact 2D
grids. A layered model / 2D grid can also be used for"regional" stations
far from the local search volume in combination with 3D models and 3D
grids for stations within the search volume. This option may introduce
some error if strong heterogeneity in the local 3D velocity structure
intersects the (usually downgoing) ray paths to the regional stations.

The Grid2Time program uses a "flat earth", rectangular, left-handed,
*x,y,z,t* coordinate system (positive X = East, positive Y = North,
positive Z = down). Distance units are kilometers, and many input/output
distance quantities can be expressed in rectangular or geographic
(latitude and logitude) coordinates.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Podvin and Lecomte, Eikonal, Finite-difference Algorithm
--------------------------------------------------------

The travel times between a station and all nodes of a 3D grid are
calculated using the Eikonal finite-difference scheme of `Podvin and
Lecomte (1991) <./references.html>`__. The algorithm is implemented in
the Grid2Time program using a C function time\_3d() due to P. Podvin
(last revision 2 January 1992). The abstract of `Podvin and Lecomte
(1991) <./references.html>`__ describes the algorithm as:

    *This method relies on a systematic application of Huygen's
    principle in the finite difference approximation. Such approximation
    explicitly takes into account the existence of different propagation
    modes (transmitted and diffracted body waves, head waves). Local
    discontinuities of the time gradient in the first arrival time field
    (e.g. caustics) are built as intersections of locally independent
    wavefronts. As a consequence, the proposed method provides accurate
    first travel times in the presence of extremely severe, arbitrarily
    shaped velocity contrasts.*

    *Associated to a simple procedure which accurately traces rays in
    the obtained time field, this method provides a very fast tool for a
    large spectrum of seismic and seismological problems.*

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Take-Off Angles Algorithm
-------------------------

The take-off angles at a node for the first-arrival ray to the source
are estimated from the gradients of travel-time at the node. Two
gradients are estimated for each axis direction *x, y, z* - one
*G\ :sub:`low`* between the node and its preceeding neighbour along the
axis, and a second *G\ :sub:`high`* between the following neighbor and
the node. The total gradient *G\ :sub:`axis`* along an axis is the mean
of these two gradients; the total gradient along the three axes
determine the vector gradient of travel-time. The ray take-off angles
*R\ :sub:`dip`* (dip, range of 0 (down) to 180 deg (up)) and
*R\ :sub:`az`* (azimuth, range of 0 to 360 deg CW from North) specify
the direction opposite to the vector gradient of travel-time.

A crude quality factor *Q\ :sub:`axis`* between 0 and 10 is determined
from the ratio

    *Q\ :sub:`axis`* = (20 *G\ :sub:`low`* *G\ :sub:`high`*) /
    (*G\ :sub:`low`*\ :sup:`2` + *G\ :sub:`high`*\ :sup:`2`)

If *Q\ :sub:`axis`* < 0 (i.e. the two gradients have opposite sign),
*Q\ :sub:`axis`* is set equal to 0. If *Q\ :sub:`axis`* = 10 then the
two gradients have the same magnitude and sign. A final quality for the
take-off angles is determined from the weighted average of the qualities
along each axis, where the weighting is given by the magnitude of the
mean gradient along each axis,

    *Q* = (\|*G\ :sub:`x`*\ \| *Q\ :sub:`x`* + \|\ *G\ :sub:`y`*\ \|
    *Q\ :sub:`y`* + \|\ *G\ :sub:`z`*\ \| *Q\ :sub:`z`*) /
    (\|*G\ :sub:`x`*\ \| + \|\ *G\ :sub:`y`*\ \| +
    \|\ *G\ :sub:`z`*\ \|).

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Running the program - Input
---------------------------

Synopsis: ``Grid2Time InputControlFile``

The Grid2Time program takes a single argument ``InputControlFile`` which
specifies the complete path and filename for an `Input Control
File <./control.html>`__ with certain required and optional statements
specifying program parameters and input/output file names and locations.
See the `Grid2Time Statements section <./control.html#_Grid2Time_>`__ of
the Input Control File for more details. Note that to run Grid2Time the
`Generic Statements section <./control.html#_generic_>`__ of the Input
Control File must contain the ``CONTROL`` and ``TRANS`` (Geographic
Transformation) statements.

In addition, the Grid2Time program requires:

#. A 2D or a 3D velocity model `3D Grid <./formats.html#_grid_>`__ file
   created using `Vel2Grid <./Vel2Grid.html>`__ or other software. One
   velocity model grid is required for each wave type (i.e. P or S).
   Note that a 3D Grid file may specify a 2D model.

The names, locations and other information for these files is specified
in the `Grid2Time Statements section <./control.html#_Grid2Time_>`__ of
the Input Control File.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Output
------

The travel-times and take-off angles throughout a grid are written to a
separate `3D Grid File <formats.html#_grid_>`__ for each phase at each
station. For a descrition of the naming convention for these grid files,
see the ```GTFILES`` <./control.html#_Grid2Time_gtfiles_>`__ statement
in the Grid2Time Statements section of the Input Control File.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Processing and Display of results
---------------------------------

The travel-time and angles grid results for a single source can be
post-processed with the program `Grid2GMT <./Grid2GMT.html>`__ to
produce a GMT command script for plotting with the `GMT plotting
package <http://gmt.soest.hawaii.edu/>`__.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Back to `the NonLinLoc site Home page <./index.html>`__.

*Anthony Lomax*
