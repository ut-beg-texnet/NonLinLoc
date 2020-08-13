| 

Loc2ssst - generate SSST corrections on a  to 3D grid for NLLoc locations
=========================================================================

Given a set of NLLoc locations, **Loc2ssst** generates station and phase specific, SSST corrections over a specified 3D grid 
and adds these corrections to the input travel-times used for NLLoc location to generate new, SSST travel-time grids.


Overview
--------

For each station and phase in the input NLLoc locations (see LSLOCFILES), 
the Loc2ssst program accumulates over a specified SSST 3D grid (see LSGRID)
a weighted sum of the residuals for the station and phase for each input NLLoc location.  Each residual is weighted based on the 
distance `D` between the grid cell and the coresponding hypocenter using the relations:

| `weight = exp(-(D^2 / char_dist^2)) + weight_floor`

where `CharDist` and `WeightFloor` are specified in LSPARAMS.


Running the program - Input
---------------------------

Synopsis: ``Loc2ssst InputControlFile``

The Loc2ssst program takes a single argument ``InputControlFile`` which
specifies the complete path and filename for an `Input Control
File` with certain required and optional statements specifying program parameters and
input/output file names and locations. See the `Loc2ssst Statements
section of the Input Control File for more details. Note that to run Loc2ssst
the `Generic Statements
section of the Input Control File must contain the ``CONTROL`` and ``TRANS``
(Geographic Transformation) statements.

In addition, the Loc2ssst program requires a set of `Loc2ssst
Statements in the Input Control File


Output
------

The calculated SSST values throughout the requested ``LSGRID`` grid are
written to a new `3D Grid`. 
The updated travel-time values throughout the requested ``LSOUTGRID`` grid are
written to a new `3D Grid`, these files can be used as travel-time files for subsequent NLLoc location. 
For a descrition of the naming convention for these grid files, see the
```LSOUT`` statement in the Loc2ssst Statements section of the Input Control File.


Processing and Display of results
---------------------------------

The 3D SSST and updated travel-time grids can be post-processed with the program
`Grid2GMT` to
produce a GMT command script for plotting with the `GMT plotting
package <http://gmt.soest.hawaii.edu/>`__.
