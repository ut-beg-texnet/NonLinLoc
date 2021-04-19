NonLinLoc Software Package Control File
=======================================

:ref:`Overview` -
:ref:`Definitions` -
:ref:`Generic Control Statements` -
:ref:`Vel2Grid Program` -
:ref:`Vel2Grid3D Program` -
:ref:`Grid2Time Program` -
:ref:`Time2EQ Program` -
:ref:`NLLoc Program` -
:ref:`Loc2ssst Program`

Overview
--------

| The various NonLinLoc programs all use the same control file syntax
  and share some "generic" control statements. The control statements
  for all the NonLinLoc programs for a project (a study with common
  "generic" control statements) may be combined into one file without
  conflicts. The basic control file statement syntax consists of a
  control *keyword* followed by one or more *parameters* on a single
  line (except when a newline is explicitely required).
 `KEYWORD parameter1 parameter2 ...`` The keyword must begin in the
  first column and is always in upper case. Keywords and parameters must
  be separated by one or more spaces or tabs. A required newline in a
  parameter list is indicated by ``[newline]``. Blank lines and lines
  with ``#`` in the first column are ignored. Use ``#`` in the first
  column for comments and to "comment out" a statement.



Definitions
-----------

| **Statement Priority**
| *required* - must be present in control file to run the coressponding
  program
| *optional* - optional in control file
| *repeatable* - may be present multiple times in control file
| *ignored* - may be present but is not used by program under certain
  conditions

| **Datatypes**
| *integer* decimal integer ( *i.e.* ``0, 5, 285``)
| *float* - decimal floating point number (*i.e.*
 `1.0, 3.68, -4.5, 5.4e6``)
| *chars* - sequence of characters without spaces (*i.e.*
 `NO_SAVE, abcdef, /data/bigevent.dat``)
| *string* - sequence of characters which is read until the end of line,
  spaces are allowed (*i.e.*
 `The biggest earthquake sequence in history``)
| *choice* - selection from a fixed list of items (*i.e.*
 `SAVE NO_SAVE``)

| **Miscellaneous**
| default: - indicates default values for parameters when an *optional*
  control statement is not present in the control file
| min: max: - indicates minimum and maximum allowed values for
  parameters
| VERY\_LARGE\_DOUBLE - a large floating point value (typically 1.0e+30)

Generic Control Statements
--------------------------

| The generic control file statements may be used by one or more of the
  programs in the NonLinLoc package. The statements ``TRANS`` and
 `MAPLINE`` and their parameters must be the same for all programs
  runs for a g*n location project.

```INCLUDE`` <#_generic_include_>`__ -
```CONTROL`` <#_generic_control_>`__ - ```TRANS`` <#_generic_trans_>`__
- ```MAPLINE`` <#_generic_mapline_>`__ - ``` `` <#_generic_maptrans_>`__
- ```MAPGRID`` <#_generic_mapgrid_>`__

| **INCLUDE - Include**
| *optional*, *repeatable*
| Syntax 1: ``INCLUDE`` ``includeFile``
| Inserts text from another file at current positon in control file.
|    ``includeFile`` (*string*) path and name of file to include
| Notes:
|    1. This statement is implemented only for NonLinLoc programs
  Vel2Grid, Grid2Time, Time2EQ and NLLoc.
|    2. The included text must contain only valid NonLinLoc control
  statements, blank lines or comment lines, but may not have another
 `INCLUDE`` statement.

| **CONTROL - Control**
| *required*, *non-repeatable*
| Syntax 1: ``CONTROL`` ``messageFlag randomNumberSeed``
| Sets various general program control parameters.
|    ``messageFlag`` (*integer*, min:\ ``-1``, default:\ ``1``) sets the
  verbosity level for messages printed to the terminal ( ``-1`` =
  completely silent, ``0`` = error messages only, ``1`` = ``0`` +
  higher-level warning and progress messages, ``2`` and higher = ``1`` +
  lower-level warning and progress messages + information messages, ...)
|    ``randomNumberSeed`` (*integer*) integer seed value for generating
  random number sequences (used by program NLLoc to generate Metropolis
  samples and by program Time2EQ to generate noisy time picks)

| **TRANS - Geographic Transformation**
| *required*, *non-repeatable*
| Syntax 1: ``TRANS`` ``GLOBAL``
| Syntax 2: ``TRANS`` ``SIMPLE latOrig longOrig rotAngle``
| Syntax 3: ``TRANS`` ``NONE``
| Syntax 4: ``TRANS`` ``SDC latOrig longOrig rotAngle``
| Syntax 5: ``TRANS`` ``LAMBERT refEllipsoid latOrig longOrig firstStdParal secondStdParal rotAngle``
| Syntax 6: ``TRANS`` ``TRANS_MERC refEllipsoid latOrig longOrig rotAngle``
| Syntax 7: ``TRANS`` ``AZIMUTHAL_EQUIDIST refEllipsoid latOrig longOrig rotAngle``
| Sets geographic to working coordinates transformation parameters. The
 `GLOBAL`` option sets spherical regional/teleseismic mode, with no
  geographic transformation - most position values are input and used
  directly as latitude and longitude in degrees. The ``SIMPLE``,
 `SDC``, ``LAMBERT``, ``TRANS_MERC`` and ``AZIMUTHAL_EQUIDIST`` options make transformations
  of geographic coordinates into a Cartesian/rectangular system. The
 `NONE`` transformation performs no geographic conversion.
|    ``latOrig`` (*float*, min:\ ``-90.0``, max:\ ``90.0``) latitude in
  decimal degrees of the rectangular coordinates origin
|    ``longOrig`` (*float*, min:\ ``-180.0``, max:\ ``180.0``) longitude
  in decimal degrees of the rectangular coordinates origin
|    ``rotAngle`` (*float*, min:\ ``-360.0``, max:\ ``360.0``) rotation
  angle of geographic north in degrees clockwise relative to the
  rectangular coordinates system Y-axis
|    ``refEllipsoid`` (*choice*:
 `WGS-84 GRS-80 WGS-72 Australian Krasovsky International Hayford-1909 Clarke-1880 Clarke-1866 Airy Bessel Hayford-1830 Sphere``)
  reference ellipsoid name
|    ``latOrig`` (*float*, min:\ ``-90.0``, max:\ ``90.0``) latitude in
  decimal degrees of the rectangular coordinates origin
|    ``longOrig`` (*float*, min:\ ``-180.0``, max:\ ``180.0``) longitude
  in decimal degrees of the rectangular coordinates origin
|    ``firstStdParal secondStdParal`` (*float*, min:\ ``-90.0``,
  max:\ ``90.0``) first and second standard parallels (meridians) in
  decimal degrees
|    ``rotAngle`` (*float*, min:\ ``-360.0``, max:\ ``360.0``) rotation
  angle of geographic north in degrees clockwise relative to the
  rectangular coordinates system Y-axis
| Notes:
|    1. `` rotAngle                ` = ``0`` gives
  North along the positive Y-axis,
 ` rotAngle                ` = ``-30`` gives
  North along the axis 30 deg counterclockwise from the positive Y-axis
  of the rotated, rectangular system.
|    2. The ``GLOBAL`` mode uses a "spherical earth",
  longitude,latitude,depth coordinate system (positive X = East,
  positive Y = North, positive Z = down). Longitude and latitude units
  are degrees, depth is in kilometers, most input/output distance
  quantities are expressed in geographic ( longitude,latitude,depth)
  coordinates, most, but not all, horizontal distances are in degrees.
|    3. The ``SIMPLE`` transformation only corrects longitudinal
  distances as a function of latitude **Algorithm:**
 ` x = (long - longOrig) * 111.111 * cos(lat_radians);                     y = (lat - latOrig) * 111.111;                     lat = latOrig + y / 111.111;                     long = longOrig + x / (111.111 * cos(lat_radians));                `
|    4. The ``NONE`` transformation performs no geographic conversion,
  The input, cartesian/rectangular XYZ coordinates are used throughout
|    5. The ``SDC`` transformation is a "Short Distance Conversion"
  projection indended for use with very small study regions. For
  algorithm details see ``MAP_TRANS_SDC`` in GridLib.c
|    6. The ``LAMBERT`` Lambert Conformal Conic projection is adapted from
  the source code of the GMT plotting package
|    7. The ``TRANS_MERC`` Transverse Mercator projection is adapted from
  the source code of the GMT plotting package
|    8. The ``AZIMUTHAL_EQUIDIST`` Transverse Mercator projection is adapted from
  the source code of the GMT plotting package

| **MAPLINE - Geographic Maplines**
| *optional*, *repeatable*
| Syntax 1: ``MAPLINE`` ``formatType name red green blue lineStyle``
| Specifies a file and drawing parameters for geographic line data.
|    ``formatType`` (*choice*:
 `GMT_LATLON GMT_LONLAT XY_LONLAT GMT_LONLATDEPTH GMT_LONLATELEV_M GMT_GRD``)
  line file format or GMT grd file format
|    ``name`` (*string*) full path and file name
|    ``red green blue`` (*float*, min:\ ``0.0``, max:\ ``1.0``) red,
  green and blue intensities (0.0-1.0) (not implemented)
|    ``lineStyle`` (*choice*: ``SOLID DASHED DOTTED DASHDOT``) line
  style (not implemented)
| Notes:
|    1. All formats except ``GMT_GRD`` specify 2D or 3D line files. Use
 `GMT_GRD`` to specify GMT grd files, these will be plotted as a
  background image with a green-greyscale by default. If a GMT cpt file
  exists with the same path and name as the GMT grd file, but ending
  with ".cpt", it will be used to determine the color scale.
|    2. A GMT grid (``GMT_GRD``) cannot be used with a rotated
  coordinate system.

| **MAPTRANS - Geographic Transformation for Grid2GMT plot output**
| *optional*, *non-repeatable*
| Syntax 1: `` `` `` ``
| Sets geographic to working coordinates transformation parameters for
  Grid2GMT plotting output. See ``TRANS`` above for syntax (replacing
 `TRANS`` by ``MAPTRANS``).
|    
| Notes:
|    1. ``MAPTRANS`` specifies the transformation for Grid2GMT output to
  GMT plotting.
|    2. ``MAPTRANS`` superseeds any other ``TRANS`` statement in the
  control file for Grid2GMT output.

| **MAPGRID - Grid Description for Grid2GMT plot output**
| *optional*, *non-repeatable*
| Syntax 1: ``MAPGRID``
 `xNum yNum zNum xOrig yOrig zOrig dx dy dz gridType``
| Specifies the size and type of the 3D GMT plotting grid.
|    ``xNum yNum zNum`` (*integer*, min:\ ``2``) number of grid nodes in
  the x, y and z directions
|    ``xOrig yOrig zOrig`` (*float*) x, y and z location of the grid
  origin in km relative to the geographic origin.
|    ``dx dy dz`` (*float*) grid node spacing in kilometers along the x,
  y and z axes
|    ``gridType`` (*choice*: ``XXX``) grid type (ignored).
| Notes:
|    1. The 3D grid dimensions are in kilometers with Z positive down
  (left-handed coordinate system).
|    2. The grid is *dx\*(xNum-1)* km long in the x direction, and
  similarly for y and z.
|    3. ``MAPGRID`` specifies the plot region for GRid2GMT output to GMT
  plotting. ``MAPGRID`` superseeds any other ``xxxGRID`` statements in
  the control file.


Vel2Grid Program
----------------


```VGOUT`` <#_Vel2Grid_vgout_>`__ - ```VGTYPE`` <#_Vel2Grid_vgtype_>`__
- ```VGGRID`` <#_Vel2Grid_vgrid_>`__ - ```LAYER`` <#_Vel2Grid_layer_>`__
- ```2DTO3DTRANS`` <#_Vel2Grid_2d3dtrans_>`__ -
```VERTEX`` <#_Vel2Grid_vertex_>`__ - ```EDGE`` <#_Vel2Grid_edge_>`__ -
```POLYGON2`` <#_Vel2Grid_polygon2_>`


| **VGOUT - Output File Root Name**
| *required*, *non-repeatable*
| Syntax 1: ``VGOUT`` ``fileRoot``
| Specifies the directory path and file *root* name (no extension) for
  the output velocity grid.
|    ``fileRoot`` (*string*) full or relative path and file *root* name
  (no extension) for output
| Notes:
|    1. The 3D velocity grid ouput files have names of the form:
|    ``fileRoot.waveType.mod``.*FileExtension*
|    where *FileExtension* is ``.buf`` or ``.hdr`` .

| **VGTYPE - Wave Type**
| *required*, *repeatable*
| Syntax 1: ``VGTYPE`` ``waveType``
| Specifies the physical wave type for a velocity grid.
|    ``waveType`` (*choice*: ``P S``) wave type

| **VGGRID - Grid Description**
| *required*, *non-repeatable*
| Syntax 1: ``VGGRID``
 `xNum yNum zNum xOrig yOrig zOrig dx dy dz gridType``
| Specifies the size and type of the 3D velocity grid.
|    ``xNum yNum zNum`` (*integer*, min:\ ``2``) number of grid nodes in
  the x, y and z directions
|    ``xOrig yOrig zOrig`` (*float*) x, y and z location of the grid
  origin in km relative to the geographic origin.
|    ``dx dy dz`` (*float*) grid node spacing in kilometers along the x,
  y and z axes
|    ``gridType`` (*choice*:
 `VELOCITY VELOCITY_METERS SLOWNESS VEL2 SLOW2 SLOW_2_METERS SLOW_LEN``)
  physical quantity to store on grid ( ``VELOCITY`` = km/s,
 `VELOCITY_METERS`` = m/s, ``SLOWNESS`` = s/km, ``VEL2`` = vel\*\*2,
 `SLOW2`` = (s/km)\*\*2, ``SLOW_2_METERS`` = slow\*\*2 ((s/m)\*\*2),
 `SLOW_LEN`` = slow\*dx (sec)).
| Notes:
|    1. The 3D grid dimensions are in kilometers with Z positive down
  (left-handed coordinate system).
|    2. The grid is *dx\*(xNum-1)* km long in the x direction, and
  similarly for y and z.
|    3. For a 2D grid xNum=2 and xOrig=yOrig=0.0 since a 2D grid
  represents a 1D model and is invariant with respect to translations in
  x or y.

| **LAYER - Velocity Model - Layer**
| *optional*, *repeatable*
| Syntax 1: ``LAYER`` ``depth VpTop VpGrad VsTop VsGrad rhoTop rhoGrad``
| Specifies a constant or gradient velocity layer.
|    ``depth`` (*float*) depth to top of layer (use negative values for
  layers above z=0)
|    ``VpTop VsTop rhoTop`` (*float*) P velocity, and S velocity in km/s
  and density in kg/m\*\*3 at the top of the layer.
|    ``VpGrad VsGrad rhoGrad`` (*float*) Linear P velocity and S
  velocity gradients in km/s/km and density gradient in kg/m\*\*3/km
  increasing directly downwards from the top of the layer.
| Notes:
|    1. Multiple layers must be specified in order of increasing depth
  of top of layer.
|    2. The layer with the deepest top extends implicitly to infinite
  depth.

| **2DTO3DTRANS - Velocity Model - 2D model to 3D model transformation**
| *optional*, *non-repeatable*
| Syntax 1: ``2DTO3DTRANS`` ``xOrig yOrig rotation``
|    ``xOrig yOrig`` (*float*) x and y coordinates in kilometers of the
  center of rotation in the 3D model.
|    ``rotation`` (*float*, min:\ ``-360.0``, max:\ ``360.0``) rotation
  angle in degreees COUNTERCLOCKWISE.
| Notes:
|    1. The 2D to 3D transformation is applied after the general
  geographic transformation specified by the Generic control statement
 `TRANS`` .
|    2. With `` rotation                ` =0 the
  2D model section will be parallel to the *x* direction in the 3D
  model, and the 2D model will be extended along the *y* direction in
  the 3D model.

| **VERTEX - Velocity Model - Vertex**
| *optional*, *repeatable*
| Syntax 1: ``VERTEX`` ``id_num zloc xloc yloc``
| Specifies a vertex in 2D or 3D space.
|    ``id_num`` (*integer*) vertex identification number (must be
  unique)
|    ``zloc xloc yloc`` (*float*) z (positive DOWN), x and y location in
  kilometers of vertex ( *yloc* ignored for 2D models)
| Notes:
|    1. A single vertex may be used in the definitions of multiple edges
  (see EDGE).

| **EDGE - Velocity Model - Edge**
| *optional*, *repeatable*
| Syntax 1: ``EDGE`` ``id_num vertex1 vertex2``
|    ``id_num`` (*integer*) edge identification number (must be unique)
| Notes:
|    1. A single edge may be used in the definitions of multiple 2D
  polygons (see POLYGON2).

| **POLYGON2 - Velocity Model - 2D polygon**
| *optional*, *repeatable*
| Syntax 1: ``POLYGON2`` ``id_num n_edges depth Vp_top Vp_grad Vs_top Vs_grad p_top p_grad   [NEW_LINE]  edge1, edge2, ...``
|    ``id_num`` (*integer*) polygon identification number (must be
  unique)
|    ``n_edges`` (*integer*, min:\ ``0``) the number of edges defining
  this polygon
|    ``depth`` (*float*) reference depth for velocity and density (use
  negative values for depths above z=0)
|    ``VpTop VsTop rhoTop`` (*float*) P velocity, and S velocity in km/s
  and density in kg/m\*\*3 at the reference depth (
 `depth                    ` ).
|    ``VpGrad VsGrad rhoGrad`` (*float*) Linear P velocity and S
  velocity gradients in km/s/km and density gradient in kg/m\*\*3/km
  increasing directly downwards from the reference depth (
 `depth                    ` ).
|    ``edge1, edge2, ...`` (*integer*) new line containing list of edge
  indexes defining this polygon
| Notes:
|    1. A 2D polygon may share edges with other 2D polygons.
|    2. The reference depth (
 ` depth                ` ) may be above,
  within, or below the polygon.

Vel2Grid3D Program
------------------


```VGINP`` <#_Vel2Grid3D_vginp_>`__ -
```VGCLIP`` <#_Vel2Grid3D_vgclip_>`


| **VGINP - Input Velocity Model File**
| *required*, *non-repeatable*
| Syntax 1: ``VGINP`` ``inputFile fileType params``
| Specifies the path/name, type and optional parameters of the input
  velocity model file.
|    ``inputFile`` (*string*) full or relative path and filename
|    ``fileType`` (*choice*: ``SIMUL2K FDTOMO``) File type of input 3D
  velocity models defined by interpolation between control point nodes.
|    ``params`` (*float*) For FDTOMO type requires: orig\_x orig\_y
  orig\_z num\_x num\_y num\_z d\_x d\_y d\_z

| **VGCLIP - Clip Limits for Output Velocity**
| *optional*, *non-repeatable*
| Syntax 1: ``VGCLIP`` ``Vmin Vmax``
| Sets minimum and maximum clip limits for the output velocity values,
  or controls sharpening of velocity difference an interface
|    ``Vmin Vmax`` (*float*) minimum and maximum clip limits for the
  output velocity values.
| Notes:
|    1. If Vmin < Vmax: sets minimum and maximum clip limits for the
  output velocity values.
|    2. IVmin > Vmax: sharpens velocity difference across an interface
  (such as the Moho): if velocity at node below current input node is >
  Vmax: set the NLL grid point velocity equal to the velocity of the
  deepest overlying input node with velocity < Vmax.


Grid2Time Program
-----------------


```GTFILES`` <#_Grid2Time_gtfiles_>`__ -
```GTMODE`` <#_Grid2Time_gtmode_>`__ -
```GTSRCE`` <#_Grid2Time_gtsrce_>`__ -
```GT_PLFD`` <#_Grid2Time_gt_plfd_>`


| **GTFILES - Input and Output File Root Name**
| *required*, *non-repeatable*
| Syntax 1: ``GTFILES``
 `ttimeFileRoot outputFileRoot waveType iSwapBytesOnInput``
| Specifies the directory path and file *root* name (no extension), and
  the wave type identifier for the input velocity grid and output time
  grids.
|    ``ttimeFileRoot`` (*string*) full or relative path and file *root*
  name (no extension) for input velocity grid (generated by program
  Vel2Grid)
|    ``outputFileRoot`` (*string*) full or relative path and file *root*
  name (no extension) for output travel-time and take-off angle grids
|    ``waveType`` (*choice*: ``P S``) wave type
|    ``iSwapBytesOnInput`` (*integer*, min:\ ``0``, max:\ ``1``,
  default:\ ``0``) flag to indicate if hi and low bytes of input
  velocity grid file should be swapped
| Notes:
|    1. The
 `ttimeFileRoot                    ` and
 `outputFileRoot                    ` are
  appended with
 `.                         waveType                    `
|    2. The 3D time grid ouput files have names of the form:

   `outputFileRoot.waveType                             .                             label                        `
    . *gridType* . *FileExtension*

where *label* is a source label ( *i.e.* a station or N_S_L_C codes code), *gridType* is
``time`` or ``angle`` , *FileExtension* is ``.buf`` or ``.hdr``.

| **GTMODE - Program Modes**
| *required*, *non-repeatable*
| Syntax 1: ``GTMODE`` ``gridMode angleMode``
| Specifies several program run modes.
|    ``gridMode`` (*choice*: ``GRID3D GRID2D``) grid type (
 `GRID3D                        ` for a
  3D, Nx\*Ny\*Nz grid or
 `GRID2D                        ` for a
  2D, 2\*Ny\*Nz grid)
|    ``angleMode`` (*choice*: ``ANGLES_YES ANGLES_NO ANGLES_INCLINATION``) sets if take-off
  angles are calculated and an angles grid is output ( ``ANGLES_YES``
  for angles calulcation or ``ANGLES_NO`` for no angles calculation,
  or ``ANGLES_INCLINATION`` for inclination angle calculation only with full precision)

| **GTSRCE - Source Description**
| *required*, *repeatable*
| Syntax 1: ``GTSRCE`` ``label XYZ xSrce ySrce zSrce elev``
| Syntax 2: ``GTSRCE`` ``label LATLON latSrce longSrce zSrce elev``
| Syntax 3: ``GTSRCE``
 `label LATLONDM latDegSrce latMinSrce latDir longDegSrce longMinSrce longDir zSrce elev``
| Syntax 4: ``GTSRCE``
 `label LATLONDS latDegSrce latMinSrce latSecSrce latDir longDegSrce longMinSrce longSecSrce longDir zSrce elev``
| Specifies a source location. One time grid and one angles grid (if
  requested) will be generated for each source. Four formats are
  supported: ``XYZ`` (rectangular grid coordinates), ``LATLON`` (decimal
  degrees for latitude/longitude), ``LATLONDM`` (degrees + decimal
  minutes for latitude/longitude) and ``LATLONDS`` (degrees + minutes +
  decimal seconds for latitude/longitude).
|    ``label`` (*string*) source label ( *i.e.* a station or N_S_L_C codes code: ``ABC``
  )
|    ``xSrce ySrce`` (*float*) x and y grid positions relative to
  geographic origin in kilometers for source
|    ``zSrce`` (*float*) z grid position (depth, positive DOWN) in
  kilometers for source
|    ``elev`` (*float*) elevation above z grid position (positive UP) in
  kilometers for source
|    ``latSrce`` (*float*, min:\ ``-90.0``, max:\ ``90.0``) latitude in
  decimal degrees for source (pos = North)
|    ``longSrce`` (*float*, min:\ ``-180.0``, max:\ ``180.0``) longitude
  in decimal degrees for source (pos = East)
|    ``latDegSrce latMinSrce latSecSrce`` (*float*) latitude degrees,
  minutes and seconds for source
|    ``longDegSrce longMinSrce longSecSrce`` (*float*) longitude
  degrees, minutes and seconds for source
|    ``latDir`` (*choice*: ``N S``) geographic direction
|    ``longDir`` (*choice*: ``W E``) geographic direction

| **GT\_PLFD - Podvin and Lecomte Finite Difference**
| *required*, *non-repeatable*, for Podvin and Lecomte finite
  difference, must not be present otherwise
| Syntax 1: ``GT_PLFD`` ``hs_eps_init message_flag``
| Selects Podvin and Lecomte finite difference method and specifies
  method parameters.
|    ``hs_eps_init`` (*float*, min:\ ``0.0``) fraction (typically
  1.0E-3) defining the tolerated model inhomogeneity for exact
  initialization. A tolerance larger than 0.01 will potentially create
  errors larger than those involved by the F.D. scheme without any exact
  initialization.
|    ``message_flag`` (*integer*, min:\ ``0``, max:\ ``2``) Message flag
  (0:silent, 1:few messages, 2:verbose) A negative value inhibits
  "clever" initialization.
| Notes:
|    1. See Podvin and Lecomte finite difference source code and Podvin
  and Lecomte, 1991 for more information.

Time2EQ Program
---------------


```EQFILES`` <#_Time2EQ_eqfiles_>`__ -
```EQEVENT`` <#_Time2EQ_eqevent_>`__ - ```EQSTA`` <#_Time2EQ_eqsta_>`__
- ```EQSRCE`` <#_Time2EQ_eqsrce_>`__ -
```EQMECH`` <#_Time2EQ_eqmech_>`__ - ```EQMODE`` <#_Time2EQ_eqmode_>`__
- ```EQQUAL2ERR`` <#_Time2EQ_eqqual2err_>`__ -
```EQVPVS`` <#_Time2EQ_eqvpvs_>`


| **EQFILES - Input and Output File Root Name**
| *required*, *non-repeatable*
| Syntax 1: ``EQFILES`` ``ttimeFileRoot outputFileName``
| Specifies the directory path and file *root* name (no extension) for
  the input time grids, and the path and filename for the output
  phase/observation file.
|    ``ttimeFileRoot`` (*string*) full or relative path and file *root*
  name (no extension) for input time grids (generated by program
  Grid2Time)
|    ``outputFileName`` (*string*) full or relative path and name for
  output phase/observation file
| Notes:
|    1. The `` ttimeFileRoot                `
  should not include the standardized phase code ( *i.e.* ``P`` or ``S``
  ).

| **EQEVENT - Hypocenter parameters**
| *optional*, *repeatable*
| Syntax 1: ``EQEVENT`` ``label xEvent yEvent zEvent originSeconds``
|    ``label`` (*string*) event identification label
|    ``xEvent yEvent zEvent`` (*float*) x, y and z grid coordinates of
  hypocenter
|    ``originSeconds`` (*float*) origin time in seconds
| Notes:
|    1. The the origin time
 ` originSeconds                ` is added to
  the travel-time read from the time grid to get the synthetic phase
  time.

| **EQSTA - Station List**
| *required*, *repeatable*
| Syntax 1: ``EQSTA``
 `label phase errorType error errorReportType errorReport probActive``
| Specifies a station or N_S_L_C code, phase and timing error to use to generate a
  synthetic phase reading.
|    ``label`` (*string*) station or N_S_L_C code label ( *i.e.* ``NN_STA``
  )
|    ``phase`` (*string*) phase type ( *i.e.* ``P`` or ``S`` )
|    ``errorType`` (*choice*: ``GAU BOX FIX NONE``) calculated random
  timing error type ( ``GAU`` for normal deviate with zero mean and
  variance = ``error                    ` ,
  or ``BOX`` for boxcar deviate with zero mean and width = 2 \*
 `error                    ` , or ``FIX``
  for time error/static =
 `error                    ` , or ``NONE``
  for time error/static =
 `0.0                    ` )
|    ``error`` (*float*) error magnitude in seconds
|    ``errorReportType`` (*choice*: ``GAU``) timing error type to write
  to output phase/observation file *Err* field (The current version of
  NLLoc recognizes only ``GAU`` )
|    ``errorReport`` (*float*) error magnitude in seconds to write to
  output phase/observation file *ErrMag* field.
|    ``probActive`` (*float*, default:\ ``1.0``) Probability (0-1) that
  a time for this station/phase should be created.
| Notes:
|    1. The `` label                ` and
 ` phase                ` when concatenated to
  the `` ttimeFileRoot                ` (i.e.
 ` ttimeFileRoot.label.phase                ` )
  should correspond to the path and root name of an existing,
  travel-time grid file.
|    2. The error is calculated stochastically and added to the
  travel-time. Use
 ` error                     = 0.0                `
  to obtain exact synthetic travel-times.

| **EQSRCE - Source Description**
| *optional*, *repeatable*
| Syntax 1: ``EQSRCE`` ``label XYZ xSrce ySrce zSrce elev``
| Syntax 2: ``EQSRCE`` ``label LATLON latSrce longSrce zSrce elev``
| Syntax 3: ``EQSRCE``
 `label LATLONDM latDegSrce latMinSrce latDir longDegSrce longMinSrce longDir zSrce elev``
| Syntax 4: ``EQSRCE``
 `label LATLONDS latDegSrce latMinSrce latSecSrce latDir longDegSrce longMinSrce longSecSrce longDir zSrce elev``
| Specifies a source location. Four formats are supported: ``XYZ``
  (rectangular grid coordinates), ``LATLON`` (decimal degrees for
  latitude/longitude), ``LATLONDM`` (degrees + decimal minutes for
  latitude/longitude) and ``LATLONDS`` (degrees + minutes + decimal
  seconds for latitude/longitude).
|    ``label`` (*string*) source label ( *i.e.* a station or N_S_L_C code ``NN_STA``
  )
|    ``xSrce ySrce`` (*float*) x and y grid positions relative to
  geographic origin in kilometers for source
|    ``zSrce`` (*float*) z grid position (depth, positive DOWN) in
  kilometers for source
|    ``elev`` (*float*) elevation above z grid position (positive UP) in
  kilometers for source
|    ``latSrce`` (*float*) latitude in decimal degrees for source (pos =
  North)
|    ``longSrce`` (*float*) longitude in decimal degrees for source (pos
  = East)
|    ``latDegSrce latMinSrce latSecSrce`` (*float*) latitude degrees,
  minutes and seconds for source
|    ``longDegSrce longMinSrce longSecSrce`` (*float*) longitude
  degrees, minutes and seconds for source
|    ``latDir`` (*choice*: ``N S``) geographic direction
|    ``longDir`` (*choice*: ``W E``) geographic direction

| **EQMECH - Event mechanism description**
| *optional*, *non-repeatable*
| Syntax 1: ``EQMECH`` ``mechType strike dip rake``
| Specifies the mechanism parameters for synthetic first motion
  calculations.
|    ``mechType`` (*choice*: ``DOUBLE ISO NONE``, default:\ ``NONE``)
  source mechanism type ( ``DOUBLE`` for double couple, or ``ISO`` for
  isotropic/explosion, or ``NONE`` for no first motion calculation)
|    ``strike`` (*float*, min:\ ``0.0``, max:\ ``360.0``) strike of
  fault plane in degrees (0,360) clockwise from North in the Geographic
  reference frame (any
 `rotAngle                    ` specified
  in the generic control statement ``GTSRCE`` will be added to
 `strike                    ` ).
|    ``dip`` (*float*, min:\ ``0.0``, max:\ ``90.0``) dip of the fault
  plane in degrees (0,90) down from the horizontal.
|    ``rake`` (*float*, min:\ ``-180.0``, max:\ ``180.0``) angle in
  degrees (-180,180) on the fault plane between the strike direction and
  the slip direction.
| Notes:
|    1. The the origin time
 ` originSeconds                ` is added to
  the travel-time read from the time grid to get the synthetic phase
  time.

| **EQMODE - Select Mode: sta->source or source->station**
| *optional*, *non-repeatable*
| Syntax 1: ``EQMODE`` ``mode``
| Selects calculation of times from single source to multiple stations,
  or from multiple sources to single station. The phase labels in the
  output phase/observation file are set to the station labels or to the
  source labels, depending on the mode.
|    ``mode`` (*choice*: ``SRCE_TO_STA STA_TO_SRCE``,
  default:\ ``SRCE_TO_STA``) ``SRCE_TO_STA`` for single sources to
  multiple stations or ``STA_TO_SRCE`` for single station to multiple
  sources.

| **EQQUAL2ERR - Quality to Error Mapping**
| *required*, *non-repeatable*
| Syntax 1: ``EQQUAL2ERR`` ``Err0 ... ... ... ...``
| Specifies the mapping of error to phase pick quality for output of
  phase/observations in HYPO71 file format (which does not include time
  uncertainties) ( *i.e.* time uncertainties in seconds ( *i.e.*
 `0.01`` or ``0.5`` ) to quality ``0,1,2,3`` or ``4`` ).
|    ``Err0 ... ErrN`` (*float*, min:\ ``0.0``) one time uncertainty
  value for each quality level that may be output to the
  phase/observation file. Synthetic errors less than or equal to the
  first value *Err0* are output with quality ``0`` , less than or equal
  to the second are output with ``1`` , etc.

| **EQVPVS - P Velocity to S Velocity Ratio**
| *optional*, *non-repeatable* (**ver 2.0**)
| Syntax 1: ``EQVPVS`` ``VpVsRatio``
| Specifies the P velocity to S velocity ratio to calculate S phase
  travel-times.
|    ``VpVsRatio`` (*float*) P velocity to S velocity ratio. If
 `VpVsRatio                    ` > 0.0 then
  only P phase travel-times grids are read and
 `VpVsRatio                    ` is used to
  calculate S phase travel-times. If
 `VpVsRatio                    ` < 0.0 then
  S phase travel-times grids are used.


NLLoc Program
-------------


```LOCSIG`` <#_NLLoc_locsig_>`__ - ```LOCCOM`` <#_NLLoc_loccom_>`__ -
```LOCSRCE`` <#_NLLoc_gtsrce_>`__ - ```LOCFILES`` <#_NLLoc_locfiles_>`__
- ```LOCHYPOUT`` <#_NLLoc_lochypout_>`__ -
```LOCSEARCH`` <#_NLLoc_locsearch_>`__ -
```LOCMETH`` <#_NLLoc_locmeth_>`__ - ```LOCGAU`` <#_NLLoc_locgau_>`__ -
```LOCGAU2`` <#_NLLoc_locgau2_>`__ -
```LOCPHASEID`` <#_NLLoc_locphaseid_>`__ -
```LOCQUAL2ERR`` <#_NLLoc_locqual2err_>`__ -
```LOCGRID`` <#_NLLoc_locgrid_>`__ -
```LOCPHSTAT`` <#_NLLoc_locphstat_>`__ -
```LOCANGLES`` <#_NLLoc_locangles_>`__ -
```LOCMAG`` <#_NLLoc_locmag_>`__ - ```LOCCMP`` <#_NLLoc_loccmp_>`__ -
```LOCALIAS`` <#_NLLoc_localias_>`__ -
```LOCEXCLUDE`` <#_NLLoc_locexclude_>`__ -
```LOCDELAY`` <#_NLLoc_locdelay_>`__ -
```LOCELEVCORR`` <#_NLLoc_elevcorr_>`__ -
```LOCTOPO_SURFACE`` <#_NLLoc_topo_surface_>`__ -
```LOCSTAWT`` <#_NLLoc_stawt_>`


| **LOCSIG - Signature text**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCSIG`` ``signature``
| Identification of an individual, institiution or other entity -
  written in some output files.
|    ``signature`` (*line*) signature text

| **LOCCOM - Comment text**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCCOM`` ``comment``
| Comment about location run - written in some output files.
|    ``comment`` (*line*) comment text

| **LOCSRCE - Source Description**
| *optional*, *repeatable* (**ver 3.0**)
| Syntax 1: ``LOCSRCE`` ``...``
| Duplicate of statement GTSRCE - Source Description. Allows
  specification of a station location when using "DEFAULT" travel-time
  grids during TRANS GLOBAL mode location. (If for a given station there
  is no travel-time file containing the station's code in its file name,
  and there is a LOCSRCE entry for this station code, then NLLoc will
  look for a travel-time file containing "DEFAULT" as station code in
  its file name to use for this station. The phase code in the
  travel-time file names must match that for the station's phase
  reading.)

| **LOCFILES - Input and Output File Root Name**
| *required*, *non-repeatable*
| Syntax 1: ``LOCFILES``
 `obsFiles obsFileType ttimeFileRoot outputFileRoot iSwapBytes``
| Specifies the directory path and filename for the phase/observation
  files, and the file *root* names (no extension) for the input time
  grids and the output files.
|    ``obsFiles`` (*string*) full or relative path and name for
  phase/observations files, mulitple files may be specified with
  standard UNIX "wild-card" characters ( ``*`` and ``?`` )
|    ``obsFileType`` (*choice*:
 `NLLOC_OBS HYPO71 HYPOELLIPSE NEIC CSEM_ALERT SIMULPS HYPOCENTER HYPODD SEISAN NORDIC NCSN_Y2K_5 NCEDC_UCB ETH_LOC RENASS_WWW RENASS_DEP INGV_BOLL INGV_BOLL_LOCAL INGV_ARCH``)
  format type for phase/observations files (see Phase File Formats)
|    ``ttimeFileRoot`` (*string*) full or relative path and file *root*
  name (no extension) for input time grids (generated by program
  Grid2Time, edu.sc.seis.TauP.TauP\_Table\_NLL, or other software.
|    ``outputFileRoot`` (*string*) full or relative path and file *root*
  name (no extension) for output files
|    ``iSwapBytes`` (*integer*, min:\ ``0``, max:\ ``1``,
  default:\ ``0``) flag to indicate if hi and low bytes of input time
  grid files should be swapped. Allows reading of travel-time grids from
  different computer architecture platforms during TRANS GLOBAL mode
  location.

| **LOCHYPOUT - Output File Types**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCHYPOUT`` ``fileType1 ... ... ... ... ...``
| Specifies the filetypes to be used for output.
|    ``fileType1 ... fileTypeN`` (*choice*:
 `SAVE_NLLOC_ALL SAVE_NLLOC_SUM NLL_FORMAT_VER_2 FILENAME_DEC_SEC SAVE_NLLOC_EXPECTATION SAVE_NLLOC_OCTREE SAVE_FMAMP SAVE_HYPOELL_ALL SAVE_HYPOELL_SUM SAVE_HYPO71_ALL SAVE_HYPO71_SUM SAVE_HYPOINV_SUM SAVE_HYPOINVERSE_Y2000_ARC SAVE_NLLOC_OCTREE``,
  default:\ ``SAVE_NLLOC_ALL SAVE_HYPOINVERSE_Y2000_ARC``) File format
  types to be output: ``SAVE_NLLOC_ALL`` = save summary and event files
  of type NLLoc Hypocenter-Phase file , Phase Statistics file , Scatter
  file and Confidence Level file ; ``SAVE_NLLOC_SUM`` = save summary
  file only of type NLLoc Hypocenter-Phase file ; ``NLL_FORMAT_VER_2`` =
  save NLLoc Hypocenter-Phase files in new format (WARNING: this new
  output format is currently under development and subject to
  modification.) NLLoc Hypocenter-Phase file , Phase Statistics file ,
  Scatter file and Confidence Level file ; ``FILENAME_DEC_SEC`` = output
  file named with 2 decimal second precision instead of default integer
  second precision - avoids overwriting of output files for multiple
  events or multiple locations with earliest observation time in same
  second ; ``SAVE_NLLOC_EXPECTATION`` = hypocenter, location statistics
  and phase statistics results are based on expectation hypocenter
  instead of maximum likelihood hypocenter (default) NLLoc
  Hypocenter-Phase file ; ``SAVE_NLLOC_OCTREE`` = saving of oct-tree
  structure to disk file when LOCSEARCH OCT used ); ``SAVE_FMAMP`` =
  saving of fmamp hypocenter-phase file for input to fmamp,
  probabilistic first-motion mechanism program ); ``SAVE_HYPOELL_ALL`` =
  save summary and event files of type Quasi-HYPOELLIPSE file ;
 `SAVE_HYPOELL_SUM`` = save summary file only of type
  Quasi-HYPOELLIPSE file ; ``SAVE_HYPO71_ALL`` = save summary and event
  files of type HYPO71 Hypocenter/Station file ; ``SAVE_HYPO71_SUM`` =
  save summary file only of type HYPO71 Hypocenter/Station file ;
 `SAVE_HYPOINV_SUM`` = save summary file only of type HypoInverse
  Archive file ; ``SAVE_HYPOINVERSE_Y2000_ARC`` = save summary file only
  of type HypoInverse Y2000 Archive file ;
| Notes:
|    1. The HypoInverse Archive format serves as input to the program
  FPFIT (Reasenberg *et al.* , 1985) for grid-search determination of
  focal mechanism solutions.

| **LOCSEARCH - Search Type**
| *required*, *non-repeatable*
| Syntax 1: ``LOCSEARCH`` ``GRID numSamplesDraw``
| Syntax 2: ``LOCSEARCH``
 `MET numSamples numLearn numEquil numBeginSave numSkip stepInit stepMin stepFact probMin``
| Syntax 3: ``LOCSEARCH``
 `OCT initNumCells_x initNumCells_y initNumCells_z minNodeSize maxNumNodes numScatter useStationsDensity stopOnMinNodeSize``
| Specifies the search type and search parameters. The possible search
  types are ``GRID`` (grid search), ``MET`` (Metropolis), and ``OCT``
  (Octtree).
|    ``numSamplesDraw`` (*integer*) specifies the number of scatter
  samples to draw from each saved PDF grid ( *i.e.* grid with
 `gridType                    ` =
 `PROB_DENSITY`` and
 `saveFlag                    ` = ``SAVE``
  ) No samples are drawn if
 `saveFlag                    ` < 0.
|    ``numSamples`` (*integer*, min:\ ``0``) total number of accepted
  samples to obtain
|    ``numLearn`` (*integer*, min:\ ``0``) number of accepted samples
  for learning stage of search
|    ``numEquil`` (*integer*, min:\ ``0``) number of accepted samples
  for equilibration stage of search
|    ``numBeginSave`` (*integer*, min:\ ``0``) number of accepted
  samples after which to begin saving stage of search, denotes end of
  equilibration stage
|    ``numSkip`` (*integer*, min:\ ``1``) number of accepted samples to
  skip between saves (
 `numSkip                    ` = ``1``
  saves every accepted sample)
|    ``stepInit`` (*float*) initial step size in km for the learning
  stage ( ``stepInit                    ` <
 `0.0`` gives automatic step size selection. If the search takes too
  long, the initial step size may be too large; this may be the case if
  the search region is very large relative to the volume of the high
  confidence region for the locations.)
|    ``stepMin`` (*float*, min:\ ``0.0``) minimum step size allowed
  during any search stage (This parameter should not be critical, set it
  to a low value.)
|    ``stepFact`` (*float*, min:\ ``0.0``) step factor for scaling step
  size during equilibration stage (Try a value of 8.0 to start.)
|    ``probMin`` (*float*) minimum value of the maximum probability
  (likelihood) that must be found by the end of learning stage, if this
  value is not reached the search is aborted (This parameters allows the
  filtering of locations outside of the search grid and locations with
  large residuals.)
|    ``initNumCells_x initNumCells_y initNumCells_z`` (*integer*)
  initial number of octtree cells in the x, y, and z directions
|    ``minNodeSize`` (*float*) smallest octtree node side length to
  process, the octree search is terminated after a node with a side
  smaller than this length is generated
|    ``maxNumNodes`` (*integer*) total number of nodes to process
|    ``numScatter`` (*integer*) the number of scatter samples to draw
  from the octtree results
|    ``useStationsDensity`` (*integer*, min:\ ``0``, max:\ ``1``,
  default:\ ``0``) flag, if 1 weights oct-tree cell probability values
  used for subdivide decision in proportion to number of stations in
  oct-tree cell; gives higher search priority to cells containing
  stations, stablises convergence to local events when global search
  used with dense cluster of local stations
|    ``stopOnMinNodeSize`` (*integer*, min:\ ``0``, max:\ ``1``,
  default:\ ``1``) flag, if 1, stop search when first min\_node\_size
  reached, if 0 stop subdividing a given cell when min\_node\_size
  reached
| Notes:
|    1. See NLLoc Program Oct-Tree Algorithm , Grid-Search Algorithm and
  Metropolis Sampling Algorithm for more information.
|    2. Samples are saved to a binary, event Scatter file (see Scatter
  file formats ). For the grid-search, because the samples are drawn
  stochastically, the number of samples actually obtained my differ
  slightly from the requested number.
|    3. If a large number of samples are saved, the spatial density of
  samples will be proportional to the PDF.
|    4. The scatter samples are useful for plotting the PDF as a
  transparent "cloud" and for relatively compact disk storage of the
  PDF.

| **LOCMETH - Location Method**
| *required*, *non-repeatable*
| Syntax 1: ``LOCMETH``
 `method maxDistStaGrid minNumberPhases maxNumberPhases minNumberSphases VpVsRatio maxNum3DGridMemory minDistStaGrid iRejectDuplicateArrivals``
| Specifies the location method (algorithm) and method parameters.
|    ``method`` (*choice*: ``GAU_ANALYTIC EDT EDT_OT_WT EDT_OT_WT_ML``)
  location method/algorithm ( ``GAU_ANALYTIC`` = the inversion approach
  of Tarantola and Valette (1982) with L2-RMS likelihood function.
 `EDT`` = Equal Differential Time likelihood function cast into the
  inversion approach of Tarantola and Valette (1982) ``EDT_OT_WT`` =
  Weights EDT-sum probabilities by the variance of origin-time estimates
  over all pairs of readings. This reduces the probability (PDF values)
  at points with inconsistent OT estimates, and leads to more compact
  location PDF's. ``EDT_OT_WT_ML`` = version of EDT\_OT\_WT with EDT
  origin-time weighting applied using a grid-search, maximum-likelihood
  estimate of the origin time. Less efficient than EDT\_OT\_WT which
  uses simple statistical estimate of the origin time.)
|    ``maxDistStaGrid`` (*float*) maximum distance in km between a
  station and the center of the initial search grid; phases from
  stations beyond this distance will not be used for event location
|    ``minNumberPhases`` (*integer*) minimum number of phases that must
  be accepted before event will be located
|    ``maxNumberPhases`` (*integer*) maximum number of accepted phases
  that will be used for event location; only the first
 `maxNumberPhases                    ` read
  from the phase/observations file are used for location
|    ``minNumberSphases`` (*integer*) minimum number of S phases that
  must be accepted before event will be located
|    ``VpVsRatio`` (*float*) P velocity to S velocity ratio. If
 `VpVsRatio                    ` > 0.0 then
  only P phase travel-times grids are read and
 `VpVsRatio                    ` is used to
  calculate S phase travel-times. If
 `VpVsRatio                    ` < 0.0 then
  S phase travel-times grids are used.
|    ``maxNum3DGridMemory`` (*integer*) maximum number of 3D travel-time
  grids to attempt to read into memory for Metropolis-Gibbs search. This
  helps to avoid time-consuming memory swapping that occurs if the total
  size of grids read exceeds the real memory of the computer. 3D grids
  not in memory are read directly from disk. If
 `maxNum3DGridMemory                    ` <
  0 then NLLoc attempts to read all grids into memory.
|    ``minDistStaGrid`` (*float*) minimum distance in km between a
  station and the center of the initial search grid; phases from
  stations closer than this distance will not be used for event location
|    ``iRejectDuplicateArrivals`` (*int*) flag indicating if duplicate
  arrivals used for location (1=reject, 0=use if time diff < sigma / 2);
  duplicate arrivals have same station label and phase name
| Notes:
|    1. See NLLoc Program Inversion Approach for more information on the
 `GAU_ANALYTIC`` method.
|    2. See NLLoc Program EDT likelihood function for more information
  on the ``EDT`` method.
|    3. Phases that are not used for location are written to output
  files and are used for calculating average residuals.

| **LOCGAU - Gaussian Model Errors**
| *required*, *non-repeatable*
| Syntax 1: ``LOCGAU`` ``SigmaTime CorrLen``
| Specifies parameters for Gaussian modelisation-error covariances
 ` Covariance                     ij                `
  between stations ``i`` and ``j`` using the relation ( Tarantola and
  Valette, 1982 ):

   `Covariance                         ij                         =                         SigmaTime                         2                         exp(-0.5(Dist                         2                         ij                         )/                         CorrLen                         2                         )                    `

| where ``Dist`` is the distance in km between stations ``i`` and ``j``
  .
|    ``SigmaTime`` (*float*, min:\ ``0.0``) typical error in seconds for
  travel-time to one station due to model errors
|    ``CorrLen`` (*float*, min:\ ``0.0``) correllaton length that
  controls covariance between stations ( *i.e.* may be related to a
  characteristic scale length of the medium if variations on this scale
  are not included in the velocity model)

| **LOCGAU2 - Travel-Time Dependent Model Errors**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCGAU2`` ``SigmaTfraction SigmaTmin SigmaTmax``
| Specifies parameters for travel-time dependent modelisation-error.
  Sets the travel-time error in proportion to the travel-time, thus
  giving effectively a station-distance weighting, which was not
  included in the standard Tarantola and Valette formulation used by
  LOCGAU. This is important with velocity model errors, because nearby
  stations would usually have less absolute error than very far
  stations, and in general it is probably more correct that travel-time
  error is a percentage of the travel-time. Preliminary results using
  LOCGAU2 indicate that this way of setting travel-time errors gives
  visible improvement in hypocenter clustering. (can currently only be
  used with the EDT location methods)
|    ``SigmaTfraction`` (*float*, min:\ ``0.0``, max:\ ``1.0``) fraction
  of of travel-time to use as error
|    ``SigmaTmin`` (*float*, min:\ ``0.0``) minimum trave-time error in
  seconds
|    ``SigmaTmax`` (*float*, min:\ ``0.0``) maximum trave-time error in
  seconds

| **LOCPHASEID - Phase Identifier Mapping**
| *optional*, *repeatable*
| Syntax 1: ``LOCPHASEID`` ``stdPhase phaseCode1 ... ... ... ... ...``
| Specifies the mapping of phase codes in the phase/observation file (
  *i.e.* ``pg`` or ``Sn`` ) to standardized phase codes ( *i.e.* ``P``
  or ``S`` ).
|    ``stdPhase`` (*string*) standardized phase code (used to generate
  time-grid file names)
|    ``phaseCode1 ... phaseCodeN`` (*string*) one or more phase codes
  that may be present in a phase/observation file that should be mapped
  to the ``stdPhase``.
| Notes:
|    1. In the current version of NLLoc, it is assumed for some
  processing (such as the calculation of average P and S station
  residuals) that the standardized phase codes are ``P`` and ``S`` .
  Thus it is important to use these codes, if possible.
|    2. A phase/observation file code will be used unchanged if no
 `LOCPHASEID`` statement is specified, or the code is not present in
  any ``LOCPHASEID`` statement.

| **LOCQUAL2ERR - Quality to Error Mapping**
| *required*, *non-repeatable*, for phase/observation file formats that
  do not include time uncertainties ; *ignored*, *non-repeatable*,
  otherwise
| Syntax 1: ``LOCQUAL2ERR`` ``Err0 ... ... ... ...``
| Specifies the mapping of phase pick qualities phase/observation file (
  *i.e.* ``0,1,2,3`` or ``4`` ) to time uncertainties in seconds (
  *i.e.* ``0.01`` or ``0.5`` ).
|    ``Err0 ... ErrN`` (*float*, min:\ ``0.0``) one time uncertainty
  value for each quality level that may be used in a phase/observation
  file. The first value *Err0* is assigned to picks with quality ``0`` ,
  the second to picks with quality ``1`` , etc.
| Notes:
|    1. NLLoc requires Gaussian timing error estimates in seconds for
  the data (phase picks), the ``LOCQUAL2ERR`` statement allows a
  conversion of commonly used integer quality codes to *float* time
  values.
|    2. Use a large, positive value ( *i.e.* ``99999.9`` ) to indicate a
  phase pick that should have zero weight (infinite uncertainty).

| **LOCGRID - Search Grid Description**
| *required*, *repeatable*
| Syntax 1: ``LOCGRID``
 `xNum yNum zNum xOrig yOrig zOrig dx dy dz gridType saveFlag``
| Specifies the size and other parameters of an initial or nested 3D
  search grid. The order of ``LOCGRID`` statements is critical (see
  Notes).
|    ``xNum yNum zNum`` (*integer*, min:\ ``2``) number of grid nodes in
  the x, y and z directions
|    ``xOrig yOrig zOrig`` (*float*) x, y and z location of the grid
  origin in km relative to the geographic origin. Use a large, negative
  value ( *i.e.* ``-1.0e30`` ) to indicate automatic positioning of grid
  along corressponding direction (valid for nested grids only, may not
  be used for initial grid).
|    ``dx dy dz`` (*float*) grid node spacing in kilometers along the x,
  y and z axes
|    ``gridType`` (*choice*: ``MISFIT PROB_DENSITY``) statistical
  quantity to calculate on grid
|    ``saveFlag`` (*choice*: ``SAVE NO_SAVE``) specifies if the results
  of the search over this grid should be saved to disk
| Notes:
|    1. The order of ``LOCGRID`` statements is critical: the first
 `LOCGRID`` is the initial search grid which may not have automatic
  positionig along any axes. The succeeding ``LOCGRID`` statements may
  specify automatic positioning along one or more axes (
 ` xOrig, yOrig, zOrig                ` =
 `-1.0e30`` ), but must all be sized ( *i.e.*
 ` dx*(xNum-1)                ` , etc.) so that
  they can be fully contained within the preceeding grid. The NLLoc
  program will attempt to translate a nested grid that intersects a
  boundary of the initial grid so that it is contained inside of the
  initial grid; if this is not possible the location will be terminated
  prematurely.
|    2. With automatic positioning (
 ` xOrig, yOrig, zOrig                ` =
 `-1.0e30`` ), a grid is shifted in x/y/z so that it is centered on
  the minimum misfit hypocenter x/y/z of the preceeding grid.
|    3. Each search over a grid with
 ` gridType                ` = ``PROB_DENSITY``
  is time consuming and should generally only be used for a nested grid
  on which the full PDF is required and will be saved to disk. Use
 ` gridType                ` = ``MISFIT`` for
  the initial grid, for larger nested grids, and for smaller nested
  grids in maximum-likelihood hypocenter searches ( *i.e.* where the PDF
  is not if interest).
|    4. The 3D grid dimensions are in kilometers with Z positive down
  (left-handed coordinate system).
|    5. The grid is *dx\*(xNum-1)* km long in the x direction, and
  similarly for y and z.
|    6. For 2D velocity and travel-time grids LOCGRID should be 3D and
  positioned absolutely in space, thus xNum >> 2 and xOrig and zOrig are
  in general != 0.0

| **LOCPHSTAT - Phase Statistics parameters**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCPHSTAT``
 `RMS_Max NRdgs_Min Gap_Max P_ResidualMax S_ResidualMax Ell_Len3_Max Hypo_Depth_Min Hypo_Depth_Max Hypo_Dist_Max``
| Specifies selection criteria for phase residuals to be included in
  calculation of average P and S station residuals. The average
  residuals are saved to a summary, phase statistics file (see Phase
  Statistics file formats ).
|    ``RMS_Max`` (*float*, default:\ ``VERY_LARGE_DOUBLE``) the maximum
  allowed hypocenter RMS in seconds
|    ``NRdgs_Min`` (*integer*, default:\ ``-1``) the minimum allowed
  hypocenter number of readings
|    ``Gap_Max`` (*float*, default:\ ``VERY_LARGE_DOUBLE``) the maximum
  allowed hypocenter gap in degrees
|    ``P_ResidualMax S_ResidualMax`` (*float*,
  default:\ ``VERY_LARGE_DOUBLE``) the maximum allowed residual in
  seconds for a P or S phase
|    ``Ell_Len3_Max`` (*float*, default:\ ``VERY_LARGE_DOUBLE``) the
  maximum allowed ellipsoid major semi-axis length (km)
|    ``Hypo_Depth_Min Hypo_Depth_Max`` (*float*,
  default:\ ``VERY_LARGE_DOUBLE``) the minimum and maximum allowed
  maximum likelihood hypocenter depth (km)
|    ``Hypo_Dist_Max`` (*float*, default:\ ``VERY_LARGE_DOUBLE``) the
  maximum allowed maximum likelihood hypocenter distance (km)
| Notes:
|    1. Because the maximum residual cutoff is abrupt, it should be
  chosen and used with care.
|    2. In the current version of NLLoc, it is assumed in the
  calculation of average P and S station residuals that the standardized
  phase codes are ``P`` and ``S`` . Thus it is important to use these
  codes, if possible.

| **LOCANGLES - Take-off Angles parameters**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCANGLES`` ``angleMode qualtiyMin``
| Specifies whether to determine take-off angles for the maximum
  likelihood hypocenter and sets minimum quality cutoff for saving
  angles and corresponding phases to the HypoInverse Archive file .
|    ``angleMode`` (*choice*: ``ANGLES_YES ANGLES_NO``,
  default:\ ``ANGLES_YES``) sets if take-off angles are read from angles
  grid files and output to locations files. ( ``ANGLES_YES`` for angles
  determination or ``ANGLES_NO`` for no angles determination)
|    ``qualtiyMin`` (*integer*, default:\ ``5``) sets the minimum
  quality (see Take-Off Angles Algorithm ) for writing take-off angles
  and corresponding phase to the HypoInverse Archive file . ( ``0`` to
 `10`` )

| **LOCMAG - Magnitude Calculation Method**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCMAG`` ``ML_HB f n K Ro Mo``
| Syntax 2: ``LOCMAG`` ``MD_FMAG  c1 c2 c3 c4 c5``
| Specifies the magnitude calculation type and parameters. The possible
  magnitude types are:
| ``ML_HB`` (Local (Richter) magnitude\ *M\ :sub:`L`*\ fromHutton and
  Boore (1987)),

    *M\ :sub:`L`* = log(\ *A f*) +\ *n*\ log(\ *r*/100)
    +\ *K*\ (*r*-100) + 3.0 +\ *S*,

|
| ``MD_FMAG`` (Duration magnitude\ *M\ :sub:`L`*\ fromLahr, J.C., (1989)
  HYPOELLIPSE),

    *MD* = *C\ :sub:`1`* + *C\ :sub:`2`*\ log(\ *Fc*) + *C\ :sub:`3`\ r*
    + *C\ :sub:`4`\ z* + *C\ :sub:`5`*\ [log(*Fc*))\ :sup:`2`,

|
|    ``f`` (*float*, min:\ ``0.0``) scaling factor to convert\ *A*\ to
  an equivalent Wood-Anderson amplitude.
|    ``n`` (*float*) *n* from Hutton and Boore (1987), related to
  geometrical spreading.
|    ``K`` (*float*) *K* from Hutton and Boore (1987).
|    ``Ro`` (*float*, default:\ ``100``) Optional *Reference distance*
  (km) from Hutton and Boore (1987).
|    ``Mo`` (*float*, default:\ ``3.0``) Optional *Reference magnitude*
  from Hutton and Boore (1987).
|    ``c1 c2 c3 c4 c5`` (*float*) *c1 c2 c3 c4 c5* from Lahr, J.C.,
  (1989) HYPOELLIPSE

| **LOCCMP - Magnitude Calculation Component**
| *optional*, *repeatable*
| Syntax 1: ``LOCCMP``
 `label inst comp ampFactor sta_corr_ml_hb sta_corr_fd_fmag``
|    ``label`` (*string*) station or N_S_L_C code label ( *i.e.* ``NN_STA``
  )
|    ``inst`` (*string*) instrument identification ( *i.e.*
 `SP, BRB, VBB`` ) If *inst* begins with ``*``, then arrival is taken
  as having no absolute timing (can currently only be used with the EDT
  location methods)
|    ``comp`` (*string*) component identification ( *i.e.*
 `Z, N, E, H`` )
|    ``ampFactor`` (*float*, min:\ ``0.0``) amplitude factor, amplitude
  read from phase file is multiplied by
 `ampFactor                    ` to obtain
  the amplitude used for magnitude calculation.
|    ``sta_corr_ml_hb`` (*float*) ``ML_HB`` station correction, from
  Hutton and Boore (1987)
|    ``sta_corr_fd_fmag`` (*float*) ``FD_FMAG`` station correction, from
  Lahr, J.C., (1989) HYPOELLIPSE
| Notes:
|    1. Component specific paramaters are applied to all phase
  observations with matching label, instrument and component. Use ``?``
  or ``*`` to disable matching of label, instrument or component.

| **LOCALIAS - Station Code Alias**
| *optional*, *repeatable*
| Syntax 1: ``LOCALIAS``
 `code alias yearStart monthStart dayStart yearEnd monthEnd dayEnd``
| Specifies (1) an alias (mapping) of station codes, and (2) start and
  end dates of validity of the alias. Allows (1) station codes that vary
  over time or in different pick files to be homogenized to match codes
  in time grid files, and (2) allows selection of station data by time.
|    ``code`` (*string*) station code (or station name or source label)
  as read from the phase/observation files, or from the result of
  another alias evaluation
|    ``alias`` (*string*) new station code which will replace
 `code                    ` if the relevant
  phase pick time falls within the start and end dates of validity of
  the alias
|    ``yearStart monthStart dayStart`` (*integer*) year (including
  century), month and day of start date of validity of the alias (
 `0 0 0`` = no start date)
|    ``yearEnd monthEnd dayEnd`` (*integer*) year (including century),
  month and day of end date of validity of the alias ( ``9999 99 99`` =
  no end date)
| Notes:
|    1. In NLLoc, the alias evaluation is applied recursively,
  regardless of the order of the ``LOCALIAS`` statements. Thus, when
  selecting and specifying alias names, beware of infinite recursion.
|    2. A trailing underscore "\_" in an alias will only be used for
  time grid identification, not for output. This allows, for example, a
  station name ``ABC`` to be aliases to the name ``ABC_`` to enforce
  certain dates of validity for the station, this requires that the time
  grids generated by Grid2Time use the station code ``ABC_`` ; in all
  NLLoc output, the code ``ABC`` will be used.

| **LOCEXCLUDE - Exclude Observations**
| *optional*, *repeatable* (**ver 2.0**)
| Syntax 1: ``LOCEXCLUDE`` ``name phase``
|    ``name`` (*string*) station or N_S_L_C code label ( *i.e.* ``NN_STA``
  ) identifier after application of any alias
|    ``phase`` (*string*) phase code beofore mapping by ``LOCPHASEID`` (
 `P`` , ``S`` , ``PN`` , etc).
| Notes:
|    1. Excluded station/phase observations are weighted to 0 and so
  will not be used for location. The residual is calculated for these
  observations and they are written to output files, if a travel-time is
  available.

| **LOCDELAY - Phase Time Delays**
| *optional*, *repeatable*
| Syntax 1: ``LOCDELAY`` ``code phase numReadings delay``
| Specifies P and S delays (station corrections) to be subtracted from
  observed P and S times.
|    ``code`` (*string*) station or N_S_L_C code code (after all alias evaluations)
|    ``phase`` (*string*) phase type ( ``P`` or ``S`` )
|    ``numReadings`` (*integer*) number of residuals used to calculate
  mean residual/delay (not used by NLLoc, included for compatibility
  with the format of a summary, phase statistics file)
|    ``delay`` (*float*) delay in seconds, subtracted from observed time
| Notes:
|    1. The body of a summary, phase statistics file (see Phase
  Statistics file formats ) can be used directly as a set of
 `LOCDELAY`` statements. Thus the average phase residuals from a run
  of NLLoc can be used as the station corrections for later runs of
  NLLoc.

| **LOCELEVCORR - Simple, vertical ray elevation correction**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCELEVCORR`` ``flag  velP  velS``
| Calculates a simple elevation correction using the travel-time of a
  vertical ray from elev 0 to the elevation of the station. This control
  statement is mean to be used in GLOBAL mode with TauP or other time
  grids which use elevation 0 for the station elevation.
|    ``flag`` (*integer*, min:\ ``0``, max:\ ``1``, default:\ ``0``)
  flag to set activation of simple elevation correction (0=NO, 1=Yes)
|    ``velP`` (*float*) sets the P velocity to use for calculation of
  the elevation correction for P type phases (last leg of phase is P or
  p)
|    ``velS`` (*float*) sets the S velocity to use for calculation of
  the elevation correction for S type phases (last leg of phase is S or
  s)

| **LOCTOPO\_SURFACE - Topographic mask for location search region**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCTOPO_SURFACE`` ``gmtGrdFile  flagDumpDecimation``
| Uses a topographic surface file in GMT grid2xyz ascii or binary format
  to mask prior search volume to the half-space below the topography.
|    ``gmtGrdFile`` (*string*) path and file name of a GMT grid2xyz
  ascii (\*.asc) or binary (\*.bin and \*.bin.hdr) file(s) defining the
  topographic surface in coordinates lat(deg)/long(deg)/elev(m)
|    ``flagDumpDecimation`` (*integer*) if
 `flagDumpDecimation                    ` >
  0 write surface data to x-y-z-elev file using decimation factor
 `flagDumpDecimation                    `.
  Output file is in NLL Scatter file format; this format can be plotted
  in SeismicityViewer.
| Notes:
|    1. Important: For binary grd file, filename must end in .bin and
  there must be the associated .bin.hdr ascii header file in the same
  directory
|    2. To convert topo.grd to GMT ascii grid format, use:
 `grdinfo topo.grd > topo.grd.asc ; grd2xyz topo.grd -Z >> topo.grd.asc``
|    3. To convert topo.grd to GMT binary grid format, use:
 `grdinfo topo.grd > topo.grd.bin.hdr ; grd2xyz topo.grd -ZTLd > topo.grd.bin``

| **LOCSTAWT - Station distribution weighting**
| *optional*, *non-repeatable*
| Syntax 1: ``LOCSTAWT`` ``flag  cutoffDist``
| Calculates a weight for each station that is a function of the average
  distance between all stations used for location. This helps to correct
  for irregular station distribution, i.e. a high density of stations in
  regions such as Europe and North America and few or no stations in
  regions such as oceans. The relative weight for station *i* is:

    *wieght\ :sub:`i`* = 1.0 / [ SUM\ :sub:`j`
    exp(-dist\ :sup:`2`/cutoffDist:sup:`2` ]

| where *j* is a station used for location and *dist* is the
  epicentral/great-circle distance between stations *i* and *j*.
|    ``flag`` (*integer*, min:\ ``0``, max:\ ``1``, default:\ ``0``)
  flag to set activation of station distribution weighting (0=NO, 1=Yes)
|    ``cutoffDist`` (*float*) sets the cutoff distance for weighting
  calculation.
 `cutoffDist                    ` < ``0.0``
  sets automatic cutoff distance: equal to the mean distance between all
  pairs of stations used for location.



Loc2ssst Program
----------------

| **LSOUT - Output File Root Name**
| *required*, *non-repeatable*
| Syntax 1: ``LSOUT outputFileRoot``
| |    ``outputFileRoot`` (*string*) full or relative path and file *root* name (no extension) for output ssst and updated travel-time grids
| Notes:
|    1. For the updated travel-time grids, the  `outputFileRoot` is appended with `_ssst_corr`
|    2. The  `outputFileRoot` is appended with `.waveType`
|    3. The 3D time grid ouput files have names of the form:
   `outputFileRoot.waveType.label`.*gridType*.*FileExtension*
where *label* is a source label ( *i.e.* a station or N_S_L_C code code), *gridType* is
``ssst`` or ``time`` , *FileExtension* is ``.buf`` or ``.hdr``.

| **LSLOCFILES - Input and Output File Root Name**
| *required*, *non-repeatable*
| Syntax 1: ``LSLOCFILES``
 ``inputFileRoot`` (*string*) full or relative path and file name (with .hyp extension and optional wild-cards) specifying input NLLoc *.hyp files

| **LOCFILES - Input and Output File Root Name**
| *required*, *non-repeatable*
| Specifies the directory path and filename for the phase/observation
  files, and the file *root* names (no extension) for the input time
  grids. These parameters should be identical to NLLoc Program->LOCFILES used to generate NLLoc *.hyp files specified in LSLOCFILES.
| See NLLoc Program->LOCFILES for syntax.

| **LOCMETH - Location Method**
| *required*, *non-repeatable*
| Specifies the location method (algorithm) and method parameters.
 These parameters should be identical to NLLoc Program->LOCMETH used to generate NLLoc *.hyp files specified in LSLOCFILES.
| See NLLoc Program->LOCMETH for syntax.

  **LOCPHASEID - Phase Identifier Mapping**
| *required*, *non-repeatable*
| Specifies the mapping of phase codes in the phase/observation file (
  *i.e.* ``pg`` or ``Sn`` ) to standardized phase codes ( *i.e.* ``P``
  or ``S`` ). These parameters should be identical to NLLoc Program->LOCPHASEID used to generate NLLoc *.hyp files specified in LSLOCFILES.
| See NLLoc Program->LOCPHASEID for syntax.

| **LSPARAMS - General parameters**
| *required*, *non-repeatable*
| Syntax 1: ``LSPARAMS`` ``CharDist WeightFloor UseRejected``
|    ``CharDist`` (*float*`) Characteristic event-station distance for weighting contribution of an event to SSST correction for a station calculation.
|    ``WeightFloor`` (*float*, min:\ ``0.0``) Small value added to events-node weights so ssst values at large event-node distance remain non-zero (station static).
|    ``UseRejected`` (*integer*, min:\ ``0``, max:\ ``1``, default:\ ``0``) flag to indicate that NLL REJECTED locations should be accepted for SSST processing.


| **LSMODE - Program Modes**
| *required*, *non-repeatable*
| Syntax 1: ``LSMODE`` ``angleMode``
| Specifies angles run modes.
|    ``angleMode`` (*choice*: ``ANGLES_YES ANGLES_NO ANGLES_INCLINATION``) sets if take-off
  angles are calculated and an angles grid is output ( ``ANGLES_YES``
  for angles calulcation or ``ANGLES_NO`` for no angles calculation,
  or ``ANGLES_INCLINATION`` for inclination angle calculation only with full precision)

| **LSGRID - ssst Grid Description**
| *required*, *non-repeatable*
| Syntax 1: ``LSGRID``
 `xNum yNum zNum xOrig yOrig zOrig dx dy dz gridType``
| Specifies the size and other parameters of the 3D grid to save ssst time corrections.
|    ``xNum yNum zNum`` (*integer*, min:\ ``2``) number of grid nodes in
  the x, y and z directions
|    ``xOrig yOrig zOrig`` (*float*) x, y and z location of the grid
  origin in km relative to the geographic origin. Use a large, negative
  value ( *i.e.* ``-1.0e30`` ) to indicate automatic positioning of grid
  along corressponding direction (valid for nested grids only, may not
  be used for initial grid).
|    ``dx dy dz`` (*float*) grid node spacing in kilometers along the x,
  y and z axes
|    ``gridType`` (*choice*: ``SSST_TIMECORR``) statistical quantity to calculate on grid
| Notes:
|    1. The 3D grid dimensions are in kilometers with Z positive down
  (left-handed coordinate system).
|    2. The grid is *dx\*(xNum-1)* km long in the x direction, and
  similarly for y and z.
|    3. For 2D velocity and travel-time grids LSGRID should be 3D and
  positioned absolutely in space, thus xNum >> 2 and xOrig and zOrig are
  in general != 0.0

| **LSOUTGRID - Output travel-time Grid Description**
| *required*, *non-repeatable*
| Syntax 1: ``LSOUTGRID``
 `xNum yNum zNum xOrig yOrig zOrig dx dy dz gridType``
| Specifies the size and other parameters of the 3D grid to save updated travel-times.
|    ``xNum yNum zNum`` (*integer*, min:\ ``2``) number of grid nodes in
  the x, y and z directions
|    ``xOrig yOrig zOrig`` (*float*) x, y and z location of the grid
  origin in km relative to the geographic origin. Use a large, negative
  value ( *i.e.* ``-1.0e30`` ) to indicate automatic positioning of grid
  along corressponding direction (valid for nested grids only, may not
  be used for initial grid).
|    ``dx dy dz`` (*float*) grid node spacing in kilometers along the x,
  y and z axes
|    ``gridType`` (*choice*: ``TIME``) statistical quantity to calculate on grid
| Notes:
|    1. The 3D grid dimensions are in kilometers with Z positive down
  (left-handed coordinate system).
|    2. The grid is *dx\*(xNum-1)* km long in the x direction, and
  similarly for y and z.
|    3. For 2D velocity and travel-time grids LSOUTGRID should be 3D and
  positioned absolutely in space, thus xNum >> 2 and xOrig and zOrig are
  in general != 0.0

| **LSPHSTAT - Phase Statistics parameters**
| *optional*, *non-repeatable*
| Syntax 1: ``LSPHSTAT``
 `RMS_Max NRdgs_Min Gap_Max P_ResidualMax S_ResidualMax Ell_Len3_Max Hypo_Depth_Min Hypo_Depth_Max Hypo_Dist_Max``
| Specifies selection criteria for phase residuals to be included in
  calculation of average P and S station residuals. The average
  residuals are saved to a summary, phase statistics file (see Phase
  Statistics file formats ).
|    ``RMS_Max`` (*float*, default:\ ``VERY_LARGE_DOUBLE``) the maximum
  allowed hypocenter RMS in seconds
|    ``NRdgs_Min`` (*integer*, default:\ ``-1``) the minimum allowed
  hypocenter number of readings
|    ``Gap_Max`` (*float*, default:\ ``VERY_LARGE_DOUBLE``) the maximum
  allowed hypocenter gap in degrees
|    ``P_ResidualMax S_ResidualMax`` (*float*,
  default:\ ``VERY_LARGE_DOUBLE``) the maximum allowed residual in
  seconds for a P or S phase
|    ``Ell_Len3_Max`` (*float*, default:\ ``VERY_LARGE_DOUBLE``) the
  maximum allowed ellipsoid major semi-axis length (km)
|    ``Hypo_Depth_Min Hypo_Depth_Max`` (*float*,
  default:\ ``VERY_LARGE_DOUBLE``) the minimum and maximum allowed
  maximum likelihood hypocenter depth (km)
|    ``Hypo_Dist_Max`` (*float*, default:\ ``VERY_LARGE_DOUBLE``) the
  maximum allowed maximum likelihood hypocenter distance (km)
| Notes:
|    1. Because the maximum residual cutoff is abrupt, it should be
  chosen and used with care.
|    2. In the current version of NLLoc, it is assumed in the
  calculation of average P and S station residuals that the standardized
  phase codes are ``P`` and ``S`` . Thus it is important to use these
  codes, if possible.

| **LSSTATIONS - Stations to process**
| *optional*, *non-repeatable*
| Syntax 1: ``LSSTATIONS sta1,sta2,...``
| Specifies a set of station or N_S_L_C codes to be included for ssst processing. If not present, all stations with travel-time grids and arrivals will be processed.
|    ``sta1,sta2,...`` (*string*) comma separated list without whitespce of stations to use for ssst processing
|



