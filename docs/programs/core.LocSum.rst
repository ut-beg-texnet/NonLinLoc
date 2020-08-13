| 

LocSum - combine location results
=================================

**LocSum** combines NLLoc location results and PDF "scatter-cloud"
samples from a number of events.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

`Overview <#_overview_>`__ - `Running the program-Input <#_running_>`__
- `Output <#_output_>`__ - `[NonLinLoc Home] <./index.html>`__

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Overview
--------

The LocSum utility combines single event NLLoc location files
(`Hypocenter-Phase files <formats.html#_location_hypphs_>`__ and `binary
Scatter files <formats.html#_location_scat_>`__) into a single set of
summary location files.

For flexibility, the LocSum utility takes most of its parameters from
the command line.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Running the program - Input
---------------------------

The LocSum utility takes several of command line arguments.

Synopsis:
``LocSum SizeGridRoot decimFactor OutRoot LocRoot [Len3Max [ProbMin [RMSMax [NRdgsMin [GapMax]]]]]``

**Parameters:**

 ``SizeGridRoot`` (*chars*)
    full or relative path and *root* name (no extension) for a `3D Grid
    Header file <formats.html#_grid_hdr_>`__. The grid dimensions in
    this header file are used to create an empty grid buffer and new
    grid header file with root name ``OutRoot``.
 ``decimFactor`` (*integer*)
    decimation factor (``decimFactor > 0``) for decimating the number of
    PDF Scatter samples. Every ``decimFactor``-th sample is saved to the
    output files.
 ``OutRoot`` (*chars*)
    full or relative path and root name for output files.
 ``LocRoot`` (*chars*)
    full or relative path and *root* name (no extension) for one or more
    `NLLoc single Event Location files <formats.html#_location_>`__.
    Mulitple root names may be specified using standard UNIX "wild-card"
    characters (``*`` and ``?``); however, if any "wild-card" characters
    are used then the path and root name must be enclosed in double
    quotes (") to prevent the shell from evaluating the "wild-card"
    characters.
 ``Len3Max`` (*float*)
    maximum length in kilometers of the longest ellipsoid semi-axis at
    maximum likelihood hypocenter.
 ``ProbMin`` (*float*)
    minimum value of probability at maximum likelihood hypocenter.
 ``RMSMax`` (*float*)
    maximum RMS in seconds at maximum likelihood hypocenter.
 ``NRdgsMin`` (*integer*)
    minimum number of readings used for location.
 ``GapMax`` (*float*)
    maximum azimuth gap in degrees at maximum likelihood hypocenter.

**Notes:**

#. See the `Definitions section <./control.html#_definitions_>`__ of the
   NonLinLoc Control File documentation for more information on
   datatypes.

**Example:**

#. ``LocSum dursum0 1 dursum "dur.*.*.grid0.loc"``

   Using an existing 3D Grid Header file ``dursum0.hdr`` to determine
   the grid size, creates a dummy grid buffer file, a grid header file,
   a set of summary Hypocenter-Phase files, binary Scatter files, and a
   set of ASCII Scatter files for each location in the current directory
   with root name ``"dur.*.*.grid0.loc"``. The output files are written
   to the root name ``dursum``. The scatter samples are not decimated
   since ``decimFactor`` = ``1``.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Output
------

The LocSum utility creates the following files:

#. A summary `Hypocenter-Phase file <formats.html#_location_hypphs_>`__
   named ``OutRoot.hyp``. This file includes ``SCATTER`` blocks.
#. A summary `binary Scatter file <formats.html#_location_scat_>`__
   named ``OutRoot.scat``.
#. A set of summary ASCII Scatter files for *x-y*, *x-y* and *z-y*
   projections, named ``OutRoot.scat.ext``, where ``ext`` =
   ``XY, XZ, ZY`` for sample locations in kilometers and ``ext`` =
   ``longlat.XY, longlat.XZ, longlat.ZY`` for sample locations in
   degrees of latitude and longitude and depth in kilometers. These
   ASCII formats are compatible with the `GMT plotting
   package <http://gmt.soest.hawaii.edu/>`__.
#. A `3D Grid Header file <formats.html#_grid_hdr_>`__ named
   ``OutRoot.hdr`` and an empty `3D Grid Buffer
   file <formats.html#_grid_buf_>`__ named ``OutRoot.buf``. These file
   are created to insure compatibility with post-processing programs and
   utilities.

.. raw:: html

   <div class="color">

 

.. raw:: html

   </div>

Back to `the NonLinLoc site Home page <./index.html>`__.

*Anthony Lomax*
