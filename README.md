Install complete NonLinLoc distribution software package
------------------------------------------------------
	Clone or download ZIP from this repository (https://github.com/alomax/NonLinLoc)

	To build:
	cd src
	rm CMakeCache.txt
	cmake .
	make
	See http://alomax.net/nlloc and http://alomax.net/nlloc -> tutorials for further information

	To bulid debug:
	# edit CMakeLists.txt to comment out Release build type and uncomment Debug bulid type:
		#set(CMAKE_BUILD_TYPE Release)
		set(CMAKE_BUILD_TYPE Debug)
	cd src
	rm CMakeCache.txt
	cmake .
	make

        Thanks to Gilles Celli (European Center for Geodynamics and Seismology) for creating and debugging the NonLinLoc CMake build system.


Documentation
------------------------------------------------------
        Future documentation: http://alomax.net/nlloc/docs
        Legacy documentation: http://alomax.net/nlloc
        Users guides and additional information: see [doc/]directory in this repository


NLLoc_func_test program demonstrating running NLLoc through a function call
------------------------------------------------------
	Install as above.

	Build with:
	TODO: how to build GNU_SOURCE with cmake?????
	GNU_SOURCE=-D _GNU_SOURCE
	and comment the line:
	#GNU_SOURCE=

	Unpack demo files: unpack: NLL[VER]_func.tgz
	To run:
		cd nll_func
		./run_func.sh
	To clean
		rm -rf out/*
		cd ..
	See nll_func/run_func.sh for more detail.


ttime_func_test program demonstrating reading values from 2D or 3D grid file through a function call
------------------------------------------------------
	Install as above.

	Unpack demo files: unpack: NLL[VER]_func.tgz
	To run:
		cd ttime_func
		./run_func.sh
		cd ..
	See ttime_func/run_func.sh for more detail.


IMPORTANT NOTE: Covariance errors (2014.10.30)
------------------------------------------------------
	Applies to: LOCSEARCH MET and LOCSEARCH OCT
	The algorithm for calculating the covariance of the PDF scatter sample used in NLLoc
	was subject to precision errors when the expectation of the scatter sample (e.g. the event location)
	was far from the origin of the NLL coordinates system relative to the extent of the PDF scatter sample.
	Several tests indicate that "Far from the origin" is of the order of 1000 times the extent of the PDF scatter sample.

	Errors in covariance will affect the ellipsoid and ellipse, standard-errors (erh, erz, etc.).

	For teleseismic locations (NLL GLOBAL mode) errors with the old algorithm occur primarily in
	longitude (X), increasing with expectation longitude and occasionally becoming large towards longitude +/-180deg.
	For local studies in rectangular coordinates where the NLL coordinates origin is far from the network
	and target sources, the errors may be large.  This may be the case, for example, for a micro-seismic study
	using a regional, metric-based cartesian grid coordinate system with origin far from the study area.
	For local studies in rectangular coordinates where the NLL coordinates origin is within or near the network
	and target sources, the errors should be very small or negligible.

	The algorithm for calculating the covariance of the PDF scatter sample used in NLLoc v.6.03.00 and later does not
	have these precision errors and is valid for locations far from the origin of the coordinate system.


