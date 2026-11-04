# NonLinLoc
The NonLinLoc (Non-Linear Location) package is a set of programs for velocity model construction, travel-time calculation and probabilistic, non-linear, global-search earthquake location in 3D structures, and for visualisation of 3D volume data and location results. Many of the programs operate with a 3D Grid structure which defines a specific, gridded, rectangular volume (Non-GLOBAL mode) or spherical section (GLOBAL mode).

## IMPORTANT: The active NonLinLoc repository has been transferred to: [https://github.com/ut-beg-texnet/NonLinLoc](https://github.com/ut-beg-texnet/NonLinLoc)
Please be aware of this change and the new URL.

## Citation
If you use the NonLinLoc package in your work, please cite the following papers:

- Lomax A., Virieux J., Volant P., Berge-Thierry C. (2000) Probabilistic Earthquake Location in 3D and Layered Models. In: Thurber C.H., Rabinowitz N. (eds) Advances in Seismic Event Location. Modern Approaches in Geophysics, vol 18. Springer, Dordrecht. [https://doi.org/10.1007/978-94-015-9536-0_5](https://doi.org/10.1007/978-94-015-9536-0_5)

- Lomax A., Michelini A., Curtis A. (2014) Earthquake Location, Direct, Global-Search Methods. In: Meyers R. (eds) Encyclopedia of Complexity and Systems Science. Springer, New York, NY. [https://doi.org/10.1007/978-3-642-27737-5_150-2](https://doi.org/10.1007/978-3-642-27737-5_150-2)

If you use NLL-SSST-coherence or any of its components, please cite:

- Lomax, A. & Savvaidis, A. (2022) High-Precision Earthquake Location Using Source-Specific Station Terms and Inter-Event Waveform Similarity. Journal of Geophysical Research: Solid Earth, 127, e2021JB023190. [https://doi.org/10.1029/2021JB023190](https://doi.org/10.1029/2021JB023190)


For other A. Lomax NonLinLoc publications, see [http://alomax.net/pub_list.html](http://alomax.net/pub_list.html)


## Install complete NonLinLoc distribution software package
Clone or download ZIP from this repository [https://github.com/alomax/NonLinLoc](https://github.com/alomax/NonLinLoc)

Set output bin/ directory, if you have not already done this:
```
cd src
rm bin
# either do this:
mkdir bin   # bin/ is a subdirectory of src/
# or this:
ln -s <my_favorite_bin_path> bin   # bin/ is somewhere else
```

To build:
```
cd src
rm CMakeCache.txt
cmake .
make
```
See [http://alomax.net/nlloc](http://alomax.net/nlloc) and [http://alomax.net/nlloc](http://alomax.net/nlloc) -> tutorials for further information

To bulid debug:
edit CMakeLists.txt to comment out Release build type and uncomment Debug bulid type:
```
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)
```
then:
```
cd src
rm CMakeCache.txt
cmake .
make
```

Thanks to Gilles Celli (European Center for Geodynamics and Seismology) for creating and debugging the NonLinLoc CMake build system.


## Documentation
- Future documentation (incomplete): [http://alomax.net/nlloc/docs](http://alomax.net/nlloc/docs)
- Legacy documentation (more extensive): [http://alomax.net/nlloc](http://alomax.net/nlloc)
- Users guides and additional information: see [doc/](doc/) directory in this repository


## Additional information


### NEW! NLLoc-SSST-coherence
NLL-SSST-coherence is an enhanced, absolute-timing earthquake location procedure, iteratively generates traveltime corrections to improve multi-scale precision and uses waveform similarity to improve fine-scale precision.

For details, see:

- Lomax, A. & Savvaidis, A. (2022) High-Precision Earthquake Location Using Source-Specific Station Terms and Inter-Event Waveform Similarity. Journal of Geophysical Research: Solid Earth, 127, e2021JB023190. [https://doi.org/10.1029/2021JB023190](https://doi.org/10.1029/2021JB023190)

- Lomax, A., & Henry, P. (2023). Major California faults are smooth across multiple scales at seismogenic depth. Seismica, 2(1). [https://doi.org/10.26443/seismica.v2i1.324] (https://doi.org/10.26443/seismica.v2i1.324)

Example/Tutorial: Files and instructions for running NLL-SSST-coherence for a subset of Parkfield events [https://zenodo.org/doi/10.5281/zenodo.4756709] (https://zenodo.org/doi/10.5281/zenodo.4756709)


### NLLoc_func_test program demonstrating running NLLoc through a function call
Install as above.
TODO: how to build GNU_SOURCE with cmake?????

Build with:
```
GNU_SOURCE=-D _GNU_SOURCE
```
and comment the line:
```
#GNU_SOURCE=
```

Unpack demo files: unpack: NLL[VER]_func.tgz
To run:
```
cd nll_func
./run_func.sh
```
To clean
```
rm -rf out/*
cd ..
```
See nll_func/run_func.sh for more detail.


### ttime_func_test program demonstrating reading values from 2D or 3D grid file through a function call
Install as above.

Unpack demo files: unpack: NLL[VER]_func.tgz
To run:
```
cd ttime_func
./run_func.sh
cd ..
```
See ttime_func/run_func.sh for more detail.


## IMPORTANT NOTES


### Grid2GMT (2023.02.19)
Grid2GMT now supports GMT 5 and 6 (thanks to https://github.com/sean0921). See src/CMakeLists.txt for details on how to set GMT version when compiling.

GMT 4 may work correctly when compile is set up for GMT 6.



### Covariance errors (2014.10.30)
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


