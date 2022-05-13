.. image:: NonLinLocLogo.gif
    :width: 100px
    :align: left
    :alt: NonLinLocLogo.gif
    :target: https://github.com/alomax/NonLinLoc

NonLinLoc
=========

A software package for Probabilistic, Non-Linear, Global-Search Earthquake Location in 3D Media.

Code is stored on |github|, the development branches are |github_dev|, or the
latest stable release can be found |releases_link|.

.. |github| raw:: html

    <a href="https://github.com/alomax/NonLinLoc" target="_blank">github</a>

.. |releases_link| raw:: html

  <a href="https://github.com/alomax/NonLinLoc/releases" target="_blank">here</a>

.. |github_dev| raw:: html

  <a href="https://github.com/alomax/NonLinLoc/tree/develop" target="_blank">here</a>

.. COMMENTED-OUT BLOCK:
   NonLinLoc uses |Obspy_link| bindings when reading and writing seismic data, and for handling most
   of the event metadata, which ensures that detections can be easily migrated between
   software.

   .. |Obspy_link| raw:: html

     <a href="https://docs.obspy.org/" target="_blank">Obspy</a>

   Also within this package are:

   * :doc:`Correlation re-picking </submodules/core.lag_calc>`;
   * :doc:`Clustering routines for seismic data </submodules/utils.clustering>`;
   * :doc:`Peak finding algorithm (basic) </submodules/utils.findpeaks>`;
   * :doc:`Stacking routines </submodules/utils.stacking>` including phase-weighted stacking based on Thurber at al. (2014);

This package is written by Anthony Lomax and the NonLinLoc developers, and is distributed under the LGPL GNU Licence,
Copyright (C) 1999-2020 Anthony Lomax and NonLinLoc developers.


Citation
--------

If you use the NonLinLoc package in your work, please cite the following papers:

* Lomax A., Virieux J., Volant P., Berge-Thierry C. (2000) Probabilistic Earthquake Location in 3D and Layered Models. In: Thurber C.H., Rabinowitz N. (eds) Advances in Seismic Event Location. Modern Approaches in Geophysics, vol 18. Springer, Dordrecht. <a href="https://doi.org/10.1007/978-94-015-9536-0_5" target="_blank">https://doi.org/10.1007/978-94-015-9536-0_5</a>

* Lomax A., Michelini A., Curtis A. (2014) Earthquake Location, Direct, Global-Search Methods. In: Meyers R. (eds) Encyclopedia of Complexity and Systems Science. Springer, New York, NY. <a href="https://doi.org/10.1007/978-3-642-27737-5_150-2" target="_blank">https://doi.org/10.1007/978-3-642-27737-5_150-2</a>

For other A. Lomax NonLinLoc publications, see http://alomax.net/pub_list.html


Contents:
---------

.. toctree::
   :numbered:
   :maxdepth: 3

   intro
   installation
   updates
   tutorial
   programs
   control file


Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
