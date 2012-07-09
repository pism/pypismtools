The ``PyPISMTools`` module
======================
 
PyPISMTools is a collection of classes and functions to evaluate PISM studies by
providing tools for binning and histogram plotting. It also includes
helper functions and wrappers for things like unit conversion,
defining figure sizes and parameters, and more.

Requirements
-------------------------

The following python modules have to be installed

- netCDF4
- gdal (with python bindings)
- pyproj
- py_udunits2 (from https://github.com/ckhroulev/py_udunits2)
- python imaging library (PIL), optional, needed for for some backgrounds)

Installation
-------------------------

To install for all users, run

``$ sudo python setup.py install``

To install for the current user, run

``$ python setup.py install --user``


Examples for basemap-plot.py
-------------------------

basemap-plot.py is a script to plot a variety of ice sheet model relevant variables from a netCDF file from Greenland and Antarctica data sets. Projection information is retrieved from the first input file, and all subsquent plots are on-the-fly reprojected, which makes the script slow but flexible. 

- Download a test data set, e.g. the SeaRISE master data set from

``$ wget -nc http://websrv.cs.umt.edu/isis/images/a/a5/Greenland_5km_v1.1.nc``

- First, plot the magnitude of horizontal surface velocities 'surfvelmag' and save as 'surfvelmag.png'.

``$ basemap-plot.py --singlerow -v surfvelmag -o surfvelmag.png Greenland_5km_v1.1.nc``

.. figure:: https://github.com/pism/PyPISMTools/raw/master/docs/surfvelmag.png
   :width: 300px
   :alt: Magnitude of surface velocities.

   Example 1: Magnitude of surface velocities.


- Now, add coastlines (intermediate resolution 'i') and plot ice thickness 'thk' over an etopo background

``$ basemap-plot.py --background etopo --coastlines --map_resolution i --singlerow -v thk -o etopothk.png Greenland_5km_v1.1.nc``

.. figure:: https://github.com/pism/PyPISMTools/raw/master/docs/etopothk.png
   :width: 300px
   :alt: Ice thickness with ETOPO background.

   Example 2: Ice thickness with ETOPO background.

- Use a geotiff file as background, plot the colorbar horizontally. In this case, projection information is taken from the geotiff:

``$ basemap-plot.py --geotiff mygeotiff.tif --singlecolumn -v surfvelmag -o geotiff.png Greenland_5km_v1.1.nc``

.. figure:: https://github.com/pism/PyPISMTools/raw/master/docs/geotiff.png
   :width: 300px
   :alt: Ice thickness with ETOPO background.

   Example 3: Magnitude of surface velocities over a MODIS mosaic of Greenland.

Examples for qgis-colorramp.py
-------------------------

qgis-colorramp-plot.py creates linear and log-scaled colorramps for QGIS_ from GMT_ colormaps.

.. _QGIS: http://www.qgis.org/ 
.. _GMT: http://gmt.soest.hawaii.edu/ 

