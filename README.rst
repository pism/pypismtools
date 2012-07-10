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
   :width: 260px
   :alt: Ice thickness with ETOPO background.

   Example 2: Ice thickness with ETOPO background.

- Use a geotiff file as background, plot the colorbar horizontally. In this case, projection information is taken from the geotiff:

``$ basemap-plot.py --geotiff mygeotiff.tif --singlecolumn -v surfvelmag -o geotiff.png Greenland_5km_v1.1.nc``

.. figure:: https://github.com/pism/PyPISMTools/raw/master/docs/geotiff.png
   :width: 260px
   :alt: Ice thickness with ETOPO background.

   Example 3: Magnitude of surface velocities over a MODIS mosaic of Greenland.

Examples for qgis-colorramp.py
-------------------------

qgis-colorramp-plot.py creates linear and log-scaled colorramps for QGIS_ from GMT_ colormaps. Many great colormap can be downloaded from http://soliton.vm.bytemark.co.uk/pub/cpt-city/.

To show the bathymetry around Greenland, you can use the IBCAO colormap. By running the following command

``qgis-colorramp.py --vmin -5000 --vmax 1400 --extend 3000 ibcao.cpt``

you get a linear colorramp from -5000m to 1400m, and the last color
will be extended to 3000m. The result should like like

.. figure:: https://github.com/pism/PyPISMTools/raw/master/docs/ibcao.png
   :width: 200px
   :alt: Linear DEM colormap IBCAO.

For a nice log-scaled colormap to show speeds, try:

``qgis-colorramp.py --a 3 --log --extend 30000 Full_saturation_spectrum_CCW.cpt``

.. figure:: https://github.com/pism/PyPISMTools/raw/master/docs/Full_saturation_spectrum_CCW.png
   :width: 200px
   :alt: Log-scaled colorramp.

To use the colorramp in QGIS, click on 'Layer Properties / Colormap'
and then click on 'Load color map from file'. Choose the txt
file. Also the colorbar is saved as a png file, and can be added in
the 'Print Composer'.

.. _QGIS: http://www.qgis.org/ 
.. _GMT: http://gmt.soest.hawaii.edu/ 

