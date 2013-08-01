The ``pypismtools`` module
======================
 
pypismtools is a collection of python classes and functions to
evaluate studies made with the Parallel Ice Sheet Model (PISM,
http://www.pism-docs.org). pypismtools provides tools for binning and
histogram plotting. It also includes helper functions and wrappers for
things like unit conversion, defining figure sizes and parameters, and
more. Additionally, some scripts are included, to plot netCDF
variables over GeoTIFF files using basemap, to generate colorramps
forQGIS, and to extract contour lines.

Requirements
-------------------------

The following python modules have to be installed

- netCDF4
- gdal (with python bindings)
- pyproj
- py_udunits2 (from https://github.com/ckhroulev/py_udunits2)
- matplotlib (http://matplotlib.org/)
- basemap (http://matplotlib.org/basemap/)
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

``$ basemap-plot.py --singlerow -v surfvelmag --colorbar_label -o surfvelmag.png Greenland_5km_v1.1.nc``

.. figure:: https://github.com/pism/pypismtools/raw/master/docs/surfvelmag.png
   :width: 300px
   :alt: Magnitude of surface velocities.

   Example 1: Magnitude of surface velocities.


- Now, add coastlines (intermediate resolution 'i') and plot ice thickness 'thk' over an etopo background

``$ basemap-plot.py --background etopo --coastlines --map_resolution i --singlerow -v thk -o etopothk.png Greenland_5km_v1.1.nc``

.. figure:: https://github.com/pism/pypismtools/raw/master/docs/etopothk.png
   :width: 260px
   :alt: Ice thickness with ETOPO background.

   Example 2: Ice thickness with ETOPO background.

- Use a geotiff file as background, plot the colorbar horizontally. In this case, projection information is taken from the geotiff:

``$ basemap-plot.py --geotiff mygeotiff.tif --singlecolumn -v
surfvelmag --colorbar_label -o geotiff.png Greenland_5km_v1.1.nc``

.. figure:: https://github.com/pism/pypismtools/raw/master/docs/geotiff.png
   :width: 260px
   :alt: Ice thickness with ETOPO background.

   Example 3: Magnitude of surface velocities over a MODIS mosaic of Greenland.

Examples for qgis-colorramp.py
-------------------------

qgis-colorramp-plot.py creates linear and log-scaled colorramps for QGIS_ from GMT_ colormaps. Many great colormap can be downloaded from http://soliton.vm.bytemark.co.uk/pub/cpt-city/.

To show the bathymetry around Greenland, you can use the IBCAO colormap. By running the following command

``qgis-colorramp.py --vmin -5000 --vmax 1400 --extend -10000 4000 ibcao.cpt``

you get a linear colorramp from -5000m to 1400m, and the last color
will be extended to 3000m. The result should like like

.. figure:: https://github.com/pism/pypismtools/raw/master/docs/ibcao.png
   :width: 200px
   :alt: Linear DEM colormap IBCAO.

For a nice log-scaled colormap to show speeds, try:

``qgis-colorramp.py --a 3 --log --extend 0 30000 Full_saturation_spectrum_CCW.cpt``

.. figure:: https://github.com/pism/pypismtools/raw/master/docs/Full_saturation_spectrum_CCW.png
   :width: 200px
   :alt: Log-scaled colorramp.

To use the colorramp in QGIS, click on 'Layer Properties / Colormap'
and then click on 'Load color map from file'. Choose the txt
file. Also the colorbar is saved as a png file, and can be added in
the 'Print Composer'.

.. _QGIS: http://www.qgis.org/ 
.. _GMT: http://gmt.soest.hawaii.edu/ 

Examples for contour2shp.py
-------------------------

contour2shp.py lets you extract a contour line from a variable in a
netCDF file, and saves it as a polygon in a shapefile. Useful to create a polygon of a drainage basin from the
mask. Or you can extract the 2000m elevation contour:

``contour2shp.py -v usrf -c 2000 -s -o poly.shp Greenland_5km_v1.1.nc``

.. figure:: https://github.com/pism/pypismtools/raw/master/docs/contour2000m.png
   :width: 200px
   :alt: 2000m contour line.

