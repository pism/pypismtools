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

```$ sudo python setup.py ```

To install for the current user, run

```$ python setup.py install --user```


Example for basemap-plot.py
-------------------------

- Download a test data set, e.g. the SeaRISE master data set from

```$ wget -nc http://websrv.cs.umt.edu/isis/images/a/a5/Greenland_5km_v1.1.nc```

- First, plot the magnitude of horizontal surface velocities 'surfvelmag' and save as 'foo.png'

```$ basemap-plot.py --singlerow -v surfvelmag -o foo.png Greenland_5km_v1.1.nc```
