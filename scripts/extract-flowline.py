#!/usr/bin/env python
# Copyright (C) 2012 Andy Aschwanden
#

from argparse import ArgumentParser
import numpy as np
from scipy.interpolate import RectBivariateSpline
from pyproj import Proj
from sys import stderr

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    import PyPISMTools.PyPISMTools as ppt
except:
    import PyPISMTools as ppt


def load_flowline(filename):
    '''
    Loads lat / lon from an ascii file.

    Paramters
    ----------
    flowline_filename: filename of ascii file.

    Returns
    -------
    lat, lon: array_like coordinates

    '''

    try:
        lat, lon = np.loadtxt(filename, usecols=(0,1), unpack=True)
    except:
        lat, lon = np.loadtxt(filename, skiprows=1, usecols=(0,1), unpack=True)

    return lat, lon


def create_flowline_axis(flowline_filename, projection):
    '''
    Create a flowline axis.

    Parameters
    -----------
    flowline_filename: filename of ascii file
    projection: proj4 projection object

    Returns
    -------
    x: array_like along-flowline axis
    lat: array_like latitudes
    lon: array_like longitudes
    '''

    fl_lat, fl_lon = load_flowline(flowline_filename)
    fl_x, fl_y = projection(fl_lon, fl_lat)

    x = np.zeros_like(fl_x)
    x[1::] = np.sqrt(np.diff(fl_x)**2 + np.diff(fl_y)**2)
    x = x.cumsum()

    return x, fl_x, fl_y, fl_lon, fl_lat


        
# Set up the option parser
parser = ArgumentParser()
parser.description = "A script to extract data along a given flowline."
parser.add_argument("FILE", nargs='*')

options = parser.parse_args()
args = options.FILE

n_args = len(args)
required_no_args = 2
max_no_args = 3
if (n_args < required_no_args):
    print(("received $i arguments, at least %i expected"
          % (n_args, required_no_args)))
    import sys.exit
    sys.exit
elif (n_args > max_no_args):
    print(("received $i arguments, no more thant %i accepted"
          % (n_args, max_no_args)))
    import sys.exit
    sys.exit
else:
    flowline_filename = args[0]
    in_filename = args[1]
    if (n_args == 2):
        out_filename = 'flowline.nc'
    else:
        out_filename = args[2]

print("  opening NetCDF file %s ..." % in_filename)
try:
    # open netCDF file in 'read' mode
    nc_in = NC(in_filename, 'r')
except:
    print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
          % in_filename))
    import sys
    sys.exit()

# get the dimensions
xdim, ydim, zdim, tdim = ppt.get_dims(nc_in)
x = nc_in.variables[xdim][:]
y = nc_in.variables[ydim][:]
# set up dimension ordering
dim_order = (tdim, xdim, ydim, zdim)
projection = ppt.get_projection_from_file(nc_in)

fl, fl_x, fl_y, fl_lon, fl_lat = create_flowline_axis(flowline_filename,
                                                      projection)
nf = len(fl)

mapplane_dim_names = (xdim, ydim)

 # create dimensions. Check for unlimited dim.
unlimdimname = False
unlimdim = None
# create global attributes.
nc = NC(out_filename, 'w')
# copy global attributes
for attname in nc_in.ncattrs():
    setattr(nc, attname,getattr(nc_in, attname))
# create dimensions
fldim = "flowline"    
nc.createDimension(fldim)
for dim_name, dim in nc_in.dimensions.iteritems():
    if dim_name not in (mapplane_dim_names or nc.dimensions):
        if dim.isunlimited():
            unlimdimname = dim_name
            unlimdim = dim
            nc.createDimension(dim_name, None)
        else:
            nc.createDimension(dim_name, len(dim))
            
try:
    time_bounds_var_name = nc_in.variables[tdim].bounds
except:
    pass

try:
    nt = len(nc_in.variables[tdim])
except:
    nt = None
try:
    nx = len(nc_in.variables[xdim])
except:
    nx = None
try:
    ny = len(nc_in.variables[ydim])
except:
    ny = None
try:
    nz = len(nc_in.variables[zdim])
except:
    nz = None

# figure out which variables not need to be copied to the new file.
# mapplane coordinate variables
vars_not_copied = ['lat', 'lon', xdim, ydim]
for var_name in nc_in.variables:
    var = nc_in.variables[var_name]
    if hasattr(var, 'grid_mapping'):
        mapping_var_name = var.grid_mapping
        vars_not_copied.append(mapping_var_name)
    if hasattr(var, 'bounds'):
        bounds_var_name = var.bounds
        vars_not_copied.append(bounds_var_name)

vars_not_copied.sort()
last = vars_not_copied[-1]
for i in range(len(vars_not_copied)-2, -1, -1):
    if last == vars_not_copied[i]:
        del vars_not_copied[i]
    else:
        last = vars_not_copied[i]
        
for var_name in nc_in.variables:
    if var_name not in vars_not_copied:
        counter = 0
        var_in = nc_in.variables[var_name]
        datatype = var_in.dtype
        in_dimensions = var_in.dimensions
        if hasattr(var_in, '_FillValue'):
            fill_value = var_in._FillValue
        else:
            fill_value = None
        if (xdim in in_dimensions and ydim in in_dimensions and zdim in in_dimensions and tdim in in_dimensions):
            in_values = ppt.permute(var_in, dim_order)
            dimensions = (tdim, fldim, zdim)
            # Create variable
            var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)            
            fl_values = np.zeros((nt, nf, nz))
            max_counter = nt * nz
            print("\nInterpolating variable %s, " % var_name)
            counter = 0
            stderr.write("percent done: ")
            stderr.write("000")
            for t in range(nt):
                for z in range(nz):
                    values = in_values[t,:,:,z]
                    f = RectBivariateSpline(x, y, values)
                    fl_values[t,:,:,z] = np.array([f(xp, yp)[0,0] for xp, yp in zip(fl_x, fl_y)])
                    stderr.write("\b\b\b%03d" % (100.0 * counter / max_counter))
                    var_out[t,:,z] = fl_values
                    counter += 1
        elif (xdim in in_dimensions and ydim in in_dimensions and tdim in in_dimensions):
            in_values = ppt.permute(var_in, dim_order)
            dimensions = (tdim, fldim)
            # Create variable
            var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)
            fl_values = np.zeros((nt, nf))
            max_counter = nt
            print("\nInterpolating variable %s, " % var_name)
            counter = 0
            stderr.write("percent done: ")
            stderr.write("000")
            for t in range(nt):
                values = in_values[t,:,:]
                f = RectBivariateSpline(x, y, values)
                fl_values[t,:] = np.array([f(xp, yp)[0,0] for xp, yp in zip(fl_x, fl_y)])
                stderr.write("\b\b\b%03d" % (100.0 * counter / max_counter))
                var_out[t,:] = fl_values
                counter += 1
        elif (xdim in in_dimensions and ydim in in_dimensions):
            in_values = np.squeeze(ppt.permute(var_in, dim_order))
            dimensions = (fldim)
            # Create variable
            var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)
            fl_values = np.zeros((nf))
            values = in_values
            f = RectBivariateSpline(x, y, values)
            fl_values[:] = np.array([f(xp, yp)[0,0] for xp, yp in zip(fl_x, fl_y)])
            var_out[:] = fl_values
        else:
            dimensions = in_dimensions
            # Create variable
            var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)
            dimensions = in_dimensions
            if (dimensions > 0):
                in_values = nc.variables[var_name][:]
                var_out[:] = in_values
        for att in var_in.ncattrs():
            if att == '_FillValue':
                continue
            else:
                setattr(var_out, att, getattr(var_in, att))
        print("Done with %s" % var_name)


#nc_in.close()
nc.close()
