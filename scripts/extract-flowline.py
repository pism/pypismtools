#!/usr/bin/env python
# Copyright (C) 2012 Andy Aschwanden
#

import time
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


def piecewise_bilinear(x, y, fl_i, fl_j, A, B, C, D):
    '''
    Returns a piece-wise bilinear interpolation.

      ^ y
      |
      |
      B-----C
      |     |
      | *   |   x
    --A-----D---->
      |

    Parameters
    ----------
    x, y: 1d coordinate arrays
    fl_i, fl_j: 1d indices arrays
    A, B, C, D: array_like containing corner values

    Returns
    -------
    pw_linear: array with shape like fl_i containing interpolated values
    '''

    delta_x = fl_x - x[fl_i]
    delta_y = fl_y - y[fl_j]

    alpha = 1. / dx * delta_x
    beta  = 1. / dy * delta_y

    pw_bilinear = ( (1 - alpha) * (1 - beta) * A +
                  (1 - alpha) *      beta  * B +
                  alpha       *      beta  * C +
                  alpha       * (1 - beta) * D )

    return pw_bilinear


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


def dim_permute(values, input_order=('time', 'z', 'zb', 'y', 'x'),
            output_order=('time', 'z', 'zb', 'y', 'x')):
    '''
    Permute dimensions of an array_like object

    Parameters
    ----------
    values : array_like
    input_order : dimension tuple
    output_order: dimension tuple (optional)
                  default ordering is ('time', 'z', 'zb', 'y', 'x')

    Returns
    -------
    values_perm : array_like
    '''

    # filter out irrelevant dimensions
    dimensions = filter(lambda(x): x in input_order,
                        output_order)

    # create the mapping
    mapping = map(lambda(x): dimensions.index(x),
                  input_order)

    if mapping:
        return np.transpose(values, mapping)
    else:
        return values  # so that it does not break processing "mapping"

        
# Set up the option parser

description = '''A script to extract data along a given flowline using
bilinear interpolation.
The flowline must be given in an ascii file with col(0)=lat, col(1)=lon.
The file may have a header in row(0).'''

parser = ArgumentParser()
parser.description = description
parser.add_argument("FILE", nargs='*')

options = parser.parse_args()
args = options.FILE
fill_value = -2e33

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

print("Opening NetCDF file %s ..." % in_filename)
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
x0 = x[0]
y0 = y[0]
dx = x[1] - x[0]
dy = y[1] - y[0]
# set up dimension ordering
dim_order = (xdim, ydim, zdim, tdim)
projection = ppt.get_projection_from_file(nc_in)

# Read in flowline data
print("Reading flowline from %s" % flowline_filename)
fl, fl_x, fl_y, fl_lon, fl_lat = create_flowline_axis(flowline_filename,
                                                      projection)

# indices (i,j)
fl_i = (np.floor((fl_x - x0) / dx)).astype('int')
fl_j = (np.floor((fl_y - y0) / dy)).astype('int')

A_i, A_j = fl_i, fl_j
B_i, B_j = fl_i, fl_j + 1
C_i, C_j = fl_i + 1, fl_j + 1
D_i, D_j = fl_i + 1, fl_j

mapplane_dim_names = (xdim, ydim)

# create dimensions. Check for unlimited dim.
print("Creating dimensions") 
unlimdimname = False
unlimdim = None
# create global attributes.
nc = NC(out_filename, 'w', format='NETCDF4')
# copy global attributes
for attname in nc_in.ncattrs():
    setattr(nc, attname,getattr(nc_in, attname))
# create dimensions
fldim = "flowline"    
nc.createDimension(fldim, len(fl))
var_out = nc.createVariable(fldim, 'f', dimensions=(fldim), fill_value=fill_value)
fldim_values = np.zeros_like(fl)
fldim_values[1::] = np.cumsum(np.sqrt(np.diff(fl_x)**2 + np.diff(fl_y)**2))
var_out[:] = fldim_values
var_out.long_name = 'distance along flowline'
var_out.units = 'm'

for dim_name, dim in nc_in.dimensions.iteritems():
    if dim_name not in (mapplane_dim_names or nc.dimensions):
        if dim.isunlimited():
            unlimdimname = dim_name
            unlimdim = dim
            nc.createDimension(dim_name, None)
        else:
            nc.createDimension(dim_name, len(dim))


# figure out which variables not need to be copied to the new file.
# mapplane coordinate variables
vars_not_copied = ['lat', 'lon', xdim, ydim, tdim]
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


var_name = tdim
try:
    var_in = nc_in.variables[tdim]
    dimensions = var_in.dimensions
    datatype = var_in.dtype
    if hasattr(var_in, 'bounds'):
        time_bounds = var_in.bounds
    var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)
    var_out[:] = var_in[:]
    for att in var_in.ncattrs():
        if att == '_FillValue':
            continue
        else:
            setattr(var_out, att, getattr(var_in, att))
except:
    pass

if time_bounds:
    var_name = time_bounds
    var_in = nc_in.variables[var_name]
    dimensions = var_in.dimensions
    datatype = var_in.dtype
    if hasattr(var, 'bounds'):
        time_bounds = var_in.bounds
    var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)
    var_out[:] = var_in[:]
    for att in var_in.ncattrs():
        if att == '_FillValue':
            continue
        else:
            setattr(var_out, att, getattr(var_in, att))

var = 'lon'
var_out = nc.createVariable(var, 'f', dimensions=(fldim))
var_out.units = "degrees_east";
var_out.valid_range = -180., 180.
var_out.standard_name = "longitude"
var_out[:] = fl_lon

var = 'lat'
var_out = nc.createVariable(var, 'f', dimensions=(fldim))
var_out.units = "degrees_north";
var_out.valid_range = -90., 90.
var_out.standard_name = "latitude"
var_out[:] = fl_lat

print("Copying variables")
for var_name in nc_in.variables:
    if var_name not in vars_not_copied:
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
            input_order = (fldim, zdim, tdim)
            # Create variable
            var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)
            A_values = dim_permute(in_values[A_i,A_j,::],input_order=input_order,
                                    output_order=dimensions)
            B_values = dim_permute(in_values[B_i,B_j,::],input_order=input_order,
                                    output_order=dimensions)
            C_values = dim_permute(in_values[C_i,C_j,::],input_order=input_order,
                                    output_order=dimensions)
            D_values = dim_permute(in_values[D_i,D_j,::],input_order=input_order,
                                    output_order=dimensions)
            var_out[:] = piecewise_bilinear(x, y, fl_i, fl_j, A_values, B_values, C_values, D_values)
        elif (xdim in in_dimensions and ydim in in_dimensions and tdim in in_dimensions):
            in_values = ppt.permute(var_in, dim_order)
            dimensions = (tdim, fldim)
            input_order = (fldim, tdim)
            # Create variable
            var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)
            A_values = dim_permute(in_values[A_i,A_j,::],input_order=input_order,
                                    output_order=dimensions)
            B_values = dim_permute(in_values[B_i,B_j,::],input_order=input_order,
                                    output_order=dimensions)
            C_values = dim_permute(in_values[C_i,C_j,::],input_order=input_order,
                                    output_order=dimensions)
            D_values = dim_permute(in_values[D_i,D_j,::],input_order=input_order,
                                    output_order=dimensions)
            var_out[:] = piecewise_bilinear(x, y, fl_i, fl_j, A_values, B_values, C_values, D_values)
        elif (xdim in in_dimensions and ydim in in_dimensions):
            in_values = np.squeeze(ppt.permute(var_in, dim_order))
            dimensions = (fldim)
            # Create variable
            var_out = nc.createVariable(var_name, datatype,
                                        dimensions=dimensions, fill_value=fill_value)
            A_values = in_values[A_i,A_j]
            B_values = in_values[B_i,B_j]
            C_values = in_values[C_i,C_j]
            D_values = in_values[D_i,D_j]
            var_out[:] = piecewise_bilinear(x, y, fl_i, fl_j, A_values, B_values, C_values, D_values)
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
        print("  - done with %s" % var_name)

# writing global attributes
script_command = ' '.join([time.ctime(), ':', __file__.split('/')[-1],
                           ' '.join([str(x) for x in args])])
if nc.history:
    history = nc.history
    nc.history = script_command + '\n ' + history
else:
    nc.history = script_command

nc_in.close()
nc.close()
print("Done")
