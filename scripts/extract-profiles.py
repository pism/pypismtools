#!/usr/bin/env python
# Copyright (C) 2012-2013 Andy Aschwanden
#

from argparse import ArgumentParser
import numpy as np

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt


## {{{ http://code.activestate.com/recipes/496938/ (r1)
"""
A module that helps to inject time profiling code
in other modules to measures actual execution times
of blocks of code.

"""

__author__ = "Anand B. Pillai"
__version__ = "0.1"

import time

def timeprofile():
    """ A factory function to return an instance of TimeProfiler """

    return TimeProfiler()

class TimeProfiler:
    """ A utility class for profiling execution time for code """
    
    def __init__(self):
        # Dictionary with times in seconds
        self.timedict = {}

    def mark(self, slot=''):
        """ Mark the current time into the slot 'slot' """

        # Note: 'slot' has to be string type
        # we are not checking it here.
        
        self.timedict[slot] = time.time()

    def unmark(self, slot=''):
        """ Unmark the slot 'slot' """
        
        # Note: 'slot' has to be string type
        # we are not checking it here.

        if self.timedict.has_key(slot):
            del self.timedict[slot]

    def lastdiff(self):
        """ Get time difference between now and the latest marked slot """

        # To get the latest slot, just get the max of values
        return time.time() - max(self.timedict.values())
    
    def elapsed(self, slot=''):
        """ Get the time difference between now and a previous
        time slot named 'slot' """

        # Note: 'slot' has to be marked previously
        return time.time() - self.timedict.get(slot)

    def diff(self, slot1, slot2):
        """ Get the time difference between two marked time
        slots 'slot1' and 'slot2' """

        return self.timedict.get(slot2) - self.timedict.get(slot1)

    def maxdiff(self):
        """ Return maximum time difference marked """

        # Difference of max time with min time
        times = self.timedict.values()
        return max(times) - min(times)
    
    def timegap(self):
        """ Return the full time-gap since we started marking """

        # Return now minus min
        times = self.timedict.values()
        return time.time() - min(times)

    def cleanup(self):
        """ Cleanup the dictionary of all marks """

        self.timedict.clear()

def create_profile_axes(filename, projection, flip):
    '''
    Create a profile axis.

    Parameters
    -----------
    filename: filename of ESRI shape file
    projection: proj4 projection object

    Returns
    -------
    list of proviles with
    x: array_like along-profile axis
    lat: array_like latitudes
    lon: array_like longitudes
    
    '''

    profiles = read_shapefile(filename)
    my_profiles = []
    for profile in range(len(profiles)):
        profile_lat, profile_lon, profile_name = profiles[profile]
        if flip:
            profile_lat = profile_lat[::-1]
            profile_lon = profile_lon[::-1]
            profile_name = profile_name[::-1]
        profile_x, profile_y = projection(profile_lon, profile_lat)
        profile = np.zeros_like(profile_x)
        profile[1::] = np.sqrt(np.diff(profile_x)**2 + np.diff(profile_y)**2)
        profile = profile.cumsum()

        my_profiles.append([profile, profile_x, profile_y, profile_lon, profile_lat, profile_name])
    return my_profiles


def read_shapefile(filename):
    '''
    Reads lat / lon from a ESRI shape file.

    Paramters
    ----------
    filename: filename of ESRI shape file.

    Returns
    -------
    lat, lon: array_like coordinates
    
    '''
    import ogr
    import osr
    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.Open(filename, 0)
    layer = data_source.GetLayer(0)
    srs=layer.GetSpatialRef()
    if not srs.IsGeographic():
        print('''Spatial Reference System in % s is not latlon. Converting.'''
              % filename)
        # Create spatialReference, EPSG 4326 (lonlat)
        srs_geo = osr.SpatialReference()
        srs_geo.ImportFromEPSG(4326)
    cnt = layer.GetFeatureCount()
    names = []
    profiles = []
    for pt in range(0, cnt):
        feature = layer.GetFeature(pt)
        try:
            name = feature.name
        except:
            name = str(pt)
        geometry = feature.GetGeometryRef()
        # Transform to latlon if needed
        if not srs.IsGeographic():
            geometry.TransformTo(srs_geo)
        lon = []
        lat = []
        for point in geometry.GetPoints():            
            lon.append(point[0])
            lat.append(point[1])
        profiles.append([lat, lon, name])
    return profiles


def get_dims_from_variable(var):
    '''
    Gets dimensions from netcdf variable

    Parameters:
    -----------
    var: netCDF variable

    Returns:
    --------
    xdim, ydim, zdim, tdim: dimensions
    '''
        
    ## a list of possible x-dimensions names
    xdims = ['x','x1']
    ## a list of possible y-dimensions names
    ydims = ['y','y1']
    ## a list of possible z-dimensions names
    zdims= ['z', 'zb']
    ## a list of possible time-dimensions names
    tdims= ['t', 'time']

    xdim = None
    ydim = None
    zdim = None
    tdim = None
    
    ## assign x dimension
    for dim in xdims:
        if dim in var.dimensions:
            xdim = dim
    ## assign y dimension
    for dim in ydims:
        if dim in var.dimensions:
            ydim = dim
    ## assign y dimension
    for dim in zdims:
        if dim in var.dimensions:
            zdim = dim
    ## assign y dimension
    for dim in tdims:
        if dim in list(nc.dimensions.keys()):
            tdim = dim
    return xdim, ydim, zdim, tdim


def piecewise_bilinear(x, y, p_x, p_y, p_i, p_j, A, B, C, D):
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
    p_i, p_j: 1d indices arrays
    A, B, C, D: array_like containing corner values

    Returns
    -------
    pw_linear: array with shape like p_i containing interpolated values
    
    '''

    delta_x = p_x - x[p_i]
    delta_y = p_y - y[p_j]

    alpha = 1./dx * delta_x
    beta  = 1./dy * delta_y

    pw_bilinear = ((1-alpha) * (1-beta) * A + (1-alpha) * beta * B +
                   alpha * beta * C + alpha * (1-beta) * D)

    return pw_bilinear


def dim_permute(
    values, input_order=('time', 'z', 'zb', 'y', 'x'),
    output_order=('time', 'z', 'zb', 'y', 'x')):
    '''
    Permute dimensions of an array_like object.

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
description = '''A script to extract data along (possibly multiple) profile using
piece-wise constant or bilinear interpolation.
The profile must be given as a ESRI shape file.'''
parser = ArgumentParser()
parser.description = description
parser.add_argument("FILE", nargs='*')
parser.add_argument(
    "-b", "--bilinear",dest="bilinear",action="store_true",
    help='''Piece-wise bilinear interpolation, Default=False''',
    default=False)
parser.add_argument(
    "-f", "--flip",dest="flip",action="store_true",
    help='''Flip profile direction, Default=False''',
    default=False)

options = parser.parse_args()
bilinear = options.bilinear
args = options.FILE
flip = options.flip
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
    p_filename = args[0]
    in_filename = args[1]
    if (n_args == 2):
        out_filename = 'profile.nc'
    else:
        out_filename = args[2]

print("-----------------------------------------------------------------")
print("Running script %s ..." % __file__.split('/')[-1])
print("-----------------------------------------------------------------")
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
x_coord = nc_in.variables[xdim][:]
y_coord = nc_in.variables[ydim][:]
x0 = x_coord[0]
y0 = y_coord[0]
dx = x_coord[1] - x_coord[0]
dy = y_coord[1] - y_coord[0]
# read projection information
projection = ppt.get_projection_from_file(nc_in)


# Read in profile data
print("  reading profile from %s" % p_filename)
profiles  = create_profile_axes(
    p_filename, projection, flip)

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
profiledim = 'profile'
nc.createDimension(profiledim)

stationdim = 'station'
nc.createDimension(stationdim)

var = 'profile_name'
var_out = nc.createVariable(var, str, dimensions=(stationdim))
var_out.cf_role = "timeseries_id"
var_out.long_name = "profile name"

var_out = nc.createVariable('profile', 'f', dimensions=(stationdim, profiledim))
var_out.long_name = 'distance along profile'
var_out.units = 'm'

var = 'lon'
var_out = nc.createVariable(var, 'f', dimensions=(stationdim, profiledim))
var_out.units = "degrees_east";
var_out.valid_range = -180., 180.
var_out.standard_name = "longitude"

var = 'lat'
var_out = nc.createVariable(var, 'f', dimensions=(stationdim, profiledim))
var_out.units = "degrees_north";
var_out.valid_range = -90., 90.
var_out.standard_name = "latitude"

for k in range(len(profiles)):
    profile = profiles[k]
    ## We have to unlimited dimensions, so we need to assign start and stop
    ## start:stop where start=0 and stop is the length of the array
    ## or netcdf4python will bail. See
    ## https://code.google.com/p/netcdf4-python/issues/detail?id=76
    pl = len(profile[0])
    nc.variables['profile'][k,0:pl] = np.squeeze(profile[0])
    nc.variables['lon'][k,0:pl] = np.squeeze(profile[3])
    nc.variables['lat'][k,0:pl] = np.squeeze(profile[4])
    nc.variables['profile_name'][k] = profile[5]


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
vars_not_copied = ['lat', 'lat_bnds', 'lat_bounds', 'lon', 'lon_bnds', 'lon_bounds', xdim, ydim, tdim]
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

if tdim is not None:
    var_name = tdim
    var_in = nc_in.variables[tdim]
    dimensions = var_in.dimensions
    datatype = var_in.dtype
    if hasattr(var_in, 'bounds'):
        time_bounds_varname = var_in.bounds
        has_time_bounds = True
    else:
        has_time_bounds = False
    var_out = nc.createVariable(
        var_name, datatype, dimensions=dimensions, fill_value=fill_value)
    var_out[:] = var_in[:]
    for att in var_in.ncattrs():
        if att == '_FillValue':
            continue
        else:
            setattr(var_out, att, getattr(var_in, att))

    has_time_bounds_var = False
    if has_time_bounds:
        try:
            var_in = nc_in.variables[var_name]
            has_time_bounds_var = True
        except:
            has_time_bounds_var = False

    if has_time_bounds_var:
        var_name = time_bounds_varname
        var_in = nc_in.variables[var_name]
        dimensions = var_in.dimensions
        datatype = var_in.dtype
        var_out = nc.createVariable(
            var_name, datatype, dimensions=dimensions, fill_value=fill_value)
        var_out[:] = var_in[:]
        for att in var_in.ncattrs():
            if att == '_FillValue':
                continue
            else:
                setattr(var_out, att, getattr(var_in, att))


print("Copying variables")
for var_name in nc_in.variables:
    profiler = timeprofile()
    if var_name not in vars_not_copied:
        print("  Reading variable %s" % var_name)
        var_in = nc_in.variables[var_name]
        xdim, ydim, zdim, tdim = get_dims_from_variable(var_in)
        in_dims = var_in.dimensions
        datatype = var_in.dtype
        if hasattr(var_in, '_FillValue'):
            fill_value = var_in._FillValue
        else:
            fill_value = None
        if in_dims:
            if len(in_dims) > 1:
                p_dims = [x for x in in_dims if x not in mapplane_dim_names]
                idx = []
                for dim in mapplane_dim_names:
                    idx.append(in_dims.index(dim))
                loc = np.min(idx)
                p_dims.insert(loc, profiledim)
                out_dims = (stationdim, profiledim, tdim, zdim)
                out_dim_order_all = (stationdim, profiledim, tdim, zdim)
                out_dim_order = [x for x in out_dim_order_all if x]
                out_dim_ordered = [x for x in in_dims if x not in mapplane_dim_names]
                out_dim_ordered.insert(0, stationdim)
                out_dim_ordered.insert(2, profiledim)
                out_dim_order = filter(lambda(x): x in out_dim_order, out_dim_ordered)

                var_out = nc.createVariable(
                    var_name, datatype, dimensions=out_dim_order,
                    fill_value=fill_value)

                for k in range(len(profiles)):
                    print("    - processing profile {0}".format(profile[5]))
                    profile = profiles[k]
                    p_x = profile[1]
                    p_y = profile[2]

                    # indices (i,j)
                    p_i = (np.floor((p_x - (x0-dx/2)) / dx)).astype('int')
                    p_j = (np.floor((p_y - (y0-dy/2)) / dy)).astype('int')

                    if bilinear:

                        A_i, A_j = p_i, p_j
                        B_i, B_j = p_i, p_j + 1
                        C_i, C_j = p_i + 1, p_j + 1
                        D_i, D_j = p_i + 1, p_j

                        profiler.mark('read')
                        dim_dict = dict([(tdim, ':'), (xdim, 'A_i'), (ydim, 'A_j'), (zdim, ':')])
                        access_str = ','.join([dim_dict[x] for x in in_dims])
                        A_values = eval('var_in[%s]' % access_str)
                        A_p_values = dim_permute(A_values,
                                                     input_order=p_dims, output_order=out_dim_order)
                        dim_dict = dict([(tdim, ':'), (xdim, 'B_i'), (ydim, 'B_j'), (zdim, ':')])
                        access_str = ','.join([dim_dict[x] for x in in_dims])
                        B_values = eval('var_in[%s]' % access_str)
                        B_p_values = dim_permute(B_values,
                                                     input_order=p_dims, output_order=out_dim_order)
                        dim_dict = dict([(tdim, ':'), (xdim, 'C_i'), (ydim, 'C_j'), (zdim, ':')])
                        access_str = ','.join([dim_dict[x] for x in in_dims])
                        C_values = eval('var_in[%s]' % access_str)
                        C_p_values = dim_permute(C_values,
                                                     input_order=p_dims, output_order=out_dim_order)
                        dim_dict = dict([(tdim, ':'), (xdim, 'D_i'), (ydim, 'D_j'), (zdim, ':')])
                        access_str = ','.join([dim_dict[x] for x in in_dims])
                        D_values = eval('var_in[%s]' % access_str)
                        D_p_values = dim_permute(D_values,
                                                     input_order=p_dims, output_order=out_dim_order)
                        p_read = profiler.elapsed('read')
                        print("    - read in %3.4f s" % p_read)
                        profiler.mark('interp')
                        p_values = piecewise_bilinear(x_coord, y_coord,
                                                            p_x, p_y,
                                                            p_i, p_j,
                                                            A_p_values,
                                                            B_p_values,
                                                            C_p_values,
                                                            D_p_values)
                        p_interp = profiler.elapsed('interp')
                        print("    - interpolated in %3.4f s" % p_interp)
                    else:
                        dim_dict = dict([(tdim, ':'), (xdim, 'p_i'), (ydim, 'p_j'), (zdim, ':')])
                        access_str = ','.join([dim_dict[x] for x in in_dims])
                        profiler.mark('read')
                        in_values = eval('var_in[%s]' % access_str)
                        p_read = profiler.elapsed('read')
                        p_values = dim_permute(in_values,
                                                     input_order=p_dims, output_order=out_dim_order)

                    profiler.mark('write')
                    access_str = 'k,' + ','.join([':'.join(['0', str(coord)]) for coord in p_values.shape])
                    exec('var_out[%s] = p_values' % access_str)
                    p_write = profiler.elapsed('write')
                    print('''    - read in %3.4f s, written in %3.4f s''' % (p_read, p_write))
        else:
            var_out = nc.createVariable(
                var_name, datatype, dimensions=var_in.dimensions,
                fill_value=fill_value)
            var_out[:] = var_in[:]

        for att in var_in.ncattrs():
            if att == '_FillValue':
                continue
            elif att == 'coordinates':
                if tdim:
                    coords = '{0} lat lon'.format(tdim)
                else:
                    coords = 'lat lon'
                setattr(var_out, 'coordinates', coords)

            else:
                setattr(var_out, att, getattr(var_in, att))
        print("  - done with %s" % var_name)


# writing global attributes
script_command = ' '.join([time.ctime(), ':', __file__.split('/')[-1],
                           ' '.join([str(l) for l in args])])
if hasattr(nc_in, 'history'):
    history = nc_in.history
    nc.history = script_command + '\n ' + history
else:
    nc.history = script_command

nc_in.close()
nc.close()
print("Extracted profiles to file %s" % out_filename)
