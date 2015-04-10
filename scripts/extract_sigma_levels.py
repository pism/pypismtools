#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden
#



"""This script containts tools for extracting sigma level.
"""

from argparse import ArgumentParser
import numpy as np
from scipy.interpolate import interp1d

from netCDF4 import Dataset as NC

try:
    import pypismtools.pypismtools as ppt
except ImportError:             # pragma: nocover
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


def output_dimensions(input_dimensions):
    """Build a list of dimension names used to define a variable in the
    output file."""
    _, _, zdim, tdim = get_dims_from_variable(input_dimensions)

    if tdim:
        result = [stationdim, tdim, profiledim]
    else:
        result = [stationdim, profiledim]

    if zdim:
        result.append(zdim)

    return result


def get_dims_from_variable(var_dimensions):
    '''
    Gets dimensions from netcdf variable

    Parameters:
    -----------
    var: netCDF variable

    Returns:
    --------
    xdim, ydim, zdim, tdim: dimensions
    '''

    def find(candidates, collection):
        """Return one of the candidates if it was found in the collection or
        None otherwise.

        """
        for name in candidates:
            if name in collection:
                return name
        return None

    # possible x-dimensions names
    xdims = ['x', 'x1']
    # possible y-dimensions names
    ydims = ['y', 'y1']
    # possible z-dimensions names
    zdims = ['z', 'zb']
    # possible time-dimensions names
    tdims = ['t', 'time']

    return [find(dim, var_dimensions) for dim in [xdims, ydims, zdims, tdims]]


def copy_dimensions(in_file, out_file, exclude_list):
    """Copy dimensions from in_file to out_file, excluding ones in
    exclude_list."""
    for name, dim in in_file.dimensions.iteritems():
        if (name not in exclude_list and
                name not in out_file.dimensions):
            if dim.isunlimited():
                out_file.createDimension(name, None)
            else:
                out_file.createDimension(name, len(dim))


def copy_time_dimension(in_file, out_file, name):
    """Copy time dimension, the corresponding coordinate variable, and the
    corresponding time bounds variable (if present) from an in_file to
    an out_file.

    """
    var_in = in_file.variables[name]
    var_out = create_variable_like(in_file, name, out_file)
    var_out[:] = in_file.variables[name][:]

    try:
        bounds_name = var_in.bounds
        var_out = create_variable_like(bounds_name, in_file, out_file)
        var_out[:] = in_file.variables[bounds_name][:]
    except AttributeError:
        # we get here if var_in does not have a bounds attribute
        pass

    
def create_variable_like(in_file, var_name, out_file, dimensions=None,
                         fill_value=-2e9):
    """Create a variable in an out_file that is the same var_name in
    in_file, except possibly depending on different dimensions,
    provided in dimensions.

    """
    var_in = in_file.variables[var_name]
    try:
        fill_value = var_in._FillValue
    except AttributeError:
        # fill_value was set elsewhere
        pass

    if dimensions is None:
        dimensions = var_in.dimensions

    dtype = var_in.dtype
    var_out = out_file.createVariable(var_name, dtype, dimensions=dimensions,
                                      fill_value=fill_value)
    copy_attributes(var_in, var_out)
    return var_out


def copy_attributes(var_in, var_out):
    """Copy attributes from var_in to var_out. Give special treatment to
    _FillValue and coordinates.

    """
    _, _, _, tdim = get_dims_from_variable(var_in.dimensions)
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

            
def copy_global_attributes(in_file, out_file):
    "Copy global attributes from in_file to out_file."
    for attribute in in_file.ncattrs():
        setattr(out_file, attribute, getattr(in_file, attribute))


def create_variable_like(in_file, var_name, out_file, dimensions=None,
                         fill_value=-2e9):
    """Create a variable in an out_file that is the same var_name in
    in_file, except possibly depending on different dimensions,
    provided in dimensions.

    """
    var_in = in_file.variables[var_name]
    try:
        fill_value = var_in._FillValue
    except AttributeError:
        # fill_value was set elsewhere
        pass

    if dimensions is None:
        dimensions = var_in.dimensions

    dtype = var_in.dtype

    var_out = out_file.createVariable(var_name, dtype, dimensions=dimensions,
                                      fill_value=fill_value)
    copy_attributes(var_in, var_out)
    return var_out


if __name__ == "__main__":
    # Set up the option parser
    description = '''A script to extract data along (possibly multiple) profile using
    piece-wise constant or bilinear interpolation.
    The profile must be given as a ESRI shape file.'''
    parser = ArgumentParser()
    parser.description = description

    parser.add_argument("INPUTFILE", nargs=1, help="input NetCDF file name")
    parser.add_argument("OUTPUTFILE", nargs=1, help="output NetCDF file name", default="out.nc")
    parser.add_argument("-n", "--n_levels", dest="n_levels",
                        help="no. of levels",
                        default=11)
    parser.add_argument("-a", "--age_iso", dest="age_iso",
                        help="list of increasing iso age levels",
                        default=[11700, 29000,115000])
    parser.add_argument("-v", "--variable", dest="variables",
                        help="comma-separated list with variables",
                        default='age')

    options = parser.parse_args()
    fill_value = -2e9
    variables = options.variables.split(',')
    n_levels = options.n_levels
    age_iso = np.array(options.age_iso)
    n_age_iso = len(age_iso)
    

    print("-----------------------------------------------------------------")
    print("Running script %s ..." % __file__.split('/')[-1])
    print("-----------------------------------------------------------------")
    print("Opening NetCDF file %s ..." % options.INPUTFILE[0])
    try:
        # open netCDF file in 'read' mode
        nc_in = NC(options.INPUTFILE[0], 'r')
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
               % options.INPUTFILE[0]))
        import sys
        sys.exit()

    # get file format
    format = nc_in.file_format
    # get the dimensions
    xdim, ydim, zdim, tdim = ppt.get_dims(nc_in)
    # read projection information
    projection = ppt.get_projection_from_file(nc_in)
    # new sigma coordinate with n_levels
    z_out = np.linspace(0, 1, n_levels)
    
    nt = len(nc_in.dimensions[tdim])
    nx = len(nc_in.dimensions[xdim])
    ny = len(nc_in.dimensions[ydim])

    # We should bail here if no z-dim is found
    if zdim is not None:
        z = nc_in.variables[zdim][:]

    # use same file format as input file
    print("Creating dimensions")
    nc_out = NC(options.OUTPUTFILE[0], 'w', format=format)
    copy_global_attributes(nc_in, nc_out)


    # re-create dimensions from an input file in an output file, but
    # skip vertical dimension
    copy_dimensions(nc_in, nc_out, zdim)
    # create new zdim
    nc_out.createDimension(zdim, n_levels)
    isodim = 'n_iso'
    nc_out.createDimension(isodim, n_age_iso) 

    out_dims = (tdim, zdim, ydim, xdim)
    
    # copy mapplane dimension variables
    for var_name in (xdim, ydim):
        var_out = create_variable_like(
            nc_in,
            var_name,
            nc_out,
            dimensions=(var_name,))
        var_out[:] = nc_in.variables[var_name][:]
        
    # create new sigma coordinate
    sigma_var = nc_out.createVariable(zdim, 'd', dimensions=(zdim,))
    sigma_var.long_name = 'Sigma-coordinate in Cartesion system'
    sigma_var.axis = 'Z'
    sigma_var.positive = 'up'
    sigma_var[:] = z_out
    
    if tdim is not None:
        copy_time_dimension(nc_in, nc_out, tdim)


    standard_name = 'land_ice_thickness'
    for name in nc_in.variables.keys():
        v = nc_in.variables[name]
        if getattr(v, "standard_name", "") == standard_name:
            print("variabe {0} found by its standard_name {1}".format(name,
                                                                      standard_name))
            myvar = name
            pass


    thickness = ppt.permute(nc_in.variables[myvar], output_order=out_dims)
    thk_min = z[2]

    iso_name = 'depth_iso'
    
    print(("    - reading variable %s" % (myvar)))
        
    print("Copying variables")

    for var_name in variables:

        print("  Reading variable %s" % var_name)
        profiler = timeprofile()
        var_in = nc_in.variables[var_name]
        in_dims = var_in.dimensions
        datatype = var_in.dtype
        profiler.mark('transpose')
        var_in_data = ppt.permute(var_in, output_order=out_dims)
        p = profiler.elapsed('transpose')
        print("    - transposed array in %3.4f s" % p)

        profiler.mark('interpolation')
        if tdim is not None:
            data = np.zeros((nt, n_levels, ny, nx))
            mask = np.ones((nt, n_levels, ny, nx))
            out_var = np.ma.array(data=data, mask=mask, fill_value=fill_value)
            out_var = create_variable_like(nc_in, var_name, nc_out, dimensions=(tdim, zdim, ydim, xdim))
            
            data = np.zeros((nt, n_age_iso, ny, nx))
            mask = np.ones((nt, n_age_iso, ny, nx))
            iso_var = np.ma.array(data=data, mask=mask, fill_value=fill_value)
            iso_var = nc_out.createVariable(iso_name, datatype='double', dimensions=(tdim, isodim, ydim, xdim))

            for t in range(nt):
                for m in range(ny):
                    for n in range(nx):
                        thk = thickness[t,m,n]
                        v = var_in_data[t,:,m,n]
                        if thk > thk_min:
                            z_surf = z[z<thk][-1]
                            z_in = z[z<thk] / z_surf
                            v_in = v[z<thk]
                            f = interp1d(z_in, v_in)
                            v_out = f(z_out)
                            v_out[np.nonzero(v_out<0)] = 0
                            out_var[t,:,m,n] = v_out
                            depth_in = z_surf - z[z<thk]
                            fi = interp1d(v_in, depth_in)
                            vi_out = fi(age_iso)
                            print vi_out
                            
        else:
            data = np.zeros((n_levels, ny, nx))
            mask = np.ones((n_levels, ny, nx))
            out_var = np.ma.array(data=data, mask=mask, fill_value=fill_value)
            out_var = create_variable_like(nc_in, var_name, nc_out, dimensions=(zdim, ydim, xdim))
            data = np.zeros((n_age_iso, ny, nx))
            mask = np.ones((n_age_iso, ny, nx))
            out_var = np.ma.array(data=data, mask=mask, fill_value=fill_value)
            out_var = nc_out.createVariable(iso_name, datatype='double', dimensions=(isodim, ydim, xdim))
            for m in range(ny):
                for n in range(nx):
                    thk = thickness[m,n]
                    v = var_in_data[:,m,n]
                    if thk > thk_min:
                        z_in = z[z<thk] / z[z<thk][-1]
                        v_in = v[z<thk]
                        f = interp1d(z_in, v_in)
                        v_out = f(z_out)
                        v_out[np.nonzero(v_out<0)] = 0
                        out_var[:,m,n] = v_out
        p = profiler.elapsed('interpolation')
        print("    - interpolated in %3.4f s" % p)                     
    
        
    # writing global attributes
    import time, sys
    script_command = ' '.join([time.ctime(), ':', __file__.split('/')[-1],
                               ' '.join([str(l) for l in sys.argv[1:]])])
    if hasattr(nc_in, 'history'):
        history = nc_in.history
        nc_out.history = script_command + '\n ' + history
    else:
        nc_out.history = script_command

    nc_in.close()
    nc_out.close()
    print("Extracted 3D variable(s) {} to file {}".format(variables, options.OUTPUTFILE[0]))