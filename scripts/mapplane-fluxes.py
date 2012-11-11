#!/usr/bin/env python

import numpy as np
import pylab as plt
from skimage import measure
from argparse import ArgumentParser

from netCDF4 import Dataset as NC

try:
    import PyPISMTools.PyPISMTools as ppt
except:
    import PyPISMTools as ppt

# Set up the argument parser
parser = ArgumentParser()
parser.description = '''A script for PISM output files to make mass time series
plots using pylab/matplotlib'''
parser.add_argument("FILE", nargs='*')

options = parser.parse_args()
args = options.FILE
nt = len(args)
var_name = 'ocean_kill_flux_cumulative'
outunit = 'Gt'

values = []
record = -1
no_cells = []
cum_lengths = []
cum_areas = []
for k in range(0, nt):

    file_name = args[k]
    print("  opening NetCDF file %s ..." % file_name)
    try:
        # open netCDF file in 'append' mode
        nc = NC(file_name, 'r')
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
              % file_name))
        import sys
        sys.exit(1)

    # get the dimensions
    xdim, ydim, zdim, tdim = ppt.get_dims(nc)
    # set up dimension ordering
    dim_order = (tdim, zdim, ydim, xdim)
    # add lat/lon values
    x = (np.squeeze(ppt.permute(nc.variables[xdim], dim_order)))
    x_units = nc.variables[xdim].units
    y = (np.squeeze(ppt.permute(nc.variables[ydim], dim_order)))
    y_units = nc.variables[ydim].units

    var = 'topg'
    print(("    - reading variable %s from file %s" % (var, file_name)))
    try:
        topg = np.squeeze(ppt.permute(nc.variables[var], dim_order))
        topg_units = nc.variables[var].units
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
              % (var, file_name)))
        import sys
        sys.exit(1)
        
    topg = ppt.unit_converter(topg, topg_units, 'm')
    mask = (topg >= 0)
    topg = np.ma.array(topg, mask = mask)

    var = 'thk'
    print(("    - reading variable %s from file %s" % (var, file_name)))
    try:
        thk = np.squeeze(ppt.permute(nc.variables[var], dim_order))
        thk_units = nc.variables[var].units
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
              % (var, file_name)))
        import sys
        sys.exit(1)
        
    thk = ppt.unit_converter(thk, topg_units, 'm')

    speed_units = 'm year-1'
    
    var = 'ubar'
    print(("    - reading variable %s" % var))
    try:
        ubar = np.squeeze(ppt.permute(nc.variables[var], dim_order))
        ubar_units = nc.variables[var].units
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
              % (var, file_name)))
        import sys
        sys.exit(1)
        
    ubar = ppt.unit_converter(ubar, ubar_units, speed_units)

    var = 'vbar'
    print(("    - reading variable %s" % var))
    try:
        vbar = np.squeeze(ppt.permute(nc.variables[var], dim_order))
        vbar_units = nc.variables[var].units
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
              % (var, file_name)))
        import sys
        sys.exit(1)
        
    vbar = ppt.unit_converter(vbar, vbar_units, speed_units)

    print(("    - reading variable %s" % var_name))
    try:
        data = np.squeeze(ppt.permute(nc.variables[var_name], dim_order))
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
              % (var_name, file_name)))
        import sys
        sys.exit(1)

    try:
        inunit = str(nc.variables[var_name].units)
    except:
        print(("ERROR:  units not found in variable '%s' in file %s ... ending ..."
              % (var_name, file_name)))
        import sys
        sys.exit(1)

    if outunit is not None:
        data = ppt.unit_converter(data, inunit, outunit)

    mask = (data >= 0)
    data = np.ma.array(data, mask = mask)

    outdimunits = 'm'
    dx = ppt.unit_converter(np.abs(x[1] - x[0]), x_units, outdimunits)
    dy = ppt.unit_converter(np.abs(y[1] - y[0]), y_units, outdimunits)

    velbar = np.sqrt(ubar ** 2 + vbar ** 2)
    
    # get number of non-zero non-masked cells
    n_cells = data[data.nonzero()].data.shape[0]
    no_cells.append(n_cells)
    # calculate cumulative length of discharge cells depending on
    # grid resolution
    cum_length = n_cells * dx  # in m
    cum_lengths.append(cum_length)
    cum_thickness = np.abs(topg[data.nonzero()].sum())
    cum_area = cum_length * cum_thickness
    cum_areas.append(cum_area)
    values.append(data)

    nc.close()
