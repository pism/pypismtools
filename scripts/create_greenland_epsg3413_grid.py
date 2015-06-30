#!/usr/bin/env python
import numpy as np
from pyproj import Proj
from argparse import ArgumentParser

from netCDF4 import Dataset as CDF

# set up the argument parser
parser = ArgumentParser()
parser.description = "Create CDO-compliant grid description"
parser.add_argument("FILE", nargs='*')
parser.add_argument("-g", "--grid_spacing", dest="grid_spacing", type=float,
                    help="use X m grid spacing", default=1800)
parser.add_argument("-f", "--format", dest="fileformat", type=str.upper,
                  choices = ['NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT'],
                  help="file format out output file", default='netcdf3_64bit')

options = parser.parse_args()
args = options.FILE
grid_spacing = options.grid_spacing  # convert

fileformat = options.fileformat.upper()

if len(args) == 0:
    nc_outfile = 'grn' + str(grid_spacing) + 'm.nc'
elif len(args) == 1:
    nc_outfile = args[0]
else:
    print('wrong number arguments, 0 or 1 arguments accepted')
    parser.print_help()
    import sys
    sys.exit(0)


if __name__ == "__main__":

    # define output grid, these are the extents of Mathieu's domain (cell
    # corners)
    e0 = -638000
    n0 = -3349600
    e1 = 864700
    n1 = -657600

    # Add a buffer on each side such that we get nice grids up to a grid spacing
    # of 36 km.

    buffer_e = 40650
    buffer_n = 22000
    e0 -= buffer_e
    n0 -= buffer_n
    e1 += buffer_e
    n1 += buffer_n

    # Shift to cell centers
    e0 += grid_spacing / 2
    n0 += grid_spacing / 2
    e1 -= grid_spacing / 2
    n1 -= grid_spacing / 2

    de = dn = grid_spacing  # m
    M = int((e1 - e0) / de) + 1
    N = int((n1 - n0) / dn) + 1

    easting = np.linspace(e0, e1, M)
    northing = np.linspace(n0, n1, N)
    ee, nn = np.meshgrid(easting, northing)

    # Set up SeaRISE Projection
    projection = "+init=epsg:3413"
    proj = Proj(projection)

    lon, lat = proj(ee, nn, inverse=True)

    nc = CDF(nc_outfile, 'w', format=fileformat)
        
    nc.createDimension("x", size=easting.shape[0])
    nc.createDimension("y", size=northing.shape[0])

    var = 'x'
    var_out = nc.createVariable(var, 'f', dimensions=("x"))
    var_out.axis = "X"
    var_out.long_name = "X-coordinate in Cartesian system"
    var_out.standard_name = "projection_x_coordinate"
    var_out.units = "meters"
    var_out[:] = easting

    var = 'y'
    var_out = nc.createVariable(var, 'f', dimensions=("y"))
    var_out.axis = "Y"
    var_out.long_name = "Y-coordinate in Cartesian system"
    var_out.standard_name = "projection_y_coordinate"
    var_out.units = "meters"
    var_out[:] = northing

    var = 'lon'
    var_out = nc.createVariable(var, 'f', dimensions=("y", "x"))
    var_out.units = "degrees_east"
    var_out.valid_range = -180., 180.
    var_out.standard_name = "longitude"
    var_out[:] = lon

    var = 'lat'
    var_out = nc.createVariable(var, 'f', dimensions=("y", "x"))
    var_out.units = "degrees_north"
    var_out.valid_range = -90., 90.
    var_out.standard_name = "latitude"
    var_out[:] = lat

    var = 'dummy'
    var_out = nc.createVariable(
        var,
        'f',
        dimensions=(
            "y",
            "x"),
        fill_value=np.nan)
    var_out.units = "meters"
    var_out.long_name = "Just A Dummy"
    var_out.comment = "This is just a dummy variable for CDO."
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    var_out[:] = np.nan

    mapping = nc.createVariable("mapping", 'c')
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.
    mapping.false_northing = 0.
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = 90.
    mapping.standard_parallel = 71.
    mapping.straight_vertical_longitude_from_pole = -39.

    from time import asctime
    historystr = 'Created ' + asctime() + '\n'
    nc.history = historystr
    nc.proj4 = projection
    nc.Conventions = 'CF-1.5'
    nc.close()
