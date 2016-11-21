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
                    choices=[
                        'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT'],
                    help="file format out output file", default='netcdf3_64bit')

options = parser.parse_args()
args = options.FILE
grid_spacing = options.grid_spacing  # convert

fileformat = options.fileformat.upper()

if len(args) == 0:
    nc_outfile = 'jif' + str(grid_spacing) + 'm.nc'
elif len(args) == 1:
    nc_outfile = args[0]
else:
    print('wrong number arguments, 0 or 1 arguments accepted')
    parser.print_help()
    import sys
    sys.exit(0)


if __name__ == "__main__":

    xdim = 'x'
    ydim = 'y'
    
    # define output grid, these are the extents of Mathieu's domain (cell
    # corners)
    e0 = 403000
    n0 = 6370000
    e1 = 683000
    n1 = 6691000

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

    # UTM 8N projection: EPSG 32608
    projection = "+init=epsg:32608"
    proj = Proj(projection)

    lon, lat = proj(ee, nn, inverse=True)

    # number of grid corners
    grid_corners = 4
    # grid corner dimension name
    grid_corner_dim_name = "nv4"

    # array holding x-component of grid corners
    gc_easting = np.zeros((M, grid_corners))
    # array holding y-component of grid corners
    gc_northing = np.zeros((N, grid_corners))
    # array holding the offsets from the cell centers
    # in x-direction (counter-clockwise)
    de_vec = np.array([-de / 2, de / 2, de / 2, -de / 2])
    # array holding the offsets from the cell centers
    # in y-direction (counter-clockwise)
    dn_vec = np.array([-dn / 2, -dn / 2, dn / 2, dn / 2])
    # array holding lat-component of grid corners
    gc_lat = np.zeros((N, M, grid_corners))
    # array holding lon-component of grid corners
    gc_lon = np.zeros((N, M, grid_corners))
    
    for corner in range(0, grid_corners):
        ## grid_corners in x-direction
        gc_easting[:, corner] = easting + de_vec[corner]
        # grid corners in y-direction
        gc_northing[:, corner] = northing + dn_vec[corner]
        # meshgrid of grid corners in x-y space
        gc_ee, gc_nn = np.meshgrid(
            gc_easting[:, corner], gc_northing[:, corner])
        # project grid corners from x-y to lat-lon space
        gc_lon[:, :, corner], gc_lat[:, :, corner] = proj(
            gc_ee, gc_nn, inverse=True)


    nc = CDF(nc_outfile, 'w', format=fileformat)

    nc.createDimension(xdim, size=easting.shape[0])
    nc.createDimension(ydim, size=northing.shape[0])
    
    var = xdim
    var_out = nc.createVariable(var, 'f', dimensions=(xdim))
    var_out.axis = xdim
    var_out.long_name = "X-coordinate in Cartesian system"
    var_out.standard_name = "projection_x_coordinate"
    var_out.units = "meters"
    var_out[:] = easting

    var = ydim
    var_out = nc.createVariable(var, 'f', dimensions=(ydim))
    var_out.axis = ydim
    var_out.long_name = "Y-coordinate in Cartesian system"
    var_out.standard_name = "projection_y_coordinate"
    var_out.units = "meters"
    var_out[:] = northing

    var = 'lon'
    var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
    var_out.units = "degrees_east"
    var_out.valid_range = -180., 180.
    var_out.standard_name = "longitude"
    var_out.bounds = "lon_bnds"
    var_out[:] = lon

    var = 'lat'
    var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
    var_out.units = "degrees_north"
    var_out.valid_range = -90., 90.
    var_out.standard_name = "latitude"
    var_out.bounds = "lat_bnds"
    var_out[:] = lat


    nc.createDimension(grid_corner_dim_name, size=grid_corners)

    var = 'lon_bnds'
    # Create variable 'lon_bnds'
    var_out = nc.createVariable(
        var, 'f', dimensions=(ydim, xdim, grid_corner_dim_name))
    # Assign units to variable 'lon_bnds'
    var_out.units = "degreesE"
    # Assign values to variable 'lon_nds'
    var_out[:] = gc_lon
        
    var = 'lat_bnds'
    # Create variable 'lat_bnds'
    var_out = nc.createVariable(
        var, 'f', dimensions=(ydim, xdim, grid_corner_dim_name))
    # Assign units to variable 'lat_bnds'
    var_out.units = "degreesN"
    # Assign values to variable 'lat_bnds'
    var_out[:] = gc_lat
    
    var = 'dummy'
    var_out = nc.createVariable(
        var,
        'f',
        dimensions=(
            ydim,
            xdim),
        fill_value=np.nan)
    var_out.units = "meters"
    var_out.long_name = "Just A Dummy"
    var_out.comment = "This is just a dummy variable for CDO."
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    var_out[:] = 1.

    mapping = nc.createVariable("mapping", 'c')
    mapping.inverse_flattening = 298.257
    mapping.utm_zone_number = 8
    mapping.semi_major_axis = 6378137
    mapping.grid_mapping_name = "universal_transverse_mercator"
    mapping._CoordinateTransformType = "Projection"
    mapping._CoordinateAxisTypes = "GeoX GeoY"

    from time import asctime
    historystr = 'Created ' + asctime() + '\n'
    nc.history = historystr
    nc.proj4 = projection
    nc.Conventions = 'CF-1.5'
    nc.close()
