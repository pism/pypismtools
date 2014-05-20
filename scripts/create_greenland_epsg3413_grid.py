#!/usr/bin/env python
import numpy as np
from pyproj import Proj
from optparse import OptionParser

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

# default values
DXY = 2000.  # km

# set up the option parser
parser = OptionParser()
parser.usage = "usage: %prog [options] FILE"
parser.description = "Create CDO-compliant grid description"
parser.add_option("-g", "--grid_spacing",dest="grid_spacing",type='float',
                  help="use X m grid spacing",
                  metavar="X",default=DXY)

(options, args) = parser.parse_args()
grid_spacing = options.grid_spacing # convert

if len(args) == 0:
    nc_outfile = 'grn' + str(grid_spacing/1e3) + 'km.nc'
elif len(args) == 1:
    nc_outfile = args[0]
else:
    print('wrong number arguments, 0 or 1 arguments accepted')
    parser.print_help()
    exit(0)


if __name__ == "__main__": 

    # define output grid

    # Create a buffer that is a multiple of the grid resolution
    buffer=15000
    cell_center_shift = 75
    e0 = -638000 - buffer - cell_center_shift
    n0 = -3349600 - buffer - cell_center_shift
    e1 = 865000 + buffer - cell_center_shift
    n1 =-657600 +  buffer - cell_center_shift

    de = dn =  grid_spacing # m
    M = int((e1 - e0)/de) + 1
    N = int((n1 - n0)/dn) + 1

    easting  = np.linspace(e0, e1, M)
    northing = np.linspace(n0, n1, N)
    ee, nn = np.meshgrid(easting,northing)


    # Set up SeaRISE Projection
    projection = "+init=epsg:3413"
    proj = Proj(projection)

    lon,lat = proj(ee,nn,inverse=True)

    nc = CDF(nc_outfile,'w',format='NETCDF3_CLASSIC')

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
    var_out.units = "meters";
    var_out[:] = northing

    var = 'lon'
    var_out = nc.createVariable(var, 'f', dimensions=("y","x"))
    var_out.units = "degrees_east";
    var_out.valid_range = -180., 180.
    var_out.standard_name = "longitude"
    var_out[:] = lon

    var = 'lat'
    var_out = nc.createVariable(var, 'f', dimensions=("y","x"))
    var_out.units = "degrees_north";
    var_out.valid_range = -90., 90.
    var_out.standard_name = "latitude"
    var_out[:] = lat
      
    var = 'dummy'
    var_out = nc.createVariable(var, 'f', dimensions=("y","x"), fill_value=np.nan)
    var_out.units = "meters";
    var_out.long_name = "Just A Dummy"
    var_out.comment = "This is just a dummy variable for CDO."
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"    
    var_out[:] = np.nan

    mapping = nc.createVariable("mapping",'c')
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
