#!/usr/bin/env python

# Copyright (C) 2013 Andy Aschwanden

from sys import stderr
from argparse import ArgumentParser
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from osgeo import ogr

# Set up the option parser
parser = ArgumentParser()
parser.description = "All values within a polygon defined by a shapefile are replaced by a scalar value."
parser.add_argument("FILE", nargs=2)
parser.add_argument("-s", "--scalar_value",dest="scalar_value", type=float,
                  help="Replace with this value",default=0.)
parser.add_argument("-v", "--variables",dest="variables",
                  help="Comma separated list of variables.",default=['bmelt'])

options = parser.parse_args()
args = options.FILE
scalar_value = options.scalar_value
variables = options.variables.split(',')

driver = ogr.GetDriverByName('ESRI Shapefile')
data_source = driver.Open(args[0], 0)
layer = data_source.GetLayer(0)
srs=layer.GetSpatialRef()
# Make sure we use lat/lon coordinates.
# Fixme: allow reprojection onto lat/lon if needed.
if not srs.IsGeographic():
    print('''Spatial Reference System in % s is not lat/lon. Exiting.'''
          % filename)
    import sys
    sys.exit(0)
feature = layer.GetFeature(0)

nc = NC(args[1], 'a')

var = 'lat'
try:
    lat = nc.variables[var]
except:
    print(("ERROR:  variable '%s' not found but needed... ending ..."
              % var))
    import sys
    sys.exit()

var = 'lon'
try:
    lon = nc.variables[var]
except:
    print(("ERROR:  variable '%s' not found but needed... ending ..."
              % var))
    import sys
    sys.exit()


for var in variables:
    try:
        data = nc.variables[var]
    except:
        print(("ERROR:  variable '%s' not found but needed... ending ..."
                  % var))
        import sys
        sys.exit()

    counter = 0
    ndim = data.ndim

    stderr.write("\n  - Processing variable %s, precent done: " % var)
    stderr.write("000")

    if (ndim==2):
        M = data.shape[0]
        N = data.shape[1]
        max_counter = M*N
        for m in range(0, M):
            for n in range(0, N):
                x = lon[m,n]
                y = lat[m,n]
                wkt = "POINT(%f %f)" % (x,y)
                point = ogr.CreateGeometryFromWkt(wkt)
                if feature.GetGeometryRef().Contains(point):
                    data[m,n] = scalar_value
                stderr.write("\b\b\b%03d" % (100.0 * counter / max_counter))
                counter += 1
        
    elif (ndim==3):
        K = data.shape[0]
        M = data.shape[1]
        N = data.shape[2]
        max_counter = K*M*N
        for k in range(0, K):
            for m in range(0, M):
                for n in range(0, N):
                    x = lon[m,n]
                    y = lat[m,n]
                    wkt = "POINT(%f %f)" % (x,y)
                    point = ogr.CreateGeometryFromWkt(wkt)
                    if feature.GetGeometryRef().Contains(point):
                        data[k,m,n] = scalar_value
                    stderr.write("\b\b\b%03d" % (100.0 * counter / max_counter))
                    counter += 1
    else:
        print(("ERROR: %i dimensions currently not supported... ending..."
               % ndim))
nc.close()
