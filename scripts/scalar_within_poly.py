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
parser.add_argument(
    "-i",
    "--invert",
    dest="invert",
    action="store_true",
    help="Replace all values outside of the polygon with this value",
    default=False,
)
parser.add_argument(
    "-s", "--scalar_value", dest="scalar_value", type=float, help="Replace with this value", default=0.0
)
parser.add_argument(
    "-v", "--variables", dest="variables", help="Comma separated list of variables.", default=["bmelt"]
)

options = parser.parse_args()
args = options.FILE
scalar_value = options.scalar_value
variables = options.variables.split(",")
invert = options.invert

driver = ogr.GetDriverByName("ESRI Shapefile")
data_source = driver.Open(args[0], 0)
if data_source is None:
    print("Couldn't open file {}.\n".format(args[0]))
    import sys

    sys.exit(1)
layer = data_source.GetLayer(0)
srs = layer.GetSpatialRef()
if not srs.IsGeographic():
    print(("""Spatial Reference System in % s is not latlon. Converting.""" % filename))
    # Create spatialReference, EPSG 4326 (lonlat)
    srs_geo = osr.SpatialReference()
    srs_geo.ImportFromEPSG(4326)

nc = NC(args[1], "a")

var = "lat"
try:
    lat = nc.variables[var]
except:
    print(("ERROR:  variable '%s' not found but needed... ending ..." % var))
    import sys

    sys.exit()

var = "lon"
try:
    lon = nc.variables[var]
except:
    print(("ERROR:  variable '%s' not found but needed... ending ..." % var))
    import sys

    sys.exit()

# get dimensions of first variables
var = variables[0]
try:
    first_var = nc.variables[var]
except:
    print(("ERROR:  variable '%s' not found but needed... ending ..." % var))

for feature in layer:
    feature = layer.GetFeature(0)
    geometry = feature.GetGeometryRef()
    # Transform to latlon if needed
    if not srs.IsGeographic():
        geometry.TransformTo(srs_geo)

    counter = 0
    ndim = first_var.ndim

    stderr.write("\n  - Processing variable %s, precent done: " % var)
    stderr.write("000")

    if ndim == 2:
        M = first_var.shape[0]
        N = first_var.shape[1]
        max_counter = M * N
        for m in range(0, M):
            for n in range(0, N):
                x = lon[m, n]
                y = lat[m, n]
                wkt = "POINT(%f %f)" % (x, y)
                point = ogr.CreateGeometryFromWkt(wkt)
                if invert:
                    if feature.GetGeometryRef().Contains(point):
                        pass
                    else:
                        for var in variables:
                            try:
                                data = nc.variables[var]
                            except:
                                print(("ERROR:  variable '%s' not found but needed... ending ..." % var))
                                import sys

                                sys.exit()
                            data[m, n] = scalar_value
                else:
                    if feature.GetGeometryRef().Contains(point):
                        for var in variables:
                            try:
                                data = nc.variables[var]
                            except:
                                print(("ERROR:  variable '%s' not found but needed... ending ..." % var))
                                import sys

                                sys.exit()
                                data[m, n] = scalar_value

                stderr.write("\b\b\b%03d" % (100.0 * counter / max_counter))
                counter += 1

    elif ndim == 3:
        K = data.shape[0]
        M = data.shape[1]
        N = data.shape[2]
        max_counter = K * M * N
        for k in range(0, K):
            for m in range(0, M):
                for n in range(0, N):
                    x = lon[m, n]
                    y = lat[m, n]
                    wkt = "POINT(%f %f)" % (x, y)
                    point = ogr.CreateGeometryFromWkt(wkt)
                    if invert:
                        if feature.GetGeometryRef().Contains(point):
                            pass
                        else:
                            for var in variables:
                                try:
                                    data = nc.variables[var]
                                except:
                                    print(("ERROR:  variable '%s' not found but needed... ending ..." % var))
                                    import sys

                                    sys.exit()
                                data[k, m, n] = scalar_value
                    else:
                        if feature.GetGeometryRef().Contains(point):
                            for var in variables:
                                try:
                                    data = nc.variables[var]
                                except:
                                    print(("ERROR:  variable '%s' not found but needed... ending ..." % var))
                                    import sys

                                    sys.exit()
                                data[k, m, n] = scalar_value

                    stderr.write("\b\b\b%03d" % (100.0 * counter / max_counter))
                    counter += 1
    else:
        print(("ERROR: %i dimensions currently not supported... ending..." % ndim))
nc.close()
