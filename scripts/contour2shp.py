#!/usr/bin/env python

import numpy as np
import pylab as plt
from skimage import measure
from argparse import ArgumentParser

from netCDF4 import Dataset as NC
from netcdftime import utime
import ogr
import osr
import os
from pyproj import Proj

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt


def validateShapePath(shapePath):
    '''Validate shapefile extension'''
    return os.path.splitext(str(shapePath))[0] + '.shp'


def validateShapeData(shapeData):
    '''Make sure we can access the shapefile'''
    # Make sure the shapefile exists
    if not shapeData:
        raise ShapeDataError('The shapefile is invalid')
    # Make sure there is exactly one layer
    if shapeData.GetLayerCount() != 1:
        raise ShapeDataError('The shapefile must have exactly one layer')


# Error
class ShapeDataError(Exception):
    pass


def get_contours(array, x, y, projection, level):
    '''
    Find contours for a given level
    '''

    # Find contours at a constant value
    contours = sorted(measure.find_contours(array, level),
                      key=lambda x: len(x))

    i = range(0, len(x))
    j = range(0, len(y))

    lon = []
    lat = []
    contour_points = []
    for k in range(0, len(contours)):
        contour = contours[k]
        contour_x = x[0] + contour[:,1] * (x[-1]-x[0])/(len(i)-1)
        contour_y = y[0] + contour[:,0] * (y[-1]-y[0])/(len(j)-1)
        # Convert to EPSG:4326
        contour_lon, contour_lat = projection(contour_x, contour_y, inverse=True)
        lon.append(contour_lon)
        lat.append(contour_lat)
        points = [(contour_lon[k], contour_lat[k]) for k in range(len(contour_lat))]
        contour_points.append(points)
    # reverse direction, last entry (longest contour) first.
    contour_points.reverse()

    if single:
        contour_points = [contour_points[0]]

    return contour_points



parser = ArgumentParser(description='''A script to extract a (closed) contour line from a variable in a netCDF file, and save it as a shapefile (polygon).''')
parser.add_argument("FILE", nargs=1)
parser.add_argument("-o", "--output_filename", dest="out_file",
                  help="Name of the output shape file", default='countour.shp')
parser.add_argument("-v", "--variable", dest="varname",
                  help='''Variable to plot, default = 'mask'.''', default='mask')
parser.add_argument("-c", "--countour_levels", nargs='*',
                    dest="contour_levels",
                    help='''Contour-levels to extract, default = 0.''', default='0')
parser.add_argument("-s","--single",dest="single", action="store_true",
                  help="save only the longest contour line, Default=False", default=False)


options = parser.parse_args()
filename = options.FILE[0]
shp_filename = options.out_file
contour_levels = options.contour_levels
varname = options.varname
single = options.single

nc = NC(filename, 'r')
nc_projection = ppt.get_projection_from_file(nc)

xdim, ydim, zdim, tdim = ppt.get_dims(nc)
var_order = (tdim, zdim, ydim, xdim)

x = np.squeeze(nc.variables[xdim])
y = np.squeeze(nc.variables[ydim])


# Get driver
driver = ogr.GetDriverByName('ESRI Shapefile')
# Create shapeData
shp_filename = validateShapePath(shp_filename)
if os.path.exists(shp_filename): 
    os.remove(shp_filename)
shapeData = driver.CreateDataSource(shp_filename)
# Create spatialReference, EPSG 4326 (lonlat)
spatialReference = osr.SpatialReference()
spatialReference.ImportFromEPSG(4326)
layerName = os.path.splitext(os.path.split(shp_filename)[1])[0]
layer = shapeData.CreateLayer(layerName, spatialReference, ogr.wkbPolygon)
layerDefinition = layer.GetLayerDefn()
field_defn = ogr.FieldDefn("level", ogr.OFTReal)
layer.CreateField(field_defn)
field_defn = ogr.FieldDefn("timestamp", ogr.OFTDateTime)
layer.CreateField(field_defn)

if tdim:
    time = nc.variables['time']
    time_units = time.units
    time_calendar = time.calendar
    cdftime = utime(time_units, time_calendar)
    for k, t in enumerate(time):
        timestamp = cdftime.num2date(t)
        print('Processing {}'.format(timestamp))
        for level in contour_levels:
            contour_var = np.array(ppt.permute(nc.variables[varname], var_order), order='C')[k, Ellipsis]
            contour_points = get_contours(contour_var, x, y, nc_projection, level)
            # For each contour
            polygon = ogr.Geometry(ogr.wkbPolygon)
            for k in range(0,len(contour_points)):
                geoLocations = contour_points[k]
                ring = ogr.Geometry(ogr.wkbLinearRing)
                # For each point,
                for pointIndex, geoLocation in enumerate(geoLocations):
                    ring.AddPoint(geoLocation[0], geoLocation[1])
                ring.CloseRings()
                polygon.AddGeometry(ring)
            # Create feature
            featureDefn = layer.GetLayerDefn()
            feature = ogr.Feature(featureDefn)
            feature.SetGeometry(polygon)
            feature.SetFID(k)
            i = feature.GetFieldIndex("level")
            feature.SetField(i, level)
            i = feature.GetFieldIndex("timestamp")
            feature.SetField(i, str(timestamp))
            polygon = None
            # Save feature
            layer.CreateFeature(feature)
            # Cleanup
            feature = None
else:
    for level in contour_levels:
        contour_var = np.array(np.squeeze(ppt.permute(nc.variables[varname], var_order)), order='C')
        contour_points = get_contours(contour_var, x, y, nc_projection, level)
        # For each contour
        polygon = ogr.Geometry(ogr.wkbPolygon)
        for k in range(0,len(contour_points)):
            geoLocations = contour_points[k]
            ring = ogr.Geometry(ogr.wkbLinearRing)
            # For each point,
            for pointIndex, geoLocation in enumerate(geoLocations):
                ring.AddPoint(geoLocation[0], geoLocation[1])
            ring.CloseRings()
            polygon.AddGeometry(ring)
        # Create feature
        featureDefn = layer.GetLayerDefn()
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(polygon)
        feature.SetFID(k)
        i = feature.GetFieldIndex("level")
        feature.SetField(i, level)
        polygon = None
        # Save feature
        layer.CreateFeature(feature)
        # Cleanup
        feature = None
# Cleanup
shapeData = None
    


# save(shp_filename, contour_points, level)


nc.close()
