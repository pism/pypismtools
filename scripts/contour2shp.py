#!/usr/bin/env python

import numpy as np
import pylab as plt
from skimage import measure
from argparse import ArgumentParser
from scipy.interpolate import RectBivariateSpline

from netCDF4 import Dataset as NC

import osgeo.ogr
import osgeo.osr
import os

try:
    import PyPISMTools.PyPISMTools as ppt
except:
    import PyPISMTools as ppt


# Shapefile related code is taken from
# http://invisibleroads.com/tutorials/gdal-shapefile-points-save.html

def save(shapePath, geoLocations, proj4):
    '''Save points in the given shapePath'''
    # Get driver
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    # Create shapeData
    shapePath = validateShapePath(shapePath)
    if os.path.exists(shapePath): 
        os.remove(shapePath)
    shapeData = driver.CreateDataSource(shapePath)
    # Create spatialReference
    spatialReference = getSpatialReferenceFromProj4(proj4)
    # Create layer
    layerName = os.path.splitext(os.path.split(shapePath)[1])[0]
    layer = shapeData.CreateLayer(layerName, spatialReference, osgeo.ogr.wkbPoint)
    layerDefinition = layer.GetLayerDefn()
    # For each point,
    for pointIndex, geoLocation in enumerate(geoLocations):
        # Create point
        geometry = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
        geometry.SetPoint(0, geoLocation[0], geoLocation[1])
        # Create feature
        feature = osgeo.ogr.Feature(layerDefinition)
        feature.SetGeometry(geometry)
        feature.SetFID(pointIndex)
        # Save feature
        layer.CreateFeature(feature)
        # Cleanup
        geometry.Destroy()
        feature.Destroy()
    # Cleanup
    shapeData.Destroy()
    # Return
    return shapePath

def load(shapePath):
    '''Given a shapePath, return a list of points in GIS coordinates'''
    # Open shapeData
    shapeData = osgeo.ogr.Open(validateShapePath(shapePath))
    # Validate shapeData
    validateShapeData(shapeData)
    # Get the first layer
    layer = shapeData.GetLayer()
    # Initialize
    points = []
    # For each point,
    for index in xrange(layer.GetFeatureCount()):
        # Get
        feature = layer.GetFeature(index)
        geometry = feature.GetGeometryRef()
        # Make sure that it is a point
        if geometry.GetGeometryType() != osgeo.ogr.wkbPoint: 
            raise ShapeDataError('This module can only load points; use geometry_store.py')
        # Get pointCoordinates
        pointCoordinates = geometry.GetX(), geometry.GetY()
        # Append
        points.append(pointCoordinates)
        # Cleanup
        feature.Destroy()
    # Get spatial reference as proj4
    proj4 = layer.GetSpatialRef().ExportToProj4()
    # Cleanup
    shapeData.Destroy()
    # Return
    return points, proj4

def getSpatialReferenceFromProj4(proj4):
    '''Return GDAL spatial reference object from proj4 string'''
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromProj4(proj4)
    return spatialReference

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


parser = ArgumentParser(description='''A script to extract a (closed) contour line from a variable in a netCDF file.''')
parser.add_argument("FILE", nargs=1)
parser.add_argument("-o", "--output_filename", dest="out_file",
                  help="Name of the output shape file", default='countour.shp')
parser.add_argument("-v", "--variable", dest="varname",
                  help='''Variable to plot, default = 'mask'.''', default='mask')
parser.add_argument("-c", "--countour_level", dest="contour_level",
                  help='''Contour-level to extract, default = 0.''', default=0.)


options = parser.parse_args()
filename = options.FILE[0]
shp_filename = options.out_file
contour_level = options.contour_level
varname = options.varname

nc = NC(filename, 'r')

xdim, ydim = ppt.get_dims(nc)[0:2]

x_var = np.squeeze(nc.variables[xdim])
y_var = np.squeeze(nc.variables[ydim])
i = range(0, len(x_var))
j = range(0, len(y_var))

xx, yy = np.meshgrid(x_var, y_var)

contour_var = np.array(np.squeeze(ppt.permute(nc.variables[varname])), order='C')

nc_projection = ppt.get_projection_from_file(nc)

contour_level = 0
# Find contours at a constant value
contours = sorted(measure.find_contours(contour_var, contour_level),
                  key=lambda x: len(x))
# We only extract the longest contigous contour
contour = contours[-1]

contour_x = x_var[0] + contour[:,1] * (x_var[-1] - x_var[0]) / (len(i) - 1)
contour_y = y_var[0] + contour[:,0] * (y_var[-1] - y_var[0]) / (len(j) - 1)
# save shapefile
print(("Saving shapefile %s" % shp_filename))
points = [(contour_x[k], contour_y[k]) for k in range(len(contour_y))]
save(shp_filename, points, nc_projection.srs)

nc.close()
