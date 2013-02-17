#!/usr/bin/env python

import numpy as np
import pylab as plt
from skimage import measure
from argparse import ArgumentParser

from netCDF4 import Dataset as NC

import osgeo.ogr
import osgeo.osr
import os
from pyproj import Proj

try:
    import PyPISMTools.PyPISMTools as ppt
except:
    import PyPISMTools as ppt


# Shapefile related code is adapted from
# http://invisibleroads.com/tutorials/gdal-shapefile-points-save.html

def save(shapePath, contour_points):
    '''Save points in the given shapePath'''
    # Get driver
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    # Create shapeData
    shapePath = validateShapePath(shapePath)
    if os.path.exists(shapePath): 
        os.remove(shapePath)
    shapeData = driver.CreateDataSource(shapePath)
    # Create spatialReference, EPSG 4326 (lonlat)
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromEPSG(4326)
    layerName = os.path.splitext(os.path.split(shapePath)[1])[0]
    layer = shapeData.CreateLayer(layerName, spatialReference, osgeo.ogr.wkbPolygon)
    layerDefinition = layer.GetLayerDefn()
    # For each contour,
    polygon = osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
    for k in range(0,len(contour_points)):
        geoLocations = contour_points[k]
        ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)
        # For each point,
        for pointIndex, geoLocation in enumerate(geoLocations):
            ring.AddPoint(geoLocation[0], geoLocation[1])
        ring.CloseRings()
        polygon.AddGeometry(ring)
    featureDefn = layer.GetLayerDefn()
    # Create feature
    feature = osgeo.ogr.Feature(featureDefn)
    feature.SetGeometry(polygon)
    feature.SetFID(k)
    # Save feature
    layer.CreateFeature(feature)
    # Cleanup
    polygon.Destroy()
    feature.Destroy()
    # Cleanup
    shapeData.Destroy()
    # Return
    return shapePath


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


parser = ArgumentParser(description='''A script to extract a (closed) contour line from a variable in a netCDF file, and save it as a shapefile (polygon).''')
parser.add_argument("FILE", nargs=1)
parser.add_argument("-o", "--output_filename", dest="out_file",
                  help="Name of the output shape file", default='countour.shp')
parser.add_argument("-v", "--variable", dest="varname",
                  help='''Variable to plot, default = 'mask'.''', default='mask')
parser.add_argument("-c", "--countour_level", dest="contour_level",
                  help='''Contour-level to extract, default = 0.''', default=0.)
parser.add_argument("-s","--single",dest="single", action="store_true",
                  help="save only the longest contour line, Default=False", default=False)


options = parser.parse_args()
filename = options.FILE[0]
shp_filename = options.out_file
contour_level = options.contour_level
varname = options.varname
single = options.single

nc = NC(filename, 'r')

xdim, ydim, zdim, tdim = ppt.get_dims(nc)
var_order = (tdim, zdim, ydim, xdim)

x_var = np.squeeze(nc.variables[xdim])
y_var = np.squeeze(nc.variables[ydim])
i = range(0, len(x_var))
j = range(0, len(y_var))

contour_var = np.array(np.squeeze(ppt.permute(nc.variables[varname], var_order)), order='C')
nc_projection = ppt.get_projection_from_file(nc)

# Find contours at a constant value
contours = sorted(measure.find_contours(contour_var, contour_level),
                  key=lambda x: len(x))

contours_x = []
contours_y = []
contours_lon = []
contours_lat = []
contour_points = []
for k in range(0,len(contours)):
    contour = contours[k]
    contour_x = x_var[0] + contour[:,1] * (x_var[-1]-x_var[0])/(len(i)-1)
    contour_y = y_var[0] + contour[:,0] * (y_var[-1]-y_var[0])/(len(j)-1)
    contour_lon, contour_lat = nc_projection(contour_x, contour_y, inverse=True)
    contours_x.append(contour_x)
    contours_y.append(contour_y)
    contours_lon.append(contour_lon)
    contours_lat.append(contour_lat)
    points = [(contour_lon[k], contour_lat[k]) for k in range(len(contour_lat))]
    contour_points.append(points)
# reverse direction, last entry (longest contour) first.
contours_x.reverse()
contours_y.reverse()
contours_lon.reverse()
contours_lat.reverse()
contour_points.reverse()

if single:
    contour_points = [contour_points[0]]
    
# save shapefile
print(("Saving shapefile %s" % shp_filename))
save(shp_filename, contour_points)


# display result
fig = plt.figure()
ax = fig.add_subplot(111)
xx, yy = np.meshgrid(x_var, y_var)
ax.pcolormesh(xx, yy, contour_var, cmap=plt.cm.Blues_r)
for k in range(0,len(contour_points)):
    plt.plot(contours_x[k], contours_y[k], linewidth=2, color='r')
plt.axis('image')
plt.xticks([])
plt.yticks([])
plt.show()

nc.close()
