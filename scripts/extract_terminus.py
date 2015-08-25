#!/usr/bin/env python

import numpy as np
import pylab as plt
from skimage import measure
from argparse import ArgumentParser

from netCDF4 import Dataset as NC
from netcdftime import utime
import gdal
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



options = parser.parse_args()
filename = options.FILE[0]
shp_filename = options.out_file

dst_fieldname = 'mask'

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

time = nc.variables['time']
time_units = time.units
time_calendar = time.calendar
cdftime = utime(time_units, time_calendar)

src_ds = gdal.Open('NETCDF:{}:{}'.format(filename, dst_fieldname))

mem_driver = ogr.GetDriverByName('Memory')
mem_driver = ogr.GetDriverByName('ESRI Shapefile')
mem_ds = mem_driver.CreateDataSource('memory_layer')

try:
    poly_layer = mem_ds.GetLayerByName(dst_layername)
except:
    poly_layer = None

if poly_layer is None:

    srs = None
    if src_ds.GetProjectionRef() != '':
        srs = osr.SpatialReference()
        srs.ImportFromWkt( src_ds.GetProjection())
        
    poly_layer = mem_ds.CreateLayer('poly', srs, ogr.wkbPolygon)

    if dst_fieldname is None:
        dst_fieldname = 'DN'
        
    fd = ogr.FieldDefn( dst_fieldname, ogr.OFTInteger )
    poly_layer.CreateField( fd )
    dst_field = 0
else:
    if dst_fieldname is not None:
        dst_field = poly_layer.GetLayerDefn().GetFieldIndex(dst_fieldname)
        if dst_field < 0:
            print("Warning: cannot find field '%s' in layer '%s'" % (dst_fieldname, dst_layername))

# try:
#     line_layer = mem_ds.GetLayerByName(dst_layername)
# except:
#     line_layer = None

# if line_layer is None:

#     srs = None
#     if src_ds.GetProjectionRef() != '':
#         srs = osr.SpatialReference()
#         srs.ImportFromWkt( src_ds.GetProjectionRef)
        
#     line_layer = mem_ds.CreateLayer('line', srs, ogr.wkbMultiLineString)

#     if dst_fieldname is None:
#         dst_fieldname = 'DN'
        
#     fd = ogr.FieldDefn( dst_fieldname, ogr.OFTInteger )
#     line_layer.CreateField( fd )
#     dst_field = 0
# else:
#     if dst_fieldname is not None:
#         dst_field = line_layer.GetLayerDefn().GetFieldIndex(dst_fieldname)
#         if dst_field < 0:
#             print("Warning: cannot find field '%s' in layer '%s'" % (dst_fieldname, dst_layername))

            
# for k, t in enumerate(time):
for k in range(src_ds.RasterCount):
    if k==0:
        #timestamp = cdftime.num2date(t)
        #print('Processing {}'.format(timestamp))
        srcband = src_ds.GetRasterBand(k+1)
        result = gdal.Polygonize(srcband, None, poly_layer, dst_field, [],
                          callback = gdal.TermProgress)
        poly_layer.SetAttributeFilter("{} = 0".format(dst_fieldname))
        # geomcol =  ogr.Geometry(ogr.wkbGeometryCollection)
        # for feature in poly_layer:
        #     geom = feature.GetGeometryRef()
        #     ring = geom.GetGeometryRef(0)
        #     geomcol.AddGeometry(ring)
        # geomcol =  ogr.Geometry(ogr.wkbGeometryCollection)
        # featureDefn = line_layer.GetLayerDefn()
        # line_feature = ogr.Feature(featureDefn)
        # line_feature.SetGeometry(geomcol)
        # line_layer.CreateFeature(line_feature)
        
        
poly_layer = None
#line_layer = None
mem_ds = None

# save(shp_filename, contour_points, level)


nc.close()
