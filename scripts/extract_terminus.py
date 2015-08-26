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
        srs.ImportFromWkt(src_ds.GetProjection())

    poly_layer = mem_ds.CreateLayer('poly', srs, ogr.wkbPolygon)

    if dst_fieldname is None:
        dst_fieldname = 'DN'
        
    fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
    poly_layer.CreateField(fd)
    poly_dst_field = 0
else:
    if dst_fieldname is not None:
        poly_dst_field = poly_layer.GetLayerDefn().GetFieldIndex(dst_fieldname)
        if poly_dst_field < 0:
            print("Warning: cannot find field '%s' in layer '%s'" % (dst_fieldname, dst_layername))

try:
    line_layer = mem_ds.GetLayerByName(dst_layername)
except:
    line_layer = None

if line_layer is None:

    srs = None
    if src_ds.GetProjectionRef() != '':
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_ds.GetProjectionRef())
        
    line_layer = mem_ds.CreateLayer('line', srs, ogr.wkbMultiLineString)

    if dst_fieldname is None:
        dst_fieldname = 'DN'
        
    fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
    line_layer.CreateField(fd)
    line_dst_field = 0
else:
    if dst_fieldname is not None:
        line_dst_field = line_layer.GetLayerDefn().GetFieldIndex(dst_fieldname)
        if line_dst_field < 0:
            print("Warning: cannot find field '%s' in layer '%s'" % (dst_fieldname, dst_layername))

ocean_value = 4
floating_value = 3
            
# for k, t in enumerate(time):
for k in range(src_ds.RasterCount):
    if k==0:
        #timestamp = cdftime.num2date(t)
        #print('Processing {}'.format(timestamp))
        srcband = src_ds.GetRasterBand(k+1)
        result = gdal.Polygonize(srcband, None, poly_layer, poly_dst_field, [],
                          callback = gdal.TermProgress)
        poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, ocean_value))
        for m, feature in enumerate(poly_layer):
            geom = feature.GetGeometryRef()
            ring = geom.GetGeometryRef(0)
            # print ring.GetPointCount()
            print feature.GetField(dst_fieldname)
            # Create feature
            featureDefn = line_layer.GetLayerDefn()
            line_feature = ogr.Feature(featureDefn)
            line_feature.SetGeometry(ring)
            line_feature.SetFID(m)
            i = line_feature.GetFieldIndex(dst_fieldname)
            line_feature.SetField(i, feature.GetField(dst_fieldname))
            line_layer.CreateFeature(line_feature)
        poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, floating_value))
        # for m, feature in enumerate(poly_layer):
        #     geom = feature.GetGeometryRef()
        #     ring = geom.GetGeometryRef(0)
        #     # print ring.GetPointCount()
        #     print feature.GetField(dst_fieldname)
        #     # Create feature
        #     featureDefn = line_layer.GetLayerDefn()
        #     line_feature = ogr.Feature(featureDefn)
        #     line_feature.SetGeometry(ring)
        #     line_feature.SetFID(m)
        #     i = line_feature.GetFieldIndex(dst_fieldname)
        #     line_feature.SetField(i, feature.GetField(dst_fieldname))
        #     line_layer.CreateFeature(line_feature)
        
        
# poly_layer = None
# line_layer = None
# mem_ds = None



nc.close()
