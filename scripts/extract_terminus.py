#!/usr/bin/env python

import numpy as np
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
ts_fieldname = 'timestamp'

nc = NC(filename, 'r')
nc_projection = ppt.get_projection_from_file(nc)

xdim, ydim, zdim, tdim = ppt.get_dims(nc)
var_order = (tdim, zdim, ydim, xdim)

x = np.squeeze(nc.variables[xdim])
y = np.squeeze(nc.variables[ydim])


time = nc.variables['time']
time_units = time.units
time_calendar = time.calendar
cdftime = utime(time_units, time_calendar)
timestamps = cdftime.num2date(time[:])
nc.close()

src_ds = gdal.Open('NETCDF:{}:{}'.format(filename, dst_fieldname))

mem_driver = ogr.GetDriverByName('Memory')
mem_ds = mem_driver.CreateDataSource('memory_layer')

def create_memory_layer():
    try:
        layer = mem_ds.GetLayerByName(dst_layername)
    except:
        layer = None

    if layer is None:

        srs = None
        if src_ds.GetProjectionRef() != '':
            srs = osr.SpatialReference()
            srs.ImportFromWkt(src_ds.GetProjection())

        layer = mem_ds.CreateLayer('poly', srs, ogr.wkbPolygon)

        if dst_fieldname is None:
            dst_fieldname = 'DN'

        fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
        layer.CreateField(fd)
        dst_field = 0
    else:
        if dst_fieldname is not None:
            dst_field = layer.GetLayerDefn().GetFieldIndex(dst_fieldname)
            if dst_field < 0:
                print("Warning: cannot find field '%s' in layer '%s'" % (dst_fieldname, dst_layername))
    return layer

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
    ocean_layer = mem_ds.GetLayerByName(dst_layername)
except:
    ocean_layer = None

if ocean_layer is None:

    srs = None
    if src_ds.GetProjectionRef() != '':
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_ds.GetProjection())

    ocean_layer = mem_ds.CreateLayer('ocean', srs, ogr.wkbPolygon)

    if dst_fieldname is None:
        dst_fieldname = 'DN'
        
    fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
    ocean_layer.CreateField(fd)
    ocean_dst_field = 0
else:
    if dst_fieldname is not None:
        ocean_dst_field = ocean_layer.GetLayerDefn().GetFieldIndex(dst_fieldname)
        if ocean_dst_field < 0:
            print("Warning: cannot find field '%s' in layer '%s'" % (dst_fieldname, dst_layername))

try:
    floating_layer = mem_ds.GetLayerByName(dst_layername)
except:
    floating_layer = None

if floating_layer is None:

    srs = None
    if src_ds.GetProjectionRef() != '':
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_ds.GetProjection())

    floating_layer = mem_ds.CreateLayer('floating', srs, ogr.wkbPolygon)

    if dst_fieldname is None:
        dst_fieldname = 'DN'
        
    fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
    floating_layer.CreateField(fd)
    floating_dst_field = 0
else:
    if dst_fieldname is not None:
        floating_dst_field = floating_layer.GetLayerDefn().GetFieldIndex(dst_fieldname)
        if floating_dst_field < 0:
            print("Warning: cannot find field '%s' in layer '%s'" % (dst_fieldname, dst_layername))

try:
    tmp_layer = mem_ds.GetLayerByName(dst_layername)
except:
    tmp_layer = None

if tmp_layer is None:

    srs = None
    if src_ds.GetProjectionRef() != '':
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_ds.GetProjection())

    tmp_layer = mem_ds.CreateLayer('floating', srs, ogr.wkbPolygon)

    if dst_fieldname is None:
        dst_fieldname = 'DN'
        
    fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
    tmp_layer.CreateField(fd)
    tmp_dst_field = 0
    fd = ogr.FieldDefn(ts_fieldname, ogr.OFTDate)
    tmp_layer.CreateField(fd)
    tmp_dst_field = 0


else:
    if dst_fieldname is not None:
        tmp_dst_field = tmp_layer.GetLayerDefn().GetFieldIndex(dst_fieldname)
        if tmp_dst_field < 0:
            print("Warning: cannot find field '%s' in layer '%s'" % (dst_fieldname, dst_layername))

# Get driver
shp_driver = ogr.GetDriverByName('ESRI Shapefile')
shp_filename = validateShapePath(shp_filename)
if os.path.exists(shp_filename): 
    os.remove(shp_filename)
shp_ds = shp_driver.CreateDataSource(shp_filename)

srs = None
if src_ds.GetProjectionRef() != '':
    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjection())


terminus_layer = shp_ds.CreateLayer('terminus', srs, ogr.wkbPolygon)

fd = ogr.FieldDefn(ts_fieldname, ogr.OFTDate)
terminus_layer.CreateField(fd)
terminus_dst_field = 0


bufferDist = 1
ocean_value = 4
floating_value = 3
            
# for k, t in enumerate(time):
for k in range(src_ds.RasterCount):
    if k<2:
        timestamp = timestamps[k]
        print('Processing {}'.format(timestamp))
        srcband = src_ds.GetRasterBand(k+1)
        poly_layer = create_memory_layer()
        result = gdal.Polygonize(srcband, None, poly_layer, poly_dst_field, [],
                          callback = gdal.TermProgress)
        poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, ocean_value))
        featureDefn = ocean_layer.GetLayerDefn()
        for m, feature in enumerate(poly_layer):
            ingeom = feature.GetGeometryRef()
            geomBuffer = ingeom.Buffer(bufferDist)

            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(geomBuffer)
            ocean_layer.CreateFeature(outFeature)

        poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, floating_value))
        featureDefn = ocean_layer.GetLayerDefn()
        for m, feature in enumerate(poly_layer):
            ingeom = feature.GetGeometryRef()
            geomBuffer = ingeom.Buffer(bufferDist)

            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(geomBuffer)
            floating_layer.CreateFeature(outFeature)

        # Now clip them
        ocean_layer.Clip(floating_layer, tmp_layer)
        
        featureDefn = terminus_layer.GetLayerDefn()
        for feature in tmp_layer:
            # create a new feature
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(feature.GetGeometryRef())
            i = outFeature.GetFieldIndex(ts_fieldname)
            outFeature.SetField(i, str(timestamp))
            # add the feature to the output layer
            terminus_layer.CreateFeature(outFeature)

        
        
poly_layer = None
terminus_layer = None
mem_ds = None
src_ds = None
