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


def create_memory_layer(dst_fieldname):
    '''
    Create a in-memory layer with 1 OFTInteger field
    '''

    srs = None
    if src_ds.GetProjectionRef() != '':
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_ds.GetProjection())

    layer = mem_ds.CreateLayer('poly', srs, ogr.wkbPolygon)

    fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
    layer.CreateField(fd)
    dst_field = 0

    return layer, dst_field

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


parser = ArgumentParser(description='''A script the terminus in a (possibly) PISM netCDF file, and save it as a shapefile (polygon).''')
parser.add_argument("FILE", nargs=1)
parser.add_argument("-o", "--output_filename", dest="out_file",
                  help="Name of the output shape file", default='terminus.shp')


options = parser.parse_args()
filename = options.FILE[0]
shp_filename = options.out_file
dst_fieldname = 'mask'
ts_fieldname = 'timestamp'

nc = NC(filename, 'r')
xdim, ydim, zdim, tdim = ppt.get_dims(nc)

if tdim:
    time = nc.variables[tdim]
    time_units = time.units
    time_calendar = time.calendar
    cdftime = utime(time_units, time_calendar)
    timestamps = cdftime.num2date(time[:])
    has_time = True
else:
    tdim = None
nc.close()

src_ds = gdal.Open('NETCDF:{}:{}'.format(filename, dst_fieldname))

# Get Memory Driver
mem_driver = ogr.GetDriverByName('Memory')
mem_ds = mem_driver.CreateDataSource('memory_layer')

# Get SHP Driver
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
            
for k in range(src_ds.RasterCount):
    if tdim is None:
        timestamp = '0-0-0'
    else:
        timestamp = timestamps[k]
    print('Processing {}'.format(timestamp))
    srcband = src_ds.GetRasterBand(k+1)
    poly_layer, dst_field = create_memory_layer(dst_fieldname)
    result = gdal.Polygonize(srcband, None, poly_layer, dst_field, [],
                      callback = gdal.TermProgress)
    poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, ocean_value))
    ocean_layer, dst_field = create_memory_layer(dst_fieldname)
    featureDefn = ocean_layer.GetLayerDefn()
    for m, feature in enumerate(poly_layer):
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(bufferDist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        ocean_layer.CreateFeature(outFeature)

    poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, floating_value))
    floating_layer, dst_field = create_memory_layer(dst_fieldname)
    featureDefn = floating_layer.GetLayerDefn()
    for m, feature in enumerate(poly_layer):
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(bufferDist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        floating_layer.CreateFeature(outFeature)

    # Now clip layers
    tmp_layer, dst_field = create_memory_layer(dst_fieldname)
    ocean_layer.Clip(floating_layer, tmp_layer)
    poly_layer = None
    ocean_layer = None
    floating_layer = None

    featureDefn = terminus_layer.GetLayerDefn()
    for feature in tmp_layer:
        # create a new feature
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(feature.GetGeometryRef())
        i = outFeature.GetFieldIndex(ts_fieldname)
        outFeature.SetField(i, str(timestamp))
        # add the feature to the output layer
        terminus_layer.CreateFeature(outFeature)

        
# Clean-up
poly_layer = None
terminus_layer = None
mem_ds = None
src_ds = None
