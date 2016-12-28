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
import logging
import logging.handlers

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
fh = logging.handlers.RotatingFileHandler('extract.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(lineno)d - %(message)s')

# add formatter to ch and fh
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)


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


parser = ArgumentParser(
    description='''A script to extract interfaces (calving front, ice-ocean, or groundling line) from a PISM netCDF file, and save it as a shapefile (polygon).''')
parser.add_argument("FILE", nargs=1)
parser.add_argument("-o", "--output_filename", dest="out_file",
                    help="Name of the output shape file", default='interface.shp')
parser.add_argument("-t", "--type" , dest="extract_type",
                    choices=['calving_front', 'grounding_line', 'ice_ocean'],
                    help="Interface to extract.", default='ice_ocean')


options = parser.parse_args()
filename = options.FILE[0]
extract_type = options.extract_type
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

interface_layer = shp_ds.CreateLayer('interface', srs, ogr.wkbPolygon)
fd = ogr.FieldDefn(ts_fieldname, ogr.OFTString)
interface_layer.CreateField(fd)
interface_dst_field = 0

bufferDist = 1
if extract_type in ('calving_front'):
    a_value = 4
    b_value = 3
elif extract_type in ('grounding_line'):
    a_value = 3
    b_value = 2
elif extract_type in ('ice_ocean'):
    a_value = [0, 4]
    b_value = [2, 3]
else:
    print('Type {} not recognized'.format(extact_type))
    import sys
    sys.exit(0)
    
for k in range(src_ds.RasterCount):
    if tdim is None:
        timestamp = '0-0-0'
    else:
        timestamp = timestamps[k]
    logger.info('Processing {}'.format(timestamp))
    srcband = src_ds.GetRasterBand(k + 1)
    poly_layer, dst_field = create_memory_layer(dst_fieldname)
    logger.info('Running gdal.Polygonize()')
    result = gdal.Polygonize(srcband, None, poly_layer, dst_field, [],
                             callback=gdal.TermProgress)
    if extract_type in ('ice_ocean'):
        poly_layer.SetAttributeFilter("{dn} = {val1} OR {dn} = {val2}".format(dn=dst_fieldname, val1=a_value[0], val2=a_value[1]))
    else:
        poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, a_value))
    logger.info('Extracting interface A')
    a_layer, dst_field = create_memory_layer(dst_fieldname)
    featureDefn = a_layer.GetLayerDefn()
    for m, feature in enumerate(poly_layer):
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(bufferDist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        a_layer.CreateFeature(outFeature)

    if extract_type in ('ice_ocean'):
        poly_layer.SetAttributeFilter("{dn} = {val1} OR {dn} = {val2}".format(dn=dst_fieldname, val1=b_value[0], val2=b_value[1]))
    else:
        poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, a_value))
    logger.info('Extracting interface B')
    b_layer, dst_field = create_memory_layer(dst_fieldname)
    featureDefn = b_layer.GetLayerDefn()
    for m, feature in enumerate(poly_layer):
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(bufferDist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        b_layer.CreateFeature(outFeature)

    # Now clip layers
    logger.info('Clipping A and B')
    tmp_layer, dst_field = create_memory_layer(dst_fieldname)
    a_layer.Clip(b_layer, tmp_layer)
    poly_layer = None
    a_layer = None
    b_layer = None

    logger.info('Saving results')
    featureDefn = interface_layer.GetLayerDefn()
    for feature in tmp_layer:
        # create a new feature
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(feature.GetGeometryRef())
        i = outFeature.GetFieldIndex(ts_fieldname)
        outFeature.SetField(i, str(timestamp))
        # add the feature to the output layer
        interface_layer.CreateFeature(outFeature)

# Clean-up
poly_layer = None
interface_layer = None
mem_ds = None
src_ds = None
    
