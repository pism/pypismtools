#!/usr/bin/env python

import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from netCDF4 import Dataset as NC
from netcdftime import utime
import gdal
import ogr
import osr
import os
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
fh = logging.handlers.RotatingFileHandler("extract.log")
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter
formatter = logging.Formatter("%(module)s:%(lineno)d - %(message)s")

# add formatter to ch and fh
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)


def create_memory_layer(dst_fieldname):
    """
    Create a in-memory layer with 1 OFTInteger field
    """

    srs = None
    if src_ds.GetProjectionRef() != "":
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_ds.GetProjection())

    layer = mem_ds.CreateLayer("contours", srs)

    fd = ogr.FieldDefn("id", ogr.OFTInteger)
    layer.CreateField(fd)
    fd = ogr.FieldDefn("level", ogr.OFTReal)
    layer.CreateField(fd)

    return layer


def validateShapePath(shapePath):
    """Validate shapefile extension"""
    return os.path.splitext(str(shapePath))[0] + ".shp"


def validateShapeData(shapeData):
    """Make sure we can access the shapefile"""
    # Make sure the shapefile exists
    if not shapeData:
        raise ShapeDataError("The shapefile is invalid")
    # Make sure there is exactly one layer
    if shapeData.GetLayerCount() != 1:
        raise ShapeDataError("The shapefile must have exactly one layer")


# Error


class ShapeDataError(Exception):
    pass


parser = ArgumentParser(
    formatter_class=ArgumentDefaultsHelpFormatter,
    description="""A script to extract contours from netCDF file, and save it as a shapefile.""",
)
parser.add_argument("FILE", nargs=1)
parser.add_argument(
    "-a",
    "--area_threshold",
    dest="area_threshold",
    type=float,
    help="Only save features with an area > area_threshold",
    default=200,
)
parser.add_argument("-e", "--epsg", dest="epsg", type=int, help="Sets EPSG code", default=None)
parser.add_argument(
    "-l", "--levels", dest="levels", help="Which contour levels to extract. Comma-separated list", default="0"
)
parser.add_argument(
    "-o", "--output_filename", dest="out_file", help="Name of the output shape file", default="interface.shp"
)
parser.add_argument("-v", "--variable", dest="dst_fieldname", help="Name of variable to use", default="usurf")

options = parser.parse_args()
filename = options.FILE[0]
area_threshold = options.area_threshold
epsg = options.epsg
levels = np.array(options.levels.split(","), dtype=float)
shp_filename = options.out_file
ts_fieldname = "timestamp"
dst_fieldname = options.dst_fieldname

nc = NC(filename, "r")
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

src_ds = gdal.Open("NETCDF:{}:{}".format(filename, dst_fieldname))

# Get Memory Driver
mem_driver = ogr.GetDriverByName("Memory")
mem_ds = mem_driver.CreateDataSource("memory_layer")

# Get SHP Driver
shp_driver = ogr.GetDriverByName("ESRI Shapefile")
shp_filename = validateShapePath(shp_filename)
if os.path.exists(shp_filename):
    os.remove(shp_filename)
dst_ds = shp_driver.CreateDataSource(shp_filename)

srs = None
if src_ds.GetProjectionRef() != "":
    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjection())

if epsg is not None:
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)


interface_layer = dst_ds.CreateLayer("interface", srs)
fd = ogr.FieldDefn("area", ogr.OFTInteger)
interface_layer.CreateField(fd)
fd = ogr.FieldDefn("level", ogr.OFTReal)
interface_layer.CreateField(fd)
fd = ogr.FieldDefn(ts_fieldname, ogr.OFTString)
interface_layer.CreateField(fd)
fd = ogr.FieldDefn("timestep", ogr.OFTInteger)
interface_layer.CreateField(fd)

time_step = 0
for k in np.arange(0, src_ds.RasterCount):

    if tdim is None:
        timestamp = "0-0-0"
    else:
        timestamp = timestamps[k]
    logger.info("Processing {}".format(timestamp))
    srcband = src_ds.GetRasterBand(int(k + 1))
    logger.debug("Running gdal.ContourGenerate()")
    tmp_layer = create_memory_layer(dst_fieldname)
    result = gdal.ContourGenerate(srcband, 0, 0, levels, 0, 0, tmp_layer, 0, 1, callback=gdal.TermProgress)

    logger.info("Saving results")
    featureDefn = interface_layer.GetLayerDefn()
    for feature in tmp_layer:
        # create a new feature
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(feature.GetGeometryRef())
        i = outFeature.GetFieldIndex("level")
        ik = feature.GetFieldIndex("level")
        outFeature.SetField(i, feature.GetField(ik))
        i = outFeature.GetFieldIndex("timestep")
        outFeature.SetField(i, int(time_step))
        i = outFeature.GetFieldIndex(ts_fieldname)
        outFeature.SetField(i, str(timestamp))
        geom = feature.GetGeometryRef()
        area = geom.GetArea()
        i = outFeature.GetFieldIndex("area")
        outFeature.SetField(i, int(area))
        # add the feature to the output layer
        if area >= area_threshold:
            interface_layer.CreateFeature(outFeature)
            print(area)
    time_step += 1

# Clean-up
interface_layer = None
mem_ds = None
src_ds = None
