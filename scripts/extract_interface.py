#!/usr/bin/env python

import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from netCDF4 import Dataset as NC
import cftime
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import os
import logging
import logging.handlers
import pandas as pd

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
formatter = logging.Formatter("%(message)s")

# add formatter to ch and fh
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)

dtype_dict = {
    np.dtype("O"): ogr.OFTString,
    np.dtype("float64"): ogr.OFTReal,
    np.dtype("int64"): ogr.OFTInteger,
    np.dtype("bool"): ogr.OFTBinary,
    str: ogr.OFTString,
}


def create_memory_layer(dst_fieldname):
    """
    Create a in-memory layer with 1 OFTInteger field
    """

    srs = None
    if src_ds.GetProjectionRef() != "":
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_ds.GetProjection())

    layer = mem_ds.CreateLayer("poly", srs, ogr.wkbMultiPolygon)

    fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
    layer.CreateField(fd)
    dst_field = 0

    return layer, dst_field


if __name__ == "__main__":

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="""A script to extract interfaces (calving front, ice-ocean, or groundling line) from a PISM netCDF file, and save it as a shapefile (polygon).""",
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
    parser.add_argument(
        "-e", "--epsg", dest="epsg", type=int, help="Sets EPSG code", default=None
    )
    parser.add_argument(
        "-l",
        "--level",
        dest="level",
        type=float,
        help="Which contour level to extract. Used in combination with 'contour'",
        default=1000,
    )
    parser.add_argument(
        "-o",
        "--output_filename",
        dest="out_file",
        help="Name of the output shape file",
        default="interface.shp",
    )
    parser.add_argument(
        "--ensemble_file",
        dest="ensemble_file",
        help="CSV file. If given, add parameter values as attributes",
        default=None,
    )
    parser.add_argument(
        "-m",
        "--mask_variable",
        dest="dst_fieldname",
        help="Name of variable to use",
        default="mask",
    )
    parser.add_argument(
        "-s",
        "--step",
        dest="step",
        type=int,
        help="Only extract every step value",
        default=1,
    )
    parser.add_argument(
        "-t",
        "--type",
        dest="extract_type",
        choices=[
            "calving_front",
            "grounded_floating",
            "ice_noice",
            "ice_ocean",
            "grounding_line",
            "ela",
            "contour",
            "sftgif",
        ],
        help="Interface to extract.",
        default="ice_ocean",
    )

    options = parser.parse_args()
    filename = options.FILE[0]
    area_threshold = options.area_threshold
    epsg = options.epsg
    ensemble_file = options.ensemble_file
    extract_type = options.extract_type
    level = options.level
    shp_filename = options.out_file
    ts_fieldname = "timestamp"
    dst_fieldname = options.dst_fieldname
    step = options.step

    if ensemble_file:
        e_df = pd.read_csv(ensemble_file)

    nc = NC(filename, "r")
    xdim = "x"
    ydim = "y"
    zdim = "z"
    tdim = "time"

    if tdim:
        time = nc.variables[tdim]
        time_units = time.units
        time_calendar = time.calendar
        timestamps = cftime.num2date(time[:], time_units, time_calendar)
        has_time = True
    else:
        tdim = None

    if ensemble_file:
        run_id = nc.id
    nc.close()

    src_ds = gdal.Open("NETCDF:{}:{}".format(filename, dst_fieldname))

    # Get Memory Driver
    mem_driver = ogr.GetDriverByName("Memory")
    mem_ds = mem_driver.CreateDataSource("memory_layer")

    # Get SHP Driver
    shp_driver = ogr.GetDriverByName("GPKG")
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

    interface_layer = dst_ds.CreateLayer("interface", srs, ogr.wkbPolygon)
    fd = ogr.FieldDefn(ts_fieldname, ogr.OFTString)
    interface_layer.CreateField(fd)
    fd = ogr.FieldDefn("area", ogr.OFTInteger)
    interface_layer.CreateField(fd)
    fd = ogr.FieldDefn("timestep", ogr.OFTDateTime)
    interface_layer.CreateField(fd)
    if ensemble_file:
        print("Creating additional fields")
        for field in e_df.keys():
            fd = ogr.FieldDefn(field, dtype_dict[e_df[field].dtype])
            interface_layer.CreateField(fd)

    interface_dst_field = 0

    bufferDist = 1
    if extract_type in ("calving_front"):
        a_value = 4
        b_value = 3
    elif extract_type in ("grounded_floating"):
        a_value = 3
        b_value = 2
    elif extract_type in ("ice_ocean"):
        a_value = 4
        b_value = [2, 3]
    elif extract_type in ("ice_noice", "sftgif"):
        a_value = 1
        b_value = 0
    elif extract_type in ("grounding_line"):
        a_value = 2
        b_value = [0, 3, 4]
    elif extract_type in ("ela"):
        a_value = 0
        b_value = 0
    elif extract_type in ("contour"):
        a_value = level
        b_value = level
    else:
        print(("Type {} not recognized".format(extact_type)))
        import sys

        sys.exit(0)

    time_step = 0
    for k in np.arange(0, src_ds.RasterCount, step):

        if tdim is None:
            timestamp = "0-0-0"
        else:
            timestamp = timestamps[k]
        logger.info("Processing {}".format(timestamp))
        srcband = src_ds.GetRasterBand(int(k + 1))
        poly_layer, dst_field = create_memory_layer(dst_fieldname)
        logger.debug("Running gdal.Polygonize()")
        result = gdal.Polygonize(
            srcband, None, poly_layer, dst_field, [], callback=gdal.TermProgress
        )
        if extract_type in ["ela", "contour"]:
            poly_layer.SetAttributeFilter("{} > {}".format(dst_fieldname, b_value))
        else:
            poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, a_value))
        logger.debug("Extracting interface A")
        a_layer, dst_field = create_memory_layer(dst_fieldname)
        featureDefn = a_layer.GetLayerDefn()
        for m, feature in enumerate(poly_layer):
            ingeom = feature.GetGeometryRef()
            geomBuffer = ingeom.Buffer(bufferDist)

            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(geomBuffer)
            a_layer.CreateFeature(outFeature)

        if extract_type in ["grounding_line"]:
            poly_layer.SetAttributeFilter(
                "{dn} = {val1} OR {dn} = {val2}  OR {dn} = {val3}".format(
                    dn=dst_fieldname, val1=b_value[0], val2=b_value[1], val3=b_value[2]
                )
            )
        elif extract_type in ["ice_ocean"]:
            poly_layer.SetAttributeFilter(
                "{dn} = {val1} OR {dn} = {val2}".format(
                    dn=dst_fieldname, val1=b_value[0], val2=b_value[1]
                )
            )
        elif extract_type in ["ela", "contour"]:
            poly_layer.SetAttributeFilter("{} < {}".format(dst_fieldname, b_value))
        else:
            poly_layer.SetAttributeFilter("{} = {}".format(dst_fieldname, b_value))
        logger.debug("Extracting interface B")
        b_layer, dst_field = create_memory_layer(dst_fieldname)
        featureDefn = b_layer.GetLayerDefn()
        for m, feature in enumerate(poly_layer):
            ingeom = feature.GetGeometryRef()
            geomBuffer = ingeom.Buffer(bufferDist)

            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(geomBuffer)
            b_layer.CreateFeature(outFeature)

        # Now clip layers
        logger.debug("Clipping A and B")
        tmp_layer, dst_field = create_memory_layer(dst_fieldname)
        a_layer.Clip(b_layer, tmp_layer)
        poly_layer = None
        a_layer = None
        b_layer = None

        logger.info("Saving results")
        featureDefn = interface_layer.GetLayerDefn()
        for feature in tmp_layer:
            # create a new feature
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(feature.GetGeometryRef())
            i = outFeature.GetFieldIndex("timestep")
            outFeature.SetField(i, int(time_step))
            i = outFeature.GetFieldIndex(ts_fieldname)
            outFeature.SetField(i, timestamp.strftime())
            geom = feature.GetGeometryRef()
            area = geom.GetArea()
            i = outFeature.GetFieldIndex("area")
            outFeature.SetField(i, int(area))
            if ensemble_file:
                df = e_df[e_df["id"] == int(run_id)]
                for d in df.items():
                    key = d[0]
                    val = d[1].values[0]
                    i = outFeature.GetFieldIndex(key)
                    print(i, key, val)
                    outFeature.SetField(i, str(val))
            # add the feature to the output layer
            if area >= area_threshold:
                interface_layer.CreateFeature(outFeature)

        time_step += 1

    # Clean-up
    poly_layer = None
    interface_layer = None
    mem_ds = None
    src_ds = None
