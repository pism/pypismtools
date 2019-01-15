#!/usr/bin/env python
# Copyright (C) 2015-17 Bob McNabb, Andy Aschwanden

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import gdal
import fiona
import numpy as np
from fiona.crs import from_epsg
from shapely.geometry import LineString, mapping
import sys
from netCDF4 import Dataset as NC
from netcdftime import utime

import logging
import logging.handlers

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
fh = logging.handlers.RotatingFileHandler("extract.log")
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(lineno)d - %(message)s"
)

# add formatter to ch and fh
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)


class GeoImgInfo:
    def __init__(self, in_filename, in_dir="."):
        self.filename = in_filename
        self.in_dir_path = in_dir
        try:
            self.gd = gdal.Open(self.filename)
        except:
            print(("could not open file %s" % self.filename))

        self.gt = self.gd.GetGeoTransform()
        self.proj = self.gd.GetProjection()
        self.npix_x = self.gd.RasterXSize
        self.npix_y = self.gd.RasterYSize
        self.xmin = self.gt[0]
        self.xmax = self.gt[0] + self.npix_x * self.gt[1] + self.npix_y * self.gt[2]
        self.ymin = self.gt[3] + self.npix_x * self.gt[4] + self.npix_y * self.gt[5]
        self.ymax = self.gt[3]
        self.dx = self.gt[1]
        self.dy = self.gt[5]
        self.nodatavalue = self.gd.GetRasterBand(1).GetNoDataValue()
        self.RasterCount = self.gd.RasterCount


def getRasterBandArray(in_filename, BandNo=1):

    gd = gdal.Open(in_filename)
    rb = gd.GetRasterBand(BandNo)
    return rb.ReadAsArray()


def get_dims(nc):
    """
    Gets dimensions from netcdf instance

    Parameters:
    -----------
    nc: netCDF instance

    Returns:
    --------
    xdim, ydim, zdim, tdim: dimensions
    """

    # a list of possible x-dimensions names
    xdims = ["x", "x1"]
    # a list of possible y-dimensions names
    ydims = ["y", "y1"]
    # a list of possible z-dimensions names
    zdims = ["z", "z1"]
    # a list of possible time-dimensions names
    tdims = ["t", "time"]

    xdim = None
    ydim = None
    zdim = None
    tdim = None

    # assign x dimension
    for dim in xdims:
        if dim in list(nc.dimensions.keys()):
            xdim = dim
    # assign y dimension
    for dim in ydims:
        if dim in list(nc.dimensions.keys()):
            ydim = dim
    # assign z dimension
    for dim in zdims:
        if dim in list(nc.dimensions.keys()):
            zdim = dim
    # assign time dimension
    for dim in tdims:
        if dim in list(nc.dimensions.keys()):
            tdim = dim
    return xdim, ydim, zdim, tdim


parser = ArgumentParser(
    formatter_class=ArgumentDefaultsHelpFormatter,
    description="Convert rasters containing (U,V) components of velocity field to vector line data.",
)
parser.add_argument("FILE", nargs=1)
parser.add_argument(
    "-U", "--Udata", dest="Udata", help="Raster containing x components of velocity"
)
parser.add_argument(
    "-V", "--Vdata", dest="Vdata", help="Raster containing y components of velocity"
)
parser.add_argument(
    "--Uerror",
    dest="Uerror",
    help="Raster containing x components of error",
    default=None,
)
parser.add_argument(
    "--Verror",
    dest="Verror",
    help="Raster containing y components of error",
    default=None,
)
parser.add_argument(
    "--epsg",
    dest="epsg",
    help="EPSG code of project. Overrides input projection",
    default=None,
)
parser.add_argument(
    "-s",
    "--scale_factor",
    type=float,
    dest="scale_factor",
    help="Scales length of line. Default=1.",
    default=1.0,
)
parser.add_argument(
    "-p",
    "--prune_factor",
    type=int,
    dest="prune_factor",
    help="Pruning. Only use every x-th value. Default=1",
    default=1,
)
parser.add_argument(
    "-t",
    "--threshold",
    type=float,
    dest="threshold",
    help="Magnitude values smaller or equal than threshold will be masked. Default=None",
    default=0.0,
)
args = parser.parse_args()
prune_factor = args.prune_factor
scale_factor = args.scale_factor
threshold = args.threshold

nc_file_u = args.Udata.split(":")[1]
nc = NC(nc_file_u, "r")
xdim, ydim, zdim, tdim = get_dims(nc)

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

URasterInfo = GeoImgInfo(args.Udata)
VRasterInfo = GeoImgInfo(args.Vdata)

URasterCount = URasterInfo.RasterCount
VRasterCount = VRasterInfo.RasterCount
assert URasterCount == VRasterCount
RasterCount = URasterCount
RasterInfo = URasterInfo
Ufill_value = URasterInfo.nodatavalue
Vfill_value = VRasterInfo.nodatavalue

if args.epsg is None:
    epsg = int(RasterInfo.proj.split(",")[-1].translate(None, '[]"' ""))
else:
    epsg = args.epsg

x = np.linspace(RasterInfo.xmin, RasterInfo.xmax, RasterInfo.npix_x)
y = np.linspace(RasterInfo.ymin, RasterInfo.ymax, RasterInfo.npix_y)

X, Y = np.meshgrid(x, np.flipud(y))
X = X[::prune_factor, ::prune_factor]
Y = Y[::prune_factor, ::prune_factor]

nx, ny = X.shape


# create the schema (fields) and get the EPSG information for the dataset
schema = {
    "properties": [
        ("ux", "float"),
        ("uy", "float"),
        ("speed", "float"),
        ("timestamp", "str"),
    ],
    "geometry": "LineString",
}


# open the shapefile
logger.info("Processing")
with fiona.open(
    args.FILE[0], "w", crs=from_epsg(epsg), driver="ESRI Shapefile", schema=schema
) as output:
    for k in range(RasterCount):
        if tdim is None:
            timestamp = "0-0-0"
        else:
            timestamp = timestamps[k]
        logger.info("Processing {}".format(timestamp))
        print(("Processing {}".format(timestamp)))

        Ux = getRasterBandArray(args.Udata, BandNo=k + 1)[
            ::prune_factor, ::prune_factor
        ]
        Uy = getRasterBandArray(args.Vdata, BandNo=k + 1)[
            ::prune_factor, ::prune_factor
        ]

        Speed = np.sqrt(Ux ** 2 + Uy ** 2)

        prop_dict = {}
        prop_dict["ux"] = Ux
        prop_dict["uy"] = Uy
        prop_dict["speed"] = Speed

        # Read and add error of U component
        if args.Uerror is not None:
            Ex = getRasterBandArray(args.Uerror, BandNo=k + 1)[
                ::prune_factor, ::prune_factor
            ]
            schema["properties"].append(("ex", "float"))
            prop_dict["ex"] = ex
        # Read and add error of V component
        if args.Verror is not None:
            Ey = getRasterBandArray(args.Verror, BandNo=k + 1)[
                ::prune_factor, ::prune_factor
            ]
            schema["properties"].append(("ey", "float"))
            prop_dict["ey"] = ey

        # create features for each x,y pair, and give them the right properties
        m = 0
        for i in range(nx):
            for j in range(ny):
                if (
                    (Ux[i, j] != Ufill_value)
                    & (Uy[i, j] != Vfill_value)
                    & (Speed[i, j] > threshold)
                ):
                    m += 1
                    sys.stdout.write("\r")
                    # Center cooridinates
                    x_c, y_c = X[i, j], Y[i, j]
                    # Start point
                    x_a, y_a = (
                        X[i, j] - scale_factor * Ux[i, j] / 2,
                        Y[i, j] - scale_factor * Uy[i, j] / 2,
                    )
                    # End point
                    x_e, y_e = (
                        X[i, j] + scale_factor * Ux[i, j] / 2,
                        Y[i, j] + scale_factor * Uy[i, j] / 2,
                    )
                    # Create LineString
                    line = LineString([[x_a, y_a], [x_c, y_c], [x_e, y_e]])
                    line_dict = dict(
                        [(k, float(v[i, j])) for (k, v) in prop_dict.items()]
                    )
                    line_dict["timestamp"] = str(timestamp)
                    output.write({"properties": line_dict, "geometry": mapping(line)})

        print("  {} points found and written".format(str(m)))

print("Done writing {}".format(args.FILE[0]))
# close the shapefile now that we're all done
output.close()
