#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import gdal
import ogr
import osr
import os
import logging
import logging.handlers
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
import fiona
import itertools


# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
fh = logging.handlers.RotatingFileHandler('extract.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter
formatter = logging.Formatter('%(module)s:%(lineno)d - %(message)s')

# add formatter to ch and fh
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)


parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
    description='''A script to extract interfaces (calving front, ice-ocean, or groundling line) from a PISM netCDF file, and save it as a shapefile (polygon).''')
parser.add_argument("FILE", nargs=1)
parser.add_argument("-a", "--attribute_field" , dest="field",
                    help="Attribute field to group by", default='timestep')
parser.add_argument("-o", "--output_filename", dest="out_file",
                    help="Name of the output shape file", default='dissolved.shp')


options = parser.parse_args()
ifile = options.FILE[0]
ofile = options.out_file
field = options.field

with fiona.open(ifile) as input:
    # preserve the schema of the original shapefile, including the crs
    meta = input.meta
    with fiona.open(ofile, 'w', **meta) as output:
        # groupby clusters consecutive elements of an iterable which have the same key so you must first sort the features by the 'STATEFP' field
        e = sorted(input, key=lambda k: k['properties'][field])
        # group by the attribute field 
        for key, group in itertools.groupby(e, key=lambda x:x['properties'][field]):
            properties, geom = zip(*[(feature['properties'], shape(feature['geometry'])) for feature in group])
            # write the feature, computing the unary_union of the elements in the group with the properties of the first element in the group
            output.write({'geometry': mapping(unary_union(geom)), 'properties': properties[0]})


# Update the area
shp_driver = ogr.GetDriverByName('ESRI Shapefile')
ds = shp_driver.Open(ofile, 1)

nl =  ds.GetLayerCount()
for k in range(nl):
    layer = ds.GetLayer(k)
    for feature in layer:
        geom = feature.GetGeometryRef()
        area = geom.GetArea()
        feature.SetField("area", int(area))
        layer.SetFeature(feature)

