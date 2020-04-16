#!/usr/bin/env python
# Copyright (C) 2020 Andy Aschwanden

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import ogr
import osr


if __name__ == "__main__":

    __spec__ = None

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="""A script to point shape file created by VectorFieldCalc to a MultiLine shape file.""",
    )
    parser.add_argument("INFILE", nargs=1)
    parser.add_argument("OUTFILE", nargs=1)
    parser.add_argument(
        "--s_srs",
        dest="s_srs",
        help="Source CRS",
        default=None,
    )
    parser.add_argument(
        "--t_srs",
        dest="t_srs",
        help="Target CRS",
        default=None,
    )

    options = parser.parse_args()
    infile = options.INFILE[0]
    outfile = options.OUTFILE[0]
    
    s_ds = ogr.Open(infile)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    t_ds = driver.CreateDataSource(outfile)

    s_layer = s_ds.GetLayer(0)
    s_spatialRef = s_layer.GetSpatialRef()

    if s_spatialRef == "":
        print(f"{infile} does not have a valid projection. Use '--s_srs projection' to set source CRS. (NOT implemented yet)")
    else:
        srs = s_spatialRef

    t_layer = t_ds.CreateLayer("profile", srs, ogr.wkbLineString)
    fieldDefn = ogr.FieldDefn('profile_id', ogr.OFTReal)
    t_layer.CreateField(fieldDefn)

    def profile(coords):
        p = ogr.Geometry(type=ogr.wkbLineString)
        for xy in coords:
            p.AddPoint_2D(xy[0],xy[1])
        return p

    profiles = {}
    for feature in s_layer:
        profile_id = feature.path_id
        point_id = feature.point_id
        geom = feature.GetGeometryRef()
        point = geom.GetPoint()
        
        if profile_id not in profiles.keys():
            profiles[profile_id] = []

        profiles[profile_id].append((point[0], point[1]))

    for coords in profiles.values():
        featureDefn = t_layer.GetLayerDefn()
        feature = ogr.Feature(featureDefn)
        p = profile(coords)
        feature.SetGeometry(p)
        t_layer.CreateFeature(feature)

        p.Destroy()
        feature.Destroy()

    t_ds.Destroy()
