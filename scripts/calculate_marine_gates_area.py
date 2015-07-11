#!/usr/bin/env python
# Copyright (C) 2014-2015 Andy Aschwanden
#
# Calculates the total aggregate cross-sectional area of marine
# gates (bed below sea level) from flux gates.
# Flux gates need to be extracted by running extract_profiles.py
# first.

import numpy as np
import pylab as plt
from argparse import ArgumentParser
from netCDF4 import Dataset as NC
from udunits2 import Converter, System, Unit
try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt

# Init Unit system
sys = System()

# Set up the option parser
parser = ArgumentParser()
parser.description = "Calculate flux gate cross sectional area."
parser.add_argument("FILE", nargs='*')

options = parser.parse_args()
args = options.FILE

varname = 'bedrock_altitude'

area_outunits = 'km2'

for k in range(len(args)):
    infile = args[k]
    print("opening file %s" % infile)
    nc = NC(infile, 'r')

    myvar = varname
    for name in nc.variables.keys():
        v = nc.variables[name]
        if getattr(v, "standard_name", "") == varname:
            print("variabe {0} found by its standard_name {1}".format(name,
                                                                      varname))
            myvar = name
    print(("    - reading variable %s from file %s" % (myvar, infile)))

    bed = nc.variables[myvar]
    bed_inunits = bed.units

    profile = nc.variables['profile']
    p = profile[:].flatten()
    dx = p[1] - p[0]
    dx_inunits = profile.units

    area_inunits = Unit(sys, bed_inunits) * Unit(sys, dx_inunits)
    area = np.sum(-bed[:] * (bed[:] < 0)) * dx
    c = Converter((area_inunits, area_outunits))
    out_area = ppt.unit_converter(area, area_inunits, area_outunits)
    print(
        'total aggregate cross-sectional area is {:4.0f} {}'.format(out_area, area_outunits))
    nc.close()
