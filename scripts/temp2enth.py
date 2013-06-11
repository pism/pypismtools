#!/usr/bin/env python
import numpy as np
from argparse import ArgumentParser

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt

import PISM
    
# Set up the option parser
parser = ArgumentParser()
parser.description = "A script to plot a variable in a netCDF file over a GeoTiff. Uses GDAL python bindings, Proj4, and Basemap. Script is fine-tuned for whole Greenland plots, but can be adapted for other needs."
parser.add_argument("FILE", nargs='*')

options = parser.parse_args()
args = options.FILE

nc = NC(args[0], 'r')

x = nc.variables['x'][:]
y = nc.variables['y'][:]
z = nc.variables['z'][:]

# The enthalpy convert expects Kelvin
inunits = nc.variables['temp'].units
outunits = 'K'
try:
    temp = ppt.unit_converter(np.squeeze(nc.variables['temp'][:]), inunits, outunits)
except:
    temp = np.squeeze(nc.variables['temp'][:])
    print('''WARNING: I don't know the units of variable temp, and assume it is Kelvin''')

enthalpy_true = np.squeeze(nc.variables['enthalpy'][:])
usurf = np.squeeze(nc.variables['usurf'][:])
topg = np.squeeze(nc.variables['topg'][:])

context = PISM.Context()
config = context.config
EC = PISM.EnthalpyConverter(config)

enthalpy = np.zeros_like(temp)

K, L, M = enthalpy.shape

p_air = config.get("surface_pressure")

for m in range(M):
    for l in range(L):
        for k in range(K):
            depth = topg[k,l] - usurf[k,l] + z[m]
            if (topg[k,l] + z[m] < usurf[k,l]):
                p = EC.getPressureFromDepth(depth)
            else:
                p = p_air
            # liquid water fraction is set to zero because
            # it is not available in this context
            e = EC.getEnth(temp[k,l,m], 0., p)
            enthalpy[k,l,m] = e

# Compare with enthalpy field in file.
# Note that difference is only zero in the absence of temperate ice
print enthalpy - enthalpy_true

nc.close()
