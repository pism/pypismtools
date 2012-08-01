#!/usr/bin/env python

import numpy as np
import pylab as plt
from argparse import ArgumentParser
from scipy.interpolate import interp1d

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

# elevation of old surface
vH_in = 1.85 * np.pi
# elevation of new surface
vH_out = 1.75 * np.pi

# grid dimensions
Mx = 5
My = 5
Mz = 401

# create some fake data
index_in = np.linspace(0, Mz-1, Mz)
z_in = np.linspace(0, 2*np.pi, Mz)
var_in = np.zeros((Mx, My, Mz))
for k in range(0, Mx):
    for l in range(0, My):
        var_in[k,l,:] = np.sin(z_in)


# do the actual interpolation
z_scale = vH_out / vH_in
index_out = index_in * z_scale
f = interp1d(index_out, var_in)
var_out = f(index_in[z_in < vH_out])

# 3d version, not yet working
# elevation of old surface
vH_in = 1.85 * np.pi*np.ones((Mx, My))
# elevation of new surface
vH_out = 1.75 * np.pi*np.ones((Mx, My))


# do the actual interpolation
z_scale = vH_out / vH_in
index_out = np.tile(index_in, (Mx, My, 1)) * np.tile(z_scale, (Mz, 1, 1)).transpose((1, 2, 0))
f = interp1d(index_out, var_in)
var_out = f(index_in[z_in < vH_out])

