#!/usr/bin/env python

import numpy as np
import pylab as plt
from argparse import ArgumentParser
from scipy.interpolate import interp1d
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC


def interp():
    # do the actual interpolation
    for k in range(0, Mx):
        for l in range(0, My):
            vH_out = usurf_out[k,l]
            vH_in = usurf_in[k,l]
            var_in_column = var_in[k,l,:]
            z_scale = vH_out / vH_in
            index_out = index_in * z_scale
            f = interp1d(index_out, var_in_column)
            index_in_belowH = index_in[z_in < vH_out]
            index_in_belowHp1 = index_in_belowH[-1] + 1
            index_out_belowH = index_in[z_in < vH_in]
            index_out_belowHp1 = index_out_belowH[-1] + 1
            # prefill with uppermost value inside
            var_out_column = np.ones_like(var_in_column) * var_in_column[index_out_belowHp1]
            var_out_column[0:index_in_belowHp1] = f(index_in_belowH)
            var_out[k,l,:] = var_out_column


if __name__ == '__main__':
    import time

    # grid dimensions
    Mx = 500
    My = 500
    Mz = 401

    # create some fake data

    # elevation of old surface
    usurf_in = 1.85 * np.pi*np.ones((Mx, My))
    # elevation of new surface
    usurf_out = 1.75 * np.pi*np.ones((Mx, My))

    index_in = np.linspace(0, Mz-1, Mz)
    z_in = np.linspace(0, 2*np.pi, Mz)
    var_in = np.zeros((Mx, My, Mz))
    for k in range(0, Mx):
        for l in range(0, My):
            var_in[k,l,:] = np.sin(z_in)


    var_out = np.zeros((Mx, My, Mz))

    t = time.time()
    interp()
    elapsed = time.time() - t
    print("time spent interpolating: %f" % elapsed)

    plt.figure()
    plt.plot(index_in, var_in[0,0,:])
    plt.plot(index_in, var_out[0,0,:])
    plt.show()
