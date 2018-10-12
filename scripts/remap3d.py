#!/usr/bin/env python

import time
import subprocess
import numpy as np
import pylab as plt
from argparse import ArgumentParser
from scipy.interpolate import interp1d
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt


def interp(z, var_old, thk_old, thk_new):
    '''

    Remaps/rescales 3D field from old surface to new surface.

    Parameters
    ----------
    z : 1d array, vertical coordinate
    var_old : (Mx,My,Mz) array, input values
    thk_old : (Mx,My) array, input upper ice surface
    thk_new : (Mx,My) array, output upper surface

    Returns
    -------

    var_new : (Mx,My,Mz) array
    '''

    indices = np.linspace(0, len(z) - 1, len(z))
    try:
        Mt, Mx, My, Mz = var_old.shape
    except:
        Mx, My, Mz = var_old.shape
        Mt = 0
    # prefill new variable (but what is a safe value?)
    var_new = np.zeros_like(var_old)
    # do the actual interpolation
    if Mt > 0:
        for j in range(0, Mt):
            for k in range(0, Mx):
                for l in range(0, My):
                    vH_new = thk_new[k, l]
                    vH_old = thk_old[k, l]
                    var_old_column = var_old[j, k, l, :]
                    if (vH_new > 0) and (vH_old > 0):
                        z_scale = vH_new / vH_old
                        index_new = indices * z_scale
                        f = interp1d(index_new, var_old_column)
                        indices_belowH = indices[z < vH_new]
                        index_belowHp1 = indices_belowH[-1] + 1
                        index_new_belowH = indices[z < vH_old]
                        index_new_belowHp1 = index_new_belowH[-1] + 1
                        # prefill with uppermost value inside
                        var_new_column = np.ones_like(
                            var_old_column) * var_old_column[index_new_belowHp1]
                        var_new_column[0:index_belowHp1] = f(indices_belowH)
                    else:
                        # fill with atmospheric value
                        var_new_column = np.ones_like(
                            var_old_column) * var_old_column[-1]
                    var_new[j, k, l, :] = var_new_column
    else:
        for k in range(0, Mx):
            for l in range(0, My):
                vH_new = thk_new[k, l]
                vH_old = thk_old[k, l]
                var_old_column = var_old[k, l, :]
                if (vH_new > 0):
                    z_scale = vH_new / vH_old
                    index_new = indices * z_scale
                    f = interp1d(index_new, var_old_column)
                    indices_belowH = indices[z < vH_new]
                    index_belowHp1 = indices_belowH[-1] + 1
                    index_new_belowH = indices[z < vH_old]
                    index_new_belowHp1 = index_new_belowH[-1] + 1
                    # prefill with uppermost value inside
                    var_new_column = np.ones_like(
                        var_old_column) * var_old_column[index_new_belowHp1]
                    var_new_column[0:index_belowHp1] = f(indices_belowH)
                else:
                    # fill with atmospheric value
                    var_new_column = np.ones_like(
                        var_old_column) * var_old_column[-1]
                var_new[k, l, :] = var_new_column
    return var_new


def create_test_data(Mx, My, Mz):
    '''
    Create test data
    '''

    # elevation of old surface
    usurf_old = 1.85 * np.pi * np.ones((Mx, My))
    # elevation of new surface
    usurf_new = 1.75 * np.pi * np.ones((Mx, My))

    z = np.linspace(0, 2 * np.pi, Mz)
    var_old = np.zeros((Mx, My, Mz))
    for k in range(0, Mx):
        for l in range(0, My):
            var_old[k, l, :] = np.sin(z)

    return z, var_old, usurf_old, usurf_new


if __name__ == '__main__':

    # Set up the argument parser
    parser = ArgumentParser()
    parser.description = '''A to remap/rescale 3D fields (e.g. enthalpy) from one ice sheet body
    to another, given their ice thicknesses. Both files need to have same (Mx,My) dimensions'''
    parser.add_argument("FILE", nargs='*')
    parser.add_argument("--test", dest="test", action='store_true',
                        help="test with some fake date", default=False)
    parser.add_argument("-c", "--copy_thickness", dest="copy_thickness", action='store_true',
                        help="copy ice thickness to new file", default=False)
    parser.add_argument("-v", "--variable", dest="varname",
                        help='''Variable used for remapping, default = "enthalpy".''', default='enthalpy')

    options = parser.parse_args()
    args = options.FILE
    copy_thickness = options.copy_thickness
    test = options.test
    interp_var_name = options.varname
    thk_var_name = 'thk'

    if not test and len(args) == 0:
        print("no arguments given, running with test data")
        test = True

    if not test:
        thk_file_name = args[0]
        from_file_name = args[1]
        to_file_name = args[2]

        subprocess.call(["ncks", "-O", from_file_name, to_file_name])
        if copy_thickness:
            print(("copy ice thickness from %s to %s." %
                   (thk_file_name, to_file_name)))
            subprocess.call(["ncks", "-A -v thk", thk_file_name, to_file_name])

        # open file in append mode
        nc_to = NC(to_file_name, 'a')
        # get dimensions from file
        xdim, ydim, zdim, tdim = ppt.get_dims(nc_to)
        # set variable order for permutation
        var_order = (tdim, xdim, ydim, zdim)
        # read in z-coordinate
        z = nc_to.variables[zdim][:]

        # read ice thickness
        print(("    - reading variable %s from file %s" %
               (thk_var_name, to_file_name)))
        try:
            thk_to = np.squeeze(
                ppt.permute(nc_to.variables[thk_var_name], var_order))
        except:
            print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
                   % (thk_var_name, to_file_name)))
            import sys
            sys.exit()

        # read interpolation variable
        print(("    - reading variable %s from file %s" %
               (interp_var_name, to_file_name)))
        try:
            # Don't know why squeezing changes dimension ordering
            ## var_old = np.squeeze(ppt.permute(nc_to.variables[interp_var_name], var_order))
            var_old = ppt.permute(nc_to.variables[interp_var_name], var_order)
        except:
            print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
                   % (interp_var_name, to_file_name)))
            import sys
            sys.exit()

        nc_thk = NC(thk_file_name, 'r')

        # read ice thickness
        print(("    - reading variable %s from file %s" %
               (thk_var_name, thk_file_name)))
        try:
            thk_from = np.squeeze(
                ppt.permute(nc_thk.variables[thk_var_name], var_order))
        except:
            print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
                   % (thk_var_name, thk_file_name)))
            import sys
            sys.exit()
    else:
        # grid dimensions
        Mx = 5
        My = 4
        Mz = 401
        z, var_old, thk_from, thk_to = create_test_data(Mx, My, Mz)

    t = time.time()
    var_new = interp(z, var_old, thk_to, thk_from)
    elapsed = time.time() - t
    print(("time spent interpolating: %f" % elapsed))

    if not test:

        interp_var = nc_to.variables[interp_var_name]
        input_dimensions = interp_var.dimensions

        # filter out irrelevant dimensions
        dimensions = [x for x in var_order if x in input_dimensions]
        # create the mapping
        mapping = [dimensions.index(x) for x in input_dimensions]

        if mapping:
            var_new = np.transpose(var_new, mapping)

        interp_var[:] = var_new

    Mt, Mx, My, Mz = var_new.shape
    i = np.floor(Mx / 2)
    j = np.floor(My / 2)

    plt.figure()
    plt.imshow(thk_to - thk_from)
    plt.colorbar()

    plt.figure()
    plt.plot(var_old[0, i, j, :], z, label='old')
    plt.plot(var_new[0, i, j, :], z, label='new')
    plt.legend()

    if not test:
        nc_thk.close()
        nc_to.close()
