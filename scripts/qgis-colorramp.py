#!/usr/bin/env python
import numpy as np
import matplotlib
import pylab as plt
import matplotlib as mpl
from argparse import ArgumentParser

try:
    from pypismtools import gmtColormap
except:
    from pypismtools.pypismtools import gmtColormap


def cmap_map(function, cmap):
    """
    Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous points
    in a colormap.

    Adapted from http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # First get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step: np.array(cmap(step)[0:3])
    old_LUT = np.array(map(reduced_cmap, step_list))
    new_LUT = np.array(map(function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(('red', 'green', 'blue')):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j, i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = sorted(map(lambda x: x + (x[1], ), this_cdict.items()))
        cdict[key] = colorvector

    return mpl.colors.LinearSegmentedColormap('colormap', cdict, 1024)

# Set up the option parser
parser = ArgumentParser()
parser.description = """
A script to convert a GMT (*.cpt) colormap or matplotlib colormap into QGIS-readable color ramp.
"""
parser.add_argument("FILE", nargs='*')
parser.add_argument("--log", dest="log", action="store_true",
                    help="make a log-normalized color scale", default=False)
parser.add_argument("--joughin_speed", dest="joughin_speed", action="store_true",
                    help='''
                  Joughin-style log''', default=False)
parser.add_argument("--joughin_speed_10k", dest="joughin_speed_10k", action="store_true",
                    help='''
                  Joughin-style log''', default=False)
parser.add_argument("--joughin_speed_1k", dest="joughin_speed_1k", action="store_true",
                    help='''
                  Joughin-style log''', default=False)
parser.add_argument("--joughin_speed_300", dest="joughin_speed_300", action="store_true",
                    help='''
                  Joughin-style log''', default=False)
parser.add_argument("--habermann_tauc", dest="habermann_tauc", action="store_true",
                    help='''
                  log tauc scaling from Habermann et al (2013)''', default=False)
parser.add_argument("--bath_topo", dest="bath_topo", action="store_true",
                    help='''
                  Scaling for bathymetry/topography for Greenland''', default=False)
parser.add_argument("--colorbar_extend", dest="cb_extend", choices=['neither', 'both', 'min', 'max'],
                    help='''Extend of colorbar. Default='both'.''', default='both')
parser.add_argument("--colorbar_label", dest="colorbar_label",
                    help='''Label for colorbar.''', default=None)
parser.add_argument("-a", "--a_log", dest="a", type=float,
                    help='''
                  a * logspace(vmin, vmax, N)''', default=1)
parser.add_argument("--vmin", dest="vmin", type=float,
                    help='''
                  a * logspace(vmin, vmax, N)''', default=-1)
parser.add_argument("--vmax", dest="vmax", type=float,
                    help='''
                  a * logspace(vmin, vmax, N)''', default=3)
parser.add_argument("--extend", dest="extend", nargs=2, type=float,
                    help='''
                  appends color ramp by repeating first and last color for value''',
                    default=None)
parser.add_argument("--N", dest="N", type=int,
                    help='''
                  a * logspace(vmin, vmax, N''', default=1022)
parser.add_argument("--topo2000m", dest="topo2000m", action="store_true",
                    help='''
                  Topo''', default=False)
parser.add_argument("-r", "--reverse", dest="reverse", action="store_true",
                    help="reverse color scale", default=False)
parser.add_argument("--orientation", dest="orientation", choices=['horizontal', 'vertical'],
                    help="Orientation, default = 'horizontal", default='horizontal')

options = parser.parse_args()
args = options.FILE
joughin_speed = options.joughin_speed
joughin_speed_10k = options.joughin_speed_10k
joughin_speed_1k = options.joughin_speed_1k
joughin_speed_300 = options.joughin_speed_300
habermann_tauc = options.habermann_tauc
bath_topo = options.bath_topo
topo2000m = options.topo2000m
a = options.a
log = options.log
extend = options.extend
N = options.N
vmin = options.vmin
vmax = options.vmax
reverse = options.reverse
colorbar_label = options.colorbar_label
cb_extend = options.cb_extend
# experimental
log_color = False
orientation = options.orientation

# read in CPT colormap
for k in range(len(args)):

    cmap_file = args[k]
    try:
        cdict = plt.cm.datad[cmap_file]
        prefix = cmap_file
    except:
        # import and convert colormap
        cdict = gmtColormap(cmap_file, log_color=log_color, reverse=reverse)
        prefix = '.'.join(cmap_file.split('.')[0:-1])
        suffix = cmap_file.split('.')[-1]

    # either log scaling or linear scaling (default)
    if joughin_speed:
        # This is a little duck-punching to get a QGIS colormap
        # similar to Joughin (2010)
        data_values = np.logspace(vmin, vmax, N)[0:889]
        data_values[-1] = 3000
        N = len(data_values)
        norm = mpl.colors.LogNorm(vmin=1, vmax=3000)
        ticks = [1, 3, 10, 30, 100, 300, 1000, 3000]
        ticks = [1, 10, 100, 1000, 10000]
        format = '%i'
        cb_extend = 'both'
        colorbar_label = 'm yr$^{\mathregular{-1}}$'
    elif joughin_speed_300:
        # This is a little duck-punching to get a QGIS colormap
        # similar to Joughin (2010)
        vmin = 1
        vmax = 300
        data_values = np.logspace(vmin, vmax, N)[0:889]
        data_values[-1] = 300
        N = len(data_values)
        norm = mpl.colors.LogNorm(vmin=1, vmax=300)
        ticks = [1, 3, 10, 30, 100, 300]
        format = '%i'
        cb_extend = 'both'
        colorbar_label = 'm yr$^{\mathregular{-1}}$'
    elif joughin_speed_1k:
        # This is a little duck-punching to get a QGIS colormap
        # similar to Joughin (2010)
        data_values = np.logspace(vmin, vmax, N)[0:889]
        data_values[-1] = 1000
        N = len(data_values)
        norm = mpl.colors.LogNorm(vmin=1, vmax=1000)
        ticks = [1, 3, 10, 30, 100, 300, 1000]
        format = '%i'
        cb_extend = 'both'
        # colorbar_label = 'm yr$^{\mathregular{-1}}$'
    elif joughin_speed_10k:
        # This is a little duck-punching to get a QGIS colormap
        # similar to Joughin (2010)
        data_values = np.logspace(vmin, vmax, N)[0:889]
        data_values[-1] = 10000
        N = len(data_values)
        norm = mpl.colors.LogNorm(vmin=1, vmax=10000)
        ticks = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
        format = '%i'
        cb_extend = 'both'
        colorbar_label = 'm yr$^{\mathregular{-1}}$'
    elif habermann_tauc:
        # This is a little duck-punching to get a QGIS colormap
        # similar to Joughin (2010)
        vmin = 1e4
        vmax = 3e5
        data_values = np.logspace(np.log10(vmin), np.log10(vmax), N)
        norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
        ticks = [vmin, vmax]
        format = '%i'
        cb_extend = 'both'
    elif bath_topo:
        # This is a little duck-punching to get a QGIS colormap
        # similar to Joughin (2010)
        vmin = -800
        vmax = 3000
        data_values = a * np.linspace(vmin, vmax, N)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        ticks = None
        format = None
        cb_extend = 'both'
        colorbar_label = 'm a.s.l.'
    elif topo2000m:
        # This is a little duck-punching to get a QGIS colormap
        # similar to Joughin (2010)
        vmin = 0
        vmax = 2000
        data_values = a * np.linspace(vmin, vmax, N)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        ticks = [0, 500, 1000, 1500, 2000]
        format = None
        cb_extend = 'both'
        colorbar_label = 'm a.s.l.'
    elif log:
        data_values = a * np.logspace(vmin, vmax, N)
        norm = mpl.colors.LogNorm(vmin=(10 ** vmin), vmax=a * (10 ** vmax))
        ticks = np.hstack(
            (np.logspace(vmin, vmax, vmax - vmin + 1), a * (10 ** vmax)))
        ticks = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
        format = '%i'
        cb_extend = 'both'
    elif log_color:
        data_values = a * np.logspace(vmin, vmax, N)
        norm = mpl.colors.LogNorm(
            vmin=(10 ** vmin) - 0.01, vmax=a * (10 ** vmax))
        ticks = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
        format = '%i'
        cb_extend = 'both'
    else:
        data_values = a * np.linspace(vmin, vmax, N)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        ticks = [-0.5, -0.1, -0.05, 0, 0.05, 0.1, 0.5]
        ticks = [-250, 0, 250]
        ticks = [0, 500, 1000, 1500, 2000, 2500, 3000]
        ticks = [-500, -250, -50, 0, 50, 250, 500]
        ticks = [-500, -250, 0, 250, 500]
        ticks = [-100, -50, 0, 50, 100]
        format = None
        cb_extend = cb_extend

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, N)

    # you could apply a function to the colormap, e.g. to desaturate the colormap:
    # cmap = cmap_map(lambda x: x/2+0.5, cmap)

    matplotlib.rc("font", **{"sans-serif": ["Helvetica"]})  # "size": fontsize}
    matplotlib.rc("font", **{"size": 20})
    # create the colorbar
    fig = plt.figure()
    if orientation == 'horizontal':
        ax1 = fig.add_axes([0.05, 0.4, 0.4, 0.04])
        # if joughin_speed:
        #     ax1 = fig.add_axes([0.05, 0.4, 0.4, 0.04])
        # else:
        #     ax1 = fig.add_axes([0.05, 0.65, 0.65, 0.04])
    else:
        ax1 = fig.add_axes([0.05, 0.05, 0.03, 0.65])
    cb1 = mpl.colorbar.ColorbarBase(ax1,
                                    cmap=cmap,
                                    norm=norm,
                                    ticks=ticks,
                                    format=format,
                                    extend=cb_extend,
                                    spacing='proportional',
                                    orientation=orientation)

    if colorbar_label:
        cb1.set_label(colorbar_label)

    # save high-res colorbar as png
    for format in ['png']:
        out_file = '.'.join([prefix, format])
        print("  writing colorbar %s ..." % out_file)
        fig.savefig(out_file, bbox_inches='tight', dpi=2400, transparent=True)

    # convert to RGBA array
    rgba = cb1.to_rgba(data_values, alpha=None)
    # QGIS wants 0..255
    rgba *= 255

    # create an output array combining data values and rgb values
    if extend:
        qgis_array = np.zeros((N + 2, 5))
        for k in range(0, N):
            qgis_array[k + 1, 0] = data_values[k]
            qgis_array[k + 1, 1:4] = rgba[k, 0:3]
            qgis_array[k + 1, 4] = 255
        # repeat first color
        qgis_array[0, 0] = extend[0]
        qgis_array[0, 1:4] = rgba[0, 0:3]
        qgis_array[0, 4] = 255
        # repeat last color
        qgis_array[-1, 0] = extend[1]
        qgis_array[-1, 1:4] = rgba[-1, 0:3]
        qgis_array[-1, 4] = 255
    else:
        qgis_array = np.zeros((N, 5))
        for k in range(N):
            qgis_array[k, 0] = data_values[k]
            qgis_array[k, 1:4] = rgba[k, 0:3]
            qgis_array[k, 4] = 255

    # save as ascii file
    out_file = '.'.join([prefix, 'txt'])
    print("  writing colorramp %s ..." % out_file)
    np.savetxt(
        out_file,
        qgis_array,
        delimiter=',',
        fmt=[
            '%10.5f',
            '%i',
            '%i',
            '%i',
            '%i,'])
