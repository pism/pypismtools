#!/usr/bin/env python
import numpy as np
import pylab as plt
from matplotlib import colors, mpl
from argparse import ArgumentParser

try:
        from PyPISMTools import gmtColormap
except:
        from PyPISMTools.PyPISMTools import gmtColormap

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
    for key in ('red' ,'green' ,'blue'):
	    step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(map(reduced_cmap, step_list))
    new_LUT = np.array(map(function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(('red', 'green', 'blue')):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector

    return colors.LinearSegmentedColormap('colormap',cdict,1024)

# Set up the option parser
parser = ArgumentParser()
parser.description = """
A script to convert a GMT (*.cpt) colormap into q QGIS-readable color ramp.
"""
parser.add_argument("FILE", nargs='*')
parser.add_argument("--log",dest="log", action="store_true",
                  help="make a log-normalized color scale", default=False)
parser.add_argument("-a", "--a_log", dest="a", type=float,
                  help='''
                  a * logspace(vmin, vmax, N''', default=1)
parser.add_argument("--vmin", dest="vmin", type=float,
                  help='''
                  a * logspace(vmin, vmax, N''', default=-1)
parser.add_argument("--vmax", dest="vmax", type=float,
                  help='''
                  a * logspace(vmin, vmax, N''', default=3)
parser.add_argument("--extend", dest="extend", nargs=2, type=float,
                  help='''
                  appends color ramp by repeating first and last color for value''',
                    default=None)
parser.add_argument("--N", dest="N", type=int,
                  help='''
                  a * logspace(vmin, vmax, N''', default=1024)
parser.add_argument("-r", "--reverse",dest="reverse", action="store_true",
                  help="reverse color scale", default=False)

options = parser.parse_args()
args = options.FILE
a = options.a
log = options.log
extend = options.extend
N = options.N
vmin = options.vmin
vmax = options.vmax
reverse = options.reverse
# experimental
log_color = False

# read in CPT colormap
cmap_file = args[0]

try:
        cdict = plt.cm.datad[cmap_file]
        prefix = cmap_file
except:
        # import and convert colormap
        cdict = gmtColormap(cmap_file, log_color=log_color, reverse=reverse)
        prefix = '.'.join(cmap_file.split('.')[0:-1])
        suffix = cmap_file.split('.')[-1]


# either log scaling or linear scaling (default)
if log:
	data_values = a * np.logspace(vmin, vmax, N)
	norm = colors.LogNorm(vmin=a * (10 ** vmin), vmax = a * (10 ** vmax))
	ticks = np.hstack((np.logspace(vmin, vmax, vmax - vmin + 1), a * (10 ** vmax)))
	format = '%i'
	cb_extend = 'both'
elif log_color:
	data_values = a * np.logspace(vmin, vmax, N)
	norm = colors.LogNorm(vmin= (10 ** vmin) - 0.01, vmax = a * (10 ** vmax))
	ticks = [1, 10, 100, 1000, 3000]
	format = '%i'
	cb_extend = 'both'
else:
	data_values = a * np.linspace(vmin, vmax, N)
	norm = colors.Normalize(vmin=vmin, vmax=vmax)
	ticks = None
	format = None
	cb_extend = 'both'

cmap = colors.LinearSegmentedColormap('my_colormap', cdict, N)

# you could apply a function to the colormap, e.g. to desaturate the colormap:
# cmap = cmap_map(lambda x: x/2+0.5, cmap)

# create the colorbar
fig = plt.figure()
ax1 = fig.add_axes([0.05, 0.65, 0.9, 0.05])
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm = norm,
				ticks=ticks,
				format=format,
				extend=cb_extend,
                                spacing='proportional',
                                orientation='horizontal')

# save high-res colorbar as png
out_file = '.'.join([prefix, 'png'])
print("  writing colorbar %s ..." % out_file)
fig.savefig(out_file, bbox_inches='tight', dpi=1200)


# convert to RGBA array
rgba = cb1.to_rgba(data_values, alpha=None)
# QGIS wants 0..255
rgba *= 255

# create an output array combining data values and rgb values
if extend:
        qgis_array = np.zeros((N + 1, 5))
        for k in range(1, N):
                qgis_array[k, 0] = data_values[k]
                qgis_array[k, 1:4] = rgba[k, 0:3]
                qgis_array[k, 4] = 255
        # repeat first color
        qgis_array[0, 0] = extend[0]
        qgis_array[0, 1:4] = rgba[1, 0:3]
        qgis_array[0, 4] = 255
        # repeat last color
        qgis_array[N, 0] = extend[1]
        qgis_array[N, 1:4] = rgba[N - 1, 0:3]
        qgis_array[N, 4] = 255
else:
        qgis_array = np.zeros((N, 5))
        for k in range(N):
                qgis_array[k, 0] = data_values[k]
                qgis_array[k, 1:4] = rgba[k, 0:3]
                qgis_array[k, 4] = 255

# save as ascii file
out_file = '.'.join([prefix, 'txt'])
print("  writing colorramp %s ..." % out_file)
np.savetxt(out_file, qgis_array, delimiter=',', fmt=['%10.5f', '%i', '%i', '%i', '%i,'])
