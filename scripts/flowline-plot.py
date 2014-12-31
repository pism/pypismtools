#!/usr/bin/env python

# Copyright (C) 2011 Andy Aschwanden

from unidecode import unidecode
import numpy as np
import pylab as plt
from argparse import ArgumentParser
import matplotlib.transforms as transforms
import matplotlib.colors as colors
import matplotlib.cm as cmx

from datetime import datetime

from netcdftime import utime
from netCDF4 import Dataset as NC

try:
    from pypismtools import unit_converter, set_mode, colorList, get_golden_mean, permute, gmtColormap
except:
    from pypismtools.pypismtools import unit_converter, set_mode, colorList, get_golden_mean, permute, gmtColormap

from udunits2 import Converter, System, Unit



# Set up the option parser
parser = ArgumentParser()
parser.description = "A script for profile plots using pylab/matplotlib."
parser.add_argument("FILE", nargs='*')
parser.add_argument("--bounds", dest="bounds", nargs=2, type=float,
                  help="lower and upper bound for ordinate, eg. -1 1", default=None)
parser.add_argument("--x_bounds", dest="x_bounds", nargs=2, type=int,
                  help="lower and upper bound for abscissa, eg. 0 200", default=None)
parser.add_argument("--colormap",dest="colormap",
                  help='''path to a cpt colormap, or a pylab colormap,
                  e.g. Blues''', default='Paired')
parser.add_argument("-l", "--labels",dest="labels",
                  help="comma-separated list with labels, put in quotes like 'label 1,label 2'",default=None)
parser.add_argument("--labelbar_title",dest="labelbar_title",
                  help='''Label bar title''',default=None)
parser.add_argument("--figure_title",dest="figure_title",
                  help='''Figure title''',default=None)
parser.add_argument("-f", "--output_format",dest="out_formats",
                      help="Comma-separated list with output graphics suffix, default = pdf",default='pdf')
parser.add_argument("-o", "--output_file",dest="outfile",
                  help="output file name without suffix, i.e. ts_control -> ts_control_variable",default='foo')
parser.add_argument("-p", "--print_size",dest="print_mode",
                    choices=['onecol','medium','twocol','height','presentation','small_font'],
                    help="sets figure size and font size. Default=medium",default="medium")
parser.add_argument("-s", "--simple", dest="simple", action="store_true",
                    help='''Just draw line''', default=False)
parser.add_argument("-r", "--output_resolution", dest="out_res",
                  help='''Resolution ofoutput graphics in dots per
                  inch (DPI), default = 300''', default=300)
parser.add_argument("-v", "--variable",dest="variables",
                  help="comma-separated list with variables",default='csurf')

options = parser.parse_args()
args = options.FILE
no_args = len(args)
if options.labels != None:
    labels = options.labels.split(',')
else:
    labels = None
bounds = options.bounds
figure_title = options.figure_title
x_bounds = options.x_bounds    
colormap = options.colormap
golden_mean = get_golden_mean()
labelbar_title = options.labelbar_title
out_res = options.out_res
outfile = options.outfile
out_formats = options.out_formats.split(',')
print_mode = options.print_mode
simple = options.simple
variables = options.variables.split(',')
dashes = ['-', '--', '-.', ':', '-', '--', '-.', ':']
output_order = ('profile', 'time')
alpha = 0.5
my_colors = colorList()

dx, dy = 4. / out_res, -4. / out_res

try:
    cdict = plt.cm.datad[colormap]
except:
    # import and convert colormap
    cdict = gmtColormap(colormap)
cmap = colors.LinearSegmentedColormap('my_colormap', cdict)

# Init Unit system
sys = System()

# Plotting styles
axisbg = '0.9'
shadow_color = '0.25'
numpoints = 1

aspect_ratio = golden_mean

# set the print mode
lw, pad_inches = set_mode(print_mode, aspect_ratio=aspect_ratio)

plt.rcParams['legend.fancybox'] = True


ne = len(args)
cNorm = colors.Normalize(vmin=0, vmax=ne)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)




profile_axis_o_units = 'm'

filename = args[0]
print("  opening NetCDF file %s ..." % filename)
try:
    nc0 = NC(filename, 'r')
except:
    print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
          % filename))
    import sys
    sys.exit(1)
profile_names = nc0.variables['profile_name'][:]
profile_axis_name = nc0.variables['profile'].long_name
nc0.close()

for in_varname in variables:

    if in_varname in ('age'):
        o_units = 'yr'
        o_units_str = 'yr'
    else:
        print("variable {} not supported".format(in_varname))


    for profile_id, profile_name in enumerate(profile_names):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        for file_no, filename in enumerate(args):

            print("  opening NetCDF file %s ..." % filename)
            try:
                nc = NC(filename, 'r')
            except:
                print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
                      % filename))
                import sys
                sys.exit(1)

            for name in nc.variables:
                v = nc.variables[name]
                if getattr(v, "standard_name", "") == in_varname:
                    print("variabe {0} found by its standard_name {1}".format(name,
                                                                              in_varname))
                    varname = name
                else:
                    varname = in_varname


            profile_axis = nc.variables['profile'][profile_id]
            profile_axis_units = nc.variables['profile'].units

            profile_axis_o_units = 'm'
            profile_axis = np.squeeze(unit_converter(profile_axis[:], profile_axis_units, profile_axis_o_units))

            z = np.squeeze(nc.variables['z'][:])
            b = np.squeeze(nc.variables['topg'][:])
            s = np.squeeze(nc.variables['usurf'][:])
            ss = np.tile(s, [len(z), 1]).transpose()
            my_var = nc.variables[varname]
            data = np.squeeze(my_var[profile_id, 0, Ellipsis])
            # time is second, for now only time=0

            ## stuff needed for contour plots
            x = profile_axis
            xx = np.squeeze(np.tile(x,[len(z),1])).transpose()
            zz = np.squeeze((np.tile(z,[len(x),1])).transpose() + b).transpose()

            mask = zz > ss
            data_ma = np.ma.array(data=data, mask=mask)
            ax.pcolor(xx, zz, data_ma)
       # plt.close()
