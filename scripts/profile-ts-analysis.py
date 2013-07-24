#!/usr/bin/env python

# Copyright (C) 2013 Andy Aschwanden

import numpy as np
import pylab as plt
from argparse import ArgumentParser
import matplotlib.transforms as transforms
import matplotlib.colors as colors
import matplotlib.cm as cmx

from datetime import datetime
import dateutil.parser

from netcdftime import utime
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    from pypismtools import unit_converter, set_mode, colorList, get_golden_mean, permute, trend_estimator
except:
    from pypismtools.pypismtools import unit_converter, set_mode, colorList, get_golden_mean, permute, trend_estimator


# Set up the option parser
parser = ArgumentParser()
parser.description = '''A script for profile time-series plots using pylab/matplotlib.'''
parser.add_argument("FILE", nargs='*')
parser.add_argument("--bounds", dest="bounds", nargs=2, type=float,
                  help="lower and upper bound for ordinate, eg. -1 1", default=None)
parser.add_argument("--x_bounds", dest="x_bounds", nargs=2,
                  help="lower and upper bound for abscissa, eg. 0 200", default=None)
parser.add_argument("-l", "--labels",dest="labels",
                  help="comma-separated list with labels, put in quotes like 'label 1,label 2'",default=None)
parser.add_argument("--index_ij", dest="index_ij", nargs=2, type=float,
                  help="i and j index for spatial fields, eg. 10 10", default=[36, 92])
parser.add_argument("-f", "--output_format",dest="out_formats",
                      help="Comma-separated list with output graphics suffix, default = pdf",default='pdf')
parser.add_argument("-n", "--normalize",dest="normalize",action="store_true",
                  help="Normalize to beginning of time series, Default=False",default=False)
parser.add_argument("-d", "--plot_detrended",dest="plot_detrended",action="store_true",
                  help="Plot detrended time series, Default=False",default=False)
parser.add_argument("-o", "--output_file",dest="outfile",
                  help="output file name without suffix, i.e. ts_control -> ts_control_variable",default='foo')
parser.add_argument("-p", "--print_size",dest="print_mode",
                  help="sets figure size and font size, available options are: \
                  'onecol','publish','medium','presentation','twocol'",default="medium")
parser.add_argument("--show",dest="show",action="store_true",
                  help="show figure (in addition to save), Default=False",default=False)
parser.add_argument("--shadow", dest="shadow", action="store_true",
                  help='''add drop shadow to line plots, Default=False''',
                  default=False)
parser.add_argument("--rotate_xticks", dest="rotate_xticks", action="store_true",
                  help="rotate x-ticks by 30 degrees, Default=False",
                  default=False)
parser.add_argument("-r", "--output_resolution", dest="out_res",
                  help='''Resolution ofoutput graphics in dots per
                  inch (DPI), default = 300''', default=300)
parser.add_argument("-t", "--twinx", dest="twinx", action="store_true",
                  help='''adds a second ordinate with units mmSLE,
                  Default=False''', default=False)
parser.add_argument("-v", "--variable",dest="variables",
                  help="comma-separated list with variables",default='csurf')

options = parser.parse_args()
args = options.FILE
if options.labels != None:
    labels = options.labels.split(',')
else:
    labels = None
bounds = options.bounds
index_i, index_j = options.index_ij[0], options.index_ij[1]
x_bounds = options.x_bounds    
golden_mean = get_golden_mean()
normalize = options.normalize
plot_detrended = options.plot_detrended
out_res = options.out_res
outfile = options.outfile
out_formats = options.out_formats.split(',')
print_mode = options.print_mode
rotate_xticks = options.rotate_xticks
shadow = options.shadow
show = options.show
twinx = options.twinx
variables = options.variables.split(',')
dashes = ['-', '--', '-.', ':', '-', '--', '-.', ':']
output_order = ('profile', 'time')
# stupid CDO changes dimension names...
output_order_cdo = ('ncells', 'time')

dx, dy = 4. / out_res, -4. / out_res

# Conversion between giga tons (Gt) and millimeter sea-level equivalent (mmSLE)
gt2mmSLE = 1. / 365


# Plotting styles
axisbg = '0.9'
shadow_color = '0.25'
numpoints = 1

aspect_ratio = golden_mean

# set the print mode
lw, pad_inches = set_mode(print_mode, aspect_ratio=aspect_ratio)

plt.rcParams['legend.fancybox'] = True

lines = []
profile = []
var_values = []
var_ylabels = []
var_longnames = []
var_units_strings = []
print("opening file %s" % args[0])
nc = NC(args[0],'r')
calendar = nc.variables["time"].calendar
time_units = nc.variables["time"].units
time_outunits = 'years since 1960-1-1'
time_units_str = 'a'
time = nc.variables['time'][:]
t = np.squeeze(unit_converter(time, time_units, time_outunits))
cdftime = utime(time_units, calendar)
date = cdftime.num2date(time)

try:
    profile_names = nc.variables['profile_name'][:]
except:
    pass

profile_var = 'profile'
if profile_var in nc.variables.keys():
    profile = nc.variables["profile"]
    profile_units = profile.units
    profile_outunits = 'km'
    profile_axis = np.squeeze(unit_converter(profile[:], profile_units, profile_outunits))

for var in variables:
    var_units = nc.variables[var].units
    var_longname = nc.variables[var].long_name
    var_longnames.append(var_longname)
    if var in ("ivol"):
        scale_exponent = 6
        scale = 10 ** scale_exponent
        out_units = "km3"
        var_units_str = ("10$^{%i}$ km$^{3}$" % scale_exponent)
        ylabel = ("volume change [%s]" % var_units_str)
    elif var in ("imass", "mass", "ocean_kill_flux_cumulative",
                 "surface_ice_flux_cumulative", "nonneg_flux_cumulative",
                 "climatic_mass_balance_cumulative",
                 "effective_climatic_mass_balance_cumulative",
                 "effective_ice_discharge_cumulative"):
        out_units = "Gt"
        var_units_str = "Gt"
        ylabel = ("mass change [%s]" % var_units_str)
    elif var in ("ocean_kill_flux"):
        out_units = "Gt year-1"
        var_units_str = "Gt a$^{-1}$"
        ylabel = ("mass change [%s]" % var_units_str)
    elif var in ("usurf", "topg"):
        out_units = "m"
        var_units_str = "m a.s.l"
        ylabel = ("elevation [%s]" % var_units_str)
    elif var in ("eigen1", "eigen2"):
        out_units = "year-1"
        var_units_str = "a$^{-1}$"
        ylabel = ("strain rate [%s]" % var_units_str)
    elif var in ("taud", "taud_mag", "taud_x", "taud_y", "bwp", "tauc"):
        out_units = "Pa"
        var_units_str = "Pa"
        ylabel = ("pressure [%s]" % var_units_str)
    elif var in ("csurf", "cbase", "cbar"):
        out_units = "m year-1"
        var_units_str = "m a$^{-1}$"
        ylabel = ("speed [%s]" % var_units_str)
    else:
        print("unit %s not recognized" % var_units)
    var_ylabels.append(ylabel)
    var_units_strings.append(var_units_str)
    try:
        var_vals = unit_converter(np.squeeze(permute(nc.variables[var],
                                                     output_order=output_order)),var_units, out_units)
    except:
        var_vals = unit_converter(np.squeeze(permute(nc.variables[var],
                                                     output_order=output_order_cdo)),var_units, out_units)
    if normalize:
        var_vals -= var_vals[0]

    var_values.append(var_vals)

nc.close()

try:
    no_profiles = len(var_values[0][:,0])
except:
    no_profiles = 1

aspect_ratio = golden_mean

# set the print mode
lw, pad_inches = set_mode(print_mode, aspect_ratio=aspect_ratio)

plt.rcParams['legend.fancybox'] = True


my_colors = colorList()
## color_converter = colors.ColorConverter()
## for c in my_colors:
##     rgb = color_converter.to_rgba_array(c)
##     hsv = colors.rgb_to_hsv(rgb[0][:3])

jet = cm = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=len(date))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

skip = 1
for k in range(len(variables)):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    lines = []
    for m in range(no_profiles):
        try:
            y = var_values[k][m,:]
        except:
            y = var_values[0][:]
        out = trend_estimator(t, y)
        colorVal = scalarMap.to_rgba(k)
        p = out[0]
        cov = out[1]
        trend = p[1]
        amplitude = np.abs(p[2])
        period = p[3]
        trend_err = np.abs(np.sqrt(cov[1][1]) * trend)
        amplitude_err = np.abs(np.sqrt(cov[2][2]) * amplitude)
        units_str = ('%s %s$^{-1}$' % (var_units_strings[k], time_units_str)) 
        if plot_detrended:
            label = ('%4.0f$\pm$%2.0f' % (amplitude, amplitude_err))
            y_detrended = y - (p[0] + p[1]*t)
            retLine, = ax.plot_date(date, y_detrended, fmt='.', color=my_colors[m], label=label)
        else:
            label = ('%4.0f$\pm$%2.0f, %4.0f$\pm$%2.0f' % (trend, trend_err, amplitude, amplitude_err))
            retLine, = ax.plot_date(date, y, fmt='.', color=my_colors[m])
            ax.plot_date(date, p[0] + p[1]*t + p[2] * np.cos(2.0 * np.pi * (t - p[3]) / 1.0),
                         fmt='-', color=my_colors[m], label=label)
            ax.plot_date(date, p[0] + p[1]*t, fmt='--', color=my_colors[m])
        lines.append(retLine)
    ax.set_xlabel("year")
    ax.set_ylabel(var_ylabels[k])
    if plot_detrended:
        plt.legend(numpoints=1, title=('amplitude\n (%s)' % (var_units_strings[k])))
    else:
        plt.legend(numpoints=1, title=('trend, amplitude\n (%s), (%s)' % (units_str, var_units_strings[k])))
        
    if x_bounds:
        x_min = dateutil.parser.parse(x_bounds[0])
        x_max = dateutil.parser.parse(x_bounds[1])
        ax.set_xlim(x_min, x_max)
    if rotate_xticks:
        ticklabels = ax.get_xticklabels()
        for tick in ticklabels:
            tick.set_rotation(30)
    else:
        ticklabels = ax.get_xticklabels()
        for tick in ticklabels:
            tick.set_rotation(0)

    for out_format in out_formats:
        out_file = outfile + '_' + variables[k] + '.' + out_format
        print "  - writing image %s ..." % out_file
        fig.savefig(out_file, bbox_inches='tight', dpi=out_res)


