#!/usr/bin/env python

# Copyright (C) 2011-2015 Andy Aschwanden

import numpy as np
import pylab as plt
from argparse import ArgumentParser
import matplotlib.transforms as transforms
import matplotlib.dates as mdates
from datetime import datetime

from netcdftime import utime
from netCDF4 import Dataset as NC

try:
    from pypismtools import unit_converter, set_mode, colorList, get_golden_mean
except:
    from pypismtools.pypismtools import unit_converter, set_mode, colorList, get_golden_mean


# Set up the option parser
parser = ArgumentParser()
parser.description = "A script for PISM output files to time series plots using pylab/matplotlib."
parser.add_argument("FILE", nargs='*')
parser.add_argument("--bounds", dest="bounds", nargs=2, type=float,
                    help="lower and upper bound for ordinate, eg. -1 1", default=None)
parser.add_argument("--time_bounds", dest="time_bounds", nargs=2, type=int,
                    help="lower and upper bound for abscissa, eg. 1990 2000", default=None)
parser.add_argument("-l", "--labels", dest="labels",
                    help="comma-separated list with labels, put in quotes like 'label 1,label 2'", default=None)
parser.add_argument("--index_ij", dest="index_ij", nargs=2, type=int,
                    help="i and j index for spatial fields, eg. 10 10", default=[0, 0])
parser.add_argument("-f", "--output_format", dest="out_formats",
                    help="Comma-separated list with output graphics suffix, default = pdf", default='pdf')
parser.add_argument("-n", "--normalize", dest="normalize", action="store_true",
                    help="Normalize to beginning of time series, Default=False", default=False)
parser.add_argument("-o", "--output_file", dest="outfile",
                    help="output file name without suffix, i.e. ts_control -> ts_control_variable", default='unnamed')
parser.add_argument("-p", "--print_size", dest="print_mode",
                    help="sets figure size and font size, available options are: \
                  'onecol','publish','medium','presentation','twocol'", default="medium")
parser.add_argument("--step", dest="step", type=int,
                    help="step for plotting values, if time-series is very long", default=1)
parser.add_argument("--show", dest="show", action="store_true",
                    help="show figure (in addition to save), Default=False", default=False)
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
parser.add_argument("-v", "--variable", dest="variables",
                    help="comma-separated list with variables", default='ivol')

options = parser.parse_args()
args = options.FILE
if options.labels != None:
    labels = options.labels.split(',')
else:
    labels = None
bounds = options.bounds
index_i, index_j = options.index_ij[0], options.index_ij[1]
time_bounds = options.time_bounds
golden_mean = get_golden_mean()
normalize = options.normalize
out_res = options.out_res
outfile = options.outfile
out_formats = options.out_formats.split(',')
print_mode = options.print_mode
rotate_xticks = options.rotate_xticks
step = options.step
shadow = options.shadow
show = options.show
twinx = options.twinx
variables = options.variables.split(',')
dashes = ['-', '--', '-.', ':', '-', '--', '-.', ':']

dx, dy = 4. / out_res, -4. / out_res

# Conversion between giga tons (Gt) and millimeter sea-level equivalent (mmSLE)
gt2mmSLE = 1. / 365


# Plotting styles
axisbg = '0.9'
shadow_color = '0.25'
numpoints = 1

colors = colorList()

aspect_ratio = golden_mean * .8

# set the print mode
lw, pad_inches = set_mode(print_mode, aspect_ratio=aspect_ratio)

plt.rcParams['legend.fancybox'] = True
no_colors = len(colors)

lines = []
var_dates = []
var_values = []
var_ylabels = []
for var in variables:
    dates = []
    values = []
    for k in range(len(args)):
        print("opening file %s" % args[k])
        nc = NC(args[k], 'r')
        try:
            t = nc.variables["time"][:]
            calendar = nc.variables["time"].calendar
        except:
            t = nc.variables["t"][:]
        try:
            units = nc.variables["time"].units
        except:
            units = nc.variables["t"].units
            calendar = nc.variables["t"].calendar
        try:
            var_units = nc.variables[var].units
        except:
            var_units = None
        # temporary fix for spin-ups with non-Gregorian calendars
        ## calendar = 'gregorian'
        cdftime = utime(units, calendar)
        try:
            date = cdftime.num2date(t[:])
            usedates = True
        except:
            date = t[:]
            usedates = False
            date = np.arange(-len(t[:]), 0) / 1e3
        dates.append(date)

        if var in ("ivol"):
            scale_exponent = 6
            scale = 10 ** scale_exponent
            out_units = "km3"
            var_unit_str = ("10$^{%i}$ km$^{3}$" % scale_exponent)
            ylabel = ("volume [%s]" % var_unit_str)
        elif var in ("imass", "mass", "ocean_kill_flux_cumulative",
                     "grounded_basal_ice_flux_cumulative", "sub_shelf_ice_flux_cumulative", "effective_discharge_flux_cumulative",
                     "surface_ice_flux_cumulative", "nonneg_flux_cumulative",
                     "climatic_mass_balance_flux_cumulative", "discharge_flux_cumulative"):
            out_units = "Gt"
            var_unit_str = "Gt"
            ylabel = ("mass change [%s]" % var_unit_str)
            sle_label = "mass change [mm SLE]"
        elif var in ("dimassdt", "ocean_kill_flux", "surface_ice_flux",
                     "grounded_basal_ice_flux", "sub_shelf_ice_flux",
                     "effective_discharge_flux",
                     "climatic_mass_balance_flux", "nonneg_flux", "discharge_flux"):
            out_units = "Gt year-1"
            var_unit_str = "Gt/yr"
            ylabel = ("mass flux [%s]" % var_unit_str)
            sle_label = "[mm SLE/yr]"
        elif var in ("usurf"):
            out_units = "m"
            var_unit_str = "m a.s.l"
            ylabel = ("elevation [%s]" % var_unit_str)
        elif var in ("iarea"):
            out_units = "km2"
            var_unit_str = "km$^2$"
            ylabel = ("ice area [%s]" % var_unit_str)
        elif var in ("thk"):
            out_units = "m"
            var_unit_str = "m"
            ylabel = ("ice thickness [%s]" % var_unit_str)
        elif var in ("h_x_i", "h_x_j", "h_y_i", "h_y_j", "grad_h"):
            out_units = " "
            var_unit_str = "-"
            ylabel = ("slope [%s]" % var_unit_str)
        elif var in ("eigen1", "eigen2"):
            out_units = "year-1"
            var_unit_str = "a$^{-1}$"
            ylabel = ("strain rate [%s]" % var_unit_str)
        elif var in ("taud", "taud_mag", "taud_x", "taud_y",
                     "bwp", "tauc"):
            out_units = "Pa"
            var_unit_str = "Pa"
            ylabel = ("pressure [%s]" % var_unit_str)
        elif var in ("csurf", "cbase", "cbar", "ubar", "vbar"):
            out_units = "m year-1"
            var_unit_str = "m a$^{-1}$"
            ylabel = ("speed [%s]" % var_unit_str)
        elif var in ("RU"):
            out_units = "km3"
            var_unit_str = "km$^3$"
            ylabel = ("runoff [%s]" % var_unit_str)
        else:
            print("unit %s not recognized" % var_units)
            ylabel = ("%s [%s]" % (var, var_units))
        var_ylabels.append(ylabel)

        if (nc.variables[var].ndim == 3):
            if var_units is not None:
                var_vals = unit_converter(np.squeeze(nc.variables[var][:, index_j, index_i]),
                                          var_units, out_units)
            else:
                var_vals = np.squeeze(nc.variables[var][:, index_j, index_i])
        else:
            if var_units is not None:
                var_vals = unit_converter(
                    np.squeeze(nc.variables[var][:]), var_units, out_units)
            else:
                var_vals = np.squeeze(nc.variables[var][:])
        if normalize:
            var_vals -= var_vals[0]

        values.append(var_vals)
        nc.close()
    var_dates.append(dates)
    var_values.append(values)

for l in range(len(variables)):
    fig = plt.figure()
    offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
    ax = fig.add_subplot(111, axisbg=axisbg)

    if usedates:

        for k in range(len(var_dates[l])):
            n = k % no_colors
            if var in ("ivol"):
                line, = ax.plot_date(
                    var_dates[l][k][::step], var_values[l][k][::step] / scale, color=colors[n])
            else:
                line, = ax.plot_date(
                    var_dates[l][k][::step], var_values[l][k][::step], '-', color=colors[n])
            lines.append(line)

            if shadow:
                shadow_transform = ax.transData + offset
                if var in ("ivol"):
                    ax.plot_date(var_dates[l][k][::step], var_values[l][k][::step] / scale, color=shadow_color, transform=shadow_transform,
                                 zorder=0.5 * line.get_zorder())
                else:
                    ax.plot_date(var_dates[l][k][::step], var_values[l][k][::step], color=shadow_color, transform=shadow_transform,
                                 zorder=0.5 * line.get_zorder())

        nd = len(var_dates[l])
        date_start = []
        date_end = []
        for m in range(0, nd):
            date_start.append(var_dates[l][m][0])
            date_end.append(var_dates[l][m][-1])

        date_start = np.min(np.array(date_start))
        date_end = np.max(np.array(date_end))

        if labels != None:
            ax.legend(lines, labels, bbox_to_anchor=(1., 1.),
                      shadow=True, numpoints=numpoints)

        if twinx:
            axSLE = ax.twinx()
            ax.set_autoscalex_on(False)
            axSLE.set_autoscalex_on(False)
        else:
            fig.autofmt_xdate()
        if bounds:
            ax.set_ylim(bounds[0], bounds[1])
        plt.hold(True)

        yearloc = mdates.YearLocator(10)
        ax.xaxis.set_major_locator(yearloc)

        if time_bounds:
            start_date = datetime(time_bounds[0], 1, 1)
            end_date = datetime(time_bounds[1], 1, 1)
            ax.set_xlim(start_date, end_date)

        ax.set_xlabel('years')
    else:
        for k in range(len(var_dates[l])):
            n = k % no_colors
            if var in ("ivol"):
                line, = ax.plot(
                    var_dates[l][k][::step], var_values[l][k][::step] / scale, color=colors[n])
            else:
                line, = ax.plot(
                    var_dates[l][k][::step], var_values[l][k][::step], '-', color=colors[n])
            lines.append(line)

            if shadow:
                shadow_transform = ax.transData + offset
                if var in ("ivol"):
                    ax.plot(var_dates[l][k][::step], var_values[l][k][::step] / scale, color=shadow_color, transform=shadow_transform,
                            zorder=0.5 * line.get_zorder())
                else:
                    ax.plot(var_dates[l][k][::step], var_values[l][k][::step], color=shadow_color, transform=shadow_transform,
                            zorder=0.5 * line.get_zorder())

        nd = len(var_dates[l])
        date_start = []
        date_end = []
        for m in range(0, nd):
            date_start.append(var_dates[l][m][0])
            date_end.append(var_dates[l][m][-1])

        date_start = np.min(np.array(date_start))
        date_end = np.max(np.array(date_end))

        if labels != None:
            ax.legend(lines, labels, loc="upper right",
                      shadow=True,
                      bbox_to_anchor=(0, 0, 1, 1),
                      bbox_transform=plt.gcf().transFigure)

        if twinx:
            axSLE = ax.twinx()
            ax.set_autoscalex_on(False)
            axSLE.set_autoscalex_on(False)

        ax.set_xlabel('ka BP')

        if bounds:
            ax.set_ylim(bounds[0], bounds[1])

    ymin, ymax = ax.get_ylim()
    if twinx:
        # Plot twin axis on the right, in mmSLE
        yminSLE = ymin * gt2mmSLE
        ymaxSLE = ymax * gt2mmSLE
        axSLE.set_xlim(date_start, date_end)
        axSLE.set_ylim(yminSLE, ymaxSLE)
        axSLE.set_ylabel(sle_label)

    ax.set_ylabel(var_ylabels[l])

    if rotate_xticks:
        ticklabels = ax.get_xticklabels()
        for tick in ticklabels:
            tick.set_rotation(30)
    else:
        ticklabels = ax.get_xticklabels()
        for tick in ticklabels:
            tick.set_rotation(0)

    for out_format in out_formats:
        out_file = outfile + '_' + variables[l] + '.' + out_format
        print "  - writing image %s ..." % out_file
        fig.savefig(out_file, bbox_inches='tight', dpi=out_res)
    if show:
        plt.show()
