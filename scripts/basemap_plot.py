#!/usr/bin/env python
# Copyright (C) 2011-2013 Andy Aschwanden
#
# Script creates a basemap plot of a variable in a netCDF file
# with a geotiff background (if given).
# Does a 1x2, 1x3, 2x2, 3x2 grid plots

from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import pylab as plt
from matplotlib import colors
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pyproj import Proj

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt


class Variable(object):

    """
    A class containing variable-specific stuff such as colorbars, tickmarks, etc
    """

    def __init__(self, var_name, kwargs):

        self.var_name = varname

        kwargsdict = {}
        expected_args = ["ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label"]
        for key in list(kwargs.keys()):
            if key in expected_args:
                kwargsdict[key] = kwargs[key]

        if "ticks" in kwargsdict:
            self.ticks = kwargsdict["ticks"]

        if "cmap" in kwargsdict:
            self.cmap = kwargsdict["cmap"]

        if "norm" in kwargsdict:
            self.norm = kwargsdict["norm"]

        if "vmin" in kwargsdict:
            self.vmin = kwargsdict["vmin"]

        if "vmax" in kwargsdict:
            self.vmax = kwargsdict["vmax"]

        if "extend" in kwargsdict:
            self.extend = kwargsdict["extend"]

        if "format" in kwargsdict:
            self.format = kwargsdict["format"]

        if "colorbar_label" in kwargsdict:
            self.colorbar_label = kwargsdict["colorbar_label"]


# Set up the option parser
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "A script to plot a variable in a netCDF file over a GeoTiff. Uses GDAL python bindings, Proj4, and Basemap. Script is fine-tuned for whole Greenland plots, but can be adapted for other needs."
parser.add_argument("FILE", nargs="*")
parser.add_argument("--alpha", dest="alpha", help="transparency of overlay", default=1.0)
parser.add_argument(
    "--background", dest="background", help="Draw a background (bluemarble, etopo, shadedrelief", default=None
)
parser.add_argument(
    "--bounds", dest="bounds", nargs=2, type=float, help="lower and upper bound for colorbar, eg. -1 1", default=None
)
parser.add_argument(
    "--boundary_tol",
    dest="boundary_tol",
    nargs=1,
    type=float,
    help="""if set, color areas brown where obs <= boundary_tol but data >= boundary_tol,
                  works for difference plots only.""",
    default=None,
)
parser.add_argument(
    "--colorbar_position",
    dest="colorbar_position",
    choices=["bottom", "right", "upper", "left"],
    help="position of the colorbar for n x m plots",
    default="bottom",
)
parser.add_argument(
    "--obs_file",
    dest="obs_file",
    help="""
                  file with observations for difference plot,
experiment - observation. Must be on same grid as experiments. Default is None""",
    default=None,
)
parser.add_argument(
    "--colormap",
    dest="colormap",
    help="""path to a cpt colormap, or a pylab colormap,
                  e.g. Blues""",
    default=None,
)
parser.add_argument("--coastlines", dest="coastlines", action="store_true", help="adds a coastlines", default=False)
parser.add_argument(
    "-c", "--colorbar", dest="colorbar", action="store_true", help="saves a colorbar seperately", default=False
)
parser.add_argument(
    "--colorbar_label", dest="colorbar_label", action="store_true", help="saves a colorbar seperately", default=False
)
parser.add_argument(
    "--drawmapscale",
    dest="drawmapscale",
    action="store_true",
    help="draws a map scale in the lower left corner",
    default=False,
)
parser.add_argument(
    "--inner_titles",
    dest="inner_titles",
    help="add an inner title, give a list like --inner_title 'a),b),c)'",
    default=None,
)
parser.add_argument(
    "--singlerow", dest="singlerow", action="store_true", help="all plots on a single row", default=False
)
parser.add_argument(
    "--singlecolumn", dest="singlecolumn", action="store_true", help="all plots on a single column", default=False
)
parser.add_argument(
    "--map_resolution",
    dest="map_res",
    choices=["l", "i", "h", "f"],
    help="Resolution of boundary database (see Basemap), default = 'l' (low)",
    default="l",
)
parser.add_argument(
    "-o",
    "--output_filename",
    dest="out_file",
    help="Name of the output file. Suffix defines output format",
    default="foo.png",
)
parser.add_argument("--geotiff_file", dest="geotiff_filename", help="GeoTIFF filename", default=None)
parser.add_argument("--shape_file", dest="shape_filename", nargs="+", help="Shapefile filename", default=None)
parser.add_argument("--out_unit", dest="outunit", help="Output unit, default is unit in file", default=None)
parser.add_argument(
    "-p",
    "--print_size",
    dest="print_mode",
    choices=["onecol", "medium", "twocol", "height", "presentation", "small_font"],
    help="sets figure size and font size, available options are: \
                    'onecol','medium','twocol','presentation'",
    default="twocol",
)
parser.add_argument(
    "-r",
    "--output_resolution",
    dest="out_res",
    help="""
                  Graphics resolution in dots per inch (DPI), default
                  = 300""",
    default=300,
)
parser.add_argument("--relative", dest="relative", action="store_true", help="do relative differences.", default=False)
parser.add_argument(
    "--no_rasterize", dest="rasterized", action="store_false", help="Don't rasterize plot. Slow.", default=True
)
parser.add_argument("--tol", dest="tol", type=float, help="tolerance", default=None)
parser.add_argument("--level", dest="level", type=int, help="level, for 3D data only. Default = 0", default=0)
parser.add_argument(
    "-s",
    "--shaded",
    dest="shaded",
    action="store_true",
    help="""Shaded topography. CAREFUL, this options is experimental.
                  It uses imshow, which does not support masked arrays,
                  and we also get the projection slighly wrong.""",
    default=False,
)
parser.add_argument(
    "-v", "--variable", dest="varname", help="""Variable to plot, default = 'csurf'.""", default="csurf"
)

options = parser.parse_args()
args = options.FILE

nt = len(args)
required_no_args = 0
max_no_args = 24
if nt < required_no_args:
    print(("received $i arguments, at least %i expected" % (nt, required_no_args)))
    import sys.exit

    sys.exit
elif nt > max_no_args:
    print(("received $i arguments, no more thant %i accepted" % (nt, max_no_args)))
    import sys.exit

    sys.exit
else:
    pass

alpha = float(options.alpha)
background = options.background
bounds = options.bounds
boundary_tol = options.boundary_tol
colormap = options.colormap
coastlines = options.coastlines
colorbar = options.colorbar
colorbar_label = options.colorbar_label
colorbar_position = options.colorbar_position
drawmapscale = options.drawmapscale
if options.inner_titles != None:
    inner_titles = options.inner_titles.split(",")
else:
    inner_titles = None
level = options.level
map_res = options.map_res
geotiff_filename = options.geotiff_filename
print_mode = options.print_mode
obs_file = options.obs_file
outunit = options.outunit
out_res = int(options.out_res)
out_file = options.out_file
shaded = options.shaded
singlerow = options.singlerow
singlecolumn = options.singlecolumn
relative = options.relative
rasterized = options.rasterized
tol = options.tol
varname = options.varname
shape_filename = options.shape_filename

cmap = None
if colormap is not None:
    try:
        cdict = plt.cm.datad[colormap]
    except:
        # import and convert colormap
        cdict = ppt.gmtColormap(colormap)
    cmap = colors.LinearSegmentedColormap("my_colormap", cdict)

# check output format
suffix = out_file.split(".")[-1]
if suffix not in ("png", "pdf", "ps", "eps", "svg"):
    print(("Requested output format %s not supported, try png, pdf, svg, ps, eps" % suffix))
    import sys.exit

    sys.exit

# set constants and other stuff
geotiff_rasterized = True

vars_speed = (
    "csurf",
    "cbase",
    "cbar",
    "magnitude",
    "balvelmag",
    "surfvelmag",
    "velbase_mag",
    "velsurf_mag",
    "velshear_mag",
)
vars_dem = ("thk", "usurf", "usrf", "surface_altitude", "surface", "land_ice_thickness")
vars_topo = ("topg", "bedrock_altitude", "bed")
vars_dh = ("dhdt", "climatic_mass_balance_cumulative")
vars_cmb = ("climatic_mass_balance", "climatic_mass_balance_original")
vars_temp = ("ice_surface_temp", "temppabase", "temppa", "temp_pa")
vars_melt = "bmelt"
vars_heat = "bheatflx"
vars_div = ("divQ", "divHU", "divUH", "divHU_umt", "divHU_cresis", "divHU_searise", "res_flux")
vars_tempice = "tempicethk_basal"
vars_stress = ("tauc", "tauc_mag", "taub_mag")
vars_hydro = "tillwat"
vars_hydro_log = "bwat"
vars_ratio_1 = "sliding_r"
vars_rel = "tau_rel"
vars_rel_log = "tau_r"

if varname in vars_speed:

    if cmap is None:
        try:
            basedir = ppt.__file__.split(ppt.__package__)
            cdict = ppt.gmtColormap(
                basedir[0] + ppt.__package__ + "/colormaps/Full_saturation_spectrum_CCW_orange.cpt"
            )
            cmap = colors.LinearSegmentedColormap("my_colormap", cdict)
        except:
            cmap = plt.cm.Blues

    vmin = 1.0
    vmax = 3e3
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = (
        [1, 3, 10, 30, 100, 300, 1000, 3000],
        cmap,
        norm,
        vmin,
        vmax,
        "both",
        "%d",
        "m yr$^{\mathregular{-1}}$",
    )
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_melt:

    if cmap is None:
        cmap = plt.cm.OrRd

    vmin = 0.001
    vmax = 1
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = ([0.001, 0.01, 0.1, 1], cmap, norm, vmin, vmax, "max", None, "m yr$^{\mathregular{-1}}$")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_heat:

    if cmap is None:
        cmap = plt.cm.jet

    vmin = 10
    vmax = 150
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = (None, cmap, norm, vmin, vmax, "both", None, "W m$^{\mathregular{-2}}$")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_stress:

    if cmap is None:
        cmap = plt.cm.jet

    vmin = 1e4
    vmax = 1.25e6
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = ([0, 1e3, 1e4, 1e5, 1e6, 1e7], cmap, norm, vmin, vmax, "both", "%1.0e", "Pa")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_tempice:

    if cmap is None:
        cmap = plt.cm.OrRd

    vmin = 0.1
    vmax = 100
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = ([25, 50, 75, 100], cmap, norm, vmin, vmax, "max", "%i", "m")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_dem:

    if cmap is None:
        cmap = plt.cm.Blues

    vmin = 0.1
    vmax = None
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = (None, cmap, norm, vmin, vmax, "max", "%d", "m")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_topo:

    if cmap is None:
        cmap = plt.cm.Blues

    ## vmin = -5000
    ## vmax = 1400
    vmin = -1000
    vmax = 2100
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = (None, cmap, norm, vmin, vmax, "both", "%d", "m a.s.l.")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_dh:

    if cmap is None:
        cmap = plt.cm.RdBu

    vmin = None
    vmax = None
    norm = None

    attr_keys = ("ticks", "vmin", "vmax", "norm", "cmap", "extend", "format", "colorbar_label")
    attr_vals = (None, vmin, vmax, norm, cmap, "both", None, "m")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_cmb:

    if cmap is None:
        cmap = plt.cm.RdBu

    vmin = None
    vmax = None
    norm = None

    attr_keys = ("ticks", "vmin", "vmax", "norm", "cmap", "extend", "format", "colorbar_label")
    attr_vals = (None, vmin, vmax, norm, cmap, "both", None, "kg m$^{\mathregular{2}}$ yr$^{\mathregular{-1}}$")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_temp:

    if cmap is None:
        cmap = plt.cm.gist_rainbow_r

    vmin = None
    vmax = None
    norm = None

    attr_keys = ("ticks", "vmin", "vmax", "norm", "cmap", "extend", "format", "colorbar_label")
    attr_vals = (None, vmin, vmax, norm, cmap, "both", None, "\u00B0C")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_div:

    if cmap is None:
        cmap = plt.cm.gist_ncar

    vmin = None
    vmax = None
    norm = None

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = (None, cmap, norm, vmin, vmax, "both", None, "m yr$^{\mathregular{-1}}$")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_hydro:

    if cmap is None:
        cmap = plt.cm.jet

    vmin = 0.001
    vmax = 2
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = (None, cmap, norm, vmin, vmax, "max", None, "m")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_hydro_log:

    if cmap is None:
        cmap = plt.cm.jet

    vmin = 0.001
    vmax = 10
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = (None, cmap, norm, vmin, vmax, "max", None, "m")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_ratio_1:

    if cmap is None:
        cmap = plt.cm.OrRd

    vmin = 0
    vmax = 1
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = ([-10, -1, 0, 1, 10], cmap, norm, vmin, vmax, "neither", "%i", "1")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_rel_log:

    if cmap is None:
        cmap = plt.cm.gist_ncar_r

    vmin = 1
    vmax = 1000
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = ([1, 10, 100, 1000], cmap, norm, vmin, vmax, "both", "%i", "1")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_rel:

    if cmap is None:
        cmap = plt.cm.PRGn

    vmin = -10
    vmax = 10
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format", "colorbar_label")
    attr_vals = ([-10, -1, 0, 1, 10], cmap, norm, vmin, vmax, "both", "%i", "1")
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

else:

    if cmap is None:
        cmap = plt.cm.gist_ncar

    vmin = None
    vmax = None
    norm = None

    attr_keys = ("ticks", "cmap", "norm", "vmin", "vmax", "extend", "format")
    attr_vals = (None, cmap, norm, vmin, vmax, "both", None)
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

bounds_min = -1
bounds_max = 1
if bounds is not None:
    bounds_min = bounds[0]
    bounds_max = bounds[1]
    variable.vmin = bounds_min
    variable.vmax = bounds_max
    variable.norm = colors.Normalize(vmin=variable.vmin, vmax=variable.vmax)

if obs_file is not None:
    variable.vmin = bounds_min
    variable.vmax = bounds_max
    variable.norm = colors.Normalize(vmin=variable.vmin, vmax=variable.vmax)
    variable.ticks = None

if geotiff_filename is not None:
    geotiff = ppt.GeoTIFF(geotiff_filename)
    width = geotiff.width
    height = geotiff.height
    lat_0 = geotiff.lat_0
    lon_0 = geotiff.lon_0
    lat = geotiff.lat
    lon = geotiff.lon
else:
    filename = args[0]
    print(("  opening NetCDF file %s ..." % filename))
    try:
        nc = NC(filename, "r")
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..." % filename))
        import sys

        sys.exit()

    xdim, ydim, zdim, tdim = ppt.get_dims(nc)

    # coordinate variable in x-direction
    x_var = np.squeeze(nc.variables[xdim][:])
    # coordinate variable in y-direction
    y_var = np.squeeze(nc.variables[ydim][:])

    center_x = (x_var[0] + x_var[-1]) / 2
    center_y = (y_var[0] + y_var[-1]) / 2
    nc_projection = ppt.get_projection_from_file(nc)
    lon_0, lat_0 = nc_projection(center_x, center_y, inverse=True)
    width = 1.2 * (np.max(x_var) - np.min(x_var))
    height = 1.0 * (np.max(y_var) - np.min(y_var))

    nc.close()


if obs_file is not None:
    print(("  opening NetCDF file %s ..." % obs_file))
    try:
        # open netCDF file in 'append' mode
        nc = NC(obs_file, "r")
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..." % obs_file))
        import sys

        sys.exit()

    # get the dimensions
    xdim, ydim, zdim, tdim = ppt.get_dims(nc)
    # set up dimension ordering
    dim_order = (tdim, zdim, ydim, xdim)

    myvar = varname
    for name in list(nc.variables.keys()):
        v = nc.variables[name]
        if getattr(v, "standard_name", "") == varname:
            print(("variabe {0} found by its standard_name {1}".format(name, varname)))
            myvar = name
    print(("    - reading variable %s from file %s" % (myvar, obs_file)))
    try:
        data = np.squeeze(ppt.permute(nc.variables[myvar], dim_order))
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..." % (variable.var_name, obs_file)))
        exit(2)

    try:
        inunit = str(nc.variables[myvar].units)
    except:
        print(("ERROR:  units not found in variable '%s' in file %s ... ending ..." % (variable.var_name, obs_file)))
        exit(2)

    if outunit is not None:
        data = ppt.unit_converter(data, inunit, outunit)

    if variable.var_name in vars_dem:
        mask = data <= variable.vmin
        obs_values = np.ma.array(data, mask=mask)
    elif variable.var_name in vars_topo:
        obs_values = data
    else:
        try:
            fill = nc.variables[var]._FillValue
            mask = data == fill
        except:
            mask = np.zeros_like(data)
            mask[data <= tol] = 1
        if tol:
            mask[data <= tol] = 1
        obs_values = np.ma.array(data, mask=mask)

    nc.close()


print("  creating Basemap ...")
m = Basemap(width=width, height=height, resolution=map_res, projection="stere", lat_0=lat_0, lon_0=lon_0)

if geotiff_filename is not None:
    xx_gtiff, yy_gtiff = m(lon, lat)

lats = []
lons = []
values = []
ocean_mask = []

for k in range(0, nt):

    filename = args[k]
    print(("  opening NetCDF file %s ..." % filename))
    try:
        # open netCDF file in 'append' mode
        nc = NC(filename, "r")
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..." % filename))
        import sys

        sys.exit(1)

    # get the dimensions
    xdim, ydim, zdim, tdim = ppt.get_dims(nc)
    # set up dimension ordering
    dim_order = (tdim, zdim, ydim, xdim)
    # add lat/lon values
    lats.append(np.squeeze(ppt.permute(nc.variables["lat"], dim_order)))
    lons.append(np.squeeze(ppt.permute(nc.variables["lon"], dim_order)))

    myvar = varname
    for name in list(nc.variables.keys()):
        v = nc.variables[name]
        if getattr(v, "standard_name", "") == varname:
            print(("variabe {0} found by its standard_name {1}".format(name, varname)))
            myvar = name
            pass
    print(("    - reading variable %s from file %s" % (myvar, filename)))
    try:
        data = np.squeeze(ppt.permute(nc.variables[myvar], dim_order))
        if data.ndim == 3:
            data = data[level, :]
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..." % (variable.var_name, filename)))
        import sys

        sys.exit(1)

    try:
        inunit = str(nc.variables[myvar].units)
    except:
        print(("ERROR:  units not found in variable '%s' in file %s ... ending ..." % (myvar, filename)))
        import sys

        sys.exit(1)

    if outunit is not None:
        data = ppt.unit_converter(data, inunit, outunit)

    if variable.var_name in vars_dem:
        mask = data <= variable.vmin
        values.append(np.ma.array(data, mask=mask))
    else:
        try:
            fill = nc.variables[var]._FillValue
            mask = data == fill
            values.append(np.ma.array(data, mask=mask))
        except:
            values.append(data)

    ocean_mask_varname = "oceanmask"
    if ocean_mask_varname in list(nc.variables.keys()):
        ocean_mask.append(np.squeeze(ppt.permute(nc.variables["mask"])))
    else:
        ocean_mask.append(np.zeros_like(data))
    nc.close()


# set the print mode
if print_mode in "height":
    ntn = nt
    if ntn == 2:
        lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=0.75)
    if ntn == 3:
        lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=0.55)
    elif ntn == 4:
        lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=0.35)
    elif ntn == 5:
        lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=0.25)
    else:
        lw, pad_inches = ppt.set_mode(print_mode)
else:
    lw, pad_inches = ppt.set_mode(print_mode)

# make a separate colorbar (if requested)
if colorbar:

    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.05, 0.05, 0.9])

    plt.matplotlib.colorbar.ColorbarBase(
        ax,
        cmap=variable.cmap,
        norm=variable.norm,
        extend=variable.extend,
        drawedges=False,
        ticks=variable.ticks,
        format=variable.format,
    )

    OUTNAME = var + "_colorbar." + suffix
    print(("  writing colorbar %s ..." % OUTNAME))
    plt.savefig(OUTNAME, bbox_inches="tight")


# create the figure
fig = plt.figure()
if singlerow:
    grid = ImageGrid(
        fig,
        111,  # similar to subplot(111)
        nrows_ncols=(1, nt),  # creates 1 x nt grid of axes
        axes_pad=0.05,  # pad between axes in inch.
        cbar_mode="single",
        cbar_size=0.115,
        cbar_location="right",
        share_all=True,
    )
elif singlecolumn:
    grid = ImageGrid(
        fig,
        111,  # similar to subplot(111)
        nrows_ncols=(nt, 1),  # creates nt x 1 grid of axes
        axes_pad=0.05,  # pad between axes in inch.
        cbar_mode="single",
        cbar_size=0.115,
        cbar_location="bottom",
        share_all=True,
    )
else:
    grid = ImageGrid(
        fig,
        111,  # similar to subplot(111)
        nrows_ncols=(nt / 3, 3),  # creates 2 x nt/2 grid of axes
        axes_pad=0.1,  # pad between axes in inch.
        cbar_mode="single",
        cbar_size=0.115,
        cbar_location=colorbar_position,
        share_all=True,
    )


if variable.var_name not in (vars_speed, vars_dem, vars_topo) and (bounds is None):
    variable.vmin = data.min()
    variable.vmax = data.max()

if bounds:
    variable.norm = colors.Normalize(vmin=variable.vmin, vmax=variable.vmax)
    variable.extend = "both"
    variable.ticks = None
    variable.format = None

for k in range(0, nt):
    ax = grid[k]
    m.ax = ax
    xx, yy = m(lons[k], lats[k])

    # Draw a background if given
    if background == "bluemable":
        m.bluemarble()
    elif background == "etopo":
        m.etopo()
    elif background == "shadedrelief":
        m.shadedrelief()
    else:
        pass

    # Plot GeoTIFF file if given
    if geotiff_filename is not None:
        if shaded:
            m.imshow(np.flipud(geotiff.RasterArray), cmap=plt.cm.gray, rasterized=geotiff_rasterized)
        else:
            m.pcolormesh(
                xx_gtiff, yy_gtiff, np.flipud(geotiff.RasterArray), cmap=plt.cm.gray, rasterized=geotiff_rasterized
            )

    # Draw a boundary mask. Areas where
    #   obs_values <= boundary_tol and values > boundary_tol
    # are colored brown.
    if boundary_tol and obs_file:
        boundary_mask = np.zeros_like(data)
        b_mask = np.ones_like(data)
        b_mask[np.logical_and((obs_values <= boundary_tol), (values[k] > boundary_tol))] = 0
        b_mask[ocean_mask[k] == 4] = 1
        boundary_mask = np.ma.array(data=boundary_mask, mask=b_mask)
        b = m.pcolormesh(xx, yy, boundary_mask, cmap=plt.cm.BrBG, alpha=alpha, rasterized=rasterized)

    # If observations are given, calculate absolute or relative differences
    if obs_file:
        if relative:
            data = (values[k] - obs_values) / obs_values
            cs = m.pcolormesh(xx, yy, data, cmap=variable.cmap, alpha=alpha, norm=variable.norm, rasterized=rasterized)
        else:
            data = values[k] - obs_values
            cs = m.pcolormesh(xx, yy, data, cmap=variable.cmap, alpha=alpha, norm=variable.norm, rasterized=rasterized)
    else:
        # otherwise just plot data
        data = values[k]
        if shaded:
            from matplotlib.colors import LightSource

            # create light source object.
            lightsource = LightSource(hsv_min_val=0.1, hsv_max_val=0.9, hsv_min_sat=0.85, hsv_max_sat=0.15)
            # convert data to rgba array including shading from light source.
            # (must specify color map)
            data = lightsource.shade(data, variable.cmap)
            cs = m.imshow(data, cmap=variable.cmap, alpha=alpha, norm=variable.norm, rasterized=rasterized)
        else:
            cs = m.pcolormesh(xx, yy, data, cmap=variable.cmap, alpha=alpha, norm=variable.norm, rasterized=rasterized)

    dlat = np.abs(lats[k][-1, -1] - lats[k][0, 0])
    dlon = np.abs(lons[k][-1, -1] - lons[k][0, 0])
    if dlat > 20:
        parallels_spacing = 5
    elif (dlat > 10) and (dlat <= 20):
        parallels_spacing = 2
    elif (dlat > 5) and (dlat <= 10):
        parallels_spacing = 1
    elif (dlat > 5) and (dlat <= 1):
        parallels_spacing = 1
    elif (dlat > 5) and (dlat <= 1):
        parallels_spacing = 1
    else:
        parallels_spacing = 0.5
    if dlon > 20:
        meridian_spacing = 10
    elif (dlon > 11) and (dlon <= 20):
        meridian_spacing = 5
    elif (dlon > 6) and (dlon <= 11):
        meridian_spacing = 2
    elif (dlon > 2) and (dlon <= 6):
        meridian_spacing = 1
    elif (dlon > 0.6) and (dlon <= 2):
        meridian_spacing = 0.5
    else:
        meridian_spacing = 0.5

    if singlerow:
        m.drawmeridians(np.arange(-175.0, 175.0, meridian_spacing), labels=[0, 0, 0, 1], linewidth=0.5)
        if k == 0:
            m.drawparallels(np.arange(-90.0, 90.0, parallels_spacing), labels=[1, 0, 0, 0], linewidth=0.5)
        else:
            m.drawparallels(np.arange(-90.0, 90.0, parallels_spacing), labels=[0, 0, 0, 0], linewidth=0.5)
    elif singlecolumn:
        m.drawparallels(np.arange(-90.0, 90.0, parallels_spacing), labels=[1, 0, 0, 0], linewidth=0.5)
        if k == nt - 1:
            m.drawmeridians(np.arange(-175.0, 175.0, meridian_spacing), labels=[0, 0, 0, 1], linewidth=0.5)
        else:
            m.drawmeridians(np.arange(-175.0, 175.0, meridian_spacing), labels=[0, 0, 0, 0], linewidth=0.5)
    else:
        if (k == 0) or (k == 3):
            m.drawparallels(np.arange(-90.0, 90.0, parallels_spacing), labels=[1, 0, 0, 0], linewidth=0.5)
        else:
            m.drawparallels(np.arange(-90.0, 90.0, parallels_spacing), labels=[0, 0, 0, 0], linewidth=0.5)
        if k >= 3:
            m.drawmeridians(np.arange(-175.0, 175.0, meridian_spacing), labels=[0, 0, 0, 1], linewidth=0.5)
        else:
            m.drawmeridians(np.arange(-90.0, 90.0, meridian_spacing), labels=[0, 0, 0, 0], linewidth=0.5)

    # add coastlines if requested (default is False)
    if coastlines:
        m.drawcoastlines(linewidth=0.25)

    if inner_titles:
        for ax in range(0, nt):
            t = ppt.add_inner_title(fig.axes[ax], inner_titles[ax], loc=2)
            t.patch.set_ec("none")

    if drawmapscale:
        x_c = m.llcrnrx + np.abs(m.urcrnrx - m.llcrnrx) * 0.15
        y_c = m.llcrnry + np.abs(m.urcrnry - m.llcrnry) * 0.075
        lon_c, lat_c = m(x_c, y_c, inverse=True)
        ms_width = np.abs(m.urcrnrx - m.llcrnrx) * 0.2 / 1e3
        m.drawmapscale(
            lon_c, lat_c, lon_0, lat_0, ms_width, units="km", fontsize=plt.rcParams["font.size"], barstyle="fancy"
        )

    contour_colors = ["white", "black"]
    if shape_filename:
        for index, shpfile in enumerate(shape_filename):
            shpfile = shpfile.split(".shp")[0]
            m.readshapefile(shpfile, "my_shapefile", linewidth=0.75, color=contour_colors[index])
            # m.readshapefile(shape_filename.split('.shp')[0],
            #             'my_shapefile', linewidth=1.1)


if singlerow:
    cbar = plt.matplotlib.colorbar.ColorbarBase(
        fig.axes[nt],
        cmap=variable.cmap,
        norm=variable.norm,
        extend=variable.extend,
        orientation="vertical",
        drawedges=False,
        ticks=variable.ticks,
        format=variable.format,
    )
elif singlecolumn:
    cbar = plt.matplotlib.colorbar.ColorbarBase(
        fig.axes[nt],
        cmap=variable.cmap,
        norm=variable.norm,
        extend=variable.extend,
        orientation="horizontal",
        drawedges=False,
        ticks=variable.ticks,
    )
else:
    if colorbar_position in ("bottom", "upper"):
        orientation = "horizontal"
    else:
        orientation = "vertical"
    cbar = plt.matplotlib.colorbar.ColorbarBase(
        fig.axes[nt],
        cmap=variable.cmap,
        norm=variable.norm,
        extend=variable.extend,
        orientation=orientation,
        drawedges=False,
        ticks=variable.ticks,
        format=variable.format,
    )

# to prevent the pdf file having white lines
cbar.solids.set_edgecolor("face")
if colorbar_label:
    cbar.set_label(variable.colorbar_label)

print(("  writing image %s ..." % out_file))
# fig.savefig(out_file, bbox_inches='tight', dpi=out_res, pad_inches=pad_inches)
fig.savefig(out_file, bbox_inches="tight", dpi=out_res)

plt.close()
del fig
