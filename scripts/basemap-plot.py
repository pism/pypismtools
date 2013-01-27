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
from argparse import ArgumentParser

from pyproj import Proj

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    import PyPISMTools.PyPISMTools as ppt
except:
    import PyPISMTools as ppt


class Variable(object):
    '''
    A class containing variable-specific stuff such as colorbars, tickmarks, etc
    '''

    def __init__(self, var_name, kwargs):

        self.var_name = varname

        kwargsdict = {}
        expected_args = ['ticks', 'cmap', 'norm', 'vmin', 'vmax',
                         'extend', 'format', 'colorbar_label']
        for key in list(kwargs.keys()):
            if key in expected_args:
                kwargsdict[key] = kwargs[key]
                
        if 'ticks' in kwargsdict:
            self.ticks = kwargsdict['ticks']

        if 'cmap' in kwargsdict:
            self.cmap = kwargsdict['cmap']

        if 'norm' in kwargsdict:
            self.norm = kwargsdict['norm']

        if 'vmin' in kwargsdict:
            self.vmin = kwargsdict['vmin']

        if 'vmax' in kwargsdict:
            self.vmax = kwargsdict['vmax']
            
        if 'extend' in kwargsdict:
            self.extend = kwargsdict['extend']

        if 'format' in kwargsdict:
            self.format = kwargsdict['format']

        if 'colorbar_label' in kwargsdict:
            self.colorbar_label = kwargsdict['colorbar_label']

# Set up the option parser
parser = ArgumentParser()
parser.description = "A script to plot a variable in a netCDF file over a GeoTiff. Uses GDAL python bindings, Proj4, and Basemap. Script is fine-tuned for whole Greenland plots, but can be adapted for other needs."
parser.add_argument("FILE", nargs='*')
parser.add_argument("--alpha",dest="alpha",
                  help="transparency of overlay", default=1.)
parser.add_argument("--background",dest="background",
                  help="Draw a background (bluemarble, etopo, shadedrelief", default=None)
parser.add_argument("--bounds", dest="bounds", nargs=2, type=float,
                  help="lower and upper bound for colorbar, eg. -1 1", default=None)
parser.add_argument("--boundary_tol", dest="boundary_tol", nargs=1, type=float,
                  help='''if set, color areas brown where obs <= boundary_tol but data >= boundary_tol,
                  works for difference plots only.''', default=None)
parser.add_argument("--obs_file",dest="obs_file",
                  help='''
                  file with observations for difference plot,
experiment - observation. Must be on same grid as experiments. Default is None''', default=None)
parser.add_argument("--colormap",dest="colormap",
                  help='''path to a cpt colormap, or a pylab colormap,
                  e.g. Blues''', default=None)
parser.add_argument("--coastlines", dest="coastlines", action="store_true",
                  help="adds a coastlines", default=False)
parser.add_argument("-c", "--colorbar", dest="colorbar", action="store_true",
                  help="saves a colorbar seperately", default=False)
parser.add_argument("--colorbar_label", dest="colorbar_label", action="store_true",
                  help="saves a colorbar seperately", default=False)
parser.add_argument("--drawmapscale", dest="drawmapscale", action="store_true",
                  help="draws a map scale in the lower left corner", default=False)
parser.add_argument("--inner_title", dest="inner_title", action="store_true",
                  help="add an inner title", default=False)
parser.add_argument("--singlerow", dest="singlerow", action="store_true",
                  help="all plots on a single row", default=False)
parser.add_argument("--singlecolumn", dest="singlecolumn", action="store_true",
                  help="all plots on a single column", default=False)
parser.add_argument("--map_resolution", dest="map_res",
                  help="Resolution of boundary database (see Basemap), default = 'l' (low)", default='l')
parser.add_argument("-o", "--output_filename", dest="out_file",
                  help="Name of the output file. Suffix defines output format", default='foo.png')
parser.add_argument("--geotiff_file", dest="geotiff_filename",
                  help="GeoTIFF filename", default=None)
parser.add_argument("--shape_file", dest="shape_filename",
                  help="Shapefile filename", default=None)
parser.add_argument("--out_unit", dest="outunit",
                  help="Output unit, default is unit in file", default=None)
parser.add_argument("-p", "--print_size", dest="print_mode",
              help="sets figure size and font size, available options are: \
              'onecol','medium','twocol','presentation'", default="twocol")
parser.add_argument("-r", "--output_resolution", dest="out_res",
                  help='''
                  Graphics resolution in dots per inch (DPI), default
                  = 300''', default=300)
parser.add_argument("--relative", dest="relative", action="store_true",
                  help="do relative differences.", default=False)
parser.add_argument("--no_rasterize", dest="rasterized", action="store_false",
                  help="Don't rasterize plot. Slow.", default=True)
parser.add_argument("--tol", dest="tol", type=float,
                  help="tolerance", default=0.)
parser.add_argument("--level", dest="level", type=int,
                  help="level, for 3D data only. Default = 0", default=0)
parser.add_argument("-s", "--shaded", dest="shaded", action="store_true",
                  help='''Shaded topography. CAREFUL, this options is experimental.
                  It uses imshow, which does not support masked arrays,
                  and we also get the projetion slighly wrong.''', default=False)
parser.add_argument("-v", "--variable", dest="varname",
                  help='''Variable to plot, default = 'csurf'.''', default='csurf')

options = parser.parse_args()
args = options.FILE

nt = len(args)
required_no_args = 0
max_no_args = 6
if (nt < required_no_args):
    print(("received $i arguments, at least %i expected"
          % (nt, required_no_args)))
    import sys.exit
    sys.exit
elif (nt > max_no_args):
    print(("received $i arguments, no more thant %i accepted"
          % (nt, max_no_args)))
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
drawmapscale = options.drawmapscale
inner_title = options.inner_title
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
    cmap = colors.LinearSegmentedColormap('my_colormap', cdict)
            
# check output format
suffix = out_file.split('.')[-1]
if suffix not in ('png', 'pdf', 'ps', 'eps', 'svg'):
    print(('Requested output format %s not supported, try png, pdf, svg, ps, eps'
          % suffix))
    import sys.exit
    sys.exit

# set constants and other stuff
meridian_spacing = 10
parallels_spacing = 5
geotiff_rasterized = True

vars_speed = ('csurf', 'cbase', 'cbar', 'magnitude', 'balvelmag', 'surfvelmag')
vars_dem = ('thk', 'usurf', 'usrf')
vars_topo = ('topg')
vars_dh = ('dhdt', 'climatic_mass_balance_cumulative')
vars_cmb = ('climatic_mass_balance')
vars_temp = ('ice_surface_temp', 'temppabase','temppa')
vars_melt = ("bmelt")
vars_div = ("divHU", "divUH", "divHU_umt", "divHU_cresis", "divHU_searise")

if varname in vars_speed:

    if cmap is None:
        try:
            basedir =  ppt.__file__.split(ppt.__package__)
            cdict = ppt.gmtColormap(basedir[0] + ppt.__package__ +
                                    '/colormaps/Full_saturation_spectrum_CCW.cpt')
            cmap = colors.LinearSegmentedColormap('my_colormap',
        cdict)
        except:
            cmap = plt.cm.Blues

    vmin = 0.3
    vmax = 3e3
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    attr_keys = ('ticks', 'cmap', 'norm',
                 'vmin', 'vmax', 'extend', 'format', 'colorbar_label')
    attr_vals = ([1, 3, 10, 30, 100, 300, 1000, 3000], cmap,
                 norm, vmin, vmax, 'both', '%d', 'm a$^{-1}$')
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_melt:

    cmap = plt.cm.OrRd

    vmin = 0.001
    vmax = 1
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    attr_keys = ('ticks', 'cmap', 'norm',
                 'vmin', 'vmax', 'extend', 'format', 'colorbar_label')
    attr_vals = ([0.001, 0.01, 0.1, 1], cmap,
                 norm, vmin, vmax, 'max', '%3.3f', 'm a$^{-1}$')
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_dem:

    if cmap is None:
        cmap = plt.cm.Blues

    vmin = 0.1
    vmax = None
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ('ticks', 'cmap', 'norm', 'vmin', 'vmax', 'extend', 'format', 'colorbar_label')
    attr_vals = (None, cmap, norm, vmin, vmax, 'max', '%d', 'm')
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

    attr_keys = ('ticks', 'cmap', 'norm', 'vmin', 'vmax', 'extend', 'format', 'colorbar_label')
    attr_vals = (None, cmap, norm, vmin, vmax, 'both', '%d', 'm a.s.l.')
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_dh:

    if cmap is None:
        cmap = plt.cm.RdBu
        
    vmin = None
    vmax = None
    norm = None

    attr_keys = ('ticks', 'vmin', 'vmax', 'norm', 'cmap', 'extend',
                 'format', 'colorbar_label')
    attr_vals = (None, vmin, vmax, norm, cmap, 'both',
                 None, 'm')
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_cmb:

    if cmap is None:
        cmap = plt.cm.RdBu
        
    vmin = None
    vmax = None
    norm = None

    attr_keys = ('ticks', 'vmin', 'vmax', 'norm', 'cmap', 'extend', 'format', 'colorbar_label')
    attr_vals = (None, vmin, vmax, norm, cmap, 'both', None, 'm a$^{-1}$')
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_temp:

    if cmap is None:
        cmap = plt.cm.gist_rainbow_r
        
    vmin = None
    vmax = None
    norm = None

    attr_keys = ('ticks', 'vmin', 'vmax', 'norm', 'cmap', 'extend', 'format', 'colorbar_label')
    attr_vals = (None, vmin, vmax, norm, cmap, 'both', None, u'\u00B0C')
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

elif varname in vars_div:

    if cmap is None:
        cmap = plt.cm.gist_ncar

    vmin = None
    vmax = None
    norm = None

    attr_keys = ('ticks', 'cmap', 'norm', 'vmin', 'vmax', 'extend', 'format', 'colorbar_label')
    attr_vals = (None, cmap, norm, vmin, vmax, 'both', None, 'm a$^{-1}$')
    var_dict = dict(list(zip(attr_keys, attr_vals)))
    variable = Variable(varname, var_dict)

else:

    if cmap is None:
        cmap = plt.cm.gist_ncar

    vmin = None
    vmax = None
    norm = None

    attr_keys = ('ticks', 'cmap', 'norm', 'vmin', 'vmax', 'extend', 'format')
    attr_vals = (None, cmap, norm, vmin, vmax, 'both', None)
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
    variable.norm = colors.Normalize(vmin=variable.vmin,
                                     vmax=variable.vmax)
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
    print("  opening NetCDF file %s ..." % filename)
    try:
        # open netCDF file in 'append' mode
        nc = NC(filename, 'r')
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
              % filename))
        import sys
        sys.exit()
        
    xdim, ydim, zdim, tdim = ppt.get_dims(nc)

    ## coordinate variable in x-direction
    x_var = np.squeeze(nc.variables[xdim][:])
    ## coordinate variable in y-direction
    y_var = np.squeeze(nc.variables[ydim][:])

    ## # Testing, doesn't work yet
    ## try:
    ##     lat1 = np.squeeze(nc.variables["lat"][:])
    ##     lon1 = np.squeeze(nc.variables["lon"][:])
    ##     lat_0 = (lat1[-1,-1] + lat1[0,0]) / 2
    ##     lon_0 = (lon1[0,0] + lon1[0,-1]) / 2
    ##     print lat_0, lon_0
    ## except:
    ##     # FIXME: We do Greenland-tuning:
    center_x = (x_var[0] + x_var[-1]) / 2
    center_y = (y_var[0] + y_var[-1]) / 2
    nc_projection = ppt.get_projection_from_file(nc)
    lon_0, lat_0 = nc_projection(center_x, center_y, inverse=True)
    ##     lon_0 = -45
    width = 1.2 * (np.max(x_var) - np.min(x_var))
    height = 1.0 * (np.max(y_var) - np.min(y_var))
    # This works for Antarctica but not for Greenland:

    ## nc_projection = ppt.get_projection_from_file(nc)
    ## srs_list = nc_projection.srs.split('+')
    ## srs_dict = {}
    ## for x in range(0, len(srs_list)):
    ##     mylist = srs_list[x].split('=')
    ##     try:
    ##         srs_dict[mylist[0]] = float(mylist[1])
    ##     except:
    ##         pass
    ## lat_0 = srs_dict['lat_0']
    ## lon_0 = srs_dict['lon_0']
    
    nc.close()


if obs_file is not None:
    print("  opening NetCDF file %s ..." % obs_file)
    try:
        # open netCDF file in 'append' mode
        nc = NC(obs_file, 'r')
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
              % obs_file))
        import sys
        sys.exit()
        
    # get the dimensions
    xdim, ydim, zdim, tdim = ppt.get_dims(nc)
    # set up dimension ordering
    dim_order = (tdim, zdim, ydim, xdim)

    if varname == 'csurf':
        if 'csurf' in list(nc.variables.keys()):
            var = 'csurf'
        else:
            var = 'magnitude'
    else:
        var = varname
    print(("    - reading variable %s from file %s" % (var, obs_file)))
    try:
        data = np.squeeze(ppt.permute(nc.variables[var], dim_order))
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
              % (variable.var_name, obs_file)))
        exit(2)

    try:
        inunit = str(nc.variables[var].units)
    except:
        print(("ERROR:  units not found in variable '%s' in file %s ... ending ..."
              % (variable.var_name, obs_file)))
        exit(2)

    if outunit is not None:
              data = ppt.unit_converter(data, inunit, outunit)

    if variable.var_name in vars_dem:
        mask = (data <= variable.vmin)
        obs_values = np.ma.array(data, mask = mask)
    elif variable.var_name in vars_topo:
        obs_values = data
    else:
        try:
            fill  = nc.variables[var]._FillValue
            mask = (data == fill)
        except:
            mask = np.zeros_like(data)
        mask[data <= tol] = 1
        obs_values = np.ma.array(data, mask = mask)
    
    nc.close()


print("  creating Basemap ...")
m = Basemap(width=width,
            height=height,
            resolution=map_res,
            projection='stere',
            lat_0=lat_0,
            lon_0=lon_0)

if geotiff_filename is not None:
    xx_gtiff, yy_gtiff = m(lon, lat)

lats = []
lons = []
values = []
ocean_mask = []

for k in range(0, nt):

    filename = args[k]
    print("  opening NetCDF file %s ..." % filename)
    try:
        # open netCDF file in 'append' mode
        nc = NC(filename, 'r')
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
              % filename))
        import sys
        sys.exit(1)

    # get the dimensions
    xdim, ydim, zdim, tdim = ppt.get_dims(nc)
    # set up dimension ordering
    dim_order = (tdim, zdim, ydim, xdim)
    # add lat/lon values
    lats.append(np.squeeze(ppt.permute(nc.variables['lat'], dim_order)))
    lons.append(np.squeeze(ppt.permute(nc.variables['lon'], dim_order)))

    if varname == 'csurf':
        if 'csurf' in list(nc.variables.keys()):
            var = 'csurf'
        else:
            var = 'magnitude'
    else:
        var = varname
    print(("    - reading variable %s from file %s" % (var, filename)))
    try:
        data = np.squeeze(ppt.permute(nc.variables[var], dim_order))
        if (data.ndim == 3):
            data = data[level,:]
    except:
        print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
              % (variable.var_name, filename)))
        import sys
        sys.exit(1)

    try:
        inunit = str(nc.variables[var].units)
    except:
        print(("ERROR:  units not found in variable '%s' in file %s ... ending ..."
              % (variable.var_name, filename)))
        import sys
        sys.exit(1)

    if outunit is not None:
              data = ppt.unit_converter(data, inunit, outunit)

    if variable.var_name in vars_dem:
        mask = (data <= variable.vmin)
        values.append(np.ma.array(data, mask = mask))
    else:
        try:
            fill  = nc.variables[var]._FillValue
            mask = (data == fill)
            values.append(np.ma.array(data, mask = mask))
        except:
            values.append(data)

    ocean_mask_varname = 'mask'
    if ocean_mask_varname in nc.variables.keys():
        ocean_mask.append(np.squeeze(ppt.permute(nc.variables['mask'])))
    else:
        ocean_mask.append(np.zeros_like(data))

    nc.close()


# set the print mode
if print_mode in 'height':
    ntn = nt
    if (ntn==2):
        lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=.75)
    if (ntn==3):
        lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=.55)
    elif (ntn==4):
        lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=.35)
    elif (ntn==5):
        lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=.25)
    else:
        lw, pad_inches = ppt.set_mode(print_mode)
else:
    lw, pad_inches = ppt.set_mode(print_mode)

# make a separate colorbar (if requested)
if colorbar:

    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.05, 0.05, 0.9])

    plt.matplotlib.colorbar.ColorbarBase(ax,
                                         cmap=variable.cmap,
                                         norm=variable.norm,
                                         extend=variable.extend,
                                         drawedges=False,
                                         ticks=variable.ticks,
                                         format=variable.format)

    OUTNAME = var + "_colorbar." + suffix
    print(("  writing colorbar %s ..." % OUTNAME))
    plt.savefig(OUTNAME, bbox_inches='tight')


# create the figure
fig = plt.figure()
if singlerow:
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (1, nt), # creates 1 x nt grid of axes
                    axes_pad=0.05, # pad between axes in inch.
                    cbar_mode='single',
                    cbar_size=0.1,
                    cbar_location='right',
                    share_all=True)
elif singlecolumn:
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (nt, 1), # creates nt x 1 grid of axes
                    axes_pad=0.05, # pad between axes in inch.
                    cbar_mode='single',
                    cbar_size=0.1,
                    cbar_location='top',
                    share_all=True)
else:
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (2, nt/2), # creates 2 x nt/2 grid of axes
                    axes_pad=0.05, # pad between axes in inch.
                    cbar_mode='single',
                    cbar_size=0.1,
                    cbar_location='top',
                    share_all=True)


if variable.var_name not in (vars_speed, vars_dem, vars_topo) and (bounds is None):
    variable.vmin = data.min()
    variable.vmax = data.max()

if bounds:
    variable.norm = colors.Normalize(vmin=variable.vmin, vmax=variable.vmax)
    variable.extend = 'both'
    variable.ticks = None
    variable.format = None

for k in range(0, nt):
    ax = grid[k]
    m.ax = ax
    xx, yy = m(lons[k], lats[k])

    # Draw a background if given
    if (background == 'bluemable'):
        m.bluemarble()
    elif (background == 'etopo'):
        m.etopo()
    elif (background == 'shadedrelief'):
        m.shadedrelief()
    else:
        pass

    # Plot GeoTIFF file if given
    if geotiff_filename is not None:
        if shaded:
            m.imshow(np.flipud(geotiff.RasterArray),
                         cmap=plt.cm.gray, rasterized=geotiff_rasterized)
        else:
            m.pcolormesh(xx_gtiff, yy_gtiff, np.flipud(geotiff.RasterArray),
                         cmap=plt.cm.gray, rasterized=geotiff_rasterized)

    # Draw a boundary mask. Areas where
    #   obs_values <= boundary_tol and values > boundary_tol
    # are colored brown.
    if boundary_tol and obs_file:
        boundary_mask = np.zeros_like(data)
        b_mask = np.ones_like(data)
        b_mask[np.logical_and((obs_values <=  boundary_tol), (values[k] > boundary_tol))] = 0
        b_mask[ocean_mask[k] == 4] = 1
        boundary_mask = np.ma.array(data=boundary_mask, mask=b_mask)
        b = m.pcolormesh(xx, yy, boundary_mask, cmap=plt.cm.BrBG,
                         alpha=alpha, rasterized=rasterized)

    # If observations are given, calculate absolute or relative differences
    if obs_file:
        if relative:
            data = (values[k] - obs_values) / obs_values
            cs = m.pcolormesh(xx, yy, data, cmap=variable.cmap,
                              alpha=alpha, norm=variable.norm,
                              rasterized=rasterized)
        else:
            data = values[k] - obs_values
            cs = m.pcolormesh(xx, yy, data, cmap=variable.cmap,
                              alpha=alpha, norm=variable.norm,
                              rasterized=rasterized)
    else:
        # otherwise just plot data
        data = values[k]
        if shaded:
            from matplotlib.colors import LightSource
            # create light source object.
            ls = LightSource(hsv_min_val=.1, hsv_max_val=.9,
                             hsv_min_sat=.85, hsv_max_sat=.15)
            # convert data to rgb array including shading from light source.
            # (must specify color map)
            data = ls.shade(data, variable.cmap)
            cs = m.imshow(data, cmap=variable.cmap, alpha=alpha,
                          norm=variable.norm, rasterized=rasterized)
        else:
            cs = m.pcolormesh(xx, yy, data, cmap=variable.cmap, alpha=alpha,
                          norm=variable.norm, rasterized=rasterized)

    if singlerow:
        m.drawmeridians(np.arange(-175., 175., meridian_spacing),
                        labels = [0, 0, 0, 1], linewidth=0.5)
        if (k==0):
            m.drawparallels(np.arange(-90., 90., parallels_spacing),
                            labels = [1, 0, 0, 0], linewidth=0.5)
        else:
            m.drawparallels(np.arange(-90., 90., parallels_spacing),
                            labels = [0, 0, 0, 0], linewidth=0.5)
    elif singlecolumn:
        m.drawparallels(np.arange(-90., 90., parallels_spacing),
                        labels = [1, 0, 0, 0], linewidth=0.5)
        if (k==nt-1):
            m.drawmeridians(np.arange(-175., 175., meridian_spacing),
                            labels = [0, 0, 0, 1], linewidth=0.5)
        else:
            m.drawmeridians(np.arange(-175., 175., meridian_spacing),
                            labels = [0, 0, 0, 0], linewidth=0.5)
    else:
        if (k==0) or (k==2):
            m.drawparallels(np.arange(-90., 90., parallels_spacing),
                            labels = [1, 0, 0, 0], linewidth=0.5)
        else:
            m.drawparallels(np.arange(-90., 90., parallels_spacing),
                            labels = [0, 0, 0, 0], linewidth=0.5)
        if (k>=2):
            m.drawmeridians(np.arange(-175., 175., meridian_spacing),
                            labels = [0, 0, 0, 1], linewidth=0.5)
        else:
            m.drawmeridians(np.arange(-90., 90., meridian_spacing),
                            labels = [0, 0, 0, 0], linewidth=0.5)


    # add coastlines if requested (default is False)
    if coastlines:
        m.drawcoastlines(linewidth=0.25)

    im_titles = ['a)','b)','c)','d)','e)','f)']
    if inner_title:
        for ax in range(0, nt):
            t = ppt.add_inner_title(fig.axes[ax], im_titles[ax], loc=2)
            t.patch.set_ec("none")

    if drawmapscale:
        m.drawmapscale(lons[0][0, 0] + 0.75, lats[0][0, 0] + 0.2, lon_0, lat_0,
                   50, fontsize=plt.rcParams['font.size'], barstyle='fancy')

    if shape_filename:
        try:
            m.readshapefile(shape_filename.split('.shp')[0],
                        'my_shapefile', linewidth=1.1)
        except:
            m.readshapefile(shape_filename,
                        'my_shapefile', linewidth=1.1)


if singlerow:
    cbar = plt.matplotlib.colorbar.ColorbarBase(fig.axes[nt],
                                         cmap=variable.cmap,
                                         norm=variable.norm,
                                         extend=variable.extend,
                                         orientation='vertical',
                                         drawedges=False,
                                         ticks=variable.ticks,
                                         format=variable.format)
elif singlecolumn:
    cbar = plt.matplotlib.colorbar.ColorbarBase(fig.axes[nt],
                                         cmap=variable.cmap,
                                         norm=variable.norm,
                                         extend=variable.extend,
                                         orientation='horizontal',
                                         drawedges=False,
                                         ticks=variable.ticks)
else:
    cbar = plt.matplotlib.colorbar.ColorbarBase(fig.axes[nt],
                                         cmap=variable.cmap,
                                         norm=variable.norm,
                                         extend=variable.extend,
                                         orientation='horizontal',
                                         drawedges=False,
                                         ticks=variable.ticks,
                                         format=variable.format)
if colorbar_label:
    cbar.set_label(variable.colorbar_label)

print("  writing image %s ..." % out_file)
# fig.savefig(out_file, bbox_inches='tight', pad_inches=pad_inches, dpi=out_res)
fig.savefig(out_file, bbox_inches='tight', dpi=out_res)

