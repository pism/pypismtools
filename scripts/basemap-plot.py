#!/usr/bin/env python

# Copyright (C) 2011 Andy Aschwanden
#
# Script creates a basemap plot of a variable in a netCDF file
# with a geotiff background (if given).
# Does a 1x2, 1x3, 2x2, 3x2 grid plots

from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import pylab as plt
from matplotlib import colors
from optparse import OptionParser

from pyproj import Proj
from osgeo import gdal
from osgeo import osr

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

from PyPISMTools import PyPISMTools as ppt

class Variable(object):
    '''
    A class containing variable-specific stuff such as colorbars, tickmarks, etc
    '''

    def __init__(self, var_name, kwargs):

        self.var_name = varname

        kwargsdict = {}
        expected_args = ['ticks', 'cmap', 'norm', 'vmin', 'vmax', 'extend']
        for key in kwargs.keys():
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

class GeoTIFF(object):
    '''
    A class to read a GeoTIFF

    Parameters
    ----------

    filename: a valid geotiff file
    '''

    def __init__(self, file_name):

        self.file_name = file_name

        try:
            print("\n  opening GeoTIFF file %s" % file_name)
            self.gtiff = gdal.Open(file_name)
        except:
            print("could not open file %s" % file_name)


        self.RasterArray = self.gtiff.ReadAsArray()
        self.gtiff_projection = self.gtiff.GetProjection()

        osr_gtiff = osr.SpatialReference()
        osr_gtiff.ImportFromWkt(self.gtiff_projection)
        self.proj4 = osr_gtiff.ExportToProj4()

        geoT = self.gtiff.GetGeoTransform()
        pxwidth = self.gtiff.RasterXSize
        pxheight = self.gtiff.RasterYSize
        ulx = geoT[0]
        uly = geoT[3]
        rezX = geoT[1]
        rezY = geoT[5]
        rx = ulx + pxwidth * rezX
        ly = uly + pxheight * rezY
        self.width = np.abs(pxwidth * rezX)
        self.height = np.abs(pxheight * rezY)
        self.center_x = ulx + pxwidth * rezX / 2
        self.center_y = uly + pxheight * rezY / 2
        self.easting = np.arange(ulx, rx + rezX, rezX)
        self.northing = np.arange(ly, uly - rezY, -rezY)
        self.X, self.Y = np.meshgrid(self.easting, self.northing)

        p_osr_gtiff = Proj(self.proj4)
        self.lon_0, self.lat_0 = p_osr_gtiff(self.center_x, self.center_y,
                                             inverse=True)
        self.lon, self.lat = p_osr_gtiff(self.X, self.Y, inverse=True)

# Set up the option parser
parser = OptionParser()
parser.usage = "%prog [options] FILE1 FILE2 ..."
parser.description = "A script to plot a variable in a netCDF file over a GeoTiff. Uses GDAL python bindings, Proj4, and Basemap. Script is fine-tuned for whole Greenland plots, but can be adapted for other needs."
parser.add_option("--alpha",dest="alpha",
                  help="transparency of overlay", default=1.)
parser.add_option("--background",dest="background", 
                  help="Draw a background (bluemarble, etopo, shadedrelief", default=None)
parser.add_option("--colormap",dest="colormap",
                  help="path to a cpt colormap", default=None)
parser.add_option("--coastlines",dest="coastlines", action="store_true",
                  help="adds a coastlines", default=False)
parser.add_option("-c", "--colorbar", dest="colorbar",action="store_true",
                  help="saves a colorbar seperately",default=False)
parser.add_option("--inner_title",dest="inner_title",action="store_true",
                  help="add an inner title",default=False)
parser.add_option("--singlerow", dest="singlerow", action="store_true",
                  help="all plots on a single row", default=False)
parser.add_option("-m", "--same_mask", dest="samemask", action="store_true",
                  help="use mask from first plot for all plots", default=False)
parser.add_option("--map_resolution", dest="map_res",
                  help="Resolution of boundary database (see Basemap), default = 'l' (low)", default='l')
parser.add_option("-o", "--output_filename", dest="out_file",
                  help="Name of the output file. Suffix defines output format", default='foo.png')
parser.add_option("--geotiff_file", dest="geotiff_filename",
                  help="GeoTIFF filename", default=None)
parser.add_option("--out_unit", dest="outunit",
                  help="Output unit, default is unit in file", default=None)
parser.add_option("-p", "--print_size", dest="print_mode",
              help="sets figure size and font size, available options are: \
              'onecol','medium','twocol','presentation'",default="twocol")
parser.add_option("-r", "--output_resolution", dest="out_res",
                  help="Graphics resolution in dots per inch (DPI), default = 300",default=300)
parser.add_option("-v", "--variable", dest="varname",
                  help='''Variable to plot, default = 'csurf'.
                  Currently supported variables are: csurf''', default='csurf')


(options, args) = parser.parse_args()
nt = len(args)
required_no_args = 1
max_no_args = 6
if (nt < required_no_args):
    print("received $i arguments, at least %i expected"
          % (nt, required_no_args))
    import sys.exit
    sys.exit
elif (nt > max_no_args):
    print("received $i arguments, no more thant %i accepted"
          % (nt, max_no_args))
    import sys.exit
    sys.exit
else:
    pass

alpha = float(options.alpha)
background = options.background
colormap = options.colormap
coastlines = options.coastlines
colorbar = options.colorbar
inner_title = options.inner_title
map_res = options.map_res
geotiff_filename = options.geotiff_filename
print_mode = options.print_mode
samemask = options.samemask
outunit = options.outunit
out_res = int(options.out_res)
out_file = options.out_file
singlerow = options.singlerow
varname = options.varname

cmap = None
if colormap is not None:
    # import and convert colormap
    cdict = ppt.gmtColormap(colormap)
    cmap = colors.LinearSegmentedColormap('my_colormap', cdict)

# check output format
pre, suffix = out_file.split('.')
if suffix not in ('png', 'pdf', 'ps', 'eps', 'svg'):
    print('Requested output format %s not supported, try png, pdf, svg, ps, eps'
          % suffix)
    import sys.exit
    sys.exit

# set constants and other stuff
meridian_spacing = 10
parallels_spacing = 5

vars_speed = ('csurf', 'cbase', 'cbar', 'magnitude', 'balvelmag', 'surfvelmag')
vars_dem = ('thk', 'usurf', 'usrf')

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
                 'vmin', 'vmax', 'extend')
    attr_vals = ([1, 3, 10, 30, 100, 300, 1000, 3000], cmap,
                 norm, vmin, vmax, 'both')
    var_dict = dict(zip(attr_keys, attr_vals))
    variable = Variable(varname, var_dict)

elif varname in vars_dem:

    if cmap is None:
        cmap = plt.cm.Blues

    vmin = 0.1
    vmax = 3e3
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    attr_keys = ('ticks', 'cmap', 'norm', 'vmin', 'vmax', 'extend')
    attr_vals = ([1000, 2000, 3000], cmap, norm, vmin, vmax, 'max')
    var_dict = dict(zip(attr_keys, attr_vals))
    variable = Variable(varname, var_dict)


if geotiff_filename is not None:
    geotiff = GeoTIFF(geotiff_filename)
    width = geotiff.width
    height = geotiff.height
    lat_0 = geotiff.lat_0
    lon_0 = geotiff.lon_0
    lat = geotiff.lat
    lon = geotiff.lon
else:
    filename = args[0]
    print "  opening NetCDF file %s ..." % filename
    try:
        # open netCDF file in 'append' mode
        nc = NC(filename, 'r')
    except:
        print("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
              % filename)
        
    ## a list of possible x-dimensions names
    xdims = ['x','x1']
    ## a list of possible y-dimensions names
    ydims = ['y','y1']

    ## assign x dimension
    for dim in xdims:
        if dim in nc.dimensions.keys():
            xdim = dim
    ## assign y dimension
    for dim in ydims:
        if dim in nc.dimensions.keys():
            ydim = dim

    ## coordinate variable in x-direction
    x_var = np.squeeze(nc.variables[xdim][:])
    ## coordinate variable in y-direction
    y_var = np.squeeze(nc.variables[ydim][:])

    center_x = (x_var[0] + x_var[-1]) / 2
    center_y = (y_var[0] + y_var[-1]) / 2
    width = 1.2 * (np.max(x_var) - np.min(x_var))
    height = 1.05 * (np.max(y_var) - np.min(y_var))

    nc_projection = ppt.get_projection_from_file(nc)

    lon_0, lat_0 = nc_projection(center_x, center_y, inverse=True)
    nc.close()

print "  creating Basemap ..."
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
var_order = ('time', 'z', ydim, xdim)

for k in range(0, nt):

    filename = args[k]
    print "  opening NetCDF file %s ..." % filename
    try:
        # open netCDF file in 'append' mode
        nc = NC(filename, 'r')
    except:
        print("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
              % filename)
        import sys.exit
        sys.exit

    lats.append(np.squeeze(ppt.permute(nc.variables['lat'], var_order)))
    lons.append(np.squeeze(ppt.permute(nc.variables['lon'], var_order)))

    if varname == 'csurf':
        if 'csurf' in nc.variables.keys():
            var = 'csurf'
        else:
            var = 'magnitude'
    else:
        var = varname
    print("    - reading variable %s from file %s" % (variable.var_name, filename))
    try:
        data = np.squeeze(ppt.permute(nc.variables[variable.var_name], var_order))
    except:
        print("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
              % (variable.var_name, filename))
        exit(2)

    try:
        inunit = str(nc.variables[var].units)
    except:
        print("ERROR:  units not found in variable '%s' in file %s ... ending ..."
              % (variable.var_name, filename))
        exit(2)

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

    nc.close()


# set the print mode
lw, pad_inches = ppt.set_mode(print_mode)


# make a separate colorbar (if requeste
if colorbar:

    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.05, 0.05, 0.9])

    plt.matplotlib.colorbar.ColorbarBase(ax,
                                         cmap=variable.cmap,
                                         norm=variable.norm,
                                         extend=variable.extend,
                                         drawedges=False,
                                         ticks=variable.ticks,
                                         format="%d")

    OUTNAME = var + "_colorbar." + suffix
    print("  writing colorbar %s ..." % OUTNAME)
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
else:
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (2, nt/2), # creates 2 x nt/2 grid of axes
                    axes_pad=0.05, # pad between axes in inch.
                    cbar_mode='single',
                    cbar_size=0.1,
                    cbar_location='top',
                    share_all=True)


for k in range(0,nt):
    ax = grid[k]
    m.ax = ax
    xx, yy = m(lons[k], lats[k])

    if (background == 'bluemable'):
        m.bluemarble()
    elif (background == 'etopo'):
        m.etopo()
    elif (background == 'shadedrelief'):
        m.shadedrelief()
    else:
        pass

    if samemask and (k != 0):
        values[k].mask = values[0].mask
        
    if geotiff_filename is not None:
        m.pcolormesh(xx_gtiff, yy_gtiff, np.flipud(geotiff.RasterArray),
                     cmap=plt.cm.gray)
        
    cs = m.pcolormesh(xx, yy, values[k], cmap=variable.cmap, alpha=alpha,
              norm=variable.norm)
    if singlerow:
        m.drawmeridians(np.arange(-175., 175., meridian_spacing),
                        labels = [0, 0, 0, 1], linewidth=0.5)
        if (k==0):
            m.drawparallels(np.arange(-90., 90., parallels_spacing),
                            labels = [1, 0, 0, 0], linewidth=0.5)
        else:
            m.drawparallels(np.arange(-90., 90., parallels_spacing),
                            labels = [0, 0, 0, 0], linewidth=0.5)
    else:
        if (k==0) or (k==3):
            m.drawparallels(np.arange(-90., 90., parallels_spacing),
                            labels = [1, 0, 0, 0], linewidth=0.5)
        else:
            m.drawparallels(np.arange(-90., 90., parallels_spacing),
                            labels = [0, 0, 0, 0], linewidth=0.5)
        if (k>=3):
            m.drawmeridians(np.arange(-175., 175., meridian_spacing),
                            labels = [0, 0, 0, 1], linewidth=0.5)
        else:
            m.drawparallels(np.arange(-90., 90., meridian_spacing),
                            labels = [0, 0, 0, 0], linewidth=0.5)


    # add coastlines if requested (default is False)
    if coastlines:
        m.drawcoastlines(linewidth=0.25)

    im_titles = ['a)','b)','c)','d)','e)','f)']
    if inner_title:
        for ax in range(0,nt):
            t = ppt.add_inner_title(fig.axes[ax], im_titles[ax], loc=2)
            t.patch.set_ec("none")

## ## Now this is a bit tricky. Without transparency (alpha) set,
## ## we could just do:
## fig.colorbar(cs,cax=grid.cbar_axes[0],
##              orientation='horizontal',
##              extend='variable.extend',
##              drawedges=False,
##              ticks=[1, 10, 100, 1000, 10000],
##              format="%d")
##
## ## With transparency, the colorbar would inhert the transparency,
## ## which we don't want. So we do instead:

if singlerow:
    plt.matplotlib.colorbar.ColorbarBase(fig.axes[nt],
                                     cmap=variable.cmap,
                                     norm=variable.norm,
                                     extend=variable.extend,
                                     orientation='vertical',
                                     drawedges=False,
                                     ticks=variable.ticks,
                                     format="%d")
else:
    plt.matplotlib.colorbar.ColorbarBase(fig.axes[nt],
                                     cmap=variable.cmap,
                                     norm=variable.norm,
                                     extend=variable.extend,
                                     orientation='horizontal',
                                     drawedges=False,
                                     ticks=variable.ticks,
                                     format="%d")

print "  writing image %s ..." % out_file
fig.savefig(out_file,bbox_inches='tight', pad_inches=pad_inches, dpi=out_res)
