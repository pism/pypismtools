'''
pypismtools: Tools to evaluate PISM parameter studies

pypismtools is a module to facilitate evaluation of PISM parameter
studies. It mainly comprises two classes, Observation and Experiment,
which act as containers for observational data and PISM model
simulations, along with helper functions. The experiment class
determines information about an experiment from the netcdf file
directly, especially from the "pism_overrides" flag. Such information
can then be used for labeling, plotting, evaluation, etc. The indend
is to provide a robust tool to evaluate data, and to avoid common mistakes
such as mis-labeling plots. Additional functions include routines to
permute (netcdf) dimension, convert units using udunits, to estimate
trends, and to import GMT colormaps.
'''

__author__ = "Andy Aschwanden"

import numpy as np
import pylab as plt

from netCDF4 import Dataset as CDF

from pyproj import Proj
from osgeo import gdal
from osgeo import osr


# FIXME: how to provide DEBUG flag to module
DEBUG = None


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


def get_dims(nc):
    '''
    Gets dimensions from netcdf instance

    Parameters:
    -----------
    nc: netCDF instance

    Returns:
    --------
    xdim, ydim, zdim, tdim: dimensions
    '''

    # a list of possible x-dimensions names
    xdims = ['x', 'x1']
    # a list of possible y-dimensions names
    ydims = ['y', 'y1']
    # a list of possible z-dimensions names
    zdims = ['z', 'z1']
    # a list of possible time-dimensions names
    tdims = ['t', 'time']

    xdim = None
    ydim = None
    zdim = None
    tdim = None

    # assign x dimension
    for dim in xdims:
        if dim in list(nc.dimensions.keys()):
            xdim = dim
    # assign y dimension
    for dim in ydims:
        if dim in list(nc.dimensions.keys()):
            ydim = dim
    # assign z dimension
    for dim in zdims:
        if dim in list(nc.dimensions.keys()):
            zdim = dim
    # assign time dimension
    for dim in tdims:
        if dim in list(nc.dimensions.keys()):
            tdim = dim
    return xdim, ydim, zdim, tdim


def get_projection_from_file(nc):
    '''
    Gets a Proj projection instance from a pointer to a netCDF file

    Parameters
    ----------
    nc : a netCDF object instance

    Returns
    -------
    p : Proj4 projection instance
    '''

    from pyproj import Proj

    # First, check if we have a global attribute 'proj4'
    # which contains a Proj4 string:
    try:
        p = Proj(str(nc.proj4))
        print(
            'Found projection information in global attribute proj4, using it')
    except:
        try:
            p = Proj(str(nc.projection))
            print(
                'Found projection information in global attribute projection, using it')
        except:
            try:
                # go through variables and look for 'grid_mapping' attribute
                for var in nc.variables.keys():
                    if hasattr(nc.variables[var], 'grid_mapping'):
                        mappingvarname = nc.variables[var].grid_mapping
                        print(
                            'Found projection information in variable "%s", using it' %
                            mappingvarname)
                        break
                var_mapping = nc.variables[mappingvarname]
                p = Proj(proj="stere",
                         ellps=var_mapping.ellipsoid,
                         datum=var_mapping.ellipsoid,
                         units="m",
                         lat_ts=var_mapping.standard_parallel,
                         lat_0=var_mapping.latitude_of_projection_origin,
                         lon_0=var_mapping.straight_vertical_longitude_from_pole,
                         x_0=var_mapping.false_easting,
                         y_0=var_mapping.false_northing)
            except:
                print('No mapping information found, exiting.')
                import sys
                sys.exit(1)

    return p


def add_inner_title(ax, title, loc, size=None, **kwargs):
    '''
    Adds an inner title to a given axis, with location loc.

    from http://matplotlib.sourceforge.net/examples/axes_grid/demo_axes_grid2.html
    '''
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    return at


def get_golden_mean():
    '''
    Returns golden mean (sqrt(5) - 1.0) / 2.0
    '''
    return (np.sqrt(5) - 1.0) / 2.0


def set_mode(mode, aspect_ratio=0.95):
    '''
    Set the print mode, i.e. document and font size. Options are:
    - onecol: width=80mm, font size=8pt. Appropriate for 1-column figures
    - twocol: width=160mm, font size=8pt. Default.
              Appropriate for 2-column figures
    - medium: width=121mm, font size=7pt.
    - small_font: width=121mm, font size=6pt.
    - height: height=2.5in.
    - small: width=80mm, font size=6pt
    - presentation: width=85mm, font size=10pt. For presentations.
    '''

    linestyle = '-'

    def set_onecol():
        '''
        Define parameters for "publish" mode and return value for pad_inches
        '''

        fontsize = 8
        lw = 1.
        markersize = 2
        fig_width = 3.15  # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.5,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'legend.fontsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.30

    def set_small():
        '''
        Define parameters for "publish" mode and return value for pad_inches
        '''

        fontsize = 6
        lw = 0.5
        markersize = 2
        fig_width = 3.15  # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.5,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'legend.fontsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'lines.markeredgewidth': 0.2,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.20

    def set_72mm():
        '''
        Define parameters for "72mm" mode and return value for pad_inches
        '''

        fontsize = 6
        markersize = 3
        lw = 0.7
        fig_width = 2.8  # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.35,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'legend.fontsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.20

    def set_50mm():
        '''
        Define parameters for "72mm" mode and return value for pad_inches
        '''

        fontsize = 5
        markersize = 2.5
        lw = 0.6
        fig_width = 2.  # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.3,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'legend.fontsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.10

    def set_medium():
        '''
        Define parameters for "medium" mode and return value for pad_inches
        '''

        fontsize = 8
        markersize = 3
        lw = 0.75
        fig_width = 3.15  # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.5,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'legend.fontsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.10

    def set_small_font():
        '''
        Define parameters for "small_font" mode and return value for pad_inches
        '''

        fontsize = 6
        markersize = 2
        lw = 0.6
        fig_width = 3.15  # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.5,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'legend.fontsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.10

    def set_large_font():
        '''
        Define parameters for "large_font" mode and return value for pad_inches
        '''

        fontsize = 10
        markersize = 9
        lw = 0.75
        fig_width = 6.2  # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.5,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'legend.fontsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.20

    def set_presentation():
        '''
        Define parameters for "presentation" mode and return value
        for pad_inches
        '''

        fontsize = 8
        lw = 1.5
        markersize = 3
        fig_width = 6.64  # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.75,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'legend.fontsize': fontsize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.2

    def set_twocol():
        '''
        Define parameters for "twocol" mode and return value for pad_inches
        '''

        fontsize = 8
        lw = 1.2
        markersize = 3
        fig_width = 6.3   # inch
        fig_height = aspect_ratio * fig_width  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.65,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'legend.fontsize': fontsize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.35

    def set_height():
        '''
        Define parameters for "twocol" mode and return value for pad_inches
        '''
        fontsize = 8
        lw = 1.1
        markersize = 1.5
        fig_height = 2.5   # inch
        fig_width = fig_height / aspect_ratio  # inch
        fig_size = [fig_width, fig_height]

        params = {'backend': 'ps',
                  'axes.linewidth': 0.65,
                  'lines.linewidth': lw,
                  'axes.labelsize': fontsize,
                  'font.size': fontsize,
                  'xtick.labelsize': fontsize,
                  'ytick.labelsize': fontsize,
                  'lines.linestyle': linestyle,
                  'lines.markersize': markersize,
                  'legend.fontsize': fontsize,
                  'font.size': fontsize,
                  'figure.figsize': fig_size}

        plt.rcParams.update(params)

        return lw, 0.025

    if (mode == "onecol"):
        return set_onecol()
    elif (mode == "small"):
        return set_small()
    elif (mode == "medium"):
        return set_medium()
    elif (mode == "72mm"):
        return set_72mm()
    elif (mode == "50mm"):
        return set_50mm()
    elif (mode == "small_font"):
        return set_small_font()
    elif (mode == "large_font"):
        return set_large_font()
    elif (mode == "presentation"):
        return set_presentation()
    elif (mode == "twocol"):
        return set_twocol()
    elif (mode == "height"):
        return set_height()
    else:
        print("%s mode not recognized, using onecol instead" % mode)
        return set_twocol()


def trend_estimator(x, y):
    '''
    Trend estimator

    Simultaneous estimation of bias, trend, annual, semi-annual and
    161-day sinusoid (alias period S2 tide errors).

    Parameters
    ----------
    x, y : array_like, x must have units "years"

    Returns
    -------
    x : ndarray
    The solution (or the result of the last iteration for an unsuccessful
    call).
    cov_x : ndarray
    Uses the fjac and ipvt optional outputs to construct an
    estimate of the jacobian around the solution.  ``None`` if a
    singular matrix encountered (indicates very flat curvature in
    some direction).  This matrix must be multiplied by the
    residual standard deviation to get the covariance of the
    parameter estimates -- see curve_fit.
    infodict : dict
    a dictionary of optional outputs with the key s::

        - 'nfev' : the number of function calls
        - 'fvec' : the function evaluated at the output
        - 'fjac' : A permutation of the R matrix of a QR
                 factorization of the final approximate
                 Jacobian matrix, stored column wise.
                 Together with ipvt, the covariance of the
                 estimate can be approximated.
        - 'ipvt' : an integer array of length N which defines
                 a permutation matrix, p, such that
                 fjac*p = q*r, where r is upper triangular
                 with diagonal elements of nonincreasing
                 magnitude. Column j of p is column ipvt(j)
                 of the identity matrix.
        - 'qtf'  : the vector (transpose(q) * fvec).

    mesg : str
    A string message giving information about the cause of failure.
    ier : int
    An integer flag.  If it is equal to 1, 2, 3 or 4, the solution was
    found.  Otherwise, the solution was not found. In either case, the
    optional output variable 'mesg' gives more information.

    Notes
    -----
    Code snipplet provided by Anthony Arendt, March 13, 2011.
    Uses scipy.optimize.leastsq, see documentation of
    scipy.optimize.leastsq for details.
    '''

    try:
        from scipy import optimize
    except:
        print("scipy.optimize not found. Please install.")
        exit(1)

    fitfunc = lambda p, x: (p[0] + p[1] * x +
                            p[2] * np.cos(2.0 * np.pi * (x - p[3]) / 1.0) +
                            p[4] * np.cos(2.0 * np.pi * (x - p[5]) / 0.5) +
                            p[6] * np.cos(2.0 * np.pi * (x - p[7]) / 0.440794))
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    p0 = [0.0, -80.0, 40.0, 0.0, 10.0, 0.0, 1.0, 0.0]

    return optimize.leastsq(errfunc, p0[:], args=(x, y), full_output=1)


def colorList():
    '''
    Returns a list with colors, e.g for line plots. etc.
    '''
    colors = ['#084594',  # dark blue
              '#FF7F00',  # orange
              '#984EA3',  # violet
              '#E41A1C',  # red
              '#4DAF4A',  # green
              '#377EB8',  # light blue
              '#FB9A99',  # light red
              '#FB9A99',  # light orange
              '#CAB2D6',  # light violet
              'brown',
              'pink']
    return colors


def gmtColormap(fileName, log_color=False, reverse=False):
    '''
    Import a CPT colormap from GMT.

    Parameters
    ----------
    fileName : a cpt file.

    Example
    -------
    >>> cdict = gmtColormap("mycolormap.cpt")
    >>> gmt_colormap = colors.LinearSegmentedColormap("my_colormap", cdict)

    Notes
    -----
    This code snipplet modified after
    http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg09547.html
    '''
    import colorsys
    import os

    try:
        try:
            f = open(fileName)
        except:
            # Check if it's a colormap provided in colormaps/
            basedir, fname = os.path.split(__file__)
            my_file = os.path.join(basedir, 'colormaps', fileName)
            f = open(my_file)
    except:
        print "file ", fileName, "not found"
        return None

    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split()
        if l[0] == "#":
            if ls[-1] == "HSV":
                colorModel = "HSV"
                continue
            else:
                continue
        if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
            pass
        else:
            x.append(float(ls[0]))
            r.append(float(ls[1]))
            g.append(float(ls[2]))
            b.append(float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    if reverse:
        r.reverse()
        g.reverse()
        b.reverse()

    x = np.array(x, np.float32)
    r = np.array(r, np.float32)
    g = np.array(g, np.float32)
    b = np.array(b, np.float32)
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
            r[i] = rr
            g[i] = gg
            b[i] = bb
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
            r[i] = rr
            g[i] = gg
            b[i] = bb
    if colorModel == "RGB":
        r = r / 255.
        g = g / 255.
        b = b / 255.

    if log_color:
        xNorm = np.zeros((len(x), ))
        xNorm[1::] = np.logspace(-1, 0, len(x) - 1)
        xNorm[1::-2] /= 4
    else:
        xNorm = (x - x[0]) / (x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i], r[i], r[i]])
        green.append([xNorm[i], g[i], g[i]])
        blue.append([xNorm[i], b[i], b[i]])
    colorDict = {'red': red, 'green': green, 'blue': blue}
    return (colorDict)


def smooth(x, window_len=11, window='hanning'):
    '''
    Smooth the data using a window with requested size (running mean,
    moving average, low pass filtering).

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Parameters
    ----------
    x : array_like, the input signal
    window_len : the dimension of the smoothing window; should be an odd integer
    window : the type of window from "flat", "hanning", "hamming",
    "bartlett", "blackman" flat window will produce a moving average smoothing.

    Returns
    -------
    y : the smoothed signal

    Example
    -------
    t = np.linspace(-2,2,0.1)
    x = np.sin(t) + np.randn(len(t))*0.1
    y = smooth(x)

    See also
    --------
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve
    scipy.signal.lfilter

    Notes
    -----
    Downloaded from http://www.scipy.org/Cookbook/SignalSmooth.

    TODO
    ----
    the window parameter could be the window itself if an array instead
    of a string
    '''

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError(
            "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[2 * x[0] - x[window_len:1:-1],
              x,
              2 * x[-1] - x[-1:-window_len:-1]]

    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='same')
    return y[window_len - 1:-window_len + 1]


def fftsmooth(x, window_len=11, window='hanning'):
    '''
    Smooth the data using a window with requested size (running mean,
    moving average, low pass filtering).

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Parameters
    ----------
    x : array_like, the input signal
    window_len : the dimension of the smoothing window; should be an odd integer
    window : the type of window from "flat", "hanning", "hamming",
    "bartlett", "blackman" flat window will produce a moving average smoothing.

    Returns
    -------
    y : the smoothed signal

    Example
    -------
    t = np.linspace(-2,2,0.1)
    x = np.sin(t) + np.randn(len(t))*0.1
    y = smooth(x)

    See also
    --------
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve
    scipy.signal.lfilter

    Notes
    -----
    Downloaded from http://www.scipy.org/Cookbook/SignalSmooth, but replaced
    np.convovle with faster scipy.signal.fftconvolve

    TODO
    ----
    the window parameter could be the window itself if an array instead
    of a string
    '''

    from scipy.signal import fftconvolve

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError(
            "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[2 * x[0] - x[window_len:1:-1],
              x,
              2 * x[-1] - x[-1:-window_len:-1]]

    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = fftconvolve(w / w.sum(), s, mode='same')
    return y[window_len - 1:-window_len + 1]


def get_rmse(a, b, N, w=None):
    '''
    Returns the (weighted) root mean square error of differences between a and b.

    Parameters
    ----------
    a, b : array_like
    N : number of valid values
    w : weights

    Returns
    -------
    rmse : scalar
    '''

    if w is None:
        w = np.ones_like(a)
    c = (a.ravel() - b.ravel()) / w.ravel()
    if isinstance(c, np.ma.MaskedArray):
        return np.sqrt(np.linalg.norm(np.ma.compressed(c), 2) ** 2.0 / N)
    else:
        return np.sqrt(np.linalg.norm(c, 2) ** 2.0 / N)


def get_avg(a, b, N, relative=False):
    '''
    Returns the average difference between a and b.

    Parameters
    ----------
    a,b : array_like
    N : number of values

    Returns
    -------
    avg : scalar

    Notes
    -----
    The average is the sum of elements of the difference (a - b)
    divided by the number of elements N.
    '''
    if relative is False:
        c = a.ravel() - b.ravel()
    else:
        c = (a.ravel() - b.ravel() / b.ravel())
    if isinstance(c, np.ma.MaskedArray):
        return np.linalg.norm(np.ma.compressed(c), 1) / N
    else:
        return np.linalg.norm(c, 1) / N


def unit_converter(data, inunit, outunit):
    '''
    Unit converter. Takes an (numpy) array, valid udunits inunits and outunits
    as strings, and returns the array in outunits.

    Parameters
    ----------
    data : array_like
    inunit : string
             unit to convert from, must be UDUNITS-compatible string
    outunit : string
              unit to conver to, must be UDUNITS-compatible string

    Returns
    -------
    out : array_like

    Example
    -------
    >>> import numpy as np
    >>> c = Converter("kg","Gt")
    >>> out = c(np.array([1,2])*1e12)
    >>> out = array([ 1.,  2.])
    '''

    inunit = str(inunit)
    outunit = str(outunit)
    data = np.array(data)
    if not (inunit == outunit):
        try:
            try:
                from cf_units import Unit
                in_unit  = Unit(inunit)
                out_unit  = Unit(outunit)
                outdata = in_unit.convert(data, out_unit)
            except:
                from udunits2 import Converter, System, Unit
                sys = System()
                c = Converter((Unit(sys, inunit), Unit(sys, outunit)))
                outdata = c(data)
        except:
            print(
                "No udunits module found, you're on your own.\n  -> I am assuming input unit is m, and will convert to km.\n  -> Installation of Constantine's awesome python wrapper for udunits is highly recommended.\n  -> Download it from https://github.com/ckhroulev/py_udunits2.")
            c = 1. / 1e3
            outdata = c * data
    else:
        outdata = data

    return outdata


def permute(variable, output_order=('time', 'z', 'zb', 'y', 'x')):
    '''
    Permute dimensions of a NetCDF variable to match the output
    storage order.

    Parameters
    ----------
    variable : a netcdf variable
               e.g. thk = nc.variables['thk']
    output_order: dimension tuple (optional)
                  default ordering is ('time', 'z', 'zb', 'y', 'x')

    Returns
    -------
    var_perm : array_like
    '''

    input_dimensions = variable.dimensions

    # filter out irrelevant dimensions
    dimensions = filter(lambda x: x in input_dimensions,
                        output_order)

    # create the mapping
    mapping = map(lambda x: dimensions.index(x),
                  input_dimensions)

    if mapping:
        return np.transpose(variable[:], mapping)
    else:
        return variable[:]  # so that it does not break processing "mapping"


def plot_mapview(var, **kwargs):
    '''Plot map view of variable var'''
    from matplotlib import colors
    kwargsdict = {}
    expected_args = ["log", "show", "title"]
    for key in kwargs.keys():
        if key in expected_args:
            kwargsdict[key] = kwargs[key]
        else:
            raise Exception("Unexpected Argument")
    log = False
    if 'log' in kwargsdict:
        log = kwargsdict['log']
    if log:
        norm = colors.LogNorm(vmin=0.3, vmax=HISTMAX)
    else:
        norm = None
    show = False
    if 'show' in kwargsdict:
        show = kwargsdict['show']
    title = None
    if 'title' in kwargsdict:
        title = kwargsdict['title']
    plt.figure()
    plt.imshow(var, origin='lower', norm=norm)
    plt.colorbar()
    plt.title(title)


def plot_histogram(data, **kwargs):
    '''Plot a simple histogram with one variable'''
    import pylab as plt

    kwargsdict = {}
    expected_args = ["log", "show"]
    for key in kwargs.keys():
        if key in expected_args:
            kwargsdict[key] = kwargs[key]
        else:
            raise Exception("Unexpected Argument")
    log = False
    if 'log' in kwargsdict:
        log = kwargsdict['log']
    show = False
    if 'show' in kwargsdict:
        show = kwargsdict['show']
    ticks = map(lambda x: "%3d" % x, bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bin, patches = ax.hist(data.values[data.values > vmin],
                              bins=bins,
                              log=log,
                              histtype='bar')
    ax.set_xlabel("ice surface velocity, m a$^{-1}$")
    ax.set_ylabel("number of grid cells")
    plt.title(data.title)
    return n, bin, patches


def plot_histogram2(data, obs, **kwargs):
    '''Plot a histogram with modeled and observed values'''
    import pylab as plt

    kwargsdict = {}
    expected_args = ["log", "show"]
    for key in kwargs.keys():
        if key in expected_args:
            kwargsdict[key] = kwargs[key]
        else:
            raise Exception("Unexpected Argument")
    log = False
    if 'log' in kwargsdict:
        log = kwargsdict['log']
    show = False
    if 'show' in kwargsdict:
        show = kwargsdict['show']
    ticks = map(lambda x: "%3d" % x, bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bin, patches = ax.hist([data.values[data.values > vmin],
                               obs.values[obs.values > vmin]],
                              bins=bins,
                              log=log,
                              histtype='bar',
                              label=['modeled', 'observed'])
    ax.set_xlabel("ice surface velocity, m a$^{-1}$")
    ax.set_ylabel("number of grid cells")
    ax.set_xlim(vmin, vmax)
    plt.title(data.title)
    ax.legend()


def print_info(experiment):
    '''Prints parameters and values of an experiment'''
    for key, val in experiment.parameter_dict.iteritems():
        if isinstance(val, basestring):
            print("      * %s = %s" % (key, val))
        else:
            print("      * %s = %3.3f" % (key, val))


def print_overall_statistics(experiments):
    '''Prints some statistics'''
    print("\nOverall statistics:\n")

    rmses = [n.rmse for n in experiments]
    rmse_min = filter(lambda x: x.rmse == min(rmses), experiments)
    rmse_max = filter(lambda x: x.rmse == max(rmses), experiments)

    print("    - smallest rmse = %3.2f" % rmse_min[0].rmse)
    print_info(rmse_min[0])
    print("    - largest rmse = %3.2f" % rmse_max[0].rmse)
    print_info(rmse_max[0])
    print("\n")


class DataObject(object):

    '''
    A base class for experiments and observations.

    Parameters
    ----------
    file_name : filename of the netcdf containing observations
    variable_name: name of the variable (e.g. csurf)
    bins : numpy array with bins, e.g.
           bins = np.linspace(0,10,11)
    '''

    def __init__(self, file_name, variable_name, bins, *args, **kwargs):
        super(DataObject, self).__init__(*args, **kwargs)
        self.file_name = file_name
        self.variable_name = variable_name
        self.bins = bins
        self.bin_bounds = [bins[0], bins[-1]]

    def __del__(self):
        # Close open file
        self.nc.close()

    def _make_histogram(self):
        '''Calculate number of grid cells n in bins bin.'''
        # np.histogram ignores masked array, so we have to exclude
        # masked values by applying the inverted mask
        n, Bins = np.histogram(self.values[~self.mask], bins=self.bins)
        self.n = n
        self.no_above_histmax = np.sum(self.values >= self.bins[-1])
        self.cells_in_histogram = np.sum(self.n)
        self.coverage = (np.float(self.cells_in_histogram) /
                         self.valid_cells * 100)
        if DEBUG:
            print("    - no of cells in histogram = %i (%3.2f %%)"
                  % (self.cells_in_histogram, self.coverage))
            print("    - no of cells above %3.2f %s = %i"
                  % (self.bin_bounds[1], self.units, self.no_above_histmax))
            print("    - no of cells in bins =\n   %s" % str(self.n))

    def add_mask(self, mask):
        '''
        Add/update a mask.

        Adds (or updates) the mask. self.mask = logical_or(self.mask,mask).

        Parameters
        ----------
        mask : array_like
               array with zeros (no mask) and ones (mask)
        '''
        self.mask = np.logical_or(self.mask, mask)
        self.values = np.ma.array(data=self.values, mask=self.mask)
        self._set_valid_cells()
        if DEBUG:
            plot_mapview(self.mask, show=True, title='updated mask')
        self._make_histogram()

    def get_total_cells(self):
        '''Return total number of grid cells.'''
        nx, ny = np.shape(self.values)
        return nx * ny

    def _set_valid_cells(self):
        '''Determine number of valid cells by summing up non-masked cells.'''
        value = np.sum(-self.mask)
        self.valid_cells = value

    def _set_grid_spacing(self):
        '''Tries to determine grid spacing in m, sets None if it fails.'''
        try:
            x_dim = self.nc.variables["x"]
            in_units = self.nc.variables["x"].units
        except:
            x_dim = self.nc.variables["x1"]
            in_units = self.nc.variables["x1"].units
        out_units = "m"
        grid_spacing = unit_converter(np.abs(x_dim[1] - x_dim[0]),
                                      in_units, out_units)
        self.grid_spacing = grid_spacing
        self.grid_spacing_units = out_units
        if DEBUG:
            print("grid spacing is %3.1f %s"
                  % (self.grid_spacing, self.grid_spacing_units))

    def get_grid_spacing(self):
        '''Return grid spacing.'''
        if self.grid_spacing is not None:
            print("grid spacing is %3.1f %s"
                  % (self.grid_spacing, self.grid_spacing_units))
            return self.grid_spacing
        else:
            return self.grid_spacing

    def get_valid_cells(self):
        '''Return number of valid cells.'''
        return self.valid_cells

    def set_rmse(self, rmse):
        '''Set root mean square error (RMSE).'''
        self.rmse = rmse

    def set_wrmse(self, wrmse):
        '''Set weighted root mean square error (WRMSE).'''
        self.wrmse = wrmse

    def set_avg(self, avg):
        '''Set average (avg).'''
        self.avg = avg

    def set_rmse_relative(self, rmse):
        '''Set root mean square error (RMSE).'''
        self.rmse_relative = rmse

    def set_avg_relative(self, avg):
        '''Set average (avg).'''
        self.avg_relative = avg


class Experiment(DataObject):

    '''
    A derived class for experiments

    A derived class for handling PISM experiments.

    Parameters
    ----------

    obs : Observation
          An instance of Observation.
    outlier: A scalar, values where (obs - experiment) > outlier will
          be ignored.
    parameter_list : list of parameters being used for evaluation
    abbr_dict : a dictionary with abbrevations for parameters
         e.g. dict(zip(['enhancement_factor','pseudo_plastic_q'],['e','q']))
    '''

    def __init__(self, obs, outlier, parameter_list, abbr_dict,
                 *args, **kwargs):
        super(Experiment, self).__init__(*args, **kwargs)

        self.valid_cells = None
        self.outliers = None
        self.no_outliers = None
        self.coverage = None
        self.avg = None
        self.rmse = None
        self.abbr_dict = abbr_dict

        def _print_info():
            '''Print infos about experiment'''
            print("  - Parameters and values for this experiment:")
            for key, val in self.parameter_dict.iteritems():
                if isinstance(val, basestring):
                    print("    * %s = %s" % (key, val))
                else:
                    print("    * %s = %3.3f" % (key, val))

        # Open netcdf file
        try:
            self.nc = CDF(self.file_name, 'r')
        except IOError:
            print("File %s does not exist or is not a valid NetCDF file"
                  % self.file_name)
        print("\nInitializing experiement from file %s" % self.file_name)
        print("  - reading variable %s" % self.variable_name)
        # Determine units from 'units' attribute, set to None if not found
        try:
            self.units = self.nc.variables[self.variable_name].units
        except:
            self.units = None
        # Read values, squeeze and permute (if appropriate)
        try:
            self.values = np.squeeze(permute(
                self.nc.variables[self.variable_name]))
        except:
            print("Could not read variable %s" % self.variable_name)
        # Determine fill value from '_FillValue' attribute, set to
        # None if not found
        try:
            self._FillValue = self.nc.variables[self.variable_name]._FillValue
        except:
            self._FillValue = None
        try:
            self.valid_min = self.nc.variables[self.variable_name].valid_min
        except:
            self.valid_min = None
        # Set grid spacing
        self._set_grid_spacing()

        nx, ny = np.shape(self.values)
        self.total_cells = nx * ny
        print("  - total cells = %i" % self.total_cells)
        # get mask from masked array "variable" if exist
        try:
            self.mask = self.values.mask
            valid_cells = np.sum(-self.mask)
        except:
            valid_cells = self.total_cells
            self.mask = np.zeros_like(self.values, dtype='bool')
        self.valid_cells = valid_cells
        print("  - valid cells = %i" % valid_cells)

        self.parameter_list = parameter_list
        pism_config = self.nc.variables["pism_config"]
        run_stats = self.nc.variables['run_stats']
        self.config = dict()
        for attr in pism_config.ncattrs():
            self.config[attr] = getattr(pism_config, attr)
        for attr in run_stats.ncattrs():
            self.config[attr] = getattr(run_stats, attr)
        self.parameter_dict = dict(
            [(param, self.config[param])
             for param in parameter_list])
        self.parameter_short_dict = dict(
            [(abbr_dict[param], self.config[param])
             for param in parameter_list])
        _print_info()

        self.title = ', '.join(["%s = %s"
                                % (key, str(val))
                                for key, val in self.parameter_short_dict.iteritems()])
        self._make_histogram()
        if DEBUG:
            plot_mapview(self.values, log=True, title='before no obs')
        print("\n  * Applying mask where no observations, updating histogram")
        self.add_mask(obs.mask)
        if DEBUG:
            plot_mapview(self.values, log=True, title='after obs')
        print("\n  * Applying mask where abs(%s - %s) > %3.2f %s, updating histogram"
              % (self.variable_name, obs.variable_name, outlier, self.units))
        outliers = (np.abs(self.values - obs.values) > outlier)
        self.add_mask(outliers)
        if DEBUG:
            plot_mapview(self.values, log=True, title='after outliers')

        print("    - no of cells in histogram = %i (%3.2f %%)"
              % (self.cells_in_histogram, self.coverage))
        print("    - no of cells above %3.2f %s = %i"
              % (self.bin_bounds[1], self.units, self.no_above_histmax))
        print("    - no of cells in bins =\n   %s" % str(self.n))

    def get_values(self):
        '''Return values'''
        return self.values

    def set_coverage(self, value):
        '''
        Set coverage defined as number of cells used for comparison divided
        by number of ice covered cells
        '''
        self.coverage = value

    def set_outliers(self, outliers):
        self.outliers = outliers
        self.no_outliers = np.sum(outliers)


class Observation(DataObject):

    """
    A derived class for observations.

    A derived class to handle observational data. The base class takes
    care of initializing.
    """

    def __init__(self, *args, **kwargs):
        super(Observation, self).__init__(*args, **kwargs)
        ## self.title = "observed " + self.variable_name
        self.title = "observed"
        # Open netcdf file
        try:
            self.nc = CDF(self.file_name, 'r')
        except IOError:
            print("File %s does not exist or is not a valid NetCDF file"
                  % self.file_name)
        print("\nInitializing '%s' observations from file %s"
              % (self.variable_name, self.file_name))
        print("  - reading variable %s" % self.variable_name)
        # Determine untis from 'units' attribute, set to None if not found
        try:
            self.units = self.nc.variables[self.variable_name].units
        except:
            self.units = None
        # Read values, squeeze and permute (if appropriate)
        try:
            self.values = np.squeeze(permute(
                self.nc.variables[self.variable_name]))
        except:
            try:
                self.values = np.squeeze(permute(
                    self.nc.variables[self.variable_name], ("time", "y1", "x1")))
            except:
                print("    could not read variable %s" % self.variable_name)
        # Determine fill value from '_FillValue' attribute, set to
        # None if not found
        try:
            self._FillValue = self.nc.variables[self.variable_name]._FillValue
        except:
            self._FillValue = None
        try:
            self.valid_min = self.nc.variables[self.variable_name].valid_min
        except:
            self.valid_min = None

        nx, ny = np.shape(self.values)
        self.total_cells = nx * ny
        print("  - total cells = %i" % self.total_cells)
        # get mask from masked array "variable" if exist
        try:
            self.mask = self.values.mask
            valid_cells = np.sum(-self.mask)
        except:
            valid_cells = self.total_cells
            self.mask = np.zeros_like(self.values, dtype='bool')
        self.valid_cells = valid_cells
        print("  - valid cells = %i" % valid_cells)
        if DEBUG:
            plot_mapview(self.mask, show=True, title='obs mask')

        self._make_histogram()
        print("  - no of cells in histogram = %i (%3.2f %%)"
              % (self.cells_in_histogram, self.coverage))
        print("  - no of cells above %3.2f %s = %i (%3.2f %%)"
              % (self.bins[1], self.units, self.valid_cells - self.n[0],
                 float(self.valid_cells - self.n[0]) / self.valid_cells * 100))
        print("  - no of cells above %3.2f %s = %i (%3.2f %%)"
              % (self.bin_bounds[1], self.units, self.no_above_histmax,
                 100 - self.coverage))
        print("  - no of cells in bins =\n   %s" % str(self.n))
