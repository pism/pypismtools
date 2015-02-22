#!/usr/bin/env python
# Copyright (C) 2012-2013, 2015 Andy Aschwanden
#

# nosetests --with-coverage --cover-branches --cover-html --cover-package=extract_profiles scripts/extract_profiles.py

from argparse import ArgumentParser
import numpy as np
import scipy.sparse

from netCDF4 import Dataset as NC

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt


## {{{ http://code.activestate.com/recipes/496938/ (r1)
"""
A module that helps to inject time profiling code
in other modules to measures actual execution times
of blocks of code.

"""

__author__ = "Anand B. Pillai"
__version__ = "0.1"

import time

def timeprofile():
    """ A factory function to return an instance of TimeProfiler """

    return TimeProfiler()

class TimeProfiler:
    """ A utility class for profiling execution time for code """

    def __init__(self):
        # Dictionary with times in seconds
        self.timedict = {}

    def mark(self, slot=''):
        """ Mark the current time into the slot 'slot' """

        # Note: 'slot' has to be string type
        # we are not checking it here.

        self.timedict[slot] = time.time()

    def unmark(self, slot=''):
        """ Unmark the slot 'slot' """

        # Note: 'slot' has to be string type
        # we are not checking it here.

        if self.timedict.has_key(slot):
            del self.timedict[slot]

    def lastdiff(self):
        """ Get time difference between now and the latest marked slot """

        # To get the latest slot, just get the max of values
        return time.time() - max(self.timedict.values())

    def elapsed(self, slot=''):
        """ Get the time difference between now and a previous
        time slot named 'slot' """

        # Note: 'slot' has to be marked previously
        return time.time() - self.timedict.get(slot)

    def diff(self, slot1, slot2):
        """ Get the time difference between two marked time
        slots 'slot1' and 'slot2' """

        return self.timedict.get(slot2) - self.timedict.get(slot1)

    def maxdiff(self):
        """ Return maximum time difference marked """

        # Difference of max time with min time
        times = self.timedict.values()
        return max(times) - min(times)

    def timegap(self):
        """ Return the full time-gap since we started marking """

        # Return now minus min
        times = self.timedict.values()
        return time.time() - min(times)

    def cleanup(self):
        """ Cleanup the dictionary of all marks """

        self.timedict.clear()

def normal(p0, p1):
    '''
    Compute the unit normal vector orthogonal to (p1-p0), pointing 'to the
    right' of (p1-p0).
    '''

    a = p0 - p1
    if a[1] != 0.0:
        n = np.array([1.0, - a[0] / a[1]])
        n = n / np.linalg.norm(n) # normalize
    else:
        n = np.array([0,1])

    # flip direction if needed:
    if np.cross(a, n) < 0:
        n = -1.0 * n

    return n

class Profile:
    def __init__(self, name, lat, lon, center_lat, center_lon, flightline, glaciertype, flowtype, projection, flip=False):
        self.name = name
        self.center_lat = center_lat
        self.center_lon = center_lon
        self.flightline = flightline
        self.glaciertype = glaciertype
        self.flowtype = flowtype
        if flip:
            self.lat = lat[::-1]
            self.lon = lon[::-1]
        else:
            self.lat = lat
            self.lon = lon
        self.x, self.y = projection(lon, lat)

        self.distance_from_start = self._distance_from_start()
        self.nx, self.ny = self._compute_normals()

    def _compute_normals(self):
        '''
        Compute normals to a flux gate described by 'p'. Normals point 'to
        the right' of the path.
        '''

        p = np.vstack((self.x, self.y)).T

        ns = np.zeros_like(p)
        ns[0] = normal(p[0], p[1])
        for j in range(1, len(p) - 1):
            ns[j] = normal(p[j-1], p[j+1])

        ns[-1] = normal(p[-2], p[-1])

        return ns[:, 0], ns[:, 1]

    def _distance_from_start(self):
        result = np.zeros_like(self.x)
        result[1::] = np.sqrt(np.diff(self.x)**2 + np.diff(self.y)**2)
        return result.cumsum()

class ProfileInterpolationMatrix:
    # sparse matrix
    A = None
    # row and column ranges for extracting array subsets
    r_min = None
    r_max = None
    c_min = None
    c_max = None
    n_rows = None
    n_cols = None

    def column(self, r, c):
        """Interpolation matrix column number corresponding to r,c of the
        array *subset*. This is the same as the linear index within
        the subset needed for interpolation.

        """
        return self.n_cols * r + c

    def grid_column(self, x, dx, X):
        "Input grid column number corresponding to X."
        return int(np.floor((X - x[0]) / dx))

    def grid_row(self, y, dy, Y):
        "Input grid row number corresponding to Y."
        return int(np.floor((Y - y[0]) / dy))

    def __init__(self, x, y, px, py, bilinear=True):
        """Interpolate values of z to points (px,py) assuming that z is on a
        regular grid defined by x and y."""

        assert len(px) == len(py)

        dx = x[1] - x[0]
        dy = y[1] - y[0]

        assert dx > 0
        assert dy > 0

        self.c_min = self.grid_column(x, dx, np.min(px))
        self.c_max = self.grid_column(x, dx, np.max(px)) + 1

        self.r_min = self.grid_row(y, dy, np.min(py))
        self.r_max = self.grid_row(y, dy, np.max(py)) + 1

        # compute the size of the subset needed for interpolation
        self.n_rows = self.r_max - self.r_min + 1
        self.n_cols = self.c_max - self.c_min + 1

        n_points = len(px)
        self.A = scipy.sparse.lil_matrix((n_points, self.n_rows * self.n_cols))

        if bilinear:
            self._compute_bilinear_matrix(x, y, dx, dy, px, py)
        else:
            raise NotImplementedError

    def _compute_bilinear_matrix(self, x, y, dx, dy, px, py):
        for k in xrange(self.A.shape[0]):
            x_k = px[k]
            y_k = py[k]

            C = self.grid_column(x, dx, x_k)
            R = self.grid_row(y, dy, y_k)

            alpha = (x_k - x[C]) / dx
            beta  = (y_k - y[R]) / dy

            # indexes within the subset needed for interpolation
            c = C - self.c_min
            r = R - self.r_min

            self.A[k, self.column(r,         c)] = (1.0 - alpha) * (1.0 - beta)
            self.A[k, self.column(r + 1,     c)] = (1.0 - alpha) * beta
            self.A[k, self.column(r,     c + 1)] = alpha * (1.0 - beta)
            self.A[k, self.column(r + 1, c + 1)] = alpha * beta

    def adjusted_matrix(self, mask):
        """Return adjusted interpolation matrix that ignores missing (masked)
        values."""

        A = self.A.tocsr()
        n_points = A.shape[0]

        output_mask = np.zeros(n_points, dtype=np.bool_)

        for r in xrange(n_points):
            # for each row, i.e. each point along the profile
            row = np.s_[A.indptr[r]:A.indptr[r+1]]
            # get the locations and values
            indexes = A.indices[row]
            values = A.data[row]

            # if a particular location is masked, set the
            # interpolation weight to zero
            for k, index in enumerate(indexes):
                if np.ravel(mask)[index]:
                    values[k] = 0.0

            # normalize so that we still have an interpolation matrix
            if values.sum() > 0:
                values = values / values.sum()
            else:
                output_mask[r] = True

            A.data[row] = values

        A.eliminate_zeros()

        return A, output_mask

    def apply(self, array):
        """Apply the interpolation to an array. Returns values at points along
        the profile."""
        subset = array[self.r_min:self.r_max+1, self.c_min:self.c_max+1]
        return self.apply_to_subset(subset)

    def apply_to_subset(self, subset):
        """Apply interpolation to an array subset."""

        if np.ma.is_masked(subset):
            A, mask = self.adjusted_matrix(subset.mask)
            data = A * np.ravel(subset)
            return np.ma.array(data, mask=mask)

        return self.A.tocsr() * np.ravel(subset)

def masked_interpolation_test():
    """Test matrix adjustment."""

    # 2x2 grid of ones
    x = [0, 1]
    y = [0, 1]
    z = np.ones((2,2))
    # set the [0,0] element to a nan and mark that value
    # as "missing" by turning it into a masked array
    z[0,0] = np.nan
    z = np.ma.array(z, mask=[[True, False],
                             [False, False]])
    # sample in the middle
    px = [0.5]
    py = [0.5]

    A = ProfileInterpolationMatrix(x, y, px, py)

    # We should get the average of the three remaining ones, i.e. 1.0.
    # (We would get a nan without adjusting the matrix.)
    assert A.apply(z)[0] == 1.0

def masked_missing_interpolation_test():
    """Test interpolation from a masked array that produces missing values
    in the output."""

    x = [-1, 0, 1]
    y = [-1, 0, 1]
    z = np.ones((3,3))

    # set the four elements in the corner to nan and mark them as
    # missing
    z[0:2,0:2] = np.nan
    mask = np.zeros_like(z, dtype=np.bool_)
    mask[0:2,0:2] = np.nan

    z = np.ma.array(z, mask=mask)

    px = [-0.5, 0.5]
    py = [-0.5, 0.5]

    A = ProfileInterpolationMatrix(x, y, px, py)

    z_interpolated = A.apply(z)

    assert z_interpolated.mask[0] == True
    assert z_interpolated[1] == 1.0

def interpolation_test():
    """Test interpolation by recovering values of a linear function."""

    Lx = 10.0                    # size of the box in the x direction
    Ly = 20.0                    # size of the box in the y direction
    P = 100                      # number of test points

    # grid size (note: it should not be a square)
    Mx = 101
    My = 201
    x = np.linspace(0, Lx, Mx)
    y = np.linspace(0, Ly, My)

    # test points
    np.random.seed([100])
    px = np.random.rand(P) * Lx
    py = np.random.rand(P) * Ly

    try:
        A = ProfileInterpolationMatrix(x, y, px, py, bilinear=False)
        raise RuntimeError("Update this test if you implemented nearest neighbor interpolation.") #pragma: nocover
    except NotImplementedError:
        pass

    # initialize the interpolation matrix
    A = ProfileInterpolationMatrix(x, y, px, py)

    # a linear function (perfectly recovered using bilinear
    # interpolation)
    def Z(x, y):
        return 0.3 * x + 0.2 * y + 0.1

    # compute values of Z on the grid
    xx, yy = np.meshgrid(x, y)
    z = Z(xx, yy)

    # interpolate
    z_interpolated = A.apply(z)

    assert np.max(np.fabs(z_interpolated - Z(px, py))) < 1e-12

def create_test_file():
    """Create an input file for testing."""

    import netCDF4
    nc = netCDF4.Dataset("test_input.nc", 'w')

    Mx = 88
    My = 152
    Mz = 11
    for name, length in [["x", Mx], ["y", My], ["z", Mz], ["time", None]]:
        nc.createDimension(name, length)
        nc.createVariable(name, "f4", (name,))

    x = np.linspace(-669650.0, 896350.0, Mx)
    y = np.linspace(-3362600.0, -644600.0, My)
    z = np.linspace(0, 4000.0, Mz)

    for name, data in [["x", x], ["y", y], ["z", z]]:
        nc.variables[name][:] = data

    nc.variables['time'][0] = 0.0

    nc.proj4 = "+init=epsg:3413"

    xx,yy = np.meshgrid(x,y)

    def F(x,y,z):
        return 0.1 * x + 0.2 * y + 0.3 + 0.4 * z

    def write(dimensions):
        name = "test_" + "_".join(dimensions)
        variable = nc.createVariable(name, "f4", dimensions)
        indexes = [Ellipsis] * len(dimensions)

        # write to the first time record only
        if "time" in dimensions:
            indexes[dimensions.index("time")] = 0

        # transpose input 2D array if needed
        if dimensions.index("y") < dimensions.index("x"):
            T = lambda x: x
        else:
            T = np.transpose

        # fill all z levels with the same
        if "z" in dimensions:
            for k in xrange(Mz):
                indexes[dimensions.index("z")] = k
                variable[indexes] = T(F(xx, yy, z[k]))
        else:
            variable[indexes] = T(F(xx, yy, 0))

    import itertools

    # 2D time-independent
    for d in itertools.permutations(["x", "y"]):
        write(d)

    # 2D time-dependent
    for d in itertools.permutations(["time", "x", "y"]):
        write(d)

    # 3D time-independent
    for d in itertools.permutations(["x", "y", "z"]):
        write(d)

    # 3D time-dependent
    for d in itertools.permutations(["time", "x", "y", "z"]):
        write(d)

    nc.close()

def profile_test():
    """Test Profile constructor."""
    import pyproj
    x = np.linspace(-1.0, 1.0, 101)
    y = np.linspace(1.0, -1.0, 101)

    projection = pyproj.Proj("+proj=latlon")

    lon,lat = projection(x, y, inverse=True)

    center_lat,center_lon = projection(0.0, 0.0, inverse=True)

    profile = Profile("test_profile", lat, lon, center_lat, center_lon, projection)

    assert profile.nx[0] == -1.0 / np.sqrt(2.0)
    assert profile.ny[0] == -1.0 / np.sqrt(2.0)

    assert np.fabs(profile.distance_from_start[1] - 0.02 * np.sqrt(2.0)) < 1e-12

    x = -1.0 * x
    lon,lat = projection(x, y, inverse=True)

    profile = Profile("flipped_profile", lat, lon, center_lat, center_lon, projection,
                      flip=True)

    assert profile.nx[0] == -1.0 / np.sqrt(2.0)
    assert profile.ny[0] == 1.0 / np.sqrt(2.0)

    x = np.linspace(-1.0, 1.0, 101)
    y = np.zeros_like(x)
    lon,lat = projection(x, y, inverse=True)

    profile = Profile("test_profile", lat, lon, center_lat, center_lon, projection)

    assert profile.nx[0] == 0.0
    assert profile.ny[0] == -1.0

def load_profiles(filename, projection, flip):
    """Load profiles from a file filename.

    Parameters
    -----------
    filename: filename of ESRI shape file

    projection: proj4 projection object. (lon,lat) coordinates of
                points along a profile are conterted to (x,y)
                coordinates in this projection. This should be the
                projection used by the dataset we're extracting
                profiles from.
    flip: boolean; set to True to flip profile directions

    Returns
    -------
    list of proviles with
    """
    profiles = []
    for lat, lon, name, clat, clon, flightline, glaciertype, flowtype in read_shapefile(filename):
        profiles.append(Profile(name, lat, lon, clat, clon, flightline, glaciertype, flowtype, projection, flip))
    return profiles

def output_dimensions(input_dimensions, stationdim, profiledim):
    """Build a list of dimension names used to define a variable in the
    output file."""
    _, _, zdim, tdim = get_dims_from_variable(input_dimensions)

    if tdim:
        result = [stationdim, tdim, profiledim]
    else:
        result = [stationdim, profiledim]

    if zdim:
        result.append(zdim)

    return result

def read_shapefile(filename):
    '''
    Reads lat / lon from a ESRI shape file.

    Paramters
    ----------
    filename: filename of ESRI shape file.

    Returns
    -------
    lat, lon: array_like coordinates

    '''
    import ogr
    import osr
    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.Open(filename, 0)
    layer = data_source.GetLayer(0)
    srs=layer.GetSpatialRef()
    if not srs.IsGeographic():
        print('''Spatial Reference System in % s is not latlon. Converting.'''
              % filename)
        # Create spatialReference, EPSG 4326 (lonlat)
        srs_geo = osr.SpatialReference()
        srs_geo.ImportFromEPSG(4326)
    cnt = layer.GetFeatureCount()
    names = []
    profiles = []
    for pt in range(0, cnt):
        feature = layer.GetFeature(pt)
        try:
            name = feature.name
        except:
            name = str(pt)
        try:
            clon = feature.clon
        except:
            clon = str(pt)
        try:
            clat = feature.clat
        except:
            clat = str(pt)
        try:
            flightline = feature.flightline
        except:
            flightline = 2
        try:
            glaciertype = feature.gtype
        except:
            glaciertype = 5
        try:
            flowtype = feature.flowtype
        except:
            flowtype = 2
        geometry = feature.GetGeometryRef()
        # Transform to latlon if needed
        if not srs.IsGeographic():
            geometry.TransformTo(srs_geo)
        lon = []
        lat = []

        # This stopped working in gdal 1.11????
        # for point in geometry.GetPoints():
        #     lon.append(point[0])
        #     lat.append(point[1])
        # So here's a bug fix??
        for i in range(0, geometry.GetPointCount()):
            # GetPoint returns a tuple not a Geometry
            pt = geometry.GetPoint(i)
            lon.append(pt[0])
            lat.append(pt[1])
        profiles.append([lat, lon, name, clat, clon, flightline, glaciertype, flowtype])
    return profiles

def get_dims_from_variable(var_dimensions):
    '''
    Gets dimensions from netcdf variable

    Parameters:
    -----------
    var: netCDF variable

    Returns:
    --------
    xdim, ydim, zdim, tdim: dimensions
    '''

    def find(candidates, list):
        """Return one of the candidates if it was found in the list or None
        otherwise."""
        for name in candidates:
            if name in list:
                return name
        return None

    ## possible x-dimensions names
    xdims = ['x','x1']
    ## possible y-dimensions names
    ydims = ['y','y1']
    ## possible z-dimensions names
    zdims = ['z', 'zb']
    ## possible time-dimensions names
    tdims = ['t', 'time']

    return [ find(dim, var_dimensions) for dim in [xdims, ydims, zdims, tdims] ]

def define_profile_variables(nc, profiledim, stationdim):
    # create dimensions
    nc.createDimension(profiledim)
    nc.createDimension(stationdim)

    variables = [("profile_name", str, (stationdim),
                  {"cf_role" : "timeseries_id",
                   "long_name" : "profile name"}),

                 ("profile", "f", (stationdim, profiledim),
                  {"long_name" : 'distance along profile',
                   "units" : "m"}),

                 ("clon", "f", (stationdim),
                  {"long_name" : "center longitude of profile",
                   "units" : "degrees_east",
                   "valid_range" : [-180.0, 180.0]}),

                 ("clat", "f", (stationdim),
                  {"long_name" : "center latitude of profile",
                   "units" : "degrees_north",
                   "valid_range" : [-90.0, 90.0]}),

                 ("lon", "f", (stationdim, profiledim),
                  {"units" : "degrees_east",
                   "valid_range" : [-180.0, 180.0],
                   "standard_name" : "longitude"}),

                 ("lat", "f", (stationdim, profiledim),
                  {"units" : "degrees_north",
                   "valid_range" : [-90.0, 90.0],
                   "standard_name" : "latitude"}),

                 ("flightline", "b", (stationdim),
                  {"long_name" : "flightline (true/false/undetermined) integer mask",
                   "flag_values": [0, 1, 2],
                   "flag_meanings": "true false undetermined",
                   "valid_range" : [0, 2]}),


                 ("flowtype", "b", (stationdim),
                  {"long_name" : "fast-flow type (isbrae/ice-stream) integer mask after Truffer and Echelmeyer (2003)",
                   "flag_values": [0, 1, 2],
                   "flag_meanings": "isbrae ice_stream undetermined",
                   "valid_range" : [0, 2]}),

                 ("glaciertype", "b", (stationdim),
                  {"long_name" : "glacier-type integer mask",
                   "comment": "glacier-type categorization after Moon et al. (2012), Science, 10.1126/science.1219985",
                   "flag_values": [0, 1, 2, 3, 4],
                   "flag_meanings": "fast_flowing_marine_terminating low_velocity_marine_terminating ice_shelf_terminating land_terminating undetermined",
                   "valid_range" : [0, 4]}),

                 ("nx", "f", (stationdim, profiledim),
                  {"long_name" : "x-component of the right-hand-pointing normal vector"}),

                 ("ny", "f", (stationdim, profiledim),
                  {"long_name" : "y-component of the right-hand-pointing normal vector"})]

    for name, type, dimensions, attributes in variables:
        variable = nc.createVariable(name, type, dimensions)
        variable.setncatts(attributes)

def copy_attributes(var_in, var_out):
    _, _, _, tdim = get_dims_from_variable(var_in.dimensions)
    for att in var_in.ncattrs():
        if att == '_FillValue':
            continue
        elif att == 'coordinates':
            if tdim:
                coords = '{0} lat lon'.format(tdim)
            else:
                coords = 'lat lon'
            setattr(var_out, 'coordinates', coords)

        else:
            setattr(var_out, att, getattr(var_in, att))

def copy_global_attributes(in_file, out_file):
    for attribute in in_file.ncattrs():
        setattr(out_file, attribute, getattr(in_file, attribute))

def extract_profile(variable, profile):
    """Extract values of variable along a profile.  """
    xdim, ydim, zdim, tdim = get_dims_from_variable(variable.dimensions)

    group = variable.group()
    x = group.variables[xdim][:]
    y = group.variables[ydim][:]

    n_points = len(profile.x)

    dim_length = dict(zip(variable.dimensions, variable.shape))

    # try to get the matrix we (possibly) pre-computed earlier:
    try:
        A = profile.A
        x_slice = profile.x_slice
        y_slice = profile.y_slice
    except AttributeError:
        # take care of the transpose the easy way (i.e. dealing with 1D objects)
        if variable.dimensions.index(ydim) < variable.dimensions.index(xdim):
            A = ProfileInterpolationMatrix(x, y, profile.x, profile.y)
            x_slice = slice(A.c_min, A.c_max+1)
            y_slice = slice(A.r_min, A.r_max+1)
        else:
            A = ProfileInterpolationMatrix(y, x, profile.y, profile.x)
            x_slice = slice(A.r_min, A.r_max+1)
            y_slice = slice(A.c_min, A.c_max+1)
        profile.A = A
        profile.x_slice = x_slice
        profile.y_slice = y_slice

    def read_subset(t=0, z=0):
        """Assemble the indexing tuple and get a sbset from a variable."""
        index = []
        indexes = {xdim : x_slice,
                   ydim : y_slice,
                   zdim : z,
                   tdim : t}
        for dim in variable.dimensions:
            try:
                index.append(indexes[dim])
            except KeyError:
                index.append(Ellipsis)
        return variable[index]

    if tdim and zdim:
        # 3D time-dependent
        result = np.zeros((dim_length[tdim], n_points, dim_length[zdim]))
        for t in xrange(dim_length[tdim]):
            for z in xrange(dim_length[zdim]):
                raise NotImplementedError
    elif tdim:
        # 2D time-dependent
        result = np.zeros((dim_length[tdim], n_points))
        for j in xrange(dim_length[tdim]):
            result[j, :] = A.apply_to_subset(read_subset(t=j))
    elif zdim:
        # 3D time-independent
        result = np.zeros((dim_length[tdim], n_points))
        for z in xrange(dim_length[zdim]):
            raise NotImplementedError
    else:
        # 2D time-independent
        result = A.apply_to_subset(read_subset())

    return result

def copy_dimensions(in_file, out_file, exclude_list):
    """Copy dimensions from in_file to out_file, excluding ones in
    exclude_list."""
    for name, dim in in_file.dimensions.iteritems():
        if (name not in exclude_list and
            name not in out_file.dimensions):
            if dim.isunlimited():
                out_file.createDimension(name, None)
            else:
                out_file.createDimension(name, len(dim))

def create_variable_like(in_file, var_name, out_file, dimensions=None,
                         fill_value=-2e9):
    """Create a variable in an out_file that is the same var_name in
    in_file, except possibly depending on different dimensions,
    provided in dimensions.

    """
    var_in = in_file.variables[var_name]
    try:
        fill_value = var_in._FillValue
    except:
        pass

    if dimensions is None:
        dimensions = var_in.dimensions

    dtype = var_in.dtype

    var_out = out_file.createVariable(var_name, dtype, dimensions=dimensions,
                                      fill_value=fill_value)
    copy_attributes(var_in, var_out)
    return var_out

def copy_time_dimension(in_file, out_file, name):
    """Copy time dimension, the corresponding coordinate variable, and the
    corresponding time bounds variable (if present) from an in_file to
    an out_file.

    """
    var_out = create_variable_like(nc_in, name, out_file)
    var_out[:] = nc_in.variables[name][:]

    try:
        bounds_name = var_in.bounds
        var_out = create_variable_like(bounds_name, in_file, out_file)
        var_out[:] = in_file.variables[bounds_name][:]

    except AttributeError:
        # we get here if var_in does not have a bounds attribute
        pass
    except Exception as e:
        print "Got an unexpected exception", e

def write_profile(out_file, index, profile):
    """Write information about a profile (name, latitude, longitude,
    center latitude, center longitude, normal x, normal y, distance
    along profile) to an output file.

    """
    ## We have two unlimited dimensions, so we need to assign start and stop
    ## start:stop where start=0 and stop is the length of the array
    ## or netcdf4python will bail. See
    ## https://code.google.com/p/netcdf4-python/issues/detail?id=76
    pl = len(profile.distance_from_start)
    out_file.variables['profile'][k,0:pl] = np.squeeze(profile.distance_from_start)
    out_file.variables['nx'][k,0:pl] = np.squeeze(profile.nx)
    out_file.variables['ny'][k,0:pl] = np.squeeze(profile.ny)
    out_file.variables['lon'][k,0:pl] = np.squeeze(profile.lon)
    out_file.variables['lat'][k,0:pl] = np.squeeze(profile.lat)
    out_file.variables['profile_name'][k] = profile.name
    out_file.variables['clat'][k] = profile.center_lat
    out_file.variables['clon'][k] = profile.center_lon
    out_file.variables['flightline'][k] = profile.flightline
    out_file.variables['glaciertype'][k] = profile.glaciertype
    out_file.variables['flowtype'][k] = profile.flowtype

if __name__ == "__main__":
    # Set up the option parser
    description = '''A script to extract data along (possibly multiple) profile using
    piece-wise constant or bilinear interpolation.
    The profile must be given as a ESRI shape file.'''
    parser = ArgumentParser()
    parser.description = description
    parser.add_argument("FILE", nargs='*')
    parser.add_argument(
        "-b", "--bilinear",dest="bilinear",action="store_true",
        help='''Piece-wise bilinear interpolation, Default=False''',
        default=False)
    parser.add_argument(
        "-f", "--flip",dest="flip",action="store_true",
        help='''Flip profile direction, Default=False''',
        default=False)
    parser.add_argument(
        "-t", "--print_timing",dest="timing",action="store_true",
        help='''Print timing information, Default=False''',
        default=False)
    parser.add_argument("-v", "--variable",dest="variables",
                        help="comma-separated list with variables",default='x,y,thk,velsurf_mag,flux_mag,uflux,vflux,pism_config,pism_overrides,run_stats,uvelsurf,vvelsurf,topg,usurf,tillphi,tauc')
    parser.add_argument(
        "-a", "--all_variables",dest="all_vars",action="store_true",
        help='''Process all variables, overwrite -v/--variable''',
        default=False)

    options = parser.parse_args()
    bilinear = options.bilinear
    args = options.FILE
    flip = options.flip
    timing = options.timing
    fill_value = -2e9
    variables = options.variables.split(',')
    all_vars = options.all_vars
    n_args = len(args)
    required_no_args = 2
    max_no_args = 3
    if (n_args < required_no_args):
        print(("received $i arguments, at least %i expected"
              % (n_args, required_no_args)))
        import sys.exit
        sys.exit
    elif (n_args > max_no_args):
        print(("received $i arguments, no more thant %i accepted"
              % (n_args, max_no_args)))
        import sys.exit
        sys.exit
    else:
        p_filename = args[0]
        in_filename = args[1]
        if (n_args == 2):
            out_filename = 'profile.nc'
        else:
            out_filename = args[2]

    print("-----------------------------------------------------------------")
    print("Running script %s ..." % __file__.split('/')[-1])
    print("-----------------------------------------------------------------")
    print("Opening NetCDF file %s ..." % in_filename)
    try:
        # open netCDF file in 'read' mode
        nc_in = NC(in_filename, 'r')
    except:
        print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
              % in_filename))
        import sys
        sys.exit()

    # get the dimensions
    xdim, ydim, zdim, tdim = ppt.get_dims(nc_in)
    # read projection information
    projection = ppt.get_projection_from_file(nc_in)

    # Read in profile data
    print("  reading profile from %s" % p_filename)
    profiles = load_profiles(p_filename, projection, flip)

    mapplane_dim_names = (xdim, ydim)

    print("Creating dimensions")
    nc = NC(out_filename, 'w', format='NETCDF4')
    copy_global_attributes(nc_in, nc)

    profiledim = 'profile'
    stationdim = 'station'

    # define variables storing profile information
    define_profile_variables(nc, profiledim, stationdim)
    # fill these variables
    for k, profile in enumerate(profiles):
        write_profile(nc, k, profile)

    # re-create dimensions from an input file in an output file, but
    # skip x and y dimensions and dimensions that are already present
    copy_dimensions(nc_in, nc, mapplane_dim_names)

    # figure out which variables do not need to be copied to the new file.
    # mapplane coordinate variables
    vars_not_copied = ['lat', 'lat_bnds', 'lat_bounds',
                       'lon', 'lon_bnds', 'lon_bounds',
                       xdim, ydim, tdim]
    for var_name in nc_in.variables:
        var = nc_in.variables[var_name]
        if hasattr(var, 'grid_mapping'):
            mapping_var_name = var.grid_mapping
            vars_not_copied.append(mapping_var_name)
        if hasattr(var, 'bounds'):
            bounds_var_name = var.bounds
            vars_not_copied.append(bounds_var_name)

    vars_not_copied.sort()
    last = vars_not_copied[-1]
    for i in range(len(vars_not_copied)-2, -1, -1):
        if last == vars_not_copied[i]:
            del vars_not_copied[i]
        else:
            last = vars_not_copied[i]

    if tdim is not None:
        copy_time_dimension(nc_in, nc, tdim)

    print("Copying variables")
    if all_vars:
        vars_list = nc_in.variables
        vars_not_found = ()
    else:
        vars_list = filter(lambda(x): x in nc_in.variables, variables)
        vars_not_found =  filter(lambda(x): x not in nc_in.variables, variables)

    for var_name in vars_list:
        profiler = timeprofile()
        if var_name in vars_not_copied:
            continue

        print("  Reading variable %s" % var_name)

        var_in = nc_in.variables[var_name]
        in_dims = var_in.dimensions
        datatype = var_in.dtype

        if in_dims and len(in_dims) > 1:
            # it is a non-scalar variable and it depends on more
            # than one dimension, so we probably need to extract profiles
            out_dims = output_dimensions(in_dims, stationdim, profiledim)
            var_out = create_variable_like(nc_in, var_name, nc, dimensions=out_dims)

            for k, profile in enumerate(profiles):
                print("    - processing profile {0}".format(profile.name))
                p_values = extract_profile(var_in, profile)

                profiler.mark('write')
                try:
                    # try without exec (should work using newer netcdf4-python)
                    indexes = np.r_[k, [np.s_[0:n] for n in p_values.shape]]
                    var_out[indexes] = p_values
                except:
                    access_str = 'k,' + ','.join([':'.join(['0', str(coord)]) for coord in p_values.shape])
                    exec('var_out[%s] = p_values' % access_str)
                p_write = profiler.elapsed('write')

                if timing:
                    print('''    - read in %3.4f s, written in %3.4f s''' % (p_read, p_write))
        else:
            # it is a scalar or a 1D variable; just copy it
            var_out = create_variable_like(nc_in, var_name, nc)
            var_out[:] = var_in[:]

        copy_attributes(var_in, var_out)
        print("  - done with %s" % var_name)

    print("The following variables were not copied because they could not be found in {}:".format(in_filename))
    print vars_not_found

    # writing global attributes
    script_command = ' '.join([time.ctime(), ':', __file__.split('/')[-1],
                               ' '.join([str(l) for l in args])])
    if hasattr(nc_in, 'history'):
        history = nc_in.history
        nc.history = script_command + '\n ' + history
    else:
        nc.history = script_command

    nc_in.close()
    nc.close()
    print("Extracted profiles to file %s" % out_filename)
