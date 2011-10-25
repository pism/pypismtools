'''
PyPISMTools: Tools to evaluate PISM parameter studies

PyPISMTools is a module to facilitate evaluation of PISM parameter
studies. It mainly comprises two classes, Observation and Experiment,
which act as containers for observational data and PISM model
simulations, along with helper functions. The experiment class
determines information about an experiment from the netcdf file
directly, especially from the 'pism_overrides' flag. Such information
can then be used for labeling, plotting, evaluation, etc. The indend
is to provide a robust tool to evaluate data, and to avoid common mistakes
such as mis-labeling plots.
'''
__author__ = "Andy Aschwanden"

__all__ = ['norm_infinity','norm_Lp','norm_2','norm_1','get_rmse','get_avg','unit_converter','permute','plot_mapview','plot_histogram','plot_histogram2','print_info','print_overall_statistics','Observation','Experiment']

import numpy as np
import pylab as plt

try:
    from netCDF3 import Dataset as CDF
except:
    from netCDF4 import Dataset as CDF

## FIXME: how to provide DEBUG flag to module
DEBUG=None

def norm_infinity(var):
    '''
    Computes the infinity norm of a numpy array.
    
    Paramters
    ---------
    var : array_like

    Returns
    -------
    norm : scalar

    Notes
    -----
    The infinity norm is the given by the absolute value of the
    maximum of a given array.
    '''
    return np.abs(var).max()

def norm_Lp(var, p, dx = 1, dy = 1):
    '''
    Computes the L^p norm of a numpy array.

    Paramters
    ---------
    var : array_like

    Returns
    -------
    norm : scalar

    Notes
    -----
    The L^p norm is the given by sum(abs(var)**p).
    '''
    return np.sum(np.abs(var)**p * dx * dy)**(1.0/p)

def norm_2(var, dx = 1, dy = 1):
    '''
    Computes the L^2 norm of a numpy array.

    Paramters
    ---------
    var : array_like

    Returns
    -------
    norm : scalar

    Notes
    -----
    The L^2 norm is the given by sum(abs(var)**2).
    '''
    return norm_Lp(var, 2, dx, dy)

def norm_1(var, dx = 1, dy = 1):
    '''
    Computes the L^1 norm of a numpy array.

    Paramters
    ---------
    var : array_like

    Returns
    -------
    norm : scalar

    Notes
    -----
    The L^1 norm is the given by sum(abs(var)).
    '''
    return norm_Lp(var, 1, dx, dy)

def get_rmse(a,b,N):
    '''
    Returns the root mean square error of differences between a and b.

    Parameters
    ----------
    a,b : array_like
    N : number of values

    Returns
    -------
    rmse : scalar

    Notes
    -----
    The average is the sum of elements of the difference (a - b)
    divided by the number of elements N.
    '''
    return np.sqrt((norm_2(a - b, dx=1, dy=1))**2.0 / N)
    
def get_avg(a,b,N):
    '''
    Returns the average difference between a and b.
    
    Parameters
    ----------
    a,b : array_like
    N : number of values

    Returns
    -------
    avg : scalar
    '''
    return norm_1(a - b, dx=1, dy=1) / N
    
def unit_converter(data,inunit,outunit):
    '''
    Unit converter. Takes an (numpy) array, valid udunits inunits and outunits as strings, and returns the array in outunits.
    
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

    if not (inunit == outunit):
        try:
            from udunits import Converter
            c = Converter(inunit,outunit)
            outdata=c(data)
        except:
            print("No udunits module found, you're on your own.\n  -> I am assuming input unit is m, and will convert to km.\n  -> Installation of Constantine's awesome python wrapper for udunits is highly recommended.\n  -> Download it from https://code.google.com/p/python-udunits/.")
            c = 1./1e3
            outdata = c*data
    else:
        outdata = data

    return outdata

def permute(variable, output_order = ('time', 'z', 'zb', 'y', 'x')):
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
    dimensions = filter(lambda(x): x in input_dimensions,
                        output_order)

    # create the mapping
    mapping = map(lambda(x): dimensions.index(x),
                  input_dimensions)

    if mapping:
        return np.transpose(variable[:], mapping)
    else:
        return variable[:]              # so that it does not break processing "mapping"

def plot_mapview(var,**kwargs):
    '''Plot map view of variable var'''
    import pylab as plt
    from matplotlib import colors
    kwargsdict = {}
    expected_args = ["log","show","title"]
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
    plt.imshow(var,origin='lower',norm=norm)
    plt.colorbar()
    plt.title(title)

def plot_histogram(data,**kwargs):
    '''Plot a simple histogram with one variable'''
    import pylab as plt

    kwargsdict = {}
    expected_args = ["log","show"]
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
    ticks = map(lambda(x): "%3d" % x, bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n,bin,patches = ax.hist(data.values[data.values>vmin],
                            bins = bins,
                            log = log,
                            histtype='bar')
    ax.set_xlabel("ice surface velocity, m a$^{-1}$")
    ax.set_ylabel("number of grid cells")
    plt.title(data.title)
    return n,bin,patches


def plot_histogram2(data,obs,**kwargs):
    '''Plot a histogram with modeled and observed values'''
    import pylab as plt

    kwargsdict = {}
    expected_args = ["log","show"]
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
    ticks = map(lambda(x): "%3d" % x, bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n,bin,patches = ax.hist([data.values[data.values>vmin],obs.values[obs.values>vmin]],
                            bins = bins,
                            log = log,
                            histtype='bar',
                            label=['modeled','observed'])
    ax.set_xlabel("ice surface velocity, m a$^{-1}$")
    ax.set_ylabel("number of grid cells")
    ax.set_xlim(vmin,vmax)
    plt.title(data.title)
    ax.legend()

def print_info(experiment):
    '''Prints parameters and values of an experiment'''
    for key,val in experiment.parameter_dict.iteritems():
        if isinstance(val,basestring):
            print("      * %s = %s" %(key,val))
        else:
            print("      * %s = %3.3f" %(key,val))


def print_overall_statistics(experiments):
    '''Prints some statistics'''
    print("\nOverall statistics:\n")

    avgs = [n.avg for n in experiments]
    avg_min = filter(lambda(x): x.avg==min(avgs),experiments)
    avg_max = filter(lambda(x): x.avg==max(avgs),experiments)

    print("    - smallest avg = %3.2f" %avg_min[0].avg)    
    print_info(avg_min[0])
    print("    - largest avg = %3.2f" %avg_max[0].avg)
    print_info(avg_max[0])

    rmses = [n.rmse for n in experiments]
    rmse_min = filter(lambda(x): x.rmse==min(rmses),experiments)
    rmse_max = filter(lambda(x): x.rmse==max(rmses),experiments)

    print("    - smallest rmse = %3.2f" %rmse_min[0].rmse)    
    print_info(rmse_min[0])
    print("    - largest rmse = %3.2f" %rmse_max[0].rmse)
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
    def __init__(self,file_name,variable_name,bins,*args, **kwargs):
        super(DataObject,self).__init__(*args, **kwargs)
        self.file_name = file_name
        self.variable_name = variable_name
        self.bins = bins
        self.bin_bounds = [bins[0],bins[-1]]

    def __del__(self):
        # Close open file
        self.nc.close()

    def _make_histogram(self):
        '''Calculate number of grid cells n in bins bin.'''
        # np.histogram ignores masked array, so we have to exclude
        # masked values by applying the inverted mask
        n,Bins = np.histogram(self.values[~self.mask],bins=self.bins)
        self.n = n
        self.no_above_histmax = np.sum(self.values>=self.bins[-1])
        self.cells_in_histogram = np.sum(self.n)
        self.coverage = np.float(self.cells_in_histogram) / self.valid_cells * 100
        if DEBUG:
            print("    - no of cells in histogram = %i (%3.2f %%)" % (self.cells_in_histogram,self.coverage))
            print("    - no of cells above %3.2f %s = %i" % (self.bin_bounds[1],self.units,self.no_above_histmax))
            print("    - no of cells in bins =\n   %s" % str(self.n))


    def add_mask(self,mask):
        '''
        Add/update a mask.

        Adds (or updates) the mask. self.mask = logical_or(self.mask,mask).

        Parameters
        ----------
        mask : array_like
               array with zeros (no mask) and ones (mask)
        '''
        self.mask = np.logical_or(self.mask,mask)
        self.values = np.ma.array(data=self.values,mask=self.mask)
        self._set_valid_cells()
        if DEBUG:
            plot_mapview(self.mask,show=True,title='updated mask')
        self._make_histogram()
        
    def get_total_cells(self):
        '''Return total number of grid cells.'''
        nx,ny = np.shape(self.values)
        return nx*ny

    def _set_valid_cells(self):
        '''Determine number of valid cells by summing up non-masked cells.'''
        value = np.sum(-self.mask)
        self.valid_cells = value

    def _set_grid_spacing(self):
        '''Tries to determine grid spacing in km, sets None if it fails.'''
        try:
            x_dim = self.nc.variables["x"]
            in_units = self.nc.variables["x"].units
        except:
            x_dim = self.nc.variables["x1"]
            in_units = self.nc.variables["x1"].units
        out_units = "km"
        grid_spacing = unit_converter(np.abs(x_dim[1]-x_dim[0]),in_units,out_units)
        self.grid_spacing = grid_spacing
        self.grid_spacing_units = out_units
        if DEBUG:
            print("grid spacing is %3.1f %s" %(self.grid_spacing,self.grid_spacing_units))

    def get_grid_spacing(self):
        '''Return grid spacing.'''
        if self.grid_spacing is not None:
            print("grid spacing is %3.1f %s" %(self.grid_spacing,self.grid_spacing_units))
            return self.grid_spacing
        else:
            return self.grid_spacing
        
    def get_valid_cells(self):
        '''Return number of valid cells.'''
        return self.valid_cells

    def set_rmse(self,rmse):
        '''Set root mean square error (RMSE).'''
        self.rmse = rmse

    def set_avg(self,avg):
        '''Set average (avg).'''
        self.avg = avg

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
    
    def __init__(self,obs,outlier,parameter_list,abbr_dict,*args, **kwargs):
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
            for key,val in self.parameter_dict.iteritems():
                if isinstance(val,basestring):
                    print("    * %s = %s" %(key,val))
                else:
                    print("    * %s = %3.3f" %(key,val))

            
        # Open netcdf file
        try:
            self.nc = CDF(self.file_name,'r')
        except IOError as e:
            print("File %s does not exist or is not a valid NetCDF file" % self.file_name)
        print("\nInitializing experiement from file %s" % self.file_name)
        print("  - reading variable %s" % self.variable_name)
        # Determine units from 'units' attribute, set to None if not found
        try:
            self.units = self.nc.variables[self.variable_name].units
        except:
            self.units = None
        # Read values, squeeze and permute (if appropriate)
        try:
            self.values = np.squeeze(permute(self.nc.variables[self.variable_name]))
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

        nx,ny = np.shape(self.values)
        self.total_cells = nx*ny
        print("  - total cells = %i" % self.total_cells)        
        # get mask from masked array "variable" if exist
        try:
            self.mask = self.values.mask
            valid_cells = np.sum(-self.mask)
        except:
            valid_cells = self.total_cells
            self.mask = np.zeros_like(self.values,dtype='bool')
        self.valid_cells = valid_cells
        print("  - valid cells = %i" % valid_cells)
            
        self.parameter_list = parameter_list
        self.pism_overrides = self.nc.variables["pism_overrides"]
        self.parameter_dict = dict([(param, self.pism_overrides.getncattr(param) ) for param in parameter_list])
        self.parameter_short_dict = dict([(abbr_dict[param], self.pism_overrides.getncattr(param) ) for param in parameter_list])
        _print_info()

        self.title = ', '.join(["%s = %s" % (key,str(val)) for key,val in self.parameter_short_dict.iteritems()])
        self._make_histogram()
        if DEBUG:
            plot_mapview(self.values,log=True,title='before no obs')
        print("\n  * Applying mask where no observations, updating histogram")
        self.add_mask(obs.mask)
        if DEBUG:
            plot_mapview(self.values,log=True,title='after obs')
        print("\n  * Applying mask where abs(%s - %s) > %3.2f %s, updating histogram" %(self.variable_name,
                                                                                       obs.variable_name,outlier,self.units))
        outliers = (np.abs(self.values - obs.values) > outlier)
        self.add_mask(outliers)
        if DEBUG:
            plot_mapview(self.values,log=True,title='after outliers')

        print("    - no of cells in histogram = %i (%3.2f %%)" % (self.cells_in_histogram,self.coverage))
        print("    - no of cells above %3.2f %s = %i" % (self.bin_bounds[1],self.units,self.no_above_histmax))
        print("    - no of cells in bins =\n   %s" % str(self.n))

 
    def get_values(self):
        '''Return values'''
        return self.values

    def set_coverage(self,value):
        '''Set coverage defined as number of cells used for comparison divided by number of ice covered cells'''
        self.coverage = value

    def set_outliers(self,outliers):
        self.outliers = outliers
        self.no_outliers = np.sum(outliers)


class Observation(DataObject):
    """
    A derived class for observations.

    A derived class to handle observational data. The base class takes
    care of initializing.
    """
    
    def __init__(self,*args, **kwargs):
        super(Observation, self).__init__(*args, **kwargs)
        ## self.title = "observed " + self.variable_name
        self.title = "observed"
        # Open netcdf file
        try:
            self.nc = CDF(self.file_name,'r')
        except IOError as e:
            print("File %s does not exist or is not a valid NetCDF file" % self.file_name)
        print("\nInitializing '%s' observations from file %s" % (self.variable_name,self.file_name))
        print("  - reading variable %s" % self.variable_name)
        # Determine untis from 'units' attribute, set to None if not found
        try:
            self.units = self.nc.variables[self.variable_name].units
        except:
            self.units = None
        # Read values, squeeze and permute (if appropriate)
        try:
            self.values = np.squeeze(permute(self.nc.variables[self.variable_name]))
        except:
            try:
                self.values = np.squeeze(permute(self.nc.variables[self.variable_name],("time","y1","x1")))
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

        nx,ny = np.shape(self.values)
        self.total_cells = nx*ny
        print("  - total cells = %i" % self.total_cells)        
        # get mask from masked array "variable" if exist
        try:
            self.mask = self.values.mask
            valid_cells = np.sum(-self.mask)
        except:
            valid_cells = self.total_cells
            self.mask = np.zeros_like(self.values,dtype='bool')
        self.valid_cells = valid_cells
        print("  - valid cells = %i" % valid_cells)
        if DEBUG:
            plot_mapview(self.mask,show=True,title='obs mask')
            
        self._make_histogram()
        print("  - no of cells in histogram = %i (%3.2f %%)" % (self.cells_in_histogram,self.coverage))
        print("  - no of cells above %3.2f %s = %i (%3.2f %%)" % (self.bins[1],self.units,self.valid_cells-self.n[0],float(self.valid_cells-self.n[0])/self.valid_cells*100))
        print("  - no of cells above %3.2f %s = %i (%3.2f %%)" % (self.bin_bounds[1],self.units,self.no_above_histmax,100-self.coverage))
        print("  - no of cells in bins =\n   %s" % str(self.n))
