#!/usr/bin/env python
# Copyright (C) 2017 Andy Aschwanden

import os
import numpy as np
import logging
import logging.handlers
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from netCDF4 import Dataset as NC

# create logger
logger = logging.getLogger('hillshade')
logger.setLevel(logging.INFO)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter
formatter = logging.Formatter('%(name)s - %(module)s - %(message)s')

# add formatter to ch and fh
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


class Hillshade(object):

    '''
    A class to add a hillshade to a netCDF time-series.

    Parameters
    ----------

    ifile: netCDF file with dimensions ('time', 'y', 'x'). Other permutations are currently not supported
    variable: variable used to create a hillshade
    threshold_masking: bool. If True, variable[threshold_masking_variable < threshold_masking_value] = fill_value is set

    kwargs
    ----------

    altitude:
    azimuth:
    fill_value:
    threshold_masking_variable: if threshold_masking is True, use this variable to mask variable
    threshold_masking_value: if threshold_masking is True, use this value to mask variable
    zf: 
    '''

    def __init__(self, ifile, variable='usurf', threshold_masking=True, variables_to_mask=None, *args, **kwargs):
        super(Hillshade, self).__init__(*args, **kwargs)

        self.threshold_masking = threshold_masking
        self.do_masking = False
        self.ifile = ifile
        self.variable = variable
        if variables_to_mask is not None:
            self.variables_to_mask = variables_to_mask.split(',')
            self.do_masking = True
        else:
            self.variables_to_mask = variables_to_mask
        self.params = {'altitude': 45,
                       'azimuth': 45,
                       'fill_value': 0,
                       'threshold_masking_variable': 'thk',
                       'threshold_masking_value': 10,
                       'zf': 5}
        for key, value in kwargs:
            if key in ('altitude', 'azimuth', 'fill_value', 'hillshade_var', 'zf'):
                self.params[key] = value

        self._check_vars()
        self.dx = self._get_dx()
        self._create_vars()

    def _check_vars(self):
       
        nc = NC(self.ifile, 'r')
        for mvar in (['time'] + [self.variable]):
            if mvar in nc.variables:
                logger.info('variable {} found'.format(mvar))
            else:
                logger.info('variable {} NOT found'.format(mvar))

        if self.do_masking:
            for mvar in (self.variables_to_mask + [self.params['threshold_masking_variable']]):
                if mvar in nc.variables:
                    logger.info('variable {} found'.format(mvar))
                else:
                    logger.info('variable {} NOT found'.format(mvar))
        nc.close()
           
    def _cart2pol(self, x, y):
        '''
        cartesian to polar coordinates
        '''
        theta = np.arctan2(y, x)
        rho = np.sqrt(x**2 + y**2)
        return (theta, rho) 

    def _create_vars(self):
        '''
        create netCDF variables if they don't exist yet
        '''

        ifile = self.ifile
        nc = NC(ifile, 'a')
        variable = self.variable
        hs_var = variable  + '_hs'
        if hs_var  not in nc.variables:
            nc.createVariable(hs_var, 'i', dimensions=('time', 'y', 'x'), fill_value=fill_value)
        nc.close()
                
    def _get_dx(self):
        
        nc = NC(self.ifile, 'r')

        x0, x1 = nc.variables['x'][0:2]
        y0, y1 = nc.variables['y'][0:2]
        
        nc.close()

        dx = x1 - x0
        dy = y1 - y0

        assert dx == dy

        return dx

    def _hillshade(self, dem):
       '''
       shaded relief using the ESRI algorithm
       '''

       # lighting azimuth
       azimuth = self.params['azimuth']
       azimuth = 360.0 - azimuth + 90 # convert to mathematic unit
       if (azimuth>360) or (azimuth==360):
          azimuth = azimuth - 360
       azimuth = azimuth * (np.pi / 180)  # convert to radians

       # lighting altitude
       altitude = self.params['altitude']
       altitude = (90 - altitude) * (np.pi / 180)  # convert to zenith angle in radians

       # calc slope and aspect (radians)
       dx = self.dx
       fx, fy = np.gradient(dem, dx)  # uses simple, unweighted gradient of immediate
       [asp, grad] = self._cart2pol(fy, fx)  # convert to carthesian coordinates

       zf = self.params['zf']
       grad = np.arctan(zf * grad)  # steepest slope
       # convert asp
       asp[asp<np.pi] = asp[asp < np.pi]+(np.pi / 2)
       asp[asp<0] = asp[asp<0] + (2 * np.pi)

       ## hillshade calculation
       h = 255.0 * ((np.cos(altitude) * np.cos(grad))
                    + (np.sin(altitude) * np.sin(grad) * np.cos(azimuth - asp)))
       h[h<0] = 0  # set hillshade values to min of 0.

       return h

    def run(self):
        logger.info('Processing file {}'.format(ifile))
        fill_value = self.params['fill_value']
        hs_var = self.variable  + '_hs'
        nc = NC(ifile, 'a')
        nt = len(nc.variables['time'][:])
        for t in range(nt):
            logger.info('Processing time {} of {}'.format(t, nt))
            dem = nc.variables['usurf'][t, Ellipsis]
            hs = self._hillshade(dem)
            hs[dem==0] = fill_value
            nc.variables[hs_var][t, Ellipsis] = hs
            if self.threshold_masking:
                m = nc.variables[self.params['threshold_masking_variable']][t, Ellipsis]
                hs[m <= self.params['threshold_masking_value']] = fill_value
            if self.do_masking:
                for mvar in self.variables_to_mask:
                    mt = nc.variables[self.params['threshold_masking_variable']][t, Ellipsis]
                    m = nc.variables[mvar][t, Ellipsis]
                    try:
                        m_fill_value = nc.variables[mvar]._FillValue
                    except:
                        m_fill_value = fill_value
                    m[mt < self.params['threshold_masking_value']] = m_fill_value
                    nc.variables[mvar][t, Ellipsis] = m
            
        nc.close()


if __name__ == "__main__":

    # set up the option parser
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.description = "Adding a hillshade to a netCDF file."
    parser.add_argument("FILE", nargs=1,
                        help="netCDF file with dimensions ('time', 'y', 'x'). Other permutations are currently not supported")
    parser.add_argument("--variable",
                        help="file", default=None)
    parser.add_argument("--variable",
                        help="variable used to create a hillshade", default='usurf')

    options = parser.parse_args()
    ifile = options.FILE[0]
    variables = options.varible

    hs = Hillshade(ifile, variable)
    hs.run() 

