#!/usr/bin/env python

import numpy as np
import pylab as plt
from skimage import measure
from argparse import ArgumentParser

from netCDF4 import Dataset as NC

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt

from udunits2 import Converter, System, Unit
sys = System()

class DischargeAnalysis(object):
    '''
    A container class for mapplane discharge objects.

    Does calculations.
    '''
    
    def __init__(self, file_names, *args, **kwargs):
        
        for ikey in kwargs.keys():
            setattr(self, ikey, kwargs[ikey])

        analysis = []
        nt = len(file_names)
        self.nt = nt
        for k in range(0, nt):
            file_name = file_names[k]
            d = Discharge(file_name)
            analysis.append(d)

        self.analysis = analysis
        self.area_units = d.area_units
        self.icethk_units = d.icethk_units
            
    def icethk_cum(self):
        nt = self.nt
        icethk = np.zeros((nt))
        for k in range(0, nt):
            icethk[k] = self.analysis[k].icethk_avg.sum()
        return icethk

    def area_cum(self):
        nt = self.nt
        area = np.zeros((nt))
        for k in range(0, nt):
            area[k] = self.analysis[k].area.sum()
        return area

    def grid_spacings(self):
        nt = self.nt
        dx = np.zeros((nt))
        for k in range(0, nt):
            dx[k] = self.analysis[k].dx
        return dx
    
    def grid_units(self):
        units = self.analysis[0].dx_units
        return units


class Discharge(DischargeAnalysis):
    '''
    A class for a mapplane discharge object
    '''
    
    def __init__(self, file_name, *args, **kwargs):

        super(DischargeAnalysis, self).__init__(**kwargs)
        
        for ikey in kwargs.keys():
            setattr(self, ikey, kwargs[ikey])

        print("  opening NetCDF file %s ..." % file_name)
        try:
            # open netCDF file
            nc = NC(file_name, 'r')
        except:
            print(("ERROR:  file '%s' not found or not NetCDF format ... ending ..."
                  % file_name))
            import sys as s
            s.exit(1)

        # get the dimensions
        xdim, ydim, zdim, tdim = ppt.get_dims(nc)
        # set up dimension ordering
        dim_order = (tdim, zdim, ydim, xdim)
        # add lat/lon values
        x = (np.squeeze(ppt.permute(nc.variables[xdim], dim_order)))
        x_units = nc.variables[xdim].units
        y = (np.squeeze(ppt.permute(nc.variables[ydim], dim_order)))
        y_units = nc.variables[ydim].units

        var = 'topg'
        print(("    - reading variable %s" % var))
        try:
            topg = np.squeeze(ppt.permute(nc.variables[var], dim_order))
            topg_units = nc.variables[var].units
        except:
            print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
                  % (var, file_name)))
            import sys as s
            s.exit(1)

        mask = (topg >= 0)
        topg = np.ma.array(topg, mask = mask)
        self.topg = topg

        var = 'thk'
        print(("    - reading variable %s" % var))
        try:
            thk = np.squeeze(ppt.permute(nc.variables[var], dim_order))
            thk_units = nc.variables[var].units
        except:
            print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
                  % (var, file_name)))
            import sys as s
            s.exit(1)

        self.thk = thk

        speed_units = 'm year-1'

        var = 'ubar'
        print(("    - reading variable %s" % var))
        try:
            ubar = np.squeeze(ppt.permute(nc.variables[var], dim_order))
            ubar_units = sys, nc.variables[var].units
        except:
            print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
                  % (var, file_name)))
            import sys as s
            s.exit(1)

        ubar = ppt.unit_converter(ubar, ubar_units, speed_units)
        self.ubar = ubar

        var = 'vbar'
        print(("    - reading variable %s" % var))
        try:
            vbar = np.squeeze(ppt.permute(nc.variables[var], dim_order))
            vbar_units = nc.variables[var].units
        except:
            print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
                  % (var, file_name)))
            import sys as s
            s.exit(1)

        vbar = ppt.unit_converter(vbar, vbar_units, speed_units)
        self.vbar = vbar

        print(("    - reading variable %s" % var_name))
        try:
            data = np.squeeze(ppt.permute(nc.variables[var_name], dim_order))
            data_units = nc.variables[var_name].units
        except:
            print(("ERROR:  unknown or not-found variable '%s' in file %s ... ending ..."
                  % (var_name, file_name)))
            import sys as s
            s.exit(1)

        try:
            inunit = str(nc.variables[var_name].units)
        except:
            print(("ERROR:  units not found in variable '%s' in file %s ... ending ..."
                  % (var_name, file_name)))
            import sys as s
            s.exit(1)

        if outunit is not None:
            data = ppt.unit_converter(data, inunit, outunit)

        mask = (data >= min_discharge)
        data = np.ma.array(data, mask = mask)
        self.data = data

        outdimunits = 'm'
        dx = ppt.unit_converter(np.abs(x[1] - x[0]), x_units, outdimunits)
        dy = ppt.unit_converter(np.abs(y[1] - y[0]), y_units, outdimunits)
        self.dx = dx
        self.dx_units = outdimunits
        self.dy = dy
        self.dy_units = outdimunits
        
        velbar = np.sqrt(ubar ** 2 + vbar ** 2)

        is_discharge = data.nonzero()
        self.is_discharge = is_discharge
        
        # get number of non-zero non-masked cells
        n_cells = data[is_discharge].data.shape[0]
        self.n_cells = n_cells

        # stencil width
        n = 3

        nx, ny = thk.shape
        ii, jj = np.indices((n,n))

        icethk_avg = np.zeros((n_cells))
        gatethk_avg = np.zeros((n_cells))
        velbar_avg = np.zeros((n_cells))

        for k in range(0, n_cells):
            r = (ii + is_discharge[0][k]) % (nx - 1)  # periodic stencil
            c = (jj + is_discharge[1][k]) % (ny - 1)  # periodic stencil
            icethk_avg[k] = thk[r,c].sum() / len(thk[r,c].nonzero())
            gatethk_avg[k] = np.abs(topg[r,c]).sum() / len(np.abs(topg[r,c]).nonzero())
            velbar_avg[k] = velbar[r,c].sum() / len(velbar[r,c].nonzero())

        self.icethk_avg = icethk_avg
        self.gatethk_avg = gatethk_avg
        self.velbar_avg = velbar_avg
        self.icethk_units = thk_units
        self.gatethk_units = topg_units
        self.velbar_units = speed_units
        
        area = icethk_avg * dx
        area_units = str(Unit(sys, thk_units) * Unit(sys, outdimunits))
        self.area = area
        self.area_units = area_units

        nc.close()

    

# Set up the argument parser
parser = ArgumentParser()
parser.description = '''Calculate mapplane fluxes, cross-sectional areas, etc.'''
parser.add_argument("FILE", nargs='*')
parser.add_argument("-p", "--print_size", dest="print_mode",
                  help='''sets figure size and font size, available options are:
                  onecol, medium, presentation, twocol''', default="medium")
parser.add_argument("-o", "--output_file", dest="out_file",
                  help='''output file name without suffix''', default='unnamed')
parser.add_argument("-r", "--output_resolution", dest="out_res",
                  help='''Resolution ofoutput graphics in dots per
                  inch (DPI), default = 300''', default=300)

options = parser.parse_args()
file_names = options.FILE
print_mode = options.print_mode
var_name = 'ocean_kill_flux'
outunit = 'Gt year-1'
out_file = options.out_file
out_formats = ['pdf']
out_res = options.out_res

min_discharges = [0,-0.01, -0.1, -0.5, -1]

d = []
for min_discharge in min_discharges:
    d_analysis = DischargeAnalysis(file_names,
                                   min_discharge=min_discharge)
    d.append(d_analysis)


# set the print mode
golden_mean = ppt.get_golden_mean()
aspect_ratio = golden_mean
lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=golden_mean)
colors = ppt.colorList()

markersize = 5
params = {'lines.markersize': markersize}
plt.rcParams.update(params)

fig = plt.figure()
ax = fig.add_subplot(111)
for k in range(0, d_analysis.nt):
    area = d_analysis.analysis[k].area
    area_units = d_analysis.area_units
    area = ppt.unit_converter(area, area_units, 'km2')
    ax.plot(area, 'd', color=colors[k])

fig = plt.figure()
ax = fig.add_subplot(111)
for k in range(0, len(min_discharges)):
    area_cum = d[k].area_cum()
    area_units = d[k].area_units
    grid_spacings = d[k].grid_spacings()
    grid_units = d[k].grid_units()
    grid_spacings = ppt.unit_converter(grid_spacings, grid_units, 'km')
    area_cum = ppt.unit_converter(area_cum, area_units, 'km2')
    ax.plot(grid_spacings, area_cum, 'd', color=colors[k], lw=lw)
plt.axhspan(140, 290, facecolor='0.5', alpha=0.5)
ax.set_xlabel('grid spacing [km]')
ax.set_ylabel('total flux gate area [km$^2$]')
ax.set_xlim([0, 25])
ax.set_yscale('log')
ax.legend(min_discharges, loc='upper right', numpoints=1,
          title='d$_{min}$ per cell [Gt a$^{-1}$]',
          bbox_to_anchor=(0, 0, 1, 1.085), borderaxespad=0.,
          bbox_transform=plt.gcf().transFigure)
for out_format in out_formats:
    out_name = out_file + '.' + out_format
    print("    - Saving %s." % out_name)
    plt.savefig(out_name, bbox_inches='tight', dpi=out_res)
