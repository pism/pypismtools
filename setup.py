from distutils.core import setup
import os
import sys

PKG_NAME = 'pypismtools'
print(('\nInstalling %s' % PKG_NAME))
print('----------------------------')

packages = ('numpy', 'osgeo', 'netCDF4', 'pyproj', 'cf_units',
            'mpl_toolkits.basemap', 'mpl_toolkits.axes_grid1')
print('\nChecking dependencies:')
not_installed = []
for package in packages:
    try:
        __import__(package)
        print(('  - % s package is installed' % package))
    except ImportError:
        print(('  - % s package NOT installed' % package))
        not_installed.append(package)
if not_installed:
    print('Installation of the following packages is optional but recommended:')
    for package in not_installed:
        if package is 'osgeo':
            print(' - GDAL python bindings')
        else:
            print(('  - %s' % package))
    print('Exiting')
    import sys
    sys.exit()

setup(name=PKG_NAME,
      version='0.15',
      description='Python tools to evaluate PISM results',
      author='Andy Aschwanden',
      author_email='aaschwanden@alaska.edu',
      url='https://github.com/pism/pypismtools',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Topic :: Utilities'
      ],
      scripts=['scripts/basemap-plot.py', 'scripts/qgis-colorramp.py',
               'scripts/contour2shp.py', 'scripts/remap3d.py',
               'scripts/flowline-plot.py',
               'scripts/flowline-2d-plot.py', 'scripts/flowline-2d-box-plot.py',
               'scripts/extract_terminus.py', 'scripts/extract_profiles.py',
               'scripts/profile-plot.py', 'scripts/profile-ts-plot.py',
               'scripts/temp2enth.py', 'scripts/im-plot.py',
               'scripts/scalar_within_poly.py', 'scripts/create_greenland_bamber_grid.py',
               'scripts/create_greenland_epsg3413_grid.py',
               'scripts/create_greenland_ext_epsg3413_grid.py',
               'scripts/create_jif_utm22n_grid.py',
               'scripts/extract_sigma_levels.py',
               'scripts/plot_timeseries.py', 'scripts/calculate_marine_gates_area.py', 'scripts/vraster2lineshapefile.py'],
      packages=[PKG_NAME],
      package_dir={PKG_NAME: '.'},
      package_data={PKG_NAME: ['colormaps/*.cpt']},
      )

print(
    '\n*********************************************************************')
print(('Make make sure you have\n %s\nin your search path!'
       % os.path.join(sys.prefix, 'bin')))
print('*********************************************************************')
