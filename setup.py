from distutils.core import setup
import os
import sys

PKG_NAME = 'PyPISMTools'
print(('\nInstalling %s' % PKG_NAME))
print('----------------------------')

packages = ('numpy', 'osgeo', 'netCDF4', 'udunits2', 'pyproj',
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
    print('Please install the following packages:')
    for package in not_installed:
        if package is 'osgeo':
            print(' - GDAL python bindings')
        elif package is 'udunits2':
            print(' - py_udunits2 from https://github.com/ckhroulev/py_udunits2')
        else:
            print(('  - %s' % package))
    print('Exiting')
    import sys
    sys.exit()

setup(name=PKG_NAME,
      version='0.12',
      description='Python tools to evaluate PISM results',
      author='Andy Aschwanden',
      author_email='aaschwanden@alaska.edu',
      url='https://github.com/pism/PyPISMTools',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Utilities'
        ],
      scripts=['scripts/basemap-plot.py','scripts/qgis-colorramp.py',
               'scripts/contour2shp.py', 'remap3d.py],
      packages=[PKG_NAME],
      package_dir={PKG_NAME: '.'},
      package_data={PKG_NAME: ['colormaps/*.cpt']},
      )

print('\n*********************************************************************')
print(('Make make sure you have\n %s\nin your search path!'
      % os.path.join(sys.prefix, 'bin')))
print('*********************************************************************')
