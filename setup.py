from distutils.core import setup

PKG_NAME = 'PyPISMTools'
print('\nInstalling %s' % PKG_NAME)
print('----------------------------')

packages = ('numpy', 'osgeo', 'netCDF4', 'udunits2', 'pyproj',
            'mpl_toolkits.basemap', 'mpl_toolkits.axes_grid1')
print('\nChecking dependencies:')
not_installed = []
for package in packages:
    try:
        __import__(package)
        print('  - % s package is installed' % package)
    except ImportError:
        print('  - % s package not installed' % package)
        not_installed.append(package)
if not_installed:
    print('Please install the following packages:')
    for package in not_installed:
        print('  - %s' % package)
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
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Utilities'
        ],
      scripts=['scripts/basemap-plot.py'],
      packages=[PKG_NAME],
      package_dir={PKG_NAME: '.'},
      package_data={PKG_NAME: ['colormaps/*.cpt']},
      )
