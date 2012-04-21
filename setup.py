from setuptools import setup

PKG_NAME = 'PyPISMTools'

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
      install_requires=['pyproj','udunits2', 'gdal', 'netCDF4'],
      scripts=['scripts/basemap-plot.py'],
      packages=[PKG_NAME],
      package_dir={PKG_NAME: '.'},
      package_data={PKG_NAME: ['colormaps/*.cpt']},
      )
