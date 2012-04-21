from distutils.core import setup

PKG_NAME = 'PyPISMTools'

setup(name=PKG_NAME,
      version='0.12',
      description='Python tools to evaluate PISM results',
      author='Andy Aschwanden',
      author_email='aaschwanden@alaska.edu',
      url='https://github.com/pism/PyPISMTools',
      scripts=['scripts/basemap-plot.py'],
      packages=[PKG_NAME],
      package_dir={PKG_NAME: '.'},
      package_data={PKG_NAME: ['colormaps/*.cpt']},
      )
