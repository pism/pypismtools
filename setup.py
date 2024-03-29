from distutils.core import setup
import os
import sys

PKG_NAME = "pypismtools"
print(("\nInstalling %s" % PKG_NAME))
print("----------------------------")

packages = ("numpy", "osgeo", "netCDF4", "pyproj")
print("\nChecking dependencies:")
not_installed = []
for package in packages:
    try:
        __import__(package)
        print(("  - % s package is installed" % package))
    except ImportError:
        print(("  - % s package NOT installed" % package))
        not_installed.append(package)
if not_installed:
    print("Installation of the following packages is optional but recommended:")
    for package in not_installed:
        if package == "osgeo":
            print(" - GDAL python bindings")
        else:
            print(("  - %s" % package))
    print("Exiting")
    import sys

    sys.exit()

setup(
    name=PKG_NAME,
    version="0.5",
    description="Python tools to evaluate PISM results",
    author="Andy Aschwanden",
    author_email="aaschwanden@alaska.edu",
    url="https://github.com/pism/pypismtools",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Utilities",
    ],
    scripts=[
        "scripts/basemap_plot.py",
        "scripts/contour2shp.py",
        "scripts/create_greenland_bamber_grid.py",
        "scripts/create_greenland_epsg3413_grid.py",
        "scripts/create_greenland_ext_epsg3413_grid.py",
        "scripts/create_jif_utm22n_grid.py",
        "scripts/dissolve_by_attribute.py",
        "scripts/extract_sigma_levels.py",
        "scripts/extract_interface.py",
        "scripts/extract_contours.py",
        "scripts/extract_profiles.py",
        "scripts/qgis_colorramp.py",
        "scripts/scalar_within_poly.py",
        "scripts/remap3d.py",
        "scripts/vraster2lineshapefile.py",
        "scripts/point2line.py",
    ],
    packages=[PKG_NAME],
    package_dir={PKG_NAME: "."},
    package_data={PKG_NAME: ["colormaps/*.cpt"]},
)

print("\n*********************************************************************")
print(("Make make sure you have\n %s\nin your search path!" % os.path.join(sys.prefix, "bin")))
print("*********************************************************************")
