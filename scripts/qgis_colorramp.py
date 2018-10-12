#!/usr/bin/env python
import numpy as np
import matplotlib
import pylab as plt
import matplotlib as mpl
from argparse import ArgumentParser

try:
    from pypismtools import gmtColormap
except:
    from pypismtools.pypismtools import gmtColormap


def cmap_map(function, cmap):
    """
    Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous points
    in a colormap.

    Adapted from http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # First get the list of points where the segments start or end
    for key in ("red", "green", "blue"):
        step_dict[key] = [x[0] for x in cdict[key]]
    step_list = sum(list(step_dict.values()), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT

    def reduced_cmap(step):
        return np.array(cmap(step)[0:3])

    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(("red", "green", "blue")):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j, i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = sorted([x + (x[1],) for x in list(this_cdict.items())])
        cdict[key] = colorvector

    return mpl.colors.LinearSegmentedColormap("colormap", cdict, 1024)


# Set up the option parser
parser = ArgumentParser()
parser.description = """
A script to convert a GMT (*.cpt) colormap or matplotlib colormap into QGIS-readable color ramp.
"""
parser.add_argument("FILE", nargs=1)
parser.add_argument("--tick_format", dest="tick_format", help="Overwrite ormat of the tick marks.", default=None)
parser.add_argument("--font_size", dest="font_size", help="Font size", default=12)
parser.add_argument(
    "--type",
    dest="colorbar_type",
    choices=[
        "linear",
        "log_speed_10_1500",
        "log_speed_10_3000",
        "log_speed_j",
        "log_speed_2",
        "log_speed_3",
        "log_speed_4",
        "gris_bath_topo",
        "gris_bath_topo_2",
        "gris_topo",
        "log_speed_m_day",
    ],
    help="Type of colorbar",
    default="linear",
)
parser.add_argument("--ticks", dest="fticks", nargs="*", type=float, help="tick marks", default=None)
parser.add_argument(
    "--colorbar_extend",
    dest="cb_extend",
    choices=["neither", "both", "min", "max"],
    help="""Extend of colorbar. Default='both'.""",
    default="both",
)
parser.add_argument("--colorbar_label", dest="colorbar_label", help="""Label for colorbar.""", default=None)
parser.add_argument(
    "--vmin",
    dest="vmin",
    type=float,
    help="""
                    Vmin""",
    default=1,
)
parser.add_argument(
    "--vmax",
    dest="vmax",
    type=float,
    help="""
                    Vmax)""",
    default=3000,
)
parser.add_argument(
    "--extend",
    dest="extend",
    nargs=2,
    type=float,
    help="""
                  appends color ramp by repeating first and last color for value""",
    default=None,
)
parser.add_argument(
    "--N",
    dest="N",
    type=int,
    help="""
                  a * logspace(vmin, vmax, N""",
    default=1022,
)
parser.add_argument("-r", "--reverse", dest="reverse", action="store_true", help="reverse color scale", default=False)
parser.add_argument(
    "--orientation",
    dest="orientation",
    choices=["horizontal", "vertical"],
    help="Orientation, default = 'horizontal",
    default="horizontal",
)

options = parser.parse_args()
args = options.FILE
colorbar_type = options.colorbar_type
font_size = options.font_size
tick_format = options.tick_format
fticks = options.fticks
extend = options.extend
N = options.N
vmin = options.vmin
vmax = options.vmax
reverse = options.reverse
colorbar_label = options.colorbar_label
cb_extend = options.cb_extend
# experimental
log_color = False
orientation = options.orientation


# read in CPT colormap
cmap_file = args[0]
try:
    cmap = getattr(plt.cm, cmap_file)
    prefix = cmap_file
except:
    # import and convert colormap
    cdict = gmtColormap(cmap_file, log_color=log_color, reverse=reverse)
    prefix = ".".join(cmap_file.split(".")[0:-1])
    suffix = cmap_file.split(".")[-1]
    cmap = mpl.colors.LinearSegmentedColormap("my_colormap", cdict, N)


class nlcmap(object):
    def __init__(self, cmap, levels):
        self.cmap = cmap
        self.levels = np.asarray(levels, dtype="float64")
        self._x = self.levels
        self.levmax = self.levels.max()
        self.transformed_levels = np.linspace(0.0, self.levmax, len(self.levels))

    def __call__(self, xi, alpha=1.0, **kw):
        yi = np.interp(xi, self._x, self.transformed_levels)
        return self.cmap(yi / self.levmax, alpha)


levels = [0, 10, 100, 250, 750, 3000]
levels.sort()


if colorbar_type in ("linear"):
    data_values = np.linspace(vmin, vmax, N)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb_extend = cb_extend
    format = "%2.0f"
    ticks = None
elif colorbar_type in ("log_speed_10_3000"):
    data_values = np.linspace(vmin, vmax, N)
    norm = mpl.colors.LogNorm(vmin=10, vmax=3000)
    cb_extend = cb_extend
    colorbar_label = "m yr$^{-1}$"
    format = "%2.0f"
    ticks = [0, 10, 100, 300, 1000, 3000]
elif colorbar_type in ("log_speed_10_1500"):
    data_values = np.linspace(vmin, vmax, N)
    norm = mpl.colors.LogNorm(vmin=10, vmax=1500)
    cb_extend = cb_extend
    if colorbar_label is None:
        colorbar_label = "m yr$^{-1}$"
    format = "%2.0f"
    ticks = [0, 10, 100, 300, 1000, 3000]
elif colorbar_type in ("gris_bath_topo"):
    vmin = -800
    vmax = 3000
    data_values = np.linspace(vmin, vmax, N)
    N = len(data_values)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb_extend = "both"
    format = "%i"
    ticks = [vmin, 0, 1000, 2000, vmax]
elif colorbar_type in ("gris_bath_topo_2"):
    vmin = -1500
    vmax = 2000
    data_values = np.linspace(vmin, vmax, N)
    N = len(data_values)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb_extend = "both"
    if colorbar_label is None:
        colorbar_label = "m a.s.l."
    format = "%i"
    ticks = [vmin, -750, 0, 1000, 2000, vmax]
elif colorbar_type in ("gris_topo"):
    vmin = 0
    vmax = 2000
    data_values = np.linspace(vmin, vmax, N)
    N = len(data_values)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb_extend = "both"
    if colorbar_label is None:
        colorbar_label = "m a.s.l."
    format = "%i"
    ticks = [vmin, 0, 1000, 2000, vmax]
    cmap.set_under("#2171b5")
elif colorbar_type in ("log_speed_j", "log_speed_3"):
    data_values = np.logspace(-1, 3, N)[0:889]
    data_values[-1] = vmax
    N = len(data_values)
    norm = mpl.colors.LogNorm(vmin=1, vmax=vmax)
    cb_extend = "both"
    format = "%i"
    ticks = [1, 3, 10, 30, 100, 300, 1000, 3000]
elif colorbar_type in ("log_speed_2"):
    data_values = np.logspace(-1, 2, N)[0:889]
    data_values[-1] = vmax
    N = len(data_values)
    norm = mpl.colors.LogNorm(vmin=1, vmax=vmax)
    cb_extend = "both"
    format = "%i"
    ticks = [1, 3, 10, 30, 100, 300]
elif colorbar_type in ("log_speed_4"):
    data_values = np.logspace(-1, 4, N)[0:889]
    data_values[-1] = vmax
    N = len(data_values)
    norm = mpl.colors.LogNorm(vmin=1, vmax=vmax)
    cb_extend = "both"
    format = "%i"
    ticks = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
elif colorbar_type in ("log_speed_m_day"):
    data_values = np.logspace(vmin, vmax, N)[0:889]
    data_values[-1] = vmax
    N = len(data_values)
    norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    cb_extend = "both"
    format = "%2.2f"
    ticks = [0.01, 0.1, 0.5, 1, 5, 10]
else:
    pass

if tick_format is not None:
    format = tick_format


# you could apply a function to the colormap, e.g. to desaturate the colormap:
# cmap = cmap_map(lambda x: x/2+0.5, cmap)
matplotlib.rc("font", **{"sans-serif": ["Helvetica"]})  # "size": fontsize}
matplotlib.rc("font", **{"size": font_size})
# create the colorbar
fig = plt.figure()
if orientation == "horizontal":
    ax1 = fig.add_axes([0.0, 0.5, 0.5, 0.03])
else:
    ax1 = fig.add_axes([0.05, 0.05, 0.03, 0.65])
if fticks is not None:
    ticks = fticks

cb1 = mpl.colorbar.ColorbarBase(
    ax1,
    cmap=cmap,
    norm=norm,
    ticks=ticks,
    format=format,
    extend=cb_extend,
    spacing="proportional",
    orientation=orientation,
)

if colorbar_label:
    cb1.set_label(colorbar_label)


# save high-res colorbar as png
for format in ["png"]:
    prefix = prefix + "_" + orientation
    out_file = ".".join([prefix, format])
    print(("  writing colorbar %s ..." % out_file))
    fig.savefig(out_file, bbox_inches="tight", dpi=1200, transparent=True)

# convert to RGBA array
rgba = cb1.to_rgba(data_values, alpha=None)
# QGIS wants 0..255
rgba *= 255

# create an output array combining data values and rgb values
if extend:
    qgis_array = np.zeros((N + 2, 5))
    for k in range(0, N):
        qgis_array[k + 1, 0] = data_values[k]
        qgis_array[k + 1, 1:4] = rgba[k, 0:3]
        qgis_array[k + 1, 4] = 255
    # repeat first color
    qgis_array[0, 0] = extend[0]
    qgis_array[0, 1:4] = rgba[0, 0:3]
    qgis_array[0, 4] = 255
    # repeat last color
    qgis_array[-1, 0] = extend[1]
    qgis_array[-1, 1:4] = rgba[-1, 0:3]
    qgis_array[-1, 4] = 255
else:
    qgis_array = np.zeros((N, 5))
    for k in range(N):
        qgis_array[k, 0] = data_values[k]
        qgis_array[k, 1:4] = rgba[k, 0:3]
        qgis_array[k, 4] = 255

# save as ascii file
out_file = ".".join([prefix, "txt"])
print(("  writing colorramp %s ..." % out_file))
np.savetxt(out_file, qgis_array, delimiter=",", fmt=["%10.5f", "%i", "%i", "%i", "%i,"])
