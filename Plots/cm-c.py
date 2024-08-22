# Imports

from pathlib import Path

import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt

# Set fixed parameters

mpl.use("pdf")

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [10.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

# Set auxiliary variables

ctcfactor = 1000

# Load data into a dataframe

simdir = Path("Simulations/30386-0.200-0.400-0.000-0.000")

cmfilepath = simdir / "contact-map.bin"
with open(cmfilepath) as cmfile:
    cm1Ddata = np.memmap(cmfile, dtype=np.dtype("f8"), mode="r")

# Make rcd-c plot

plotsdir = Path("Plots")

fig, ax = plt.subplots()

px_side = int(np.sqrt(2 * cm1Ddata.size + 1 / 4))
cm2Ddata = np.zeros((px_side, px_side))
cm2Ddata[np.tril_indices(px_side)] = cm1Ddata / ctcfactor
cm2Ddata = cm2Ddata + np.tril(cm2Ddata, -1).T
map = ax.imshow(cm2Ddata, norm=mpl.colors.LogNorm(vmax=1.0), cmap="OrRd")
map.get_cmap().set_bad("white")
ax.autoscale(tight=True)
cbar = fig.colorbar(map, ax=ax, aspect=64, pad=1 / 64)
cbar.ax.set_ylabel("$P(i,j)$")

fig.savefig(plotsdir / "cm-c.pdf", dpi=1600)
