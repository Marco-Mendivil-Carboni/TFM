#!/usr/bin/python3

# Imports

from sys import argv

from pathlib import Path

import pandas as pd
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt

from io import StringIO

# Set fixed parameters

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [40.00 * cm, 18.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

mpl.rcParams["figure.constrained_layout.h_pad"] = cm / 4
mpl.rcParams["figure.constrained_layout.w_pad"] = cm / 4

mpl.rcParams["figure.constrained_layout.hspace"] = 0.0
mpl.rcParams["figure.constrained_layout.wspace"] = 0.0

# Set simulation directory

if len(argv) == 2:
    simdir = Path(argv[1])
else:
    print("One command line argument required")
    exit()
if not simdir.exists():
    print("Simulation not found")
    exit()

# Load data into dataframes

anafilepath = simdir / "analysis-fin.dat"
with open(anafilepath) as anafile:
    blocklist = anafile.read().split("\n\n\n")

df_s = pd.read_csv(StringIO(blocklist[0]), sep="\\s+")
df_s.insert(loc=0, column="simobs", value=["dcm", "rg2", "nop", "ncf"])
print(df_s.to_string(index=False))

df_rcd = list()
for i_t in range(3):
    df_rcd.append(pd.read_csv(StringIO(blocklist[i_t + 1]), sep="\\s+"))

df_sd = list()
for i_c in range(6):
    df_sd.append(pd.read_csv(StringIO(blocklist[i_c + 4]), sep="\\s+"))

df_cp = list()
for i_c in range(6):
    df_cp.append(pd.read_csv(StringIO(blocklist[i_c + 10]), sep="\\s+"))

cmfilepath = simdir / "contact-map.bin"
with open(cmfilepath) as cmfile:
    cm1Ddata = np.memmap(cmfile, dtype=np.dtype("f8"), mode="r")

# Make analysis plots

fig = plt.figure()
grid = fig.add_gridspec(2, 4)
ax0 = fig.add_subplot(grid[0, 0:2])
ax1 = fig.add_subplot(grid[1, 0])
ax2 = fig.add_subplot(grid[1, 1])
ax3 = fig.add_subplot(grid[0:2, 2:4])

fig.canvas.manager.set_window_title(str(simdir) + " analysis")

colorlist = ["#d81e2c", "#a31cc5", "#194bb2"]

ax0.set_xlabel("r_b")
ax0.set_ylabel("rcd")
ax0.autoscale(tight=True)
for i_t in range(3):
    x = df_rcd[i_t]["r_b"]
    y = df_rcd[i_t]["avg"]
    y_min = df_rcd[i_t]["avg"] - df_rcd[i_t]["sem"]
    y_max = df_rcd[i_t]["avg"] + df_rcd[i_t]["sem"]
    ax0.step(x, y, color=colorlist[i_t])
    ax0.fill_between(
        x, y_min, y_max, step="pre", color=colorlist[i_t], linewidth=0.0, alpha=0.50
    )

colorlist = ["#221ab9", "#194bb2", "#1880ac", "#17a69b", "#169f62", "#15992c"]

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("s")
ax1.set_ylabel("sd")
ax1.autoscale(tight=True)
for i_c in range(6):
    if not df_sd[i_c].empty:
        x = df_sd[i_c]["s"]
        y = df_sd[i_c]["avg"]
        y_min = df_sd[i_c]["avg"] - df_sd[i_c]["sem"]
        y_max = df_sd[i_c]["avg"] + df_sd[i_c]["sem"]
        ax1.plot(x, y, color=colorlist[i_c])
        ax1.fill_between(
            x, y_min, y_max, color=colorlist[i_c], linewidth=0.0, alpha=0.50
        )

ctcfactor = 1000

ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("s")
ax2.set_ylabel("cp")
ax2.autoscale(tight=True)
for i_c in range(6):
    if not df_cp[i_c].empty:
        df_cp[i_c] = df_cp[i_c].loc[df_cp[i_c]["avg"] != 0]
        x = df_cp[i_c]["s"]
        y = df_cp[i_c]["avg"] / ctcfactor
        e = df_cp[i_c]["sem"] / ctcfactor
        ax2.scatter(x, y, color=colorlist[i_c])
        ax2.errorbar(x, y, yerr=e, color=colorlist[i_c], linestyle="None", alpha=0.50)

ax3.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax3.autoscale(tight=True)
px_side = int(np.sqrt(2 * cm1Ddata.size + 1 / 4))
cm2Ddata = np.zeros((px_side, px_side))
cm2Ddata[np.tril_indices(px_side)] = cm1Ddata / ctcfactor
cm2Ddata = cm2Ddata + np.tril(cm2Ddata, -1).T
map = ax3.imshow(
    cm2Ddata,
    norm=mpl.colors.LogNorm(vmin=np.min(cm2Ddata[np.nonzero(cm2Ddata)]), clip=True),
    cmap="OrRd",
)

fig.colorbar(map, ax=ax3, aspect=64, pad=1 / 64)

# View analysis

plt.show()
