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

mpl.rcParams["figure.constrained_layout.h_pad"] = cm / 4
mpl.rcParams["figure.constrained_layout.w_pad"] = cm / 4

mpl.rcParams["figure.constrained_layout.hspace"] = 0.0
mpl.rcParams["figure.constrained_layout.wspace"] = 0.0

# Set auxiliary variables

colorlist_3 = ["#d81e2c", "#a31cc5", "#194bb2"]
colorlist_6 = ["#221ab9", "#194bb2", "#1880ac", "#17a69b", "#169f62", "#15992c"]
ctcfactor = 1000
px_sz = 4

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
df_s.loc[1] = np.sqrt(df_s.loc[1])
df_s.insert(loc=0, column="simobs", value=["dcm", "rog", "nop", "ncf"])
df_s.set_index("simobs", inplace=True)
print(df_s.to_string())

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

for i_t in range(3):
    x = df_rcd[i_t]["r_b"]
    y = df_rcd[i_t]["avg"]
    e = df_rcd[i_t]["sem"]
    ax0.step(x, y, color=colorlist_3[i_t])
    ax0.fill_between(
        x, y - e, y + e, step="pre", color=colorlist_3[i_t], linewidth=0.0, alpha=0.50
    )
ax0.autoscale(tight=True)
ax0.set_ylim(bottom=0.0)

ax1.set_xscale("log")
ax1.set_yscale("log")
for i_c in range(6):
    if not df_sd[i_c].empty:
        x = df_sd[i_c]["s"]
        y = df_sd[i_c]["avg"]
        e = df_sd[i_c]["sem"]
        ax1.scatter(x, y, s=8, color=colorlist_6[i_c])
        ax1.errorbar(x, y, yerr=e, color=colorlist_6[i_c], linestyle="None", alpha=0.50)
ax1.autoscale(tight=True)
ax1.set_ylim(bottom=1.0)

ax2.set_xscale("log")
ax2.set_yscale("log")
for i_c in range(6):
    if not df_cp[i_c].empty:
        df_cp[i_c] = df_cp[i_c].loc[df_cp[i_c]["avg"] != 0]
        x = df_cp[i_c]["s"]
        y = df_cp[i_c]["avg"] / ctcfactor
        e = df_cp[i_c]["sem"] / ctcfactor
        ax2.scatter(x, y, s=8, color=colorlist_6[i_c])
        ax2.errorbar(x, y, yerr=e, color=colorlist_6[i_c], linestyle="None", alpha=0.50)
ax2.autoscale(tight=True)

ax3.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
px_side = int(np.sqrt(2 * cm1Ddata.size + 1 / 4))
cm2Ddata = np.zeros((px_side, px_side))
cm2Ddata[np.tril_indices(px_side)] = cm1Ddata / ctcfactor
cm2Ddata = cm2Ddata + np.tril(cm2Ddata, -1).T
len_side = px_sz * px_side
map = ax3.imshow(
    cm2Ddata,
    cmap="OrRd",
    norm=mpl.colors.LogNorm(vmax=1.0),
    extent=[0, len_side, len_side, 0],
)
map.get_cmap().set_bad("white")
ax3.autoscale(tight=True)
cbar = fig.colorbar(map, ax=ax3, aspect=64, pad=1 / 64)

# View analysis

plt.show()
