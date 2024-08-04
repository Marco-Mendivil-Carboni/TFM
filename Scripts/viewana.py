#!/usr/bin/python3

# Imports

from sys import argv

from pathlib import Path

import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt

from io import StringIO

# Set fixed parameters

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [32.00 * cm, 16.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

mpl.rcParams["toolbar"] = "None"
mpl.rcParams["figure.constrained_layout.h_pad"] = 0.25
mpl.rcParams["figure.constrained_layout.w_pad"] = 0.25
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

parfilepath = simdir / "adjustable-parameters.dat"
df_par = pd.read_csv(parfilepath, sep="\\s+", header=None)

anafilepath = simdir / "analysis-fin.dat"
with open(anafilepath) as anafile:
    blocklist = anafile.read().split("\n\n\n")

df_s = pd.read_csv(StringIO(blocklist[0]), sep="\\s+")
df_s.insert(loc=0, column="simobs", value=["dcm", "rg2", "nop", "ncf"])

df_rcd = list()
for i_t in range(3):
    df_rcd.append(pd.read_csv(StringIO(blocklist[i_t + 1]), sep="\\s+"))

df_msd = list()
for i_c in range(6):
    df_msd.append(pd.read_csv(StringIO(blocklist[i_c + 4]), sep="\\s+"))

# Make analysis plots

fig, ax = plt.subplots(2, 2, height_ratios=[0.25, 1])

fig.canvas.manager.set_window_title(str(simdir) + " analysis")

ax[0, 0].axis("off")
ax[0, 0].table(
    cellText=df_par.values, loc="center", cellLoc="left", edges="horizontal"
).scale(1, 1.5)

ax[0, 1].axis("off")
ax[0, 1].table(
    cellText=df_s.values,
    colLabels=df_s.columns,
    loc="center",
    colLoc="right",
    edges="horizontal",
).scale(1, 1.5)

colorlist = ["#d81e2c", "#a31cc5", "#194bb2"]
ax[1, 0].set_xlabel("r_b")
ax[1, 0].set_ylabel("rcd")
ax[1, 0].autoscale(tight=True)
for i_t in range(3):
    x = df_rcd[i_t]["r_b"]
    y = df_rcd[i_t]["avg"]
    y_min = df_rcd[i_t]["avg"] - df_rcd[i_t]["sem"]
    y_max = df_rcd[i_t]["avg"] + df_rcd[i_t]["sem"]
    ax[1, 0].step(x, y, color=colorlist[i_t])
    ax[1, 0].fill_between(
        x, y_min, y_max, step="pre", color=colorlist[i_t], linewidth=0.0, alpha=0.50
    )

colorlist = ["#1875ad", "#189da8", "#17a384", "#169e57", "#15992c", "#259415"]
ax[1, 1].set_xscale("log")
ax[1, 1].set_yscale("log")
ax[1, 1].set_xlabel("s")
ax[1, 1].set_ylabel("msd")
ax[1, 1].autoscale(tight=True)
for i_c in range(6):
    if not df_msd[i_c].empty:
        x = df_msd[i_c]["s"]
        y = df_msd[i_c]["avg"]
        y_min = df_msd[i_c]["avg"] - df_msd[i_c]["sem"]
        y_max = df_msd[i_c]["avg"] + df_msd[i_c]["sem"]
        ax[1, 1].plot(x, y, color=colorlist[i_c])
        ax[1, 1].fill_between(
            x, y_min, y_max, color=colorlist[i_c], linewidth=0.0, alpha=0.50
        )

# View analysis

plt.show()
