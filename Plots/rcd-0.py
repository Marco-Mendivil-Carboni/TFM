# Imports

from pathlib import Path

import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt

from io import StringIO

# Set fixed parameters

mpl.use("pdf")

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [14.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

# Load data into a dataframe

simdir = Path("Simulations/40.91-18240-07284-00.00-00.00")

anafilepath = simdir / "analysis-fin.dat"
with open(anafilepath) as anafile:
    blocklist = anafile.read().split("\n\n\n")

df_rcd = list()
for i_t in range(3):
    df_rcd.append(pd.read_csv(StringIO(blocklist[i_t + 1]), sep="\\s+"))

# Make rcd-0 plot

plotsdir = Path("Plots")

fig, ax = plt.subplots()

colorlist = ["#d81e2c", "#a31cc5", "#194bb2"]

ax.set_xlabel("$r$")
ax.set_ylabel("$\\rho(r)$")
ax.autoscale(tight=True)
for i_t in range(3):
    x = df_rcd[i_t]["r_b"]
    y = df_rcd[i_t]["avg"]
    y_min = df_rcd[i_t]["avg"] - df_rcd[i_t]["sem"]
    y_max = df_rcd[i_t]["avg"] + df_rcd[i_t]["sem"]
    ax.step(x, y, color=colorlist[i_t])
    ax.fill_between(
        x, y_min, y_max, step="pre", color=colorlist[i_t], linewidth=0.0, alpha=0.50
    )

fig.savefig(plotsdir / "rcd-0.pdf")
