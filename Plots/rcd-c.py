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

# Set auxiliary variables

lenfactor = 1000 / 33

# Load data into a dataframe

simdir = Path("Simulations/30386-0.200-0.400-0.000-0.000")

anafilepath = simdir / "analysis-fin.dat"
with open(anafilepath) as anafile:
    blocklist = anafile.read().split("\n\n\n")

df_rcd = list()
for i_t in range(3):
    df_rcd.append(pd.read_csv(StringIO(blocklist[i_t + 1]), sep="\\s+"))

# Make rcd-c plot

plotsdir = Path("Plots")

fig, ax = plt.subplots()

colorlist = ["#d81e2c", "#a31cc5", "#194bb2"]
labellist = ["LADh", "LNDe", "total"]
ax.set_xlabel("$r$ ($\\mu$m)")
ax.set_ylabel("$\\rho(r)$")
for i_t in range(3):
    x = df_rcd[i_t]["r_b"] / lenfactor
    y = df_rcd[i_t]["avg"]
    y_min = df_rcd[i_t]["avg"] - df_rcd[i_t]["sem"]
    y_max = df_rcd[i_t]["avg"] + df_rcd[i_t]["sem"]
    ax.step(x, y, color=colorlist[i_t], label=labellist[i_t])
    ax.fill_between(
        x, y_min, y_max, step="pre", color=colorlist[i_t], linewidth=0.0, alpha=0.50
    )
ax.autoscale(tight=True)
ax.set_ylim(bottom=0.0)
ax.legend(loc="upper left")

fig.savefig(plotsdir / "rcd-c.pdf")
