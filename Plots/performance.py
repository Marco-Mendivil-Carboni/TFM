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
mpl.rcParams["figure.figsize"] = [12.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

# Load data into dataframes

simrootdir = Path("Simulations")
datafilepath = simrootdir / "performance.txt"

with open(datafilepath) as datafile:
    blocklist = datafile.read().split("\n\n\n")

df_1 = pd.read_csv(StringIO(blocklist[0]), sep="\\s+", header=None)
df_1.columns = ["N", "t_e"]

df_2 = pd.read_csv(StringIO(blocklist[1]), sep="\\s+", header=None)
df_2.columns = ["N", "t_e"]

# Make performance plot

plotsrootdir = Path("Plots")

fig, ax = plt.subplots()

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("$N$")
ax.set_ylabel("$t_e$ (ms)")

ax.plot(df_1["N"], df_1["t_e"], marker="o", color="#d81e2c", label="GeForce RTX 3050")
ax.plot(df_2["N"], df_2["t_e"], marker="o", color="#a31cc5", label="RTX A4000")

ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

ax.legend(loc="upper left")

fig.savefig(plotsrootdir / "performance.pdf")
