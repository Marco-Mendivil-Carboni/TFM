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

mpl.rcParams["legend.frameon"] = False

# Load data into dataframes

simrootdir = Path("Simulations")
datafilepath = simrootdir / "performance.dat"

with open(datafilepath) as datafile:
    blocklist = datafile.read().split("\n\n\n")

df_1 = pd.read_csv(
    StringIO(blocklist[0]), delim_whitespace=True, comment="#", header=None
)
df_1.columns = ["N", "t_e"]

df_2 = pd.read_csv(
    StringIO(blocklist[1]), delim_whitespace=True, comment="#", header=None
)
df_2.columns = ["N", "t_e"]

df_3 = pd.read_csv(
    StringIO(blocklist[2]), delim_whitespace=True, comment="#", header=None
)
df_3.columns = ["N", "t_e"]

# Make performance plot

plotsrootdir = Path("Plots")

fig, ax = plt.subplots()

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("$N$")
ax.set_ylabel("$t_e$ (ms)")

ax.plot(df_1["N"], df_1["t_e"], marker="o", color="#d81e2c", label="GeForce 920M")
ax.plot(df_2["N"], df_2["t_e"], marker="o", color="#a31cc5", label="GeForce RTX 3050")
ax.plot(df_3["N"], df_3["t_e"], marker="o", color="#194bb2", label="RTX A4000")

ax.legend(loc="upper left")

fig.savefig(plotsrootdir / "performance.pdf")
