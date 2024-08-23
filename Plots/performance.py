# Imports

from pathlib import Path

import pandas as pd
from scipy.optimize import curve_fit

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

# Define fitting function


def scaling_law(x: float, a: float, p: float, c: float) -> float:
    return a * x**p + c


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

df = [df_1, df_2]

color = ["#d81e2c", "#a31cc5"]
label = ["GeForce RTX 3050", "RTX A4000"]

for i in range(2):
    ax.plot(df[i]["N"], df[i]["t_e"], marker="o", color=color[i], label=label[i])

for i in range(2):
    popt, _ = curve_fit(scaling_law, df[i]["N"], df[i]["t_e"], p0=[0.0001, 1.0, 0.1])
    ax.plot(
        df[i]["N"],
        scaling_law(df[i]["N"], *popt),
        color=color[i],
        alpha=0.5,
        linestyle="dashed",
        label="${:.1f}\\cdot 10^{{-5}}\\cdot N^{{{:.1f}}}+{:.1f}$".format(
            popt[0] * 1e5, popt[1], popt[2]
        ),
    )

ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

ax.legend(loc="upper left")

fig.savefig(plotsrootdir / "performance.pdf")
