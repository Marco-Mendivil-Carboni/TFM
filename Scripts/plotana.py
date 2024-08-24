#!/usr/bin/python3

# Imports

from sys import argv

from pathlib import Path

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import matplotlib as mpl
from matplotlib import pyplot as plt

from io import StringIO

# Set fixed parameters

mpl.use("pdf")

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"

cm = 1 / 2.54

mpl.rcParams["figure.constrained_layout.use"] = True

# Set auxiliary variables

lenfactor = 1000 / 33
bp_part = 33 * 200
ctcfactor = 1000
px_sz = 4
colorlist_rcd = ["#d81e2c", "#a31cc5", "#194bb2"]
colorlist_sd_cp = ["#221ab9", "#194bb2", "#1880ac", "#17a69b", "#169f62", "#15992c"]
labellist_rcd = ["LADh", "LNDe", "total"]
labellist_sd_cp = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"]
fitcolorlist_sd_cp = ["#df591f", "#d81e2c"]
fitrangelist_sd = [[1, 4], [16, 256]]
fitrangelist_cp = [[2, 6], [16, 256]]
fitpopts_sd = [list() for _ in range(len(fitrangelist_sd))]
fitpopts_cp = [list() for _ in range(len(fitrangelist_cp))]
fitxrangelist_sd = [list() for _ in range(len(fitrangelist_sd))]
fitxrangelist_cp = [list() for _ in range(len(fitrangelist_cp))]
for i_f in range(len(fitrangelist_sd)):
    fitxrangelist_sd[i_f] = [lim / (1e6 / bp_part) for lim in fitrangelist_sd[i_f]]
for i_f in range(len(fitrangelist_cp)):
    fitxrangelist_cp[i_f] = [lim / (1e6 / bp_part) for lim in fitrangelist_cp[i_f]]

# Define fitting function


def scaling_law(x: float, a: float, p: float) -> float:
    return a * x**p


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
df_s.loc[0:1] = df_s.loc[0:1] / lenfactor
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

plotsdir = Path("Plots") / simdir
plotsdir.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(14.00 * cm, 8.00 * cm))
ax.set_xlabel("$r$ ($\\mu$m)")
ax.set_ylabel("$\\rho(r)$")
for i_t in range(3):
    x = df_rcd[i_t]["r_b"] / lenfactor
    y = df_rcd[i_t]["avg"]
    e = df_rcd[i_t]["sem"]
    ax.step(x, y, color=colorlist_rcd[i_t], label=labellist_rcd[i_t])
    ax.fill_between(
        x, y - e, y + e, step="pre", color=colorlist_rcd[i_t], linewidth=0.0, alpha=0.50
    )
ax.autoscale(tight=True)
ax.set_ylim(bottom=0.0)
ax.legend(loc="upper left")
fig.savefig(plotsdir / "rcd.pdf")

fig, ax = plt.subplots(figsize=(12.00 * cm, 8.00 * cm))
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$s$ (Mb)")
ax.set_ylabel("$d(s)$ ($\\mu$m)")
for i_c in range(6):
    if not df_sd[i_c].empty:
        x = df_sd[i_c]["s"] / (1e6 / bp_part)
        y = df_sd[i_c]["avg"] / lenfactor
        e = df_sd[i_c]["sem"] / lenfactor
        ax.scatter(x, y, s=8, color=colorlist_sd_cp[i_c], label=labellist_sd_cp[i_c])
        ax.errorbar(
            x, y, yerr=e, color=colorlist_sd_cp[i_c], linestyle="None", alpha=0.50
        )
        if i_c != 5:
            for i_f in range(len(fitrangelist_sd)):
                inxrange = x.between(*fitxrangelist_sd[i_f])
                popt, _ = curve_fit(
                    scaling_law, x.loc[inxrange], y.loc[inxrange], p0=[1.0, 1.0]
                )
                fitpopts_sd[i_f].append(popt)
for i_f in range(len(fitrangelist_sd)):
    x_f = np.linspace(*fitxrangelist_sd[i_f], 2)
    popt = np.mean(fitpopts_sd[i_f], axis=0)
    ax.plot(
        x_f,
        scaling_law(x_f, *popt),
        color=fitcolorlist_sd_cp[i_f],
        alpha=0.50,
        linestyle="dashed",
        label="$\\propto s^{{{:.1f}}}$".format(popt[1]),
    )
ax.autoscale(tight=True)
ax.set_ylim(bottom=1.0 / lenfactor)
ax.legend(loc="lower right")
fig.savefig(plotsdir / "sd.pdf")

fig, ax = plt.subplots(figsize=(12.00 * cm, 8.00 * cm))
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$s$ (Mb)")
ax.set_ylabel("$P(s)$")
for i_c in range(6):
    if not df_cp[i_c].empty:
        df_cp[i_c] = df_cp[i_c].loc[df_cp[i_c]["avg"] != 0]
        x = df_cp[i_c]["s"] / (1e6 / bp_part)
        y = df_cp[i_c]["avg"] / ctcfactor
        e = df_cp[i_c]["sem"] / ctcfactor
        ax.scatter(x, y, s=8, color=colorlist_sd_cp[i_c], label=labellist_sd_cp[i_c])
        ax.errorbar(
            x, y, yerr=e, color=colorlist_sd_cp[i_c], linestyle="None", alpha=0.5
        )
        if i_c != 5:
            for i_f in range(len(fitrangelist_cp)):
                inxrange = x.between(*fitxrangelist_cp[i_f])
                popt, _ = curve_fit(
                    scaling_law, x.loc[inxrange], y.loc[inxrange], p0=[1e-4, -1.0]
                )
                fitpopts_cp[i_f].append(popt)
for i_f in range(len(fitrangelist_cp)):
    x_f = np.linspace(*fitxrangelist_cp[i_f], 2)
    popt = np.mean(fitpopts_cp[i_f], axis=0)
    ax.plot(
        x_f,
        scaling_law(x_f, *popt),
        color=fitcolorlist_sd_cp[i_f],
        alpha=0.50,
        linestyle="dashed",
        label="$\\propto s^{{{:.1f}}}$".format(popt[1]),
    )
ax.autoscale(tight=True)
ax.legend(loc="lower left")
fig.savefig(plotsdir / "cp.pdf")

fig, ax = plt.subplots(figsize=(14.00 * cm, 12.00 * cm))
ax.set_ylabel("$i$ (Mb)")
ax.set_xlabel("$j$ (Mb)")
ax.xaxis.set_label_position("top")
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
px_side = int(np.sqrt(2 * cm1Ddata.size + 1 / 4))
cm2Ddata = np.zeros((px_side, px_side))
cm2Ddata[np.tril_indices(px_side)] = cm1Ddata / ctcfactor
cm2Ddata = cm2Ddata + np.tril(cm2Ddata, -1).T
len_side = px_sz * px_side / (1e6 / bp_part)
map = ax.imshow(
    cm2Ddata,
    cmap="OrRd",
    norm=mpl.colors.LogNorm(vmax=1.0, clip=True),
    extent=[0, len_side, len_side, 0],
)
ax.autoscale(tight=True)
cbar = fig.colorbar(map, ax=ax, aspect=64, pad=1 / 64)
cbar.ax.set_ylabel("$P(i,j)$")
fig.savefig(plotsdir / "cm.pdf", dpi=1600)
