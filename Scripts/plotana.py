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

mpl.use("pdf")

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"

cm = 1 / 2.54

mpl.rcParams["figure.constrained_layout.use"] = True

# Set auxiliary variables

lenfactor = 1000 / 33
bp_part = 33 * 200

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
labellist = ["LADh", "LNDe", "total"]
for i_t in range(3):
    x = df_rcd[i_t]["r_b"] / lenfactor
    y = df_rcd[i_t]["avg"]
    y_min = df_rcd[i_t]["avg"] - df_rcd[i_t]["sem"]
    y_max = df_rcd[i_t]["avg"] + df_rcd[i_t]["sem"]
    ax.step(x, y, color=colorlist_3[i_t], label=labellist[i_t])
    ax.fill_between(
        x, y_min, y_max, step="pre", color=colorlist_3[i_t], linewidth=0.0, alpha=0.50
    )
ax.autoscale(tight=True)
ax.set_ylim(bottom=0.0)
ax.legend(loc="upper left")
fig.savefig(plotsdir / "rcd.pdf")

# ax1.set_xscale("log")
# ax1.set_yscale("log")
# for i_c in range(6):
#     if not df_sd[i_c].empty:
#         x = df_sd[i_c]["s"]
#         y = df_sd[i_c]["avg"]
#         y_min = df_sd[i_c]["avg"] - df_sd[i_c]["sem"]
#         y_max = df_sd[i_c]["avg"] + df_sd[i_c]["sem"]
#         ax1.plot(x, y, color=colorlist_6[i_c])
#         ax1.fill_between(
#             x, y_min, y_max, color=colorlist_6[i_c], linewidth=0.0, alpha=0.50
#         )
# ax1.autoscale(tight=True)
# ax1.set_ylim(bottom=1.0)

# ax2.set_xscale("log")
# ax2.set_yscale("log")
# for i_c in range(6):
#     if not df_cp[i_c].empty:
#         df_cp[i_c] = df_cp[i_c].loc[df_cp[i_c]["avg"] != 0]
#         x = df_cp[i_c]["s"]
#         y = df_cp[i_c]["avg"] / ctcfactor
#         e = df_cp[i_c]["sem"] / ctcfactor
#         ax2.scatter(x, y, color=colorlist_6[i_c])
#         ax2.errorbar(x, y, yerr=e, color=colorlist_6[i_c], linestyle="None", alpha=0.50)
# ax2.autoscale(tight=True)

fig, ax = plt.subplots(figsize=(10.00 * cm, 8.00 * cm))
ax.set_ylabel("$i$ (Mb)")
ax.set_xlabel("$j$ (Mb)")
ax.xaxis.set_label_position("top")
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
px_side = int(np.sqrt(2 * cm1Ddata.size + 1 / 4))
cm2Ddata = np.zeros((px_side, px_side))
cm2Ddata[np.tril_indices(px_side)] = cm1Ddata / ctcfactor
cm2Ddata = cm2Ddata + np.tril(cm2Ddata, -1).T
len_side = bp_part * px_sz * px_side / 1e6
map = ax.imshow(
    cm2Ddata,
    cmap="OrRd",
    norm=mpl.colors.LogNorm(vmax=1.0),
    extent=[0, len_side, len_side, 0],
)
map.get_cmap().set_bad("white")
ax.autoscale(tight=True)
cbar = fig.colorbar(map, ax=ax, aspect=64, pad=1 / 64)
cbar.ax.set_ylabel("$P(i,j)$")
fig.savefig(plotsdir / "cm.pdf", dpi=1600)
