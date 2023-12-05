#!/usr/bin/python

#Imports

from pathlib import Path

import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt

from io import StringIO

#Set fixed parameters

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [32.00 * cm, 16.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

#Select simulation

simrootdir = Path("Simulations")
while True:
    simname = input("Simulation name? ")
    simdir = simrootdir/simname
    if simdir.exists():
        break
    else:
        print("Simulation not found")

#Load data into dataframes

parfilepath = simdir/"adjustable-parameters.dat"
df_par = pd.read_csv(parfilepath,
    delim_whitespace=True,comment="#",header=None)

anafilepath = simdir/"analysis-fin.dat"
with open(anafilepath) as anafile:
    blocklist = anafile.read().split("\n\n\n")

df_s = pd.read_csv(StringIO(blocklist[0]),
    delim_whitespace=True,comment="#",header=None)
df_s.columns = ["avg","sqrt(var)","sem"]
df_s.insert(loc=0,column="simobs",
    value=["dcm","rg2","nop","ncf"])

df_rcd = list()
for i_t in range(3):
    df_rcd.append(pd.read_csv(StringIO(blocklist[i_t+1]),
        delim_whitespace=True,comment="#",header=None))
    df_rcd[i_t].columns = ["r_b","avg","sqrt(var)","sem"]

df_msd = pd.read_csv(StringIO(blocklist[4]),
    delim_whitespace=True,comment="#",header=None)
df_msd.columns = ["s","avg","sqrt(var)","sem"]

#Make analysis plots

fig,ax = plt.subplots(2,2,height_ratios=[0.25,1])

ax[0,0].axis("off")
ax[0,0].table(cellText=df_par.values,
    loc="center",cellLoc="left",edges="horizontal")

ax[0,1].axis("off")
ax[0,1].table(cellText=df_s.values,colLabels=df_s.columns,
    loc="center",colLoc="right",edges="horizontal")

colorlist = ["#d81e2c","#a31cc5","#194bb2"]
ax[1,0].set_xlabel("r_b")
ax[1,0].set_ylabel("rcd")
ax[1,0].autoscale(tight=True)
for i_t in range(3):
    x = df_rcd[i_t]["r_b"]
    y = df_rcd[i_t]["avg"]
    y_min = df_rcd[i_t]["avg"]-df_rcd[i_t]["sem"]
    y_max = df_rcd[i_t]["avg"]+df_rcd[i_t]["sem"]
    ax[1,0].step(x,y,color=colorlist[i_t])
    ax[1,0].fill_between(x,y_min,y_max,step="pre",
        color=colorlist[i_t],linewidth=0.0,alpha=0.50)

ax[1,1].set_xscale("log")
ax[1,1].set_yscale("log")
ax[1,1].set_xlabel("s")
ax[1,1].set_ylabel("msd")
ax[1,1].autoscale(tight=True)
x = df_msd["s"]
y = df_msd["avg"]
y_min = df_msd["avg"]-df_msd["sem"]
y_max = df_msd["avg"]+df_msd["sem"]
ax[1,1].plot(x,y,color="#169f62")
ax[1,1].fill_between(x,y_min,y_max,
    color="#169f62",linewidth=0.0,alpha=0.50)

#View analysis

man = plt.get_current_fig_manager()
man.set_window_title(simname+" analysis")

plt.show()
