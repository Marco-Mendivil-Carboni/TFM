#Imports

from pathlib import Path

import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt

from io import StringIO

#Set fixed parameters

mpl.rcParams["figure.figsize"] = [12.0,6.0]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

#Define readblock function

def readblock(block):
    blockio = StringIO(block)
    return pd.read_csv(StringIO(block),delim_whitespace=True,comment="#",header=None)

#Read analysis data

simrootdir = Path("Simulations/test")
filepath = simrootdir/"analysis-fin.dat"

with open(filepath) as file:
    blocklist = file.read().split("\n\n\n")

df_s = readblock(blocklist[0])
df_s.columns = ["avg","sqrt(var)","sem"]
df_s.index = ["dcm","rg2","nop","ncf"]

df_rcd = list()
for i_t in range(3):
    df_rcd.append(readblock(blocklist[i_t+1]))
    df_rcd[i_t].columns = ["r_b","avg","sqrt(var)","sem"]

df_msd = readblock(blocklist[4])
df_msd.columns = ["s","avg","sqrt(var)","sem"]

#Make plots

fig, ax = plt.subplots(2,2,height_ratios=[0.25,1])

ax[0,0].table(cellText=df_s.values,rowLabels=df_s.index,colLabels=df_s.columns,loc="center")
ax[0,0].axis("off")
ax[0,1].axis("off")

for i_t in range(3):

    x = df_rcd[i_t]["r_b"]
    y = df_rcd[i_t]["avg"]
    y_min = df_rcd[i_t]["avg"]-df_rcd[i_t]["sem"]
    y_max = df_rcd[i_t]["avg"]+df_rcd[i_t]["sem"]

    ax[1,0].step(x,y)
    ax[1,0].fill_between(x,y_min,y_max,step="pre",color="k",alpha=0.25)

x = df_msd["s"]
y = df_msd["avg"]
y_min = df_msd["avg"]-df_msd["sem"]
y_max = df_msd["avg"]+df_msd["sem"]

ax[1,1].set_xscale("log")
ax[1,1].set_yscale("log")
ax[1,1].plot(x,y)
ax[1,1].fill_between(x,y_min,y_max,color="k",alpha=0.25)

plt.show()
