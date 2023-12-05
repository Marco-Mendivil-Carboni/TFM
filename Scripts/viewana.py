#Imports

from pathlib import Path

import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt

from io import StringIO

#Set fixed parameters

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [12.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

#Define readblock function

def readblock(block):
    blockio = StringIO(block)
    return pd.read_csv(blockio,delim_whitespace=True,header=None,comment="#")

#Read analysis data

simrootdir = Path("Simulations/test")
filepath = simrootdir/"analysis-fin.dat"

with open(filepath) as file:
    blocklist = file.read().split("\n\n\n")

df_s = readblock(blocklist[0])
df_s.columns = ["avg","sqrt(var)","sem"]
df_s.index = ["dcm","rg2","nop","ncf"]
print(df_s)

df_rcd = list()
for i_t in range(3):
    df_rcd.append(readblock(blocklist[i_t+1]))
    df_rcd[i_t].columns = ["r_b","avg","sqrt(var)","sem"]
    print(df_rcd[i_t])

df_msd = readblock(blocklist[4])
df_msd.columns = ["s","avg","sqrt(var)","sem"]
print(df_msd)

#Make plots

for i_t in range(3):

    x = df_rcd[i_t]["r_b"]
    y = df_rcd[i_t]["avg"]
    y_min = df_rcd[i_t]["avg"]-df_rcd[i_t]["sem"]
    y_max = df_rcd[i_t]["avg"]+df_rcd[i_t]["sem"]

    plt.step(x,y)
    plt.fill_between(x,y_min,y_max,step="pre",color="k",alpha=0.25)

plt.show()

x = df_msd["s"]
y = df_msd["avg"]
y_min = df_msd["avg"]-df_msd["sem"]
y_max = df_msd["avg"]+df_msd["sem"]

plt.xscale("log")
plt.yscale("log")
plt.plot(x,y)
plt.fill_between(x,y_min,y_max,color="k",alpha=0.25)

plt.show()
