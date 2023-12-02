#Imports

from pathlib import Path

import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt

#Set fixed parameters

mpl.use("pdf")

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [12.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

#Load data into dataframe

plotsrootdir = Path("Plots")
datafilepath = plotsrootdir/"performance.dat"

df_all = pd.read_csv(datafilepath,sep=" ",names=["GPU","N","t_e"])

df_1 = df_all.loc[df_all["GPU"] == "GF-920M"]
df_2 = df_all.loc[df_all["GPU"] == "GF-RTX-3050"]
df_3 = df_all.loc[df_all["GPU"] == "RTX-A4000"]

#Make performance plot

plt.xscale("log")
plt.yscale("log")

plt.xlabel("$N$")
plt.ylabel("$t_e$ (ms)")

plt.plot(df_1.N,df_1.t_e,marker="o",color="#d81e2c",label="GeForce-920M")
plt.plot(df_2.N,df_2.t_e,marker="o",color="#a31cc5",label="GeForce-RTX-3050")
plt.plot(df_3.N,df_3.t_e,marker="o",color="#194bb2",label="RTX-A4000")

plt.legend(loc="upper left")

plt.savefig(plotsrootdir/"performance.pdf")
