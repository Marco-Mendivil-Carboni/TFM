#Imports

import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

#Test

cm = 1 / 2.54

mpl.use("pgf")
mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["pgf.rcfonts"] = False#neccesary?
mpl.rcParams["pgf.texsystem"] = "pdflatex"#neccesary?
plt.rcParams["figure.figsize"] = [12.00*cm,8.00*cm]
plt.rcParams["figure.autolayout"] = True
plt.rcParams["legend.frameon"] = False

columns = ["N","t_e"]

df_A4000 = pd.read_csv("Plots/RTX-A4000.dat",sep=" ",header=None,names=columns)
df_3050 = pd.read_csv("Plots/GeForce-RTX-3050.dat",sep=" ",header=None, 
    names=columns)
df_920M = pd.read_csv("Plots/GeForce-920M.dat",sep=" ",header=None,
    names=columns)

plt.xscale("log")
plt.yscale("log")
plt.xlim([128,65536])
plt.ylim([0.1,100.0])
plt.xlabel("$N$")
plt.ylabel("$t_e$ (ms)")
plt.plot(df_920M.N,df_920M.t_e,marker="o",color="#d81e2c",label="GeForce-920M")
plt.plot(df_3050.N,df_3050.t_e,marker="o",color="#a31cc5",
         label="GeForce-RTX-3050")
plt.plot(df_A4000.N,df_A4000.t_e,marker="o",color="#194bb2",label="RTX-A4000")
plt.legend(loc="upper left")

plt.savefig("Plots/performance.pdf")
