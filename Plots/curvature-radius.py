#Imports

from pathlib import Path

import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt

#Set fixed parameters

mpl.use("pdf")

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [12.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

#Define curvature radius function

R_o = 2
R_n = 2 * np.sqrt(2)

def Volume(x):
    T_1 = 4 * R_n ** 3 * x / (1 - x)
    T_a = 2 + np.sqrt(1 - R_o ** 2 / R_n ** 2)
    T_b = 1 - np.sqrt(1 - R_o ** 2 / R_n ** 2)
    T_2 = R_n ** 3 * T_a * T_b **2
    return T_1 + T_2

def Gamma(v):
    T_1 = R_o ** 12
    T_2 = 8 * R_o ** 6 * v ** 2
    T_a = R_o ** 18 * v ** 2
    T_b = 5 * R_o ** 12 * v ** 4
    T_c = 8 * R_o ** 6 * v ** 6
    T_d = 4 * v ** 8
    T_3 = 4 * np.sqrt(T_a + T_b + T_c + T_d)
    T_4 = 8 * v ** 4
    return T_1 + T_2 + T_3 + T_4

def Radius(v):
    T_1 = Gamma(v) ** (1 / 3)
    T_2 = R_o ** 4
    T_3 = R_o ** 8 * Gamma(v) ** (- 1 / 3)
    return (T_1 + T_2 + T_3)/(4 * v)

#Make curvature radius plot

plotsrootdir = Path("Plots")

fig,ax = plt.subplots()

ax.set_xlabel("$N'/N$")
ax.set_ylabel("$R_c$")

x = np.arange(0.0,0.5,0.001)
y = Radius(Volume(x))
x = x[y < R_n]
y = y[y < R_n]

ax.plot(x,y,color="#d81e2c")

fig.savefig(plotsrootdir/"curvature-radius.pdf")
