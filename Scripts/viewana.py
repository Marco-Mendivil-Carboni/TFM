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
simrootdir = Path("Simulations/without-bleb/04096-0.100-0.050")
datafilepath = simrootdir/"analysis-fin.dat"

with open(datafilepath) as datafile:
    datablocks = datafile.read().split("\n\n")
    print(datablocks[2])
