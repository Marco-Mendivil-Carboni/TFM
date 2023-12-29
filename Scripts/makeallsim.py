#!/usr/bin/python

#Imports

from math import cos, asin

from subprocess import run

from pathlib import Path

from itertools import product

from signal import signal
from signal import SIGUSR1

#Define setstop function

stop = False

def setstop(sig,stack):
    global stop
    stop = True

signal(SIGUSR1,setstop)

#Define getsimname function

def getsimname(N,R_n,R_o,R_b,n_l):
    simname = "{:05.0f}".format(N)
    simname += "-{:05.2f}".format(R_n)
    simname += "-{:05.2f}".format(R_o)
    simname += "-{:05.2f}".format(R_b)
    simname += "-{:05.0f}".format(n_l)
    return simname

#Define writeparam function

def writeparam(simdir,N,R_n,R_o,R_b,n_l):
    with open(simdir/"adjustable-parameters.dat","w") as parfile:
        parfile.write("number_of_particles {:05.0f}\n".format(N))
        parfile.write("nucleus_radius {:05.2f}\n".format(R_n))
        parfile.write("opening_radius {:05.2f}\n".format(R_o))
        parfile.write("bleb_radius {:05.2f}\n".format(R_b))
        parfile.write("number_of_lbs {:05.0f}\n".format(n_l))

#Define makesim function

numberofsim = 2
filespersim = 8

def makesim(simdir):

    newsim = False
    pattern = "initial-condition-*"
    while len(list(simdir.glob(pattern))) < numberofsim:
        run(["./Program/bin/ccp-perform",str(simdir)])
        newsim = True

    for simidx in range(numberofsim):
        pattern = "trajectory-{:03}-*".format(simidx)
        while len(list(simdir.glob(pattern))) < filespersim:
            run(["./Program/bin/ccp-perform",str(simdir),str(simidx)])
            newsim = True

    pattern = "analysis-*"
    if len(list(simdir.glob(pattern))) == 0 or newsim:
        run(["./Program/bin/ccp-analyze",str(simdir)])

    if stop:
        exit()

#Make simulations

simrootdir = Path("Simulations")

N = 32768

for i, j in product(range(3),range(1)):

    cvf = 0.2 + 0.1 * i
    laf = 0.00 + 0.25 * j

    R_n = 0.5 + 0.5 * ((N / cvf) ** (1 / 3))

    for k, l in product(range(1),range(2)):

        R_b = R_n * ((1 + k) / 3)  ** (1 / 3)
        R_o = R_b * ((1 + l) / 3)  ** (1 / 2)

        noacf = 2.0 / (1.0 + cos(asin(R_o / R_n)))

        n_l = laf * 4.0 / (((0.5 / (R_n - 1.154701)) ** 2) * noacf)

        simdir = simrootdir/getsimname(N,R_n,R_o,R_b,n_l)
        simdir.mkdir(exist_ok=True)

        writeparam(simdir,N,R_n,R_o,R_b,n_l)
        makesim(simdir)
