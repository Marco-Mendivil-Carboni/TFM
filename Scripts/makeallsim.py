#!/usr/bin/python3

# Imports

from subprocess import run

from pathlib import Path

from signal import signal
from signal import SIGUSR1

# Define setstop function

stop = False


def setstop(sig, stack):
    global stop
    stop = True


signal(SIGUSR1, setstop)

# Define simname function


def simname(N, R_n, n_l, R_o, R_b):
    simname = "{:05.0f}".format(N)
    simname += "-{:05.2f}".format(R_n)
    simname += "-{:05.0f}".format(n_l)
    simname += "-{:05.2f}".format(R_o)
    simname += "-{:05.2f}".format(R_b)
    return simname


# Define writeparam function


def writeparam(simdir, N, R_n, n_l, R_o=0.0, R_b=0.0):
    with open(simdir / "adjustable-parameters.dat", "w") as parfile:
        parfile.write("number_of_particles {:05.0f}\n".format(N))
        parfile.write("nucleus_radius {:05.2f}\n".format(R_n))
        parfile.write("number_of_lbs {:05.0f}\n".format(n_l))
        parfile.write("opening_radius {:05.2f}\n".format(R_o))
        parfile.write("bleb_radius {:05.2f}\n".format(R_b))


# Define makesim function

simrootdir = Path("Simulations")

numberofsim = 4
filespersim = 8


def makesim(N, R_n, n_l, R_o, R_b):
    simdir = simrootdir / simname(N, R_n, n_l, R_o, R_b)
    simdir.mkdir(exist_ok=True)

    writeparam(simdir, N, R_n, n_l, R_o, R_b)

    newsim = False
    pattern = "initial-condition-*"
    while len(list(simdir.glob(pattern))) < numberofsim:
        run(["./Program/bin/ccp-perform", str(simdir)])
        newsim = True
        if stop:
            exit()

    for simidx in range(numberofsim):
        pattern = "trajectory-{:03}-*".format(simidx)
        while len(list(simdir.glob(pattern))) < filespersim:
            run(["./Program/bin/ccp-perform", str(simdir), str(simidx)])
            newsim = True
            if stop:
                exit()

    pattern = "analysis-*"
    if len(list(simdir.glob(pattern))) == 0 or newsim:
        run(["./Program/bin/ccp-analyze", str(simdir)])


# Make simulations

N = 18240
R_n = 40.91
n_l = 0

R_o = 0.0
R_b = 0.0

makesim(N, R_n, n_l, R_o, R_b)

n_l = 7284

makesim(N, R_n, n_l, R_o, R_b)
