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

# Define check function


def check(returnc: int) -> None:
    if returnc != 0:
        exit()


# Define simparam class


class simparam:
    N_def = 18240
    n_l_def = 0
    R_o_def = 0.0
    R_b_def = 0.0

    def __init__(
        self, R_n: float, N=N_def, n_l=n_l_def, R_o=R_o_def, R_b=R_b_def
    ) -> None:
        self.R_n = R_n
        self.N = N
        self.n_l = n_l
        self.R_o = R_o
        self.R_b = R_b


# Define simname function


def simname(sp: simparam) -> str:
    simname = "{:05.2f}".format(sp.R_n)
    simname += "-{:05.0f}".format(sp.N)
    simname += "-{:05.0f}".format(sp.n_l)
    simname += "-{:05.2f}".format(sp.R_o)
    simname += "-{:05.2f}".format(sp.R_b)
    return simname


# Define writeparam function


def writeparam(simdir: Path, sp: simparam) -> None:
    with open(simdir / "adjustable-parameters.dat", "w") as parfile:
        parfile.write("nucleus_radius {:05.2f}\n".format(sp.R_n))
        if sp.N != sp.N_def:
            parfile.write("number_of_particles {:05.0f}\n".format(sp.N))
        if sp.n_l != sp.n_l_def:
            parfile.write("number_of_lbs {:05.0f}\n".format(sp.n_l))
        if sp.R_o != sp.R_o_def:
            parfile.write("opening_radius {:05.2f}\n".format(sp.R_o))
        if sp.R_o != sp.R_b_def:
            parfile.write("bleb_radius {:05.2f}\n".format(sp.R_b))


# Define makesim function

simrootdir = Path("Simulations")

numberofsim = 8
filespersim = 16


def makesim(sp: simparam) -> None:
    simdir = simrootdir / simname(sp)
    simdir.mkdir(exist_ok=True)

    writeparam(simdir, sp)

    newsim = False
    pattern = "initial-condition-*"
    while len(list(simdir.glob(pattern))) < numberofsim:
        check(run(["./Program/bin/ccp-perform", str(simdir)]).returncode)
        newsim = True
        if stop:
            exit()

    for simidx in range(numberofsim):
        pattern = "trajectory-{:03}-*".format(simidx)
        while len(list(simdir.glob(pattern))) < filespersim:
            check(
                run(["./Program/bin/ccp-perform", str(simdir), str(simidx)]).returncode
            )
            newsim = True
            if stop:
                exit()

    pattern = "analysis-*"
    if len(list(simdir.glob(pattern))) == 0 or newsim:
        check(run(["./Program/bin/ccp-analyze", str(simdir)]).returncode)


# Make simulations

makesim(simparam(R_n=40.91))

makesim(simparam(R_n=40.91, n_l=7284))
