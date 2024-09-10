#!/usr/bin/python3

# Imports

from subprocess import run

from pathlib import Path

from signal import signal
from signal import SIGUSR1

from math import cos, asin, floor

# Define setstop function

stop = False


def setstop(sig, stack):
    global stop
    stop = True


signal(SIGUSR1, setstop)

# Define check function


def check(returnc: int) -> None:
    if returnc != 0 or stop:
        exit()


# Define simparam class


class simparam:
    N_def = 30386
    R_n_def = 0.0
    R_o_def = 0.0
    R_b_def = 0.0
    n_l_def = 0

    def __init__(self, *, N=N_def, cvf, laf, brv=0.0, ora=0.0) -> None:
        self.N = N
        self.cvf = cvf
        self.laf = laf
        self.brv = brv
        self.ora = ora

        rco = 1.154701
        self.R_n = (rco / 2) + (rco / 2) * ((N / cvf) ** (1 / 3))
        self.R_b = self.R_n * brv ** (1 / 3)
        self.R_o = self.R_b * ora ** (1 / 2)

        noacf = 2.0 / (1.0 + cos(asin(self.R_o / self.R_n)))
        self.n_l = floor(laf * 4 / ((0.5 / (self.R_n - rco)) ** 2) * noacf)


# Define simname function


def simname(sp: simparam) -> str:
    simname = "{:05.0f}".format(sp.N)
    simname += "-{:03.1f}".format(sp.cvf)
    simname += "-{:03.1f}".format(sp.laf)
    simname += "-{:03.1f}".format(sp.brv)
    simname += "-{:03.1f}".format(sp.ora)
    return simname


# Define writeparam function


def writeparam(simdir: Path, sp: simparam) -> None:
    with open(simdir / "adjustable-parameters.dat", "w") as parfile:
        if sp.N != sp.N_def:
            parfile.write("number_of_particles {:05.0f}\n".format(sp.N))
        if sp.R_n != sp.R_n_def:
            parfile.write("nucleus_radius {:09.6f}\n".format(sp.R_n))
        if sp.R_o != sp.R_o_def:
            parfile.write("opening_radius {:09.6f}\n".format(sp.R_o))
        if sp.R_b != sp.R_b_def:
            parfile.write("bleb_radius {:09.6f}\n".format(sp.R_b))
        if sp.n_l != sp.n_l_def:
            parfile.write("number_of_lbs {:05.0f}\n".format(sp.n_l))


# Define makesim function

simrootdir = Path("Simulations")

numberofsim = 8
filespersim = 128


def makesim(sp: simparam) -> None:
    simdir = simrootdir / simname(sp)
    simdir.mkdir(exist_ok=True)

    writeparam(simdir, sp)
    newsim = False

    pattern = "initial-condition-*"
    while len(list(simdir.glob(pattern))) < numberofsim:
        check(
            run(
                ["./Program/bin/ccp-perform", str(simdir)],
            ).returncode
        )
        newsim = True

    for simidx in range(numberofsim):
        pattern = "trajectory-{:03}-*".format(simidx)
        while len(list(simdir.glob(pattern))) < filespersim:
            check(
                run(
                    ["./Program/bin/ccp-perform", str(simdir), str(simidx)],
                ).returncode
            )
            newsim = True

    pattern = "analysis-*"
    if len(list(simdir.glob(pattern))) == 0 or newsim:
        check(
            run(
                ["./Program/bin/ccp-analyze", str(simdir)],
            ).returncode
        )


# Make simulations

makesim(simparam(cvf=0.200, laf=0.400))

makesim(simparam(cvf=0.200, laf=0.000))
