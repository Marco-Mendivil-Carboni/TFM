#Imports

from subprocess import run

from pathlib import Path

from itertools import product

#Define makesim function

numberofsim = 4
filespersim = 4

def makesim(simdir):

    newsim = False
    pattern = "initial-condition-*"
    while len(list(simdir.glob(pattern))) < numberofsim:
        run(["./Program/bin/ccp-perform",str(simdir)])
        newsim = True

    for simidx in range(numberofsim):
        pattern = "trajectory-{:03}-*".format(simidx)
        while len(list(simdir.glob(pattern))) < filespersim:
            run(["./Program/bin/ccp-perform",str(simdir)],str(simidx))
            newsim = True

    pattern = "analysis-*"
    if len(list(simdir.glob(pattern))) == 0 or newsim:
        run(["./Program/bin/ccp-analyze",str(simdir)])

#Make simulations without bleb

simrootdir = Path("Simulations/without-bleb")

for i, j, k in product(range(4),range(4),range(4)):

    N = 4096 * 2 ** i
    cvf = 0.1 + 0.1 * j
    laf = 0.05 + 0.15 * k
    simdir = simrootdir/"{:05}-{:5.3f}-{:5.3f}".format(N, cvf, laf)
    simdir.mkdir(exist_ok=True)

    R_n = 0.5 + 0.5 * ((N / cvf) ** (1 / 3))
    n_l = laf * 4.0 / ((0.5 / (R_n - 1.154701)) ** 2)

    with open(simdir/"adjustable-parameters.dat","w") as parfile:
        parfile.write("number_of_particles {:05.0f}\n".format(N))
        parfile.write("nucleus_radius {:5.2f}\n".format(R_n))
        parfile.write("number_of_lbs {:05.0f}\n".format(n_l))

    makesim(simdir)

#Make simulations with bleb
