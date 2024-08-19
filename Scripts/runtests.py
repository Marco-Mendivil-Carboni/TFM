#!/usr/bin/python3

# Imports

from subprocess import run
from subprocess import DEVNULL

from pathlib import Path
from shutil import rmtree

from math import floor

# Create test directory

testdir = Path("Simulations/test")

if testdir.exists():
    rmtree(testdir)

testdir.mkdir()

# Define check function

testidx = 0  # test index


def check(returnc: int) -> None:
    global testidx
    testidx += 1
    if returnc == 0:
        print("Test " + str(testidx) + ": ✅")
    else:
        print("Test " + str(testidx) + ": ❌")
    print("---")


# Write parameter file

parfilepath = testdir / "adjustable-parameters.dat"

N = 30386
rco = 1.154701

cvf = 0.200
laf = 0.400

R_n = (rco / 2) + (rco / 2) * ((N / cvf) ** (1 / 3))
n_l = floor(laf * 4 / ((0.5 / (R_n - rco)) ** 2))

with open(parfilepath, "w") as parfile:
    parfile.write("nucleus_radius {:09.6f}\n".format(R_n))
    parfile.write("number_of_lbs {:05.0f}\n".format(n_l))

# Run tests

print("---")

check(run(["./Program/bin/ccp-perform", str(testdir)]).returncode)
check(run(["./Program/bin/ccp-perform", str(testdir), "0"]).returncode)

check(run(["./Program/bin/ccp-perform", str(testdir)]).returncode)
check(run(["./Program/bin/ccp-perform", str(testdir), "1"]).returncode)

check(run(["./Program/bin/ccp-analyze", str(testdir)]).returncode)

check(
    run(
        ["vmd", "-e", "./Scripts/viewsim.tcl", "-args", str(testdir), "0"],
        stdout=DEVNULL,
        stderr=DEVNULL,
    ).returncode
)

check(run(["./Scripts/viewana.py", str(testdir)]).returncode)

# Delete test directory

deletetest = input("Delete test? ").lower()
if deletetest.startswith("y"):
    rmtree(testdir)
