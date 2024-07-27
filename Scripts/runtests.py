#!/usr/bin/python3

# Imports

from subprocess import run
from subprocess import DEVNULL

from pathlib import Path
from shutil import rmtree

# Create test directory

testdir = Path("Simulations/test")

if testdir.exists():
    rmtree(testdir)

testdir.mkdir()

# Define check function

testidx = 0  # test index


def check(returnc):
    global testidx
    testidx += 1
    if returnc == 0:
        print("Test " + str(testidx) + ": ✅")
    else:
        print("Test " + str(testidx) + ": ❌")
    print("---")


# Write parameter file

parfilepath = testdir / "adjustable-parameters.dat"

with open(parfilepath, "w") as parfile:
    parfile.write("number_of_particles 8192\n")
    parfile.write("nucleus_radius 20.0\n")
    parfile.write("number_of_lbs 1024\n")

# Run tests

print("---")

check(run(["./Program/bin/ccp-perform", str(testdir)]).returncode)

check(run(["./Program/bin/ccp-perform", str(testdir)]).returncode)

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
