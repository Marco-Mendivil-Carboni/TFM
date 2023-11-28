#Imports

from subprocess import run
from subprocess import DEVNULL

from pathlib import Path
from shutil import rmtree

#Create test directory

testdir = Path("Simulations/test")

if testdir.exists():
    rmtree(testdir)

testdir.mkdir()

#Define check function

testidx = 0 #test index

def check(returnc):
    global testidx
    testidx += 1
    if returnc == 0:
        print("Test "+str(testidx)+": ✅")
    else:
        print("Test "+str(testidx)+": ❌")
    print("---")

#Write parameter file

parfilepath = testdir/"adjustable-parameters.dat"

with open(parfilepath,"w") as parfile:
    parfile.write("number_of_particles 8192\n")
    parfile.write("nucleus_radius 14.0\n")
    parfile.write("opening_radius 8.0\n")
    parfile.write("bleb_radius 10.0\n")
    parfile.write("number_of_lbs 512\n")

#Run tests

print("---")

check(run(["./Program/bin/performsim",str(testdir)]).returncode)

check(run(["./Program/bin/performsim",str(testdir)]).returncode)

check(run(["./Program/bin/performsim",str(testdir),"0"]).returncode)

check(run(["./Program/bin/analyzesim",str(testdir)]).returncode)

check(run(["vmd","-e","./Program/visualize.tcl","-args",str(testdir),"0"],
    stdout=DEVNULL,stderr=DEVNULL).returncode)

check(run(["vmd","-e","./Program/visualize.tcl","-args",str(testdir),"1"],
    stdout=DEVNULL,stderr=DEVNULL).returncode)

#Remove test directory

removetest = input("Remove test? (yes/no): ")
if removetest == "yes":
    rmtree(testdir)
