#!/bin/bash

testdir="Simulations/test"
testidx=1

check () {
    if [[ $1 -eq 0 ]]
    then
        echo "Test $testidx: ✅"
    else
        echo "Test $testidx: ❌"
    fi
    ((testidx=testidx+1))
    echo "---"
}

rm -rf $testdir
mkdir $testdir

echo "---"

./Program/bin/performsim | grep "no arguments"
check $?

./Program/bin/performsim 1 2 3 | grep "extra arguments"
check $?

./Program/bin/performsim wrong-dir | grep "unable to open wrong-dir"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "number_of_particles -10"; echo "confinement_radius 10.00";
} >> "${testdir}/adjustable-parameters.dat"
./Program/bin/performsim $testdir | grep "number_of_particles out of range"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "number_of_particles 512";
} >> "${testdir}/adjustable-parameters.dat"
./Program/bin/performsim $testdir | grep "confinement_radius out of range"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "number_of_particles 512"; echo "confinement_radius 1.00";
} >> "${testdir}/adjustable-parameters.dat"
./Program/bin/performsim $testdir | grep "chromatin volume fraction above 0.5"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "number_of_particles 1024"; echo "confinement_radius 8.0";
  echo "steps_per_frame 1024";
} >> "${testdir}/adjustable-parameters.dat"

./Program/bin/performsim $testdir
[[ -f "${testdir}/all-messages.log" && \
-f "${testdir}/initial-condition-000.gro" && \
-f "${testdir}/trajectory-000-000.trr" && \
-f "${testdir}/checkpoint-000.bin" ]]
check $?

./Program/bin/performsim $testdir
[[ -f "${testdir}/initial-condition-001.gro" && \
-f "${testdir}/trajectory-001-000.trr" && \
-f "${testdir}/checkpoint-001.bin" ]]
check $?

nvprof ./Program/bin/performsim $testdir 0
[[ -f "${testdir}/trajectory-000-001.trr" ]]
check $?

vmd -e ./Program/visualize.tcl -args $testdir 0 > /dev/null
vmd -e ./Program/visualize.tcl -args $testdir 1 > /dev/null

rm -rI $testdir
