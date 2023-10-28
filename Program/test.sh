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

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "number_of_particles 2048"; echo "confinement_radius 10.0";
  echo "number_of_lbs 128"; echo "frames_per_file 256";
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

./Program/bin/performsim $testdir 0
[[ -f "${testdir}/trajectory-000-001.trr" ]]
check $?

./Program/bin/analyzesim $testdir
[[ -f "${testdir}/analysis-000.dat" && \
-f "${testdir}/analysis-001.dat" && \
-f "${testdir}/analysis-fin.dat" ]]
check $?

vmd -e ./Program/visualize.tcl -args $testdir 0 > /dev/null
check $?

vmd -e ./Program/visualize.tcl -args $testdir 1 > /dev/null
check $?

rm -rI $testdir
