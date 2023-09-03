#!/bin/bash

testdir="Simulations/bash-script-tests"
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

./Program/bin/simulate | grep "no arguments"
check $?

./Program/bin/simulate 1 2 3 | grep "extra arguments"
check $?

./Program/bin/simulate wrong-dir | grep "unable to open wrong-dir"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "N   512"; echo "T   298.0"; echo "R   10.00"; echo "F   100";
} >> "${testdir}/adjustable-parameters.dat"
./Program/bin/simulate $testdir | grep "error reading T"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "T   298.0"; echo "N   000"; echo "R   10.00"; echo "F   100";
} >> "${testdir}/adjustable-parameters.dat"
./Program/bin/simulate $testdir | grep "error reading N"
check $?

#check chromatin volume fraction too high error

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "T   298.0"; echo "N   512"; echo "R   10.00"; echo "F   100";
} >> "${testdir}/adjustable-parameters.dat"

./Program/bin/simulate $testdir
[[ -f "${testdir}/complete-history.log" && \
-f "${testdir}/initial-configuration-000.gro" && \
-f "${testdir}/trajectory-positions-000-000.trr" ]]
check $?

./Program/bin/simulate $testdir
[[ -f "${testdir}/initial-configuration-001.gro" && \
-f "${testdir}/trajectory-positions-001-000.trr" ]]
check $?

./Program/bin/simulate $testdir 0
[[ -f "${testdir}/trajectory-positions-000-001.trr" ]]
check $?

# vmd -e ./Program/visualize-chromatin.tcl -args $testdir 0 > /dev/null

# vmd -e ./Program/visualize-chromatin.tcl -args $testdir 1 > /dev/null

rm -rI $testdir
