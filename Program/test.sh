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

./Program/bin/simulate | grep "no arguments\|$"
check $?

./Program/bin/simulate 1 2 3 | grep "extra arguments\|$"
check $?

./Program/bin/simulate wrong-dir | grep "unable to open wrong-dir\|$"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "N   512"; echo "T   298.0"; echo "R   10.00"; echo "F   100";
} >> "${testdir}/adjustable-parameters.dat"
./Program/bin/simulate $testdir | grep "error reading T\|$"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "T   298.0"; echo "N   000"; echo "R   10.00"; echo "F   100";
} >> "${testdir}/adjustable-parameters.dat"
./Program/bin/simulate $testdir | grep "error reading N\|$"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "T   298.0"; echo "N   512"; echo "R   01.00"; echo "F   100";
} >> "${testdir}/adjustable-parameters.dat"
./Program/bin/simulate $testdir | grep "chromatin volume fraction above 0.5\|$"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
{ echo "T   298.0"; echo "N   512"; echo "R   10.00"; echo "F   100";
} >> "${testdir}/adjustable-parameters.dat"

./Program/bin/simulate $testdir
[[ -f "${testdir}/all-messages.log" && \
-f "${testdir}/initial-condition-000.gro" && \
-f "${testdir}/trajectory-000-000.trr" && \
-f "${testdir}/checkpoint-000.bin" ]]
check $?

./Program/bin/simulate $testdir
[[ -f "${testdir}/initial-condition-001.gro" && \
-f "${testdir}/trajectory-001-000.trr" && \
-f "${testdir}/checkpoint-001.bin" ]]
check $?

./Program/bin/simulate $testdir 0
[[ -f "${testdir}/trajectory-000-001.trr" ]]
check $?

vmd -e ./Program/visualize.tcl -args $testdir 0 > /dev/null

vmd -e ./Program/visualize.tcl -args $testdir 1 > /dev/null

rm -rI $testdir
