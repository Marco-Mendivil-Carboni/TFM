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
echo "N   512" >> "${testdir}/adjustable-parameters.dat"
echo "T   298.0" >> "${testdir}/adjustable-parameters.dat"
./Program/bin/simulate $testdir | grep "error reading T"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
echo "T   298.0" >> "${testdir}/adjustable-parameters.dat"
echo "N   0" >> "${testdir}/adjustable-parameters.dat"
./Program/bin/simulate $testdir | grep "error reading N"
check $?

echo -n > "${testdir}/adjustable-parameters.dat"
echo "T   298.0" >> "${testdir}/adjustable-parameters.dat"
echo "N   512" >> "${testdir}/adjustable-parameters.dat"
echo "R   10.00" >> "${testdir}/adjustable-parameters.dat"
echo "F   100" >> "${testdir}/adjustable-parameters.dat"

./Program/bin/simulate $testdir
[[ $(ls $testdir | grep -c -e "complete-history.log" \
    -e "initial-configuration-000.gro" \
    -e "trajectory-positions-000-000.trr") -eq 3 ]]
check $?

./Program/bin/simulate $testdir
[[ $(ls $testdir | grep -c -e "initial-configuration-001.gro" \
    -e "trajectory-positions-001-000.trr") -eq 2 ]]
check $?

./Program/bin/simulate $testdir 0
[[ $(ls $testdir | grep -c "trajectory-positions-000-001.trr") -eq 1 ]]
check $?

# vmd -e ./Program/visualize-chromatin.tcl -args $testdir 0 > /dev/null

rm -rI $testdir
