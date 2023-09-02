#!/bin/bash

mkdir tmp

cp Simulations/test/adjustable-parameters.dat tmp

#test 1

if ./Program/bin/simulate | grep -q 'no arguments';
then
    echo "test 1: ✅"
else
    echo "test 1: ❌"
fi

#test 2

if ./Program/bin/simulate 1 2 3 | grep -q 'extra arguments';
then
    echo "test 2: ✅"
else
    echo "test 2: ❌"
fi

#test 3

if ./Program/bin/simulate wrong-dir | grep -q 'unable to open wrong-dir';
then
    echo "test 3: ✅"
else
    echo "test 3: ❌"
fi

#test 4

./Program/bin/simulate tmp > /dev/null

if ls tmp | grep -q -e 'complete-history.log' \
                    -e 'initial-configration-000.gro' \
                    -e 'trajectory-positions-000-000.trr';
then
    echo "test 4: ✅"
else
    echo "test 4: ❌"
fi

#test 5

./Program/bin/simulate tmp > /dev/null

if ls tmp | grep -q -e 'initial-configration-001.gro' \
                    -e 'trajectory-positions-001-000.trr';
then
    echo "test 5: ✅"
else
    echo "test 5: ❌"
fi

#test 6

./Program/bin/simulate tmp 0 > /dev/null

if ls tmp | grep -q 'trajectory-positions-000-001.trr';
then
    echo "test 6: ✅"
else
    echo "test 6: ❌"
fi

rm -r tmp
