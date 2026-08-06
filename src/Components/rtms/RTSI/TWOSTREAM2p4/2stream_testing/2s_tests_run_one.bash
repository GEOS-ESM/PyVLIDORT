#!/bin/bash

# Sample commands:
# 2s_tests_run_one.bash ONLY gfortran
# 2s_tests_run_one.bash LPCS ifort

echo
echo "making test_2s_$1_2p4.exe ..."
echo
make test_2s_$1_2p4.exe FC=$2 $3
echo
echo "running test_2s_$1_2p4.exe ..."
echo
./test_2s_$1_2p4.exe
#make clean

echo
echo 'done'
echo

exit 0
