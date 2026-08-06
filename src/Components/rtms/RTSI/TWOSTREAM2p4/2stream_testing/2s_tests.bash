#!/bin/bash

# Sample commands:
# 2s_tests.bash gfortran
# 2s_tests.bash ifort

echo
echo 'making test executables ...'
echo
make main FC=$1 $2

echo
echo 'running 2stream test: test_2s_ONLY_2p4.exe'
echo
./test_2s_ONLY_2p4.exe

echo
echo 'running 2stream test: test_2s_LPCS_2p4.exe'
echo
./test_2s_LPCS_2p4.exe

#make clean

echo
echo 'done'
echo

exit 0
