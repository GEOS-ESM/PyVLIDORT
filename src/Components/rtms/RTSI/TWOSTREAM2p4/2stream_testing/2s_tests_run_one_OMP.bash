#!/bin/bash

# This script runs the two 2S tests which have been set up to run using OpenMP
# (see http://openmp.org for further information about OpenMP)

# Sample commands:
# 2s_tests_run_one_OMP.bash ONLY gfortran OPENMP=t
# 2s_tests_run_one_OMP.bash LPCS ifort OPENMP=t (note: ifort not active yet)

# Do some set up

ulimit -s unlimited         # Maximum stack size of the main thread
ulimit -c unlimited         # Maximum size of core file created if a problem occurs

# OMP_STACKSIZE -->         # Maximum stack size of OpenMP-spawned threads

# gfortran requirements for:
#export OMP_STACKSIZE=512M
#export OMP_STACKSIZE=750M
#export OMP_STACKSIZE=1G
#export OMP_STACKSIZE=16M
export OMP_STACKSIZE=50M     # 2S ONLY

#export GOMP_STACKSIZE=16M
#export GOMP_STACKSIZE=10M
#export GOMP_STACKSIZE=50M
#export GOMP_STACKSIZE=250M
#export GOMP_STACKSIZE=500M



#echo
#echo 'OMP_STACKSIZE = '$OMP_STACKSIZE
#echo

#ulimit -a # -a All current limits are reported
#exit 0

echo
echo "making test_2s_$1_2p4_OMP.exe ..."
echo
make test_2s_$1_2p4_OMP.exe FC=$2 $3 $4
echo
echo "running test_2s_$1_2p4_OMP.exe ..."
echo
./test_2s_$1_2p4_OMP.exe
#valgrind -v --tool=memcheck ./test_2s_$1_2p4_OMP.exe
#valgrind --tool=memcheck --leak-check=full ./test_2s_$1_2p4_OMP.exe

echo
echo 'done'
echo

exit 0
