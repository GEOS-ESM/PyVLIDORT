#!/bin/bash

# Sample commands:
# 2s_check2.bash gfortran
# 2s_check2.bash ifort


# Run 2stream checks

file_count=0
echo
echo 'checking results of tests...'
echo

# for serial tests
echo
echo 'doing serial test loop'
for filename in results_2s_*.* ; do
  ref_filename=saved_results/$1/$filename
  if [ -e $ref_filename ]; then
    2s_diff $ref_filename $filename
    echo $((++file_count)) files processed
  fi
done

# for OpenMP tests (if the result files exist)
echo
echo
exist_flag="dud$(ls results_2s_*.*_nt*)"
#echo "exist_flag:0:6 = |${exist_flag:0:6}|"
if [ ${exist_flag:0:6} = 'dudres' ]; then
  #have OpenMP results - do OpenMP test loop
  skip='0'
elif [ ${exist_flag:0:6} = 'dud' ]; then
  #do not have OpenMP results - skip OpenMP test loop
  skip='1'
fi

if [ $skip = '0' ]; then
  #echo
  echo 'doing OpenMP test loop'
  for filename in results_2s_*.*_nt* ; do
    #chop off filename suffix
    fn_prefix=${filename%_nt*}
    #echo 'fn_prefix = '$fn_prefix
    mod_fn=${fn_prefix}
    #echo 'mod_fn    = '$mod_fn
    ref_filename=saved_results/$1/$mod_fn
    #echo
    #echo 'ref_filename = '$ref_filename
    #echo 'filename     = '$filename
    if [ -e $ref_filename ]; then
      2s_diff $ref_filename $filename
      echo $((++file_count)) files processed
    fi
  done
else
  echo 'no OpenMP test results to check'
fi

echo
echo 'done'
echo

exit 0
