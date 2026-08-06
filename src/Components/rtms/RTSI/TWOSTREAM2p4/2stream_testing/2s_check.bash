#!/bin/bash

# Sample commands:
# 2s_check.bash gfortran
# 2s_check.bash ifort


# Run 2stream checks

file_count=0
echo
echo 'checking results of tests...'
echo
for filename in results_2s_*.* ; do
  ref_filename=saved_results/$1/$filename
  if [ -e $ref_filename ]; then
    diff -b $ref_filename $filename > diff_${filename}
    echo $((++file_count)) files processed
  fi
done

echo
echo 'done'
echo

exit 0
