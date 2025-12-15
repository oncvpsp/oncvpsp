#!/bin/bash
#runs all test data and compares with reference out files

sr_tests="03_Li 07_N 07_N_dl 08_O 14_Si 14_Si_st 14_Si_UPF 17_Cl 19_K 19_K_st 20_Ca 22_Ti 29_Cu 32_Ge 56_Ba 57_La 60_Nd_GHOST 83_Bi"
r_tests="52_Te 74_W 80_Hg"

for test in $sr_tests; do
  echo "[oncvpsp.x | $test]"
  ../run.sh $test -np
  ../compare.sh $test;
  summary=$(grep "Summary" $test.diff);
  echo "    $summary";
done

for test in $r_tests; do
  echo "[oncvpspr.x | $test]"
  ../run_r.sh $test -np
  ../compare.sh ${test}_r;
  summary=$(grep "Summary" ${test}_r.diff);
  echo "    $summary";
done

grep Summary *.diff | sed -e /diff/s/^/\\\n/ >TEST.report
