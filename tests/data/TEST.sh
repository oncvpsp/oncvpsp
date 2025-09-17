#!/bin/bash
#runs all test data and compares with reference out files

../run.sh 03_Li -np	; 	../compare.sh 03_Li
../run.sh 07_N -np	; 	../compare.sh 07_N
../run.sh 07_N_dl -np	; 	../compare.sh 07_N_dl
../run.sh 08_O -np	; 	../compare.sh 08_O
../run.sh 14_Si -np	; 	../compare.sh 14_Si
../run.sh 14_Si_st -np	; 	../compare.sh 14_Si_st
../run.sh 14_Si_UPF -np	; 	../compare.sh 14_Si_UPF
../run.sh 17_Cl -np	; 	../compare.sh 17_Cl
../run.sh 19_K -np	; 	../compare.sh 19_K
../run.sh 19_K_st -np	; 	../compare.sh 19_K_st
../run.sh 20_Ca -np	; 	../compare.sh 20_Ca
../run.sh 22_Ti -np	; 	../compare.sh 22_Ti
../run.sh 29_Cu -np	; 	../compare.sh 29_Cu
../run.sh 32_Ge -np	; 	../compare.sh 32_Ge
../run.sh 34_Se -np	; 	../compare.sh 34_Se
../run.sh 38_Sr_lxc -np	; 	../compare.sh 38_Sr_lxc
../run.sh 40_Zr -np	; 	../compare.sh 40_Zr
../run_r.sh 52_Te -np	; 	../compare.sh 52_Te_r
../run.sh 56_Ba -np	; 	../compare.sh 56_Ba
../run.sh 57_La -np	; 	../compare.sh 57_La
../run.sh 60_Nd_GHOST -np	; 	../compare.sh 60_Nd_GHOST
../run.sh 73_Ta -np	; 	../compare.sh 73_Ta
../run_r.sh 74_W -np	; 	../compare.sh 74_W_r
../run.sh 79_Au_lxc -np	; 	../compare.sh 79_Au_lxc
# this also produces a psml example file for Hg
../run_r.sh 80_Hg -np	; 	../compare.sh 80_Hg_r; ../extract.sh 80_Hg_r $(pwd); ../compare_psml.sh 80_Hg_r $(pwd)
../run.sh 83_Bi -np	; 	../compare.sh 83_Bi

grep Summary *.diff | sed -e /diff/s/^/\\\n/ >TEST.report
