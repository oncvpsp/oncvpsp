enable_testing()

set(oncvpspsr_tests
    "03_Li" "07_N" "07_N_dl" "08_O" "14_Si"
    "14_Si_st" "14_Si_UPF" "17_Cl" "19_K"
    "19_K_st" "20_Ca" "22_Ti" "29_Cu" "32_Ge"
    "34_Se" "38_Sr_lxc" "40_Zr" "56_Ba" "57_La"
    "60_Nd_GHOST" "73_Ta" "79_Au_lxc" "83_Bi"
)
foreach(test_name ${oncvpspsr_tests})
    set(test_dat "${CMAKE_SOURCE_DIR}/tests/data/${test_name}.dat")
    set(test_out "${CMAKE_SOURCE_DIR}/tests/data/${test_name}.out")
    set(ref_out "${CMAKE_SOURCE_DIR}/tests/refs/${test_name}.out")
    set(diff_out "${CMAKE_SOURCE_DIR}/tests/data/${test_name}.diff")
    add_test(
        NAME ${test_name}
        COMMAND bash -c "
            TEMP1=${test_name}.tmp1
            TEMP2=${test_name}.tmp2
            $<TARGET_FILE:oncvpspsr> < ${test_dat} > ${test_out}
            if [ \$? -ne 0 ]; then
                echo 'oncvpspsr execution failed'
                exit 1
            fi
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${ref_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ >$TEMP1
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${test_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ >$TEMP2
            perl ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl -easy $TEMP1 $TEMP2 > ${diff_out}
            rm -f $TEMP1 $TEMP2
            grep -E '^Summary  : no significant difference|^Summary  : different lines=\\s*\\d+ , max abs_diff= -?0\\.000e[\\+\\-]00 \\(.*\\), max rel_diff= -?0\\.000e[\\+\\-]00 \\(.*\\)' ${diff_out}
            success=\$?
            cat ${diff_out}
            exit \$success
        "
    )
    set_property(TEST ${test_name} PROPERTY LABEL "sr")
endforeach(test_name)

set(oncvpspr_tests
    "52_Te" "74_W" "80_Hg"
)
foreach(test_name ${oncvpspr_tests})
    set(test_dat "${CMAKE_SOURCE_DIR}/tests/data/${test_name}.dat")
    set(test_out "${CMAKE_SOURCE_DIR}/tests/data/${test_name}_r.out")
    set(ref_out "${CMAKE_SOURCE_DIR}/tests/refs/${test_name}_r.out")
    set(diff_out "${CMAKE_SOURCE_DIR}/tests/data/${test_name}_r.diff")
    add_test(
        NAME ${test_name}
        COMMAND bash -c "
            TEMP1=${test_name}.tmp1
            TEMP2=${test_name}.tmp2
            $<TARGET_FILE:oncvpspr> < ${test_dat} > ${test_out}
            if [ \$? -ne 0 ]; then
                echo 'oncvpspr execution failed'
                exit 1
            fi
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${ref_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ >$TEMP1
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${test_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ >$TEMP2
            perl ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl -easy $TEMP1 $TEMP2 > ${diff_out}
            rm -f $TEMP1 $TEMP2
            grep -E '^Summary  : no significant difference|^Summary  : different lines=\\s*[0-9]+ , max abs_diff= -?0\\.000e[\\+\\-]00 \\(.*\\), max rel_diff= -?0\\.000e[\\+\\-]00 \\(.*\\)' ${diff_out}
            success=\$?
            cat ${diff_out}
            exit \$success
        "
    )
    set_property(TEST ${test_name} PROPERTY LABEL "r")
endforeach(test_name)
