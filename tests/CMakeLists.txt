enable_testing()

# Configure the interactive helper scripts (run.sh, run_r.sh, run_nr.sh,
# compare.sh) so their @PROJECT_BINARY_DIR@ placeholder resolves to this
# build tree, and they can be used standalone from the build directory.
configure_file(${CMAKE_SOURCE_DIR}/scripts/run.sh
                ${PROJECT_BINARY_DIR}/scripts/run.sh @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/run_r.sh
                ${PROJECT_BINARY_DIR}/scripts/run_r.sh @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/run_nr.sh
                ${PROJECT_BINARY_DIR}/scripts/run_nr.sh @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/scripts/compare.sh
                ${PROJECT_BINARY_DIR}/scripts/compare.sh @ONLY)
file(COPY ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl
     DESTINATION ${PROJECT_BINARY_DIR}/scripts)
file(COPY ${CMAKE_SOURCE_DIR}/scripts/extract.sh
     DESTINATION ${PROJECT_BINARY_DIR}/scripts)

# Each test is run in its own working directory: the Fortran binaries
# unconditionally write fixed relative filenames (INPUT_FILE, INP_COPY,
# ONCVPSPPSML) in the current directory (see m_psmlout.f90 and the
# oncvpsp*.F90 main programs), so tests running in a shared directory
# clobber each other's files under parallel ctest (-j).
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
    set(test_run_dir "${CMAKE_CURRENT_BINARY_DIR}/run/${test_name}")
    file(MAKE_DIRECTORY ${test_run_dir})
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
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ | grep -v 'lmax,lloc,mxprj' >$TEMP1
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${test_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ | grep -v 'lmax,lloc,mxprj' >$TEMP2
            perl ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl -easy $TEMP1 $TEMP2 > ${diff_out}
            rm -f $TEMP1 $TEMP2
            grep -E 'Summary  : no significant difference|Summary  : different lines=\\s*[0-9]+ , max abs_diff= -?0\\.000e[\\+\\-]00 \\(.*\\), max rel_diff= -?0\\.000e[\\+\\-]00 \\(.*\\)' ${diff_out}
            success=\$?
            cat ${diff_out}
            exit \$success
        "
        WORKING_DIRECTORY ${test_run_dir}
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
    set(test_run_dir "${CMAKE_CURRENT_BINARY_DIR}/run/${test_name}_r")
    file(MAKE_DIRECTORY ${test_run_dir})

    # 80_Hg additionally gets its PSML output checked against a reference
    # PSML file, exercising the relativistic PSML code path end to end.
    # fldiff.pl -easy (same tool and pass criteria as the .out comparison
    # above) tolerates the uuid/date attributes that differ on every run,
    # while still catching real content/structural differences -- this is
    # the same approach scripts/compare_psml.sh already uses.
    set(psml_check "")
    if(test_name STREQUAL "80_Hg")
        set(psml_ref "${CMAKE_SOURCE_DIR}/tests/refs/80_Hg_r.oncvpsp.psml")
        set(psml_out "${CMAKE_SOURCE_DIR}/tests/data/80_Hg_r.oncvpsp.psml")
        set(psml_diff_out "${CMAKE_SOURCE_DIR}/tests/data/80_Hg_r.oncvpsp.diff")
        set(psml_check "
            mv ONCVPSPPSML ${psml_out}
            if [ \$? -ne 0 ]; then
                echo '80_Hg_r PSML output missing'
                exit 1
            fi
            perl ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl -easy ${psml_ref} ${psml_out} > ${psml_diff_out}
            grep -E 'Summary  : no significant difference|Summary  : different lines=\\s*[0-9]+ , max abs_diff= -?0\\.000e[\\+\\-]00 \\(.*\\), max rel_diff= -?0\\.000e[\\+\\-]00 \\(.*\\)' ${psml_diff_out}
            psml_success=\$?
            cat ${psml_diff_out}
            exit \$psml_success
        ")
    endif()

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
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ | grep -v 'lmax,lloc,mxprj' >$TEMP1
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${test_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ | grep -v 'lmax,lloc,mxprj' >$TEMP2
            perl ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl -easy $TEMP1 $TEMP2 > ${diff_out}
            rm -f $TEMP1 $TEMP2
            grep -E 'Summary  : no significant difference|Summary  : different lines=\\s*[0-9]+ , max abs_diff= -?0\\.000e[\\+\\-]00 \\(.*\\), max rel_diff= -?0\\.000e[\\+\\-]00 \\(.*\\)' ${diff_out}
            success=\$?
            cat ${diff_out}
            if [ \$success -ne 0 ]; then
                exit \$success
            fi
            ${psml_check}
        "
        WORKING_DIRECTORY ${test_run_dir}
    )
    set_property(TEST ${test_name} PROPERTY LABEL "r")
endforeach(test_name)
