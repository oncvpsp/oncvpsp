###########################################################
# ONCVPSP library
###########################################################
set(src_liboncvpsp
    aeo.f90
    check_data.f90
    const_basis.f90
    der2exc.f90
    dpnint.f90
    eresid.f90
    exc_off_pbe.f
    exchdl.f90
    excpzca.f90
    excwig.f90
    gg1cc.f90
    gp1cc.f90
    gpp1cc.f90
    ldiracfb.f90
    ldiracfs.f90
    lschfb.f90
    lschfs.f90
    lschkb.f90
    lschpb.f90
    lschps.f90
    lschpsbar.f90
    lschpse.f90
    lschvkbb.f90
    lschvkbbe.f90
    lschvkbs.f90
    m_libxc_list.F90
    m_psmlout.f90
    modcore1_m.f90
    modcore2_m.f90
    modcore3_m.F90
    optimize.f90
    pspot.f90
    qroots.f90
    run_ghosts.f90
    run_optimize.f90
    sbf_basis_con.f90
    sbf_basis.f90
    sbf_rc_der.f90
    sbf8.f90
    tfapot.f90
    vkboutwf.f90
    vkbphsft.f90
    vout.f90
    vpinteg.f90
    vploc.f90
    vrel.f90
    wf_rc_der.f90
    output_text_m.f90
    input_text_m.f90
)
if (Libxc_FOUND)
    list(APPEND src_liboncvpsp functionals.F90 exc_libxc.f90)
else()
    list(APPEND src_liboncvpsp exc_libxc_stub.f90)
endif()
if (WITH_HDF5 AND HDF5_FOUND)
    add_compile_definitions(WITH_HDF5)
    list(APPEND src_liboncvpsp hdf5_utils_m.f90 output_hdf5_m.f90)
    find_package(HDF5 COMPONENTS Fortran HL REQUIRED)
endif()
if (WITH_TOML)
    list(APPEND src_liboncvpsp input_toml_m.F90)
    add_compile_definitions(WITH_TOML)
endif()
add_library(liboncvpsp STATIC ${src_liboncvpsp})
add_dependencies(liboncvpsp wxml)
target_link_libraries(liboncvpsp wxml ${COMMON_LIBRARIES})

###########################################################
# ONCVPSP non- and scalar-relativistic library
###########################################################
set(src_liboncvpsp_nr_sr
    fphsft.f90
    fpovlp.f90
    gnu_script.f90
    linout.f90
    psatom.f90
    relatom.f90
    run_config.f90
    run_diag.f90
    run_phsft.f90
    run_vkb.f90
    sratom.f90
    upfout.f90
    wellstate.f90
    postprocess_m.f90
)
add_library(liboncvpsp_nr_sr STATIC ${src_liboncvpsp_nr_sr})
target_compile_definitions(liboncvpsp_nr_sr PUBLIC RELATIVISTIC=1)
add_dependencies(liboncvpsp_nr_sr liboncvpsp wxml)
target_link_libraries(liboncvpsp_nr_sr liboncvpsp wxml ${COMMON_LIBRARIES})

###########################################################
# oncvpsp.x (scalar-relativistic)
###########################################################
add_executable(oncvpspsr oncvpsp.F90)
target_compile_definitions(oncvpspsr PUBLIC RELATIVISTIC=1)
add_dependencies(oncvpspsr liboncvpsp_nr_sr)
target_link_libraries(oncvpspsr liboncvpsp_nr_sr)
set_target_properties(oncvpspsr
    PROPERTIES
        OUTPUT_NAME oncvpsp.x
)

##########################################################
# oncvpspnr.x (non-relativistic)
##########################################################
add_executable(oncvpspnr oncvpsp.F90)
target_compile_definitions(oncvpspnr PUBLIC RELATIVISTIC=0)
add_dependencies(oncvpspnr liboncvpsp_nr_sr)
target_link_libraries(oncvpspnr liboncvpsp_nr_sr)
set_target_properties(oncvpspnr
    PROPERTIES
        OUTPUT_NAME oncvpspnr.x
)

###########################################################
# ONCVPSP relativistic library
###########################################################
set(src_liboncvpsp_r
    fphsft_r.f90
    fpovlp.f90
    gnu_script_r.f90
    linout_r.f90
    psatom_r.f90
    relatom.f90
    renorm_r.f90
    run_config_r.f90
    run_diag_r.f90
    run_diag_sr_so_r.f90
    run_phsft_r.f90
    run_plot_r.f90
    run_vkb_r.f90
    sr_so_r.f90
    upfout_r.f90
    wellstate_r.f90
)
add_library(liboncvpsp_r STATIC ${src_liboncvpsp_r})
target_compile_definitions(liboncvpsp_r PUBLIC RELATIVISTIC=2)
add_dependencies(liboncvpsp_r liboncvpsp wxml)
target_link_libraries(liboncvpsp_r liboncvpsp wxml ${COMMON_LIBRARIES})

##########################################################
# oncvpspr.x (fully-relativistic)
##########################################################
add_executable(oncvpspr oncvpsp_r.F90)
target_compile_definitions(oncvpspr PUBLIC RELATIVISTIC=2)
add_dependencies(oncvpspr liboncvpsp_r)
target_link_libraries(oncvpspr liboncvpsp_r)
set_target_properties(oncvpspr
    PROPERTIES
        OUTPUT_NAME oncvpspr.x
)

###########################################################
# Subdirectories
###########################################################
add_subdirectory(xmlf90-wxml)
