###########################################################
# ONCVPSP library
###########################################################
set(src_oncvpsp
    # Common files
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
    modcore.f90
    modcore2.f90
    modcore3.f90
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
    input_text_m.f90
    # Non-/scalar-relativistic files
    fphsft.f90
    fpovlp.f90
    gnu_script.f90
    linout.f90
    psatom.f90
    relatom.f90
    run_config.f90
    run_diag.f90
    run_phsft.f90
    run_plot.f90
    run_vkb.f90
    sratom.f90
    upfout.f90
    wellstate.f90
    # Fully-relativistic files
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
# Append Libxc interface or stub
if (Libxc_FOUND AND ONCVPSP_WITH_Libxc)
    list(APPEND src_oncvpsp functionals.F90 exc_libxc.f90)
else ()
    list(APPEND src_oncvpsp exc_libxc_stub.f90)
endif ()
# Append TOML input parser if enabled
if (ONCVPSP_WITH_TOML)
    list(APPEND src_oncvpsp input_toml_m.f90)
endif ()

add_library(oncvpsp STATIC ${src_oncvpsp})
add_dependencies(oncvpsp wxml)
target_link_libraries(oncvpsp wxml ${ONCVPSP_DEPENDENCY_LIBRARIES})
if (ONCVPSP_WITH_TOML)
    target_compile_definitions(oncvpsp PRIVATE WITH_TOML)
endif ()

###########################################################
# oncvpsp.x (scalar-relativistic)
###########################################################
add_executable(oncvpspsr oncvpsp.F90)
target_compile_definitions(oncvpspsr PUBLIC RELATIVISTIC=1)
add_dependencies(oncvpspsr oncvpsp)
target_link_libraries(oncvpspsr PRIVATE oncvpsp)
set_target_properties(oncvpspsr
    PROPERTIES
        OUTPUT_NAME oncvpsp.x
)
if (ONCVPSP_WITH_TOML)
    target_compile_definitions(oncvpspsr PRIVATE WITH_TOML)
endif ()
if (HDF5_FOUND AND ONCVPSP_WITH_HDF5)
    target_compile_definitions(oncvpspsr PRIVATE WITH_HDF5)
endif ()

##########################################################
# oncvpspnr.x (non-relativistic)
##########################################################
add_executable(oncvpspnr oncvpsp.F90)
target_compile_definitions(oncvpspnr PUBLIC RELATIVISTIC=0)
add_dependencies(oncvpspnr oncvpsp)
target_link_libraries(oncvpspnr PRIVATE oncvpsp)
set_target_properties(oncvpspnr
    PROPERTIES
        OUTPUT_NAME oncvpspnr.x
)
if (ONCVPSP_WITH_TOML)
    target_compile_definitions(oncvpspnr PRIVATE WITH_TOML)
endif ()
if (HDF5_FOUND AND ONCVPSP_WITH_HDF5)
    target_compile_definitions(oncvpspnr PRIVATE WITH_HDF5)
endif ()

##########################################################
# oncvpspr.x (fully-relativistic)
##########################################################
add_executable(oncvpspr oncvpsp_r.F90)
target_compile_definitions(oncvpspr PUBLIC RELATIVISTIC=2)
add_dependencies(oncvpspr oncvpsp)
target_link_libraries(oncvpspr PRIVATE oncvpsp)
set_target_properties(oncvpspr
    PROPERTIES
        OUTPUT_NAME oncvpspr.x
)
if (ONCVPSP_WITH_TOML)
    target_compile_definitions(oncvpspr PRIVATE WITH_TOML)
endif ()
if (HDF5_FOUND AND ONCVPSP_WITH_HDF5)
    target_compile_definitions(oncvpspr PRIVATE WITH_HDF5)
endif ()

###########################################################
# Subdirectories
###########################################################
add_subdirectory(xmlf90-wxml)
