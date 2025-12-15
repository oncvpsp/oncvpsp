###########################################################
# ONCVPSP library
###########################################################
set(src_oncvpsp
    aeo.f90
    check_data.f90
    const_basis.f90
    dp3int.f90
    eresid.f90
    exc_off_pbe.f
    exchdl.f90
    excpwca.f90
    excwig.f90
    fphsft.f90
    fphsft_r.f90
    fpovlp.f90
    fpovlp_r.f90
    gnu_script.f90
    gnu_script_r.f90
    ldiracfb.f90
    ldiracfs.f90
    linout.f90
    linout_r.f90
    lschfb.f90
    lschfs.f90
    lschkb.f90
    lschpb.f90
    lschps.f90
    lschpse.f90
    lschvkbb.f90
    lschvkbbe.f90
    lschvkbs.f90
    modcore.f90
    optimize.f90
    psatom.f90
    psatom_r.f90
    pspot.f90
    qroots.f90
    relatom.f90
    renorm_r.f90
    run_config.f90
    run_config_r.f90
    run_diag.f90
    run_diag_r.f90
    run_diag_sr_so_r.f90
    run_optimize.f90
    run_phsft.f90
    run_phsft_r.f90
    run_plot.f90
    run_plot_r.f90
    run_vkb.f90
    run_vkb_r.f90
    sbf8.f90
    sbf_basis.f90
    sbf_rc_der.f90
    sr_so_r.f90
    sratom.f90
    tfapot.f90
    upfout.f90
    upfout_r.f90
    vkboutwf.f90
    vkbphsft.f90
    vout.f90
    vpinteg.f90
    vploc.f90
    wellstate.f90
    wellstate_r.f90
    wf_rc_der.f90
)
add_library(oncvpsp STATIC ${src_oncvpsp})
target_link_libraries(oncvpsp ${ONCVPSP_DEPENDENCY_LIBRARIES})

###########################################################
# oncvpsp.x
###########################################################
add_executable(oncvpspsr oncvpsp.f90)
target_link_libraries(oncvpspsr PRIVATE oncvpsp)
set_target_properties(oncvpspsr
    PROPERTIES
        OUTPUT_NAME oncvpsp.x
)

##########################################################
# oncvpspnr.x
##########################################################
add_executable(oncvpspnr oncvpsp_nr.f90)
target_link_libraries(oncvpspnr PRIVATE oncvpsp)
set_target_properties(oncvpspnr
    PROPERTIES
        OUTPUT_NAME oncvpspnr.x
)

##########################################################
# oncvpspr.x
##########################################################
add_executable(oncvpspfr oncvpsp_r.f90)
target_link_libraries(oncvpspfr PRIVATE oncvpsp)
set_target_properties(oncvpspfr
    PROPERTIES
        OUTPUT_NAME oncvpspr.x
)
