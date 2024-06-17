# Makefile for ONCVPSP
#
# Copyright (c) 1989-2016 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
#

include ../make.inc

OBJS_ALL =	aeo.o check_data.o const_basis.o dpnint.o eresid.o exchdl.o \
		exc_off_pbe.o excpzca.o excwig.o ldiracfb.o ldiracfs.o \
		lschfb.o lschfs.o lschkb.o lschpb.o lschpse.o lschps.o \
		lschvkbbe.o lschvkbb.o lschvkbs.o modcore.o optimize.o \
		pspot.o qroots.o run_optimize.o sbf8.o sbf_basis.o \
		sbf_rc_der.o tfapot.o vkboutwf.o vkbphsft.o vout.o vpinteg.o \
		vploc.o wf_rc_der.o modcore2.o gg1cc.o gp1cc.o gpp1cc.o \
		der2exc.o modcore3.o sbf_basis_con.o run_ghosts.o \
		lschpsbar.o vrel.o

OBJS_NRSR =	fphsft.o fpovlp.o gnu_script.o linout.o psatom.o \
		relatom.o run_config.o run_diag.o run_phsft.o run_plot.o \
		run_vkb.o sratom.o upfout.o wellstate.o

OBJS_SR =	m_psmlout.o oncvpsp.o

OBJS_NR =	m_psmlout.o oncvpsp_nr.o

OBJS_R =	m_psmlout.o oncvpsp_r.o fphsft_r.o fpovlp.o gnu_script_r.o linout_r.o \
		psatom_r.o relatom.o renorm_r.o run_config_r.o run_diag_r.o \
		run_diag_sr_so_r.o run_phsft_r.o run_plot_r.o run_vkb_r.o \
		sr_so_r.o upfout_r.o wellstate_r.o

OBJS_OPSP =     $(OBJS_SR) $(OBJS_NRSR) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_OPSPNR = 	$(OBJS_NR) $(OBJS_NRSR) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_OPSPR = 	$(OBJS_R) $(OBJS_ALL) $(OBJS_LIBXC)

all:	oncvpsp oncvpspnr oncvpspr

WXML_LIB=xmlf90-wxml/libwxml.a
wxml:
	(cd xmlf90-wxml; $(MAKE) module)
#
m_psmlout.o: wxml m_libxc_list.o
INCLUDE_WXML=-I./xmlf90-wxml



oncvpsp:	$(OBJS_OPSP)
		$(FLINKER) -o oncvpsp.x $(OBJS_OPSP) $(LIBS) $(WXML_LIB)

oncvpspnr:	$(OBJS_OPSPNR)
		$(FLINKER) -o oncvpspnr.x $(OBJS_OPSPNR) $(LIBS) $(WXML_LIB)

oncvpspr:	$(OBJS_OPSPR)
		$(FLINKER) -o oncvpspr.x $(OBJS_OPSPR) $(LIBS) $(WXML_LIB)

.SUFFIXES: 
.SUFFIXES: .c .o .f .f90 .F90

.F90.o:
	$(FCCPP) $(FCCPPFLAGS) $*.F90  > $*_cpp.f90
	$(F90) $(FFLAGS) $(INC) $(INCLUDE_WXML) -o $*.o -c $*_cpp.f90

.f90.o:
	$(F90) $(FFLAGS) $(INC) $(INCLUDE_WXML) -c $*.f90

.f.o:
	$(F90) $(FFLAGS) $(INC) $(INCLUDE_WXML) -c $*.f

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $*.c

clean:
	 /bin/rm -f *.o  *.d  *.mod *.log *_cpp.f90

