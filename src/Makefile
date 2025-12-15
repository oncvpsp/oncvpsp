# Makefile for ONCVPSP
#
# Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
#

F77        = gfortran	#EDIT
F90        = gfortran	#EDIT
CC         = gcc	#EDIT


FLINKER     = $(F90)

LIBS = -L/usr/local/atlas/lib -llapack -lf77blas -lcblas -latlas	#EDIT



FFLAGS     = -O3 -msse3 -m32	#EDIT
CFLAGS     = -O3		#EDIT

OBJS_OPSP   =	oncvpsp.o vkboutwf.o vkbphsft.o \
		const_basis.o dp3int.o eresid.o exchdl.o exc_off_pbe.o \
		excpwca.o excwig.o fphsft.o vpinteg.o \
		lschvkbb.o lschfb.o lschfs.o lschkb.o lschpb.o lschps.o \
		lschpse.o lschvkbs.o modcore.o optimize.o psatom.o \
		pspot.o qroots.o run_optimize.o sbf8.o sbf_basis.o \
		sbf_rc_der.o sratom.o vout.o wellstate.o wf_rc_der.o aeo.o \
		run_vkb.o run_phsft.o run_diag.o run_plot.o \
		linout.o run_config.o vploc.o gnu_script.o \
		check_data.o lschvkbbe.o tfapot.o fpovlp.o upfout.o

OBJS_OPSPNR   =	oncvpsp_nr.o vkboutwf.o vkbphsft.o \
		const_basis.o dp3int.o eresid.o exchdl.o exc_off_pbe.o \
		excpwca.o excwig.o fphsft.o vpinteg.o \
		lschvkbb.o lschfb.o lschfs.o lschkb.o lschpb.o lschps.o \
		lschpse.o lschvkbs.o modcore.o optimize.o psatom.o \
		pspot.o qroots.o run_optimize.o sbf8.o sbf_basis.o \
		sbf_rc_der.o sratom.o vout.o wellstate.o wf_rc_der.o aeo.o \
		run_vkb.o run_phsft.o run_diag.o run_plot.o \
		linout.o run_config.o vploc.o gnu_script.o \
		check_data.o lschvkbbe.o tfapot.o fpovlp.o upfout.o

OBJS_OPSPR   =	oncvpsp_r.o vkboutwf.o vkbphsft.o \
		const_basis.o dp3int.o eresid.o exchdl.o exc_off_pbe.o \
		excpwca.o excwig.o fphsft_r.o vpinteg.o \
		lschvkbb.o lschfb.o lschfs.o lschkb.o lschpb.o lschps.o \
		lschpse.o lschvkbs.o modcore.o optimize.o psatom_r.o \
		pspot.o qroots.o run_optimize.o sbf8.o sbf_basis.o \
		sbf_rc_der.o relatom.o vout.o wellstate_r.o wf_rc_der.o aeo.o \
		run_vkb_r.o run_config_r.o run_phsft_r.o run_diag_r.o \
		run_plot_r.o linout_r.o vploc.o gnu_script_r.o \
		check_data.o lschvkbbe.o tfapot.o upfout_r.o \
		fpovlp_r.o ldiracfb.o ldiracfs.o renorm_r.o sr_so_r.o \
		run_diag_sr_so_r.o

all:	oncvpsp oncvpspnr oncvpspr

oncvpsp:	$(OBJS_OPSP)
		$(FLINKER) -o oncvpsp.x $(OBJS_OPSP) $(LIBS)

oncvpspnr:	$(OBJS_OPSPNR)
		$(FLINKER) -o oncvpspnr.x $(OBJS_OPSPNR) $(LIBS)

oncvpspr:	$(OBJS_OPSPR)
		$(FLINKER) -o oncvpspr.x $(OBJS_OPSPR) $(LIBS)
.f90.o:
	$(F90) $(FFLAGS) $(INC) -c $*.f90

.f.o:
	$(F90) $(FFLAGS) $(INC) -c $*.f

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $*.c

.SUFFIXES: .c .o .f .f90

clean:
	 /bin/rm -f *.o  *.d  *.mod *.log

