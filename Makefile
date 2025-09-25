# Makefile for ONCVPSP
#
# Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
#

MAKE = make

include make.inc

all:
	cd src ; $(MAKE) all
	./set_path
	cd tests/data ; ./TEST.sh

test:
	cd tests/data ; ./TEST.sh

clean:
	cd src ; $(MAKE) clean ; rm *.x; cd xmlf90-wxml/ ; $(MAKE) clean
	cd tests/data ; /bin/rm -f *.out *.diff TEST.report
