#!/bin/sh
#invoking this script with <prefix> as the command-line argument in a
#directory with a corresponding <prefix>.out file lets you review the
#graphics that you saw when you ran the calculation.

INFILE=$1.dat

OUTFILE=$1.out

GNUFILE=$$.scr

PLOTFILE=$1.plot

TEMP=$$.tmp

awk 'BEGIN{out=0};/GNUSCRIPT/{out=0}; {if(out == 1) {print}};\
	/DATA FOR PLOTTING/{out=1}' $OUTFILE >$PLOTFILE

awk 'BEGIN{out=0};/END_GNU/{out=0}; {if(out == 1) {print}};\
	/GNUSCRIPT/{out=1}' $OUTFILE >$TEMP

sed -e s/t1/$PLOTFILE/ $TEMP | sed -e s/t2/$1/ >$GNUFILE

gnuplot $GNUFILE

rm  $GNUFILE $TEMP $PLOTFILE
