#!/bin/bash
#runs ONCVPSP with the command-line argument <prefix> and the graphics
#which review the results
#uses the non-relativistic all-electron atom calculation

PREFIX=@PROJECT_BINARY_DIR@
BIN_DIR=@PROJECT_BINARY_DIR@/bin

INFILE=$1.dat

OUTFILE=$1_nr.out

GNUFILE=$$.scr

PLOTFILE=$1_nr.plot

TEMP=$$.tmp

# PSML output: capture ONCVPSPPSML and rename to a per-element filename.
# $3 $4 $5 $6 $7 pass through PSML CLI options (-w, -c, -r ...)
PSMLFILE=$1_nr.psml

$BIN_DIR/oncvpspnr.x $3 $4 $5 $6 $7 <$INFILE >$OUTFILE

grep GHOST $OUTFILE

awk 'BEGIN{out=0};/GNUSCRIPT/{out=0}; {if(out == 1) {print}};\
	/DATA FOR PLOTTING/{out=1}' $OUTFILE >$PLOTFILE

awk 'BEGIN{out=0};/END_GNU/{out=0}; {if(out == 1) {print}};\
	/GNUSCRIPT/{out=1}' $OUTFILE >$TEMP

sed -e s/t1/$PLOTFILE/ $TEMP | sed -e s/t2/$1_nr/ >$GNUFILE

mv ONCVPSPPSML $PSMLFILE

if [ "$2" != "-np" ]
then
	gnuplot $GNUFILE
fi

rm  $GNUFILE $TEMP $PLOTFILE
