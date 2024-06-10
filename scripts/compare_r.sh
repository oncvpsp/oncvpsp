#!/bin/bash
#compares <prefix>.out in the tests/data directory to the corresponing
#file in test/refs
#command line argument 1 is <prefix>

#comparison starts after the "ATOM AND REFERENCE CONFIGURATION" line
#so that output with built-in exc functions and corresponding libxc 
#functions can be compared (libxc produces two extra lines which 
#stops fldiff)

PREFIX=/Users/mverstra/CODES/ONCVPSP/GITHUB_VERSION/oncvpsp

#HEAD=`echo $1 | sed s/.dat//`
HEAD=`echo $1 | sed s/.dat/_r/`


OUTFILE1=$PREFIX/tests/refs/$HEAD.out
#OUTFILE1=$PREFIX/tests/refs/$1.out

OUTFILE2=$PREFIX/tests/data/$HEAD.out
#OUTFILE2=$PREFIX/tests/data/$1.out

DIFFFILE=$HEAD.diff

TEMP1=$$.tmp1

TEMP2=$$.tmp2

awk 'BEGIN{out=0}; {if(out == 1) {print}};\
	/ATOM AND REFERENCE CONFIGURATION/{out=1}' $OUTFILE1 | \
	sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ >$TEMP1

awk 'BEGIN{out=0}; {if(out == 1) {print}};\
	/ATOM AND REFERENCE CONFIGURATION/{out=1}' $OUTFILE2 | \
	sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ >$TEMP2

#$PREFIX/scripts/fldiff.pl -easy $TEMP1 $TEMP2 >& $DIFFFILE
$PREFIX/scripts/fldiff.pl -medium $TEMP1 $TEMP2 >& $DIFFFILE

rm $TEMP1 $TEMP2
