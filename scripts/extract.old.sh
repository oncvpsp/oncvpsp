#!/bin/bash

#extracts the psp file from <prefix>.out file
#command line argument 1 is <prefix>
#commmand line argument 2 is the subdirectory of the user's home
#directory where the psp is to be written as <prefix>.oncvpsp.<type>
#whree <type> is psp8 or upf.

OUTFILE=$1.out

STR1=`grep PSPCODE8 $OUTFILE`
STR2=`grep PSP_UPF $OUTFILE`

if [ "$STR1" ]
	then
	PSPFILE=~/$2/$1.oncvpsp.psp8

	awk 'BEGIN{out=0};{if(out == 1) {print}};/PSPCODE8/{out=1}' \
		$OUTFILE >$PSPFILE

	echo "$PSPFILE written"
fi

if [ "$STR2" ]
	then
	PSPFILE=~/$2/$1.oncvpsp.upf

	awk 'BEGIN{out=0};{if(out == 1) {print}};/PSP_UPF/{out=1}' \
		$OUTFILE >$PSPFILE

	echo "$PSPFILE written"
fi
