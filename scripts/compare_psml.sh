#!/bin/bash
#compares <prefix>.out in the tests/data directory to the corresponing
#file in test/refs
#command line argument 1 is <prefix>

#comparison starts after the "ATOM AND REFERENCE CONFIGURATION" line
#so that output with built-in exc functions and corresponding libxc 
#functions can be compared (libxc produces two extra lines which 
#stops fldiff)

PREFIX=/Users/mverstra/CODES/ONCVPSP/GITHUB_VERSION/oncvpsp

OUTFILE1=$PREFIX/tests/refs/$1.oncvpsp.psml

OUTFILE2=$PREFIX/tests/data/$1.oncvpsp.psml

DIFFFILE=$1.oncvpsp.diff

$PREFIX/scripts/fldiff.pl -easy $OUTFILE1 $OUTFILE2 >& $DIFFFILE

