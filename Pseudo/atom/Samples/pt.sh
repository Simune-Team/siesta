#!/bin/sh -f
#
# pt.sh -- Script to run pseudopotential test calculations
#
# Usage: pt.sh <ptname.inp> <psname.vps>
#
#
prog="../../atm"
#
if [ "$#" != 2 ] 
then
	echo "Usage: $0 <ptname.inp> <psname.vps>"
	exit
fi
#
file=$1
psfile=$2
ptname=`basename $file .inp`
psname=`basename $psfile .vps`
name="$ptname-$psname"
#
#
if [ -d $name ] 
then
	echo "Directory $name exists. Please delete it first"
	exit
fi
#
mkdir $name ; cd $name
cp ../$file ./INP
cp ../$psfile ./VPSIN
$prog
#
echo "Calculation for $name completed. Output data in directory $name..."
#


