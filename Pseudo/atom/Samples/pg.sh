#!/bin/sh -f
#
# pg.sh -- Script to run pseudopotential generation calculations
#
# Usage: pg.sh <name.inp>
#
####set -x
#
prog="../../atm"
#
if [ "$#" != 1 ] 
then
	echo "Usage: $0 <name.inp>"
	exit
fi
#
file=$1
name=`basename $1 .inp`
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
$prog
cp VPSOUT ../$name.vps
#
echo "Calculation for $name completed. Output data in directory $name..."
echo ".... and pseudopotential in $name.vps"
#


