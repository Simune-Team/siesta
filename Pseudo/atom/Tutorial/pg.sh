#!/bin/sh -f
#
# pg.sh -- Script to run pseudopotential generation calculations
#
# Usage: pg.sh <name.inp>
#
DEFAULT_DIR=../../
ATOM_UTILS_DIR=${ATOM_UTILS_DIR:-${DEFAULT_DIR}}
#
default="../../../atm"
prog=${ATOM_PROGRAM:-$default}
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
#
$prog
#
cp VPSOUT ../$name.vps
cp VPSFMT ../$name.psf
#
echo "==> Output data in directory $name"
echo "==> Pseudopotential in $name.vps and $name.psf"
#
#  Copy plotting scripts
#
for i in charge vcharge vspin coreq pots pseudo scrpots subps ; do
            cp -f ${ATOM_UTILS_DIR}/$i.gp* .
done



