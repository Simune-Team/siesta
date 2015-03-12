#!/bin/sh -f
#
# pg.sh -- Script to run pseudopotential generation calculations
#          followed by KB generation and generation of
#          a PSML file
#
# Usage: pgvnl.sh [psop-options] <name.inp>
#
prog=${ATOM_PROGRAM}
psop_prog=${PSOP_PROGRAM}
#
if [ -z "$prog" ]
then
   echo "Need to define ATOM_PROGRAM"
   exit
fi
if [ -z "${psop_prog}" ]
then
   echo "Need to define PSOP_PROGRAM"
   exit
fi
#
if [ "$#" == "0" ] 
then
	echo "Usage: $0 <name.inp> [psop_options]"
	echo "This is a (partial) list of options for psop:"
	${psop_prog} -h
	exit
fi
#
file=$1
echo "Processing file: $file"
shift   # Rest of arguments are for psop
#
name=`basename $file .inp`
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
cp -p VPSFMT $name.psf
${psop_prog} $@ $name.psf 2>&1 >| $name.psop.log
#
# This will remove any existing VNL section in the file,
# remove the markup closing the psml element,
# insert the contents of VNL, and re-add the closing
# markup
#
rm -f _tmp
sed '/<pseudopotential-operator /,/<\/pseudopotential-operator>/d' PSML | \
                               grep -v "<\/psml>" | cat - VNL > _tmp
echo "</psml>" | cat _tmp - > PSML_VNL
cp -p PSML_VNL ../${name}-vnl.psml
#
echo "==> Output data in directory $name"
echo "==> Pseudopotential with VNL in ${name}-vnl.psml"



