#!/bin/sh -f
#
# add-siesta-vnl.sh -- Script to insert a Siesta-style 
#                      "pseudopotential-operator" section in a psml file
#
# Usage: add-siesta-vnl.sh  <name.psml> [psop options]
#
psop_prog=${PSOP_PROGRAM}
#
if [ -z "${psop_prog}" ]
then
   echo "Need to define PSOP_PROGRAM"
   exit
fi
#
if [ "$#" -lt 1 ] 
then
	echo "Usage: $0 <name.psml> [psop options]"
	echo "This is a (partial) list of options for psop:"
	${psop_prog} -h
	exit
fi
#
file=$1
echo "Processing file: $file"
shift   # Rest of arguments are for psop
#
name=`basename $file .psml`
#
if [ "${name}.psml" != "$file" ]
then
	echo "The file must be a .psml file!"
	exit
fi
#
if [ -d $name ] 
then
	echo "Directory $name exists. Please delete it first"
	exit
fi
#
mkdir $name ; cd $name
cp ../$file .
#
# psop will generate PSML_BASE (with a new provenance
# record and without any local and nl parts) and VNL,
# with the new local and nl sections.
#
${psop_prog} $@ $file 2>&1 >| $name.psop.log
#
# This will remove the markup closing the psml element,
# insert the contents of VNL into the base psml
# file, and re-add the closing markup, also filtering
# the artificial wrapper in VNL.
#
rm -f _tmp
grep -v "<\/psml>" PSML_BASE | cat - VNL > _tmp
echo "</psml>" | cat _tmp - | grep -v 'tmp-wrapper' > PSML_VNL
cp -p PSML_VNL ../${name}-siesta-vnl.psml
#
echo "==> Output data in directory $name"
echo "==> Pseudopotential with Siesta-style VNL in ${name}-siesta-vnl.psml"



