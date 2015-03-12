#!/bin/sh -f
#
# add-siesta-vnl.sh -- Script to insert a Siesta-style 
#                      "pseudopotential-operator" section in a psml file
#
# Usage: add-siesta-vnl.sh [psop options] <name.psml>
#
psop_prog=${PSOP_PROGRAM}
#
if [ -z "${psop_prog}" ]
then
   echo "Need to define PSOP_PROGRAM"
   exit
fi
#
if [ "$#" != 1 ] 
then
	echo "Usage: $0 [psop options] <name.psml>"
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
${psop_prog} $@ $file 2>&1 >| $name.psop.log
#
# This will remove any existing VNL section in the file,
# remove the markup closing the psml element,
# insert the contents of VNL, and re-add the closing
# markup
#
rm -f _tmp
sed '/<pseudopotential-operator /,/<\/pseudopotential-operator>/d' $file | \
                               grep -v "<\/psml>" | cat - VNL > _tmp
echo "</psml>" | cat _tmp - > PSML_VNL
cp -p PSML_VNL ../${name}-siesta-vnl.psml
#
echo "==> Output data in directory $name"
echo "==> Pseudopotential with Siesta-style VNL in ${name}-siesta-vnl.psml"



