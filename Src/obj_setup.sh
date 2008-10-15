#!/bin/sh
#
##set -x
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
user_specified_dir=$(dirname $0)
testdir=$(dirname $srcdir)/Tests
#
destdir=$(pwd)
#
# Replicate the hierarchy of makefiles
#
for i in wxml xmlparser MPI Libs fdf ; do
    mkdir $i
    cp ${srcdir}/$i/*akefile ${destdir}/$i
done
#
sed "s#VPATH=\.#VPATH=${srcdir}#g" ${srcdir}/Makefile > ${destdir}/Makefile

#
# Tests directory
# Create a list of files and use tar to process the list and copy the files
# to the destination directory
#
( cd ${testdir} ; cd .. ; find Tests  \
              -path *Reference -prune -o  \
              -path *Reference-xml -prune -o  \
              -path *work -prune      -o  \
              -path *.arch-ids  -prune -o -print \
              | tar -cf - --no-recursion -T- )   | ( cd ${destdir} ; tar xf -)
#
# Now make a symbolic link in the destination directory
#
ln -s ${destdir} ${destdir}/Src
#
echo " *** Compilation setup done. "
echo " *** Remember to copy an arch.make file or run configure as:"
echo "    ${user_specified_dir}/configure [configure_options]"
