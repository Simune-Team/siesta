#!/bin/sh
#
# build all Utils
#
# Usage: [ OBJDIR=Objdir] sh build_all
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
topdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
#---------------------------------------------
#
if [ -z "$OBJDIR" ] ; then
    OBJDIR=Obj
fi
echo ${OBJDIR}
set -x
for i in $(find . -name \[mM\]akefile | grep -v \\./Makefile); do
      relpath=${i%/*}
      cd $relpath
      make OBJDIR=${OBJDIR} clean
      make OBJDIR=${OBJDIR} 
      cd $topdir
done
