#!/bin/sh
#
# Creates the 'parallel' distribution of siesta.
#
#  This script should be run on an "EXPORTED" directory:
#
#  cvs export { -r <appropriate tag>, -d <date> } [ -d directory ] siesta
#

dir=siesta-1.3p

if [ -d $dir ] ; then
      echo "Directory exists"
      exit
fi

mkdir $dir
echo "Copying Pseudo Tutorials Examples Docs Util Src ..."
cp -rp Pseudo Tutorials Examples Docs Util Src $dir
cp -rp README $dir
rm -rf $dir/Docs/Tech
rm -rf $dir/Src/Include

echo ""
echo "**** Remember to include the .ps.gz User guide in 2up form"





