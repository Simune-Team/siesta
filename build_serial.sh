#!/bin/sh
#
# Creates the 'serial' distribution of siesta. Needs python.
#
dir=siesta_serial

if [ -d $dir ] ; then
      echo "Directory exists"
      exit
fi

mkdir $dir
cp -rp Pseudo Examples Docs Util $dir
rm -f $dir/Docs/CHANGES
mkdir $dir/Src

cd Src
cp -rp fdf Sys NetCDF Libs ../$dir/Src
cp -rp makefile.serial ../$dir/Src/Makefile

for i in *.F *.F90; do
  python ../dempi.py $i >| ../$dir/Src/$i
done
for i in *.f ; do
  cp -p $i ../$dir/Src
done
