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
echo "Copying Pseudo Examples Docs Util..."
cp -rp Pseudo Examples Docs Util $dir
rm -f $dir/Docs/CHANGES
mkdir $dir/Src

cd Src
echo "Copying fdf Sys NetCDF and Libs to Src..."
cp -rp fdf Sys NetCDF Libs ../$dir/Src
cp -rp Makefile ../$dir/Src

echo "De-MPI'ng .F and .F90 files..."
for i in *.F *.F90; do
  python ../dempi.py $i >| ../$dir/Src/$i
done
echo "Copying serial files..."
for i in *.f *.f90 ; do
  cp -p $i ../$dir/Src
done
