#!/bin/sh
#
#
# Use template to generate appropriate headers
#
rm -f copyright.bang copyright.sharp

sed 's/%%/!/g' copyright > copyright.bang
sed 's/%%/\#/g' copyright > copyright.sharp

#
# This is all a bit kludgy:
# 
# 1. We do not process .f90 files, since those are currently only
#    in NetCDF, which is copyright by another author.
#
# 2. We make a copy of mpi.F and restore it later for the same reasons.
# 
# 3. We do not process Pseudo, only Src and Util.
#
#
rm -f MPISAVED
cp Src/MPI/mpi.F  MPISAVED
#
rm -f tmp.tmp

header=copyright.bang
#
list=`find Src Util -name '*.f'`
for i in $list ; do
    cat $header $i > tmp.tmp
    mv tmp.tmp $i
done
list=`find Src Util -name '*.F'`
for i in $list ; do
    cat $header $i > tmp.tmp
    mv tmp.tmp $i
done
list=`find Src Util -name '*.F90'`
for i in $list ; do
    cat $header $i > tmp.tmp
    mv tmp.tmp $i
done
list=`find Src Util -name '*.h'`
for i in $list ; do
    cat $header $i > tmp.tmp
    mv tmp.tmp $i
done
#
# restore mpi.F
#
mv -f MPISAVED Src/MPI/mpi.F

#
# Now change the initial character of the header...
#
header=copyright.sharp
#
list=`find Src Util -name '*akefile'`
for i in $list ; do
    cat $header $i > tmp.tmp
    mv tmp.tmp $i
done

list=`find Src Util -name '*.py'`
for i in $list ; do
    cat $header $i > tmp.tmp
    mv tmp.tmp $i
done










