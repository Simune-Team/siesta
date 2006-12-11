#!/bin/sh

#
# Select the appropriate run below (at the end).
#
# Take care to use the appropriate (e.g., parallel or serial) copy
# of siesta. To run parallel jobs, you have to figure out how to
# set the MPI environment. For a simple single-node interactive calculation
# (where possible, such as in a Rocks cluster), you can use the parallel.sh
# script.
#

ROOT="../../.."
PSEUDOS=${ROOT}/Tests/Pseudos
SRC=${ROOT}/Src
#

if [ -d work ] ; then
   echo "Work directory exists. Please delete it"
   exit
else
   mkdir work
fi

cp -p h2o.fast.fdf h2o.conv.fdf driver.dat work
#
cp ${PSEUDOS}/H.psf  work
cp ${PSEUDOS}/O.psf work
#
cd work
ln -s ../${SRC}/siesta ./siesta

../Src/driver < driver.dat | tee driver.out
##../Src/simple  | tee simple.out
##../Src/para  | tee para.out

