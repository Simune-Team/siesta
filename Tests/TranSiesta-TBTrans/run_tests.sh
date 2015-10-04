#!/bin/sh
#
# Run TS tests
#
#  To run in serial mode, replace the 'mpirun' line
#  by the appropriate incantation.
#
for d in ts_*; do
 cd $d
 make clean
 make TS="mpirun -np 4 ../../../../transiesta"
 cd ..
done
