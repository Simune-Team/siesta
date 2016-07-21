#!/bin/sh

# Run TS tests

# To run in serial mode, replace the 'mpirun' line
# by the appropriate incantation.

# Here we run the tests 
for d in ts_au \
	     ts_au_100_repetition_0.25V \
	     ts_graphene \
	     ts_term3 \
	     ts_term4 \
	     ts_au_repetition
do
    cd $d
    #make clean
    make TS="mpirun -np 4 ../../../../transiesta"
    cd ..
done
