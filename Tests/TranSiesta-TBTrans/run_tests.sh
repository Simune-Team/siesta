#!/bin/bash

# Run TS tests
MPI=${MPI:-mpirun -np 4}

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
    make MPI="$MPI"
    cd ..
done
