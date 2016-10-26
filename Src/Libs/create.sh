#!/bin/bash

cur_dir=$(pwd)

# First gather the lapack sources
python linalg2file.py -d $cur_dir/.. $cur_dir/../../Util/TS/TBtrans \
       -l ~/LA/lapack/SRC ~/LA/lapack/SRC/DEPRECATED ~/LA/lapack/INSTALL \
       --list-add-file lapack_add.files \
       --list-add-routine lapack_add.routines \
       --list-remove-file lapack_remove.files \
       --list-remove-routine lapack_remove.routines \
       -o lapack.f


# Now we need to locate the lapack sources as well
python linalg2file.py -d $cur_dir/.. $cur_dir/../../Util/TS/TBtrans \
       -f lapack.f \
       -l ~/LA/lapack/BLAS/SRC \
       --list-add-file blas_add.files -o blas.f

