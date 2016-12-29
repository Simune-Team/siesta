#!/bin/bash

# This should be altered in case your
# lapack source is somewhere else
la_dir=$HOME/LA/lapack

cur_dir=$(pwd)

# First gather the lapack sources
python linalg2file.py -d $cur_dir/.. $cur_dir/../../Util/TS/TBtrans \
       -l $la_dir/SRC $la_dir/SRC/DEPRECATED $la_dir/INSTALL \
       --list-add-file lapack_add.files \
       --list-add-routine lapack_add.routines \
       --list-remove-file lapack_remove.files \
       --list-remove-routine lapack_remove.routines \
       -o lapack_tmp.f

# Append the license to the sources
cat $la_dir/LICENSE | sed -e 's/^/!/' > lapack_license_tmp
cat lapack_license_tmp lapack_tmp.f > lapack.F
rm lapack_license_tmp lapack_tmp.f


# Now we need to locate the lapack sources as well
python linalg2file.py -d $cur_dir/.. $cur_dir/../../Util/TS/TBtrans \
       -f lapack.F \
       -l $la_dir/BLAS/SRC \
       --list-add-file blas_add.files -o blas_tmp.f

# Append the license to the sources
cat $la_dir/LICENSE | sed -e 's/^/!/' > blas_license_tmp
cat blas_license_tmp blas_tmp.f > blas.F
rm blas_license_tmp blas_tmp.f
