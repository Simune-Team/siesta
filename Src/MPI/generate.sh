#!/bin/sh
#
echo " ===> Generating module files from templates..."

for tag in v s sv vs ; do

sed -e "/_type/s//_integer/" -e "/type/s//integer/" \
    mpi__type_${tag}.f90 > mpi__integer_${tag}.f90

sed -e "/_type/s//_single/" -e "/type/s//real*4/" \
    mpi__type_${tag}.f90 > mpi__single_${tag}.f90
sed -e "/_type/s//_double/" -e "/type/s//real*8/" \
    mpi__type_${tag}.f90 > mpi__double_${tag}.f90

sed -e "/_type/s//_complex/" -e "/type/s//complex*8/" \
    mpi__type_${tag}.f90 > mpi__complex_${tag}.f90
sed -e "/_type/s//_dcomplex/" -e "/type/s//complex*16/" \
    mpi__type_${tag}.f90 > mpi__dcomplex_${tag}.f90

sed -e "/_type/s//_logical/" -e "/type/s//logical/" \
    mpi__type_${tag}.f90 > mpi__logical_${tag}.f90

sed -e "/_type/s//_character/" -e "/type/s//character(1)/" \
    mpi__type_${tag}.f90 > mpi__character_${tag}.f90

done











