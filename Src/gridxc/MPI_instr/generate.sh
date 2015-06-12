#!/bin/sh
#
### set -x
dir=$(dirname $0)
#
echo " ===> Generating module files from templates..."

if [ -z "$@" ] ; then
    KINDS=`./kind_explorer`
else
    KINDS=$@
fi

echo $KINDS
rm -f *.uses Interfaces.f90

for kind in ${KINDS} ; do

echo "         USE MPI__GXC_r${kind}_V      ;  USE MPI__GXC_r${kind}_S" >> V_S.uses
echo "         USE MPI__GXC_c${kind}_V      ;  USE MPI__GXC_c${kind}_S" >> V_S.uses

echo "         USE MPI__GXC_r${kind}_VS      ;  USE MPI__GXC_r${kind}_SV" >> VS.uses
echo "         USE MPI__GXC_c${kind}_VS     ;  USE MPI__GXC_c${kind}_SV" >> VS.uses

done

#
#
for tag in v s sv vs ; do

echo "         USE MPI__GXC_integer_${tag}" >> int_logical_char.uses
echo "         USE MPI__GXC_integer8_${tag}" >> int_logical_char.uses
echo "         USE MPI__GXC_logical_${tag}" >> int_logical_char.uses
echo "         USE MPI__GXC_character_${tag}" >> int_logical_char.uses

done

#
#
for tag in v s sv vs ; do

sed -e "/_type/s//_GXC_integer/" -e "/type/s//integer/" \
    ${dir}/mpi__type_${tag}.f90 >> Interfaces.f90

sed -e "/_type/s//_GXC_integer8/" -e "/type/s//integer\*8/" \
    ${dir}/mpi__type_${tag}.f90 >> Interfaces.f90

for kind in ${KINDS} ; do
sed -e "/_type/s//_GXC_r${kind}/" -e "/type/s//real(${kind})/" \
    ${dir}/mpi__type_${tag}.f90 >> Interfaces.f90
sed -e "/_type/s//_GXC_c${kind}/" -e "/type/s//complex(${kind})/" \
    ${dir}/mpi__type_${tag}.f90 >> Interfaces.f90

done

sed -e "/_type/s//_GXC_logical/" -e "/type/s//logical/" \
    ${dir}/mpi__type_${tag}.f90 >> Interfaces.f90

sed -e "/_type/s//_GXC_character/" -e "/type/s//character(*)/" \
    ${dir}/mpi__type_${tag}.f90  >> Interfaces.f90

done
