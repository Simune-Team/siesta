#!/bin/sh
#
echo " ===> Generating module files from templates..."

KINDS=$1

rm -f *.uses Interfaces.f90

for kind in ${KINDS} ; do

echo "         USE MPI__r${kind}_V      ;  USE MPI__r${kind}_S" >> V_S.uses
echo "         USE MPI__c${kind}_V      ;  USE MPI__c${kind}_S" >> V_S.uses

echo "         USE MPI__r${kind}_VS      ;  USE MPI__r${kind}_SV" >> VS.uses
echo "         USE MPI__c${kind}_VS     ;  USE MPI__c${kind}_SV" >> VS.uses

done


for tag in v s sv vs ; do

sed -e "/_type/s//_integer/" -e "/type/s//integer/" \
    mpi__type_${tag}.f90 >> Interfaces.f90

for kind in ${KINDS} ; do
sed -e "/_type/s//_r${kind}/" -e "/type/s//real(${kind})/" \
    mpi__type_${tag}.f90 >> Interfaces.f90
sed -e "/_type/s//_c${kind}/" -e "/type/s//complex(${kind})/" \
    mpi__type_${tag}.f90 >> Interfaces.f90

done

sed -e "/_type/s//_logical/" -e "/type/s//logical/" \
    mpi__type_${tag}.f90 >> Interfaces.f90

sed -e "/_type/s//_character/" -e "/type/s//character(1)/" \
    mpi__type_${tag}.f90  >> Interfaces.f90

done











