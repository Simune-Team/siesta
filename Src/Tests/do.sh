#!/bin/sh
#
name=$1
cd $name
#
if [ -f completed ] ; then echo " ---> skipping test $name"; exit; fi
if [ -d work ] ; then rm -rf work ; fi; mkdir work
#
for i in `cat $name.pseudos` ; do
     echo "    ==> Copying pseudopotential file for $i..."
      cp ../Pseudos/$i.psf work/$i.psf
done
#
echo ; echo "    ==> Running SIESTA as ${SIESTA}"
#
cd work
#
if (${SIESTA} 2>&1 > $name.out < ../$name.fdf) ; then
       cp $name.out .. ;
       cp SIG  ../$name.sig;
       touch ../completed;
       echo "    ===> SIESTA finished successfully";
else
      echo " **** Test $name did not complete successfully";
fi
      echo "    -------------------------------------"

   

