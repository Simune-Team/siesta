#!/bin/sh
#
# Sets up plotting dir after a gen-basis o siesta basis generation
# calculation
#
# Set this or use environment variable...
#
DEFAULT_DIR=${HOME}/beta/siesta/Util/Basis
#
SIESTA_BASIS_UTILS_DIR=${SIESTA_BASIS_UTILS_DIR:-${DEFAULT_DIR}}
#
if [ "$#" != 1 ] 
then
      echo "Usage: sh setup_plot.sh atom_label"
      exit
fi
#
label=$1
#
dirname=$label.plotdir
#
if [ -d $dirname ] 
then
        echo "Directory $dirname exists. Please delete it first"
        exit
fi
#
mkdir $dirname 
#
#--------- to be redesigned...
#
for i in ORB.S?.?.$label KB.L?.?.$label \
         CHLOCAL.$label RED_VLOCAL.$label  VNA.$label ; do
name=`basename $i .$label`
cp $i $dirname/$name
done

if [ -f CHCORE.$label ]
then
      cp CHCORE.$label $dirname/CHCORE
fi

#
#  Copy plotting scripts
#
echo "Trying to copy scripts from  ${SIESTA_BASIS_UTILS_DIR} ..."
cp -f ${SIESTA_BASIS_UTILS_DIR}/*.gplot $dirname
#
#
#
echo " "
echo "The files necessary to plot the Basis information for $label"
echo "are now in directory $dirname. Go there now and use the Gnuplot"
echo "scripts orbs.gplot, kbs.gplot, and vdens.gplot to visualize"
echo "the PAOs, KB projectors, and other interesting functions."
echo " "
echo " *** (Be sure to use """gnuplot -persist""" --- an alias might be handy)"


