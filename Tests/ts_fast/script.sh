#!/bin/sh
#
BASE_DIR=`pwd`
TS="${BASE_DIR}"/../../../"$1"/transiesta
echo "Running script with TranSIESTA=$TS"


#
# Start with the electrode calculation
#
echo "Electrode Calculation"
mkdir Elec
cd Elec
cp ../../H.psf .
cp ../../elec.fast.fdf .
$TS < elec.fast.fdf > elec.fast.out 2>../../err_elec.out
cp elec.fast.out ../../ts_fast_elec.out
#
# Check if no error occured
#
if (( `wc -l < ../../err_elec.out` > "0" ))
then
   echo "The electrode calculation did not go well ..."
   exit
fi
rm ../../err_elec.out
#
# Go back to base directory
#
cd ..

#
# Scattering region calculation
#
echo "Scattering Region Calculation"
mkdir Scat
cd Scat
cp ../../scat.fast.fdf .
cp ../../H.psf .
# Copy the electrode's .TSHS
cp ../Elec/elec.fast.TSHS .
$TS < scat.fast.fdf > scat.fast.out 2>../../err_scat.out
cp scat.fast.out ../../ts_fast_scat.out
#
# Check if no error occured
#
if (( `wc -l < ../../err_scat.out` > "0" ))
then
   echo "The scattering region calculation did not go well ..."
   exit
fi
rm ../../err_scat.out
#
# Go back to base directory
#
cd ..

#
# TBTrans calculation
#
echo "TBTrans Calculation"
#
TBT="${BASE_DIR}"/../../../Util/TBTrans/tbtrans
if [ ! -x $TBT ] ; then
  echo "Compiling $TBT..."
  (cd "${BASE_DIR}"/../../../Util/TBTrans ; make OBJDIR="$1")
fi
echo "Running script with tbtrans=$TBT"
mkdir TBT
cd TBT
#Copy input files
cp ../Elec/elec.fast.TSHS .
cp ../Scat/scat.fast.TSHS .
cp ../../scat.fast.fdf .
$TBT < scat.fast.fdf > tbt.out 2>../../err_tbt.out
cp tbt.out ../../ts_fast_tbt.out
#
# Check if no error occured
#
if (( `wc -l < ../../err_tbt.out` > "0" ))
then
   echo "The tbtrans calculation did not go well ..."
   exit
fi
rm ../../err_tbt.out
#
# Go back to base directory
#
cd ..



# If it gets here it's because it finished without error
touch ../completed
