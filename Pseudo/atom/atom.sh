#!/bin/sh
#
# atom -- Runs atom in a scratch area and copies back a tar file with
#         all the results.
#
# Usage:  atom input_file [ vps_file ]
#
dir=$PWD
echo $dir
#
input_file=$1
#
SCR="$TMP/$input_file"
mkdir $SCR
#
cp $input_file  $SCR/INP
#
if [ ! -z "$2" ]
then
	cp $2 $SCR/VPS
fi
#
cd $SCR
atm
tar cf $dir/$input_file.tar .
#
