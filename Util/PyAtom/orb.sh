#!/bin/sh
#
# Driver for GNUPLOT plotting of orbitals
#
GFILE=orb.gplot

rm -f $GFILE

for i in *.ORB.?? ; do

line=`head -1 $i`

#
#  read does not seem to work right after a pipe, hence
#  the dummy_file kludge
#
rm -f dummy_file
echo $line > dummy_file
read dum L N Z Pol Pop < dummy_file
rm -f dummy_file
#
#
#
legend="Name L=$L N=$N Z=$Z Pop=$Pop"
cat >> $GFILE << EOF

plot "$i" title "$legend" with lines 
pause -1 "L, N, Z, is_pol, Pop: $line ...| Press return"

EOF
#
done

gnuplot orb.gplot
