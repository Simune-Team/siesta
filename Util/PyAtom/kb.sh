#!/bin/sh
#
# Driver for GNUPLOT plotting of orbitals
#
GFILE=kb.gplot

rm -f $GFILE

for i in *.KB.?? ; do

line=`head -1 $i`

#
#  read does not seem to work right after a pipe, hence
#  the dummy_file kludge
#
rm -f dummy_file
echo $line > dummy_file
read dum L N Ref < dummy_file
rm -f dummy_file
#
#
#
legend="Name L=$L N=$N Ref e=$Ref"
cat >> $GFILE << EOF

plot "$i" title "$legend" with lines 
pause -1 "L, N, Ref energy: $line ...| Press return"

EOF
#
done

gnuplot kb.gplot
