#!/bin/sh
#
gnubands < Fe.bands > bands.data

echo "Now go to the SGI window and type"
echo "gnuplot"
echo "gnuplot> plot 'bands.data'"
echo 
