#!/usr/bin/python

#
from Numeric import *
from Scientific.IO.NetCDF import *
import Gnuplot      
import sys


fname=sys.argv[1]

of = NetCDFFile(fname)
print of.dimensions
nkbs=of.dimensions['nkbs']
ntb=of.dimensions['ntb']
g=Gnuplot.Gnuplot()
g('set data style lines')
for i in range(nkbs):
        L = of.variables['pjnl_l'][i]
        N = of.variables['pjnl_n'][i]
        print "KB proj number" ,  i+1, "l=",L, "n=",N
	delta=of.variables['kbdelta'][i]
	r = arange(1,ntb+1)*float(delta)
	d = Gnuplot.Data(r,of.variables['proj'][i,:])
	g.plot(d)
	raw_input('Please press return to continue...\n')
#
