#!/usr/bin/python

#
from Numeric import *
from Scientific.IO.NetCDF import *
import Gnuplot
import sys


fname=sys.argv[1]

of = NetCDFFile(fname)
print of.dimensions
norbs=of.dimensions['norbs']
ntb=of.dimensions['ntb']
g=Gnuplot.Gnuplot()
g('set data style lines')
for i in range(norbs):
        L = of.variables['orbnl_l'][i]
        Z = of.variables['orbnl_z'][i]
        print "Orb number" ,  i+1, "l=",L, "z=",Z
        delta=of.variables['delta'][i]
        r = arange(1,ntb+1)*float(delta)
        d = Gnuplot.Data(r,of.variables['orb'][i,:])
        g.plot(d)
        raw_input('Please press return to continue...\n')
#
