#!/usr/bin/python

#
from Numeric import *
from Scientific.IO.NetCDF import *
import biggles
import sys


fname=sys.argv[1]

of = NetCDFFile(fname)
ntb=of.dimensions['ntb']

nrows = 2
if of.variables.has_key('chcore'): nrows = nrows+1

t=biggles.Table(nrows,1)

delta=of.variables['vlocal'].Vlocal_delta[0]
r = arange(1,ntb+1)*float(delta)
p=biggles.FramedPlot()
p.title="Vlocal"
p.add(biggles.Curve(r,of.variables['vlocal'][:]))
t.set(0,0,p)

delta=of.variables['chlocal'].Chlocal_delta[0]
r = arange(1,ntb+1)*float(delta)
p=biggles.FramedPlot()
p.title="Chlocal"
p.add(biggles.Curve(r,of.variables['chlocal'][:]))
t.set(1,0,p)

if of.variables.has_key('chcore'):
	delta=of.variables['chcore'].Chcore_delta[0]
	r = arange(1,ntb+1)*float(delta)
	p=biggles.FramedPlot()
	p.title="Chcore"
	p.add(biggles.Curve(r,of.variables['chcore'][:]))
	t.set(2,0,p)

t.title="OTHER STUFF from " + fname
t.show()

t.save_as_img( "svg", 400, 400, fname+"kbs.svg" )
t.save_as_eps( fname+".kbs.eps" )
