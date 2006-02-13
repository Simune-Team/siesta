#!/usr/bin/python

#
from Numeric import *
from Scientific.IO.NetCDF import *

example="Ba.ion.nc"

class rfunc:
	def __init__(self,f,delta,cutoff):
		self.delta = delta
		self.cutoff = cutoff
		self.f = f

class orbital(rfunc):
  def __init__(self,f,delta,cutoff,n,l,z):
     rfunc.__init__(self,f,delta,cutoff)
     self.l = l
     self.z = z
     self.n = n
  def __str__(self):
     return ("orb%s_%s" % (self.n,self.l))
     
	
class atom:
  def __init__(self,filename):
    f=NetCDFFile(filename)
    self.element=f.Element
    vars = f.variables
    p = vars['vlocal']
    self.vlocal=rfunc(p[:],p.Vlocal_delta[0], p.Vlocal_cutoff[0])
    norbs = f.dimensions['norbs']
    self.base= []
    for i in range(norbs):
      p = vars['orb'][i]
      l = vars['orbnl_l'][i]
      n = vars['orbnl_n'][i]
      z = vars['orbnl_z'][i]
      delta = vars['delta'][i]
      cutoff = vars['cutoff'][i]
      orb=orbital(p[:],delta,cutoff,n,l,z)
      self.base.append(orb)
    del p, l, n, z, delta, cutoff, norbs, vars
    del f

a=atom(example)

# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
