# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=intel-nolibs
#
# Serial compilation without the need of any installed libraries.
# (That would be too simple: ifc needs the -Vaxlib option at link time)
#
# Following line for Pentium 4 , with static alloc
#FC=ifc  -tpp7 -O2  -Vaxlib -static -xW
#
# More conservative
FC=ifc  -tpp5 -O2 -w -mp -Vaxlib 
#
FC_ASIS=$(FC)
#
FFLAGS= -O
FFLAGS_DEBUG= -g 
LDFLAGS=-Vaxlib 
RANLIB=echo
LIBS=  
SYS=bsd
FPPFLAGS=
COMP_LIBS=linalg.a
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#

