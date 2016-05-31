# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=nolibs
#
# Serial compilation without the need of any installed libraries.
# You still need to change the name of the compiler, etc.
#
FC=gfortran
FC_ASIS=$(FC)
#
FFLAGS= -O
FFLAGS_DEBUG= -g -O0
RANLIB=echo "do we need to ranlib this? : "
LIBS=  
SYS=nag
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
