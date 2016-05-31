# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=itanium2
#
# Serial compilation
#
FC=efc -Vaxlib
#
FC_ASIS=$(FC)
#
FFLAGS= -O2 -tpp2 -W0
FFLAGS_DEBUG= -g -O0
LDFLAG=#-Vaxlib 
RANLIB=echo
LAPACK=-L/opt/intel/mkl72/lib/64 -lmkl_lapack64 -lmkl -lguide -lpthread
LIBS=$(LAPACK)
SYS=bsd
FPPFLAGS= -DWXML_INIT_FIX -DALLOC_SAVE_BUG
COMP_LIBS=
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



