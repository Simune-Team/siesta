# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile include file for IBM SP2 with Power3 processors
# No use of ESSL (it lacks zhegv for complex generalized diagonalization).
#
# Emilio Artacho, October 2000

SIESTA_ARCH=ibmp3

FC=xlf 
#FC=xlf -bloadmap:MAP
FC_ASIS=$(FC)

FFLAGS= -O3 -qarch=auto -qtune=auto -qcache=auto -qnolm
FFLAGS_DEBUG= -g -C -qinitauto -qsave -qmaxmem=16000 -qnolm
# -qipa gives speed optimization of a few percent, but takes long to link

LIBS= -Wl,-bD:2000000000 -qnolm 
COMP_LIBS=linalg.a

MPILIB=
SYS=ibm
RANLIB=ranlib
MPI_INCLUDE=/usr/local/include
FPPFLAGS=
FREE_F90=-qsuffix=f=f90 -qfree=f90
FREE_F90_CPP=-qsuffix=cpp=F90 -qfree=f90
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FREE_F90) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FREE_F90_CPP) $<
