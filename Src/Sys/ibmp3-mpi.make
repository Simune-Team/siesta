# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile include file for IBM SP2 with Power3 procs. MPI.
# No use of ESSL (it lacks zhegv for complex generalized diagonalization).
#
# Emilio Artacho, October 2000

SIESTA_ARCH=ibmp3-mpi

FC=mpxlf
FC_ASIS=$(FC)

FFLAGS= -O3 -qarch=auto -qtune=auto -qcache=auto -qmaxmem=16000 -qnolm
FFLAGS_DEBUG= -g -C -qflag=W:W -qinitauto -qsave -qmaxmem=16000 -qnolm
# -qipa gives speed optimization of a few percent, but takes long to link

LIBS= -Wl,-bD:2000000000 -qsave \
      /home/root/SCA/SCALAPACK/scalapack_RS6K.a \
      /home/root/SCA/SCALAPACK/pblas_RS6K.a \
      /home/root/SCA/SCALAPACK/tools_RS6K.a \
      -lblacs \
      /usr/lpp/ppe.poe/lib/libmpi.a \
      -qnolm
# essl lacks zhegv of lapack used for scalar diagonalization.
# It is not called in mpi runs, but it is within a not preprocessed
# if, and requires to be resolved. The following line is thus added
# to link the scalar lapack and blas routines provided within the siesta
# package, even though it has no influence on the run
COMP_LIBS=linalg.a
#
SYS=ibm
RANLIB=ranlib
MPI_INCLUDE=/usr/lpp/ppe.poe/include
MPI_INTERFACE=libmpi_f90.a
#
FPPFLAGS=-WF,-DMPI
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
