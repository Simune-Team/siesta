# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=itanium2-mpi
#
#  itanium2 machine at Leioa with myrinet
#
FC=efc -Vaxlib
#
FC_ASIS=$(FC)
#
FFLAGS= -O2 -tpp2 -W0 
FFLAGS_DEBUG= -g -O0
LDFLAG=#-Vaxlib 
RANLIB=echo
LIBS=  
SYS=scalapack_extra
FPPFLAGS=-DMPI -DNODAT -DWXML_INIT_FIX -DHAS_DLAMC3 -DALLOC_SAVE_BUG
MPI_INTERFACE=libmpi_f90.a
#
MPIROOT=/usr/local/mpich
MPI_INCLUDE=$(MPIROOT)/include
#
#MPI_LIBS= $(MPIROOT)/lib/libmpich.a  $(MPIROOT)/lib/libpmpich.a  \
#         /usr/local/gm-2.0.6_Linux/binary/lib/libgm.a \
#          -L/opt/intel/mkl61/lib/64 -lpthread
#
MPI_LIBS=-L$(MPIROOT)/lib \
         -lmpichf90 -lpmpich -lmpich -lpmpich -lmpich \
         -L/usr/local/gm/binary/lib/ -L/usr/local/gm/lib/ \
         -lgm -lpthread -lPEPCF90
#
LAPACK=-L/opt/intel/mkl72/lib/64 -lmkl_lapack64 -lmkl -lguide -lpthread
#
BLACS_LIBS=$(HOME)/lib/libblacs.a
SCALAPACK_LIBS=-L$(HOME)/lib -lscalapack-arina.GM  \
                                    -lpblas-arina.GM \
                                    -ltools-arina.GM \
			            -lredist-arina.GM
#
LIBS=  $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK) \
       $(MPI_LIBS)  $(NETCDF_LIBS)
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



