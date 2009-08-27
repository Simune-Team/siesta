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

.SUFFIXES: .o .F90 .f90 .F .f .a 
#
default: module_built 
tests:   vpath_warning  pi3 blacs_prb pblas_prb
#
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
#
TEMPLATES= mpi__type_s.f90 mpi__type_sv.f90 mpi__type_v.f90 mpi__type_vs.f90
#
INCFLAGS=-I$(MPI_INCLUDE)
#
Interfaces.f90 V_S.uses VS.uses: $(TEMPLATES)
	if [ -z "$(KINDS)" ] ; then  $(MAKE) kind_explorer ; fi
	@echo "The kind numbers for single and double precision reals follow"
	sh ${<D}/generate.sh "$(KINDS)"
#
kind_explorer: kind_explorer.o
	$(FC) -o kind_explorer $(LDFLAGS) kind_explorer.o
mpi.o: mpi__include.o Interfaces.o V_S.uses VS.uses
#
module_built: libmpi_f90.a
	@mv libmpi_f90.a ..
	@cp *.mod ..
	@touch module_built
libmpi_f90.a: mpi.o mpi__include.o Interfaces.o
	@ar $(ARFLAGS_EXTRA) cru libmpi_f90.a mpi.o mpi__include.o Interfaces.o
	-$(RANLIB) libmpi_f90.a
#
#-----------------------------------------------------------------------
# Tests of the MPI installation
#
PARENT_VPATH=`grep '^VPATH' ../Makefile`
vpath_warning:
	@echo "** If you are using VPATH, compilations launched from this directory"
	@echo "** need to provide an explicit VPATH to the 'make' command."
	@echo "** Use the VPATH in Makefile in the top directory, followed by /MPI"
	@echo "** For example, to compile the tests:"
	@echo
	@echo   make $(PARENT_VPATH)/MPI tests
	@echo
#
# Straight MPI test
#
pi3:    vpath_warning pi3.o
	$(FC) -o pi3 $(LDFLAGS) pi3.o ../$(MPI_INTERFACE) $(MPI_LIBS) 
#
# Blacs test
#
blacs_prb:  vpath_warning  blacs_prb.o
	$(FC) -o blacs_prb $(LDFLAGS) blacs_prb.o ../$(MPI_INTERFACE) $(LIBS) 
#
# Pblas test. You might need to uncomment the following line
# in those systems where pdlaprnt is not in the standard library.
#
#EXTRA_ROUTINE=pdlaprnt.o
#
pblas_prb:   vpath_warning pblas_prb.o $(EXTRA_ROUTINE)
	@echo
	@echo --Attempting to build Pblas test. 
	@echo --You might need to enable compilation of pdlaprnt.f in the Makefile
	@echo --if that routine is not in the standard library on your system...
	@echo
	$(FC) -o pblas_prb $(LDFLAGS) pblas_prb.o $(EXTRA_ROUTINE) \
              ../$(MPI_INTERFACE) $(LIBS) 
#
clean:
	@rm -f Interfaces.f90 *.o *.mod
	@rm -f module_built *.uses kind_explorer
	@rm -f pi3 blacs_prb pblas_prb




