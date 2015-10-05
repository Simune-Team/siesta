# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#
# Makefile for PLSTM 
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90

OBJDIR=Obj

all: plstm

include ../../../$(OBJDIR)/arch.make

#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

plstm: plstm.o
	$(FC) $(LDFLAGS) -o $@  plstm.o 

clean: 
	rm -f *.o plstm *.*d


dep:
	-sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
		*.f *.f90 *.F *.F90
	@sed '/^siesta_cmlsubs.o:/s/version.o//'  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile

# DO NOT DELETE THIS LINE - used by make depend
