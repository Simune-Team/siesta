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
#
# Makefile for PDOS XML processing
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: pdos
include ../../Src/arch.make
#
FFLAGS=$(FFLAGS_DEBUG) -I../../Src
#
LIBS= ../../Src/libxmlparser.a
#
OBJS=  m_orbital_chooser.o m_pdos.o f2kcli.o pdos.o
     
pdos: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)
#
clean: 
	rm -f *.o pdos *.o *.*d
#










