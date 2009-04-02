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
# Makefile for PDOS XML processing
#
# The VPATH directive below helps to re-use f2kcli.F90 from the top Src.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Src, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#  
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: pdosxml
OBJDIR=Src
VPATH=../../Src
#
include ../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
# This should be INCFLAGS...
#
FFLAGS=$(FFLAGS_DEBUG) -I../../$(OBJDIR)
#
LIBS= ../../$(OBJDIR)/libxmlparser.a
#
OBJS=  m_orbital_chooser.o m_pdos.o f2kcli.o pdosxml.o
     
pdosxml: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)
#
clean: 
	rm -f *.o pdosxml *.o *.*d
#










