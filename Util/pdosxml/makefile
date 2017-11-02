# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for PDOS XML processing
#
# The VPATH directive below helps to re-use f2kcli.F90 from the top Src.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Obj, but if you are using architecture-dependent build directories
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
.PHONY: xmlparser clean
#
default: pdosxml
OBJDIR=Obj
VPATH:=$(shell pwd)/../../Src
#
MAIN_OBJDIR:=$(shell cd -P -- ../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
XMLF90_INCFLAGS= -I $(XMLF90_ROOT)/include 
XMLF90_LIBS=$(XMLF90_ROOT)/lib/libxmlf90.a
INCFLAGS:=$(INCFLAGS) $(XMLF90_INCFLAGS)
#
LIBS= $(XMLF90_LIBS)
#
OBJS=  m_orbital_chooser.o m_pdos.o f2kcli.o pdosxml.o
#     
pdosxml:$(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)
#
clean: 
	rm -f *.o pdosxml *.o *.*d
#










