# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for PLSTM 
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90

OBJDIR=Obj

all: plstm plsts

include ../../../$(OBJDIR)/arch.make
vpath %.F90 ../../Grid

#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
# However, it breaks compatibility since we can't be sure
# that users define NETCDF_LIBS.
# For now we will remove this but let users know that
# they can revert changes by uncommenting out the following
# two lines:
#INCFLAGS += $(NETCDF_INCFLAGS)
#LIBS = $(NETCDF_LIBS)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

plstm.o: m_gridfunc.o options.o
plsts.o: m_gridfunc.o
plstm: options.o plstm.o m_gridfunc.o
	$(FC) $(LDFLAGS) -o $@  plstm.o m_gridfunc.o options.o $(LIBS)


plsts: plsts.o m_gridfunc.o 
	$(FC) $(LDFLAGS) -o $@  plsts.o m_gridfunc.o $(LIBS)


clean: 
	rm -f *.o plstm plsts

