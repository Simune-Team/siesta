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
VPATH=../../Grid

INCFLAGS += $(NETCDF_INCFLAGS)
LIBS = $(NETCDF_LIBS)
#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

plstm: plstm.o m_gridfunc.o
	$(FC) $(LDFLAGS) -o $@  plstm.o m_gridfunc.o $(LIBS)

plstm.o: m_gridfunc.o

plsts: plsts.o m_gridfunc.o
	$(FC) $(LDFLAGS) -o $@  plsts.o m_gridfunc.o $(LIBS)

plsts.o: m_gridfunc.o

clean: 
	rm -f *.o plstm plsts

