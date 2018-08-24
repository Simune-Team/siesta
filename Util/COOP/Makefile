# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
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
.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o

VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj

default: mprop fat

all: mprop dm_creator fat

include ../../$(OBJDIR)/arch.make

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI 
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI 
INCFLAGS:=$(INCFLAGS) $(NETCDF_INCFLAGS)

SIESTA_OBJS = precision.o alloc.o parallel.o m_io.o pxf.o moreParallelSubs.o

DEF_OBJS = m_getopts.o io.o subs.o local_sys.o units.o f2kcli.o
DEF_OBJS += io_hs.o
DEF_OBJS += main_vars.o

# Create objects for each routine
MPROP_OBJS = $(SIESTA_OBJS) $(DEF_OBJS) read_curves.o orbital_set.o
FAT_OBJS   = $(SIESTA_OBJS) $(DEF_OBJS) read_curves.o orbital_set.o
DM_OBJS    = $(SIESTA_OBJS) $(DEF_OBJS) iodm_netcdf.o write_dm.o

dep: 
	sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
		*.f *.f90 *.F *.F90

mprop: $(MPROP_OBJS) mprop.o
	$(FC) -o $@ $(LDFLAGS) $(MPROP_OBJS) mprop.o
fat: $(FAT_OBJS) fat.o
	$(FC) -o $@ $(LDFLAGS) $(FAT_OBJS) fat.o
dm_creator: $(DM_OBJS) dm_creator.o
	$(FC) -o $@ $(LDFLAGS) $(DM_OBJS) dm_creator.o $(NETCDF_LIBS)

clean:
	rm -f *.o *.mod mprop dm_creator fat

# DO NOT DELETE THIS LINE - used by make depend
