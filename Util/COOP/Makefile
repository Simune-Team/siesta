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

MAIN_OBJDIR:=$(shell cd -P -- ../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI 
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI 
INCFLAGS:=$(INCFLAGS) $(NETCDF_INCFLAGS)

# Two ways of splitting the list of object files:
# 1. By the needs of the three different programs
SIESTA_OBJS = precision.o alloc.o parallel.o m_io.o pxf.o moreParallelSubs.o
DEF_OBJS = m_getopts.o io.o subs.o local_sys.o units.o f2kcli.o
DEF_OBJS += io_hs.o
DEF_OBJS += main_vars.o

# Create objects for each routine
MPROP_OBJS = $(SIESTA_OBJS) $(DEF_OBJS) read_curves.o orbital_set.o handlers.o
FAT_OBJS   = $(SIESTA_OBJS) $(DEF_OBJS) read_curves.o orbital_set.o handlers.o
DM_OBJS    = $(SIESTA_OBJS) $(DEF_OBJS) iodm_netcdf.o write_dm.o handlers.o

# 2. By whether they are in VPATH or in the local directory (for deps processing)
DEP_OBJS = $(SIESTA_OBJS)
DEP_OBJS += m_getopts.o io.o  units.o f2kcli.o iodm_netcdf.o write_dm.o
LOCAL_OBJS= subs.o io_hs.o local_sys.o main_vars.o read_curves.o orbital_set.o handlers.o

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true


mprop: $(MPROP_OBJS) mprop.o
	$(FC) -o $@ $(LDFLAGS) $(MPROP_OBJS) mprop.o
fat: $(FAT_OBJS) fat.o
	$(FC) -o $@ $(LDFLAGS) $(FAT_OBJS) fat.o
dm_creator: $(DM_OBJS) dm_creator.o
	$(FC) -o $@ $(LDFLAGS) $(DM_OBJS) dm_creator.o $(NETCDF_LIBS)

clean:
	rm -f *.o *.mod mprop dm_creator fat

# DO NOT DELETE THIS LINE - used by make depend
io.o: local_sys.o
io.o: local_sys.o
iodm_netcdf.o: alloc.o local_sys.o parallel.o precision.o
iodm_netcdf.o: alloc.o local_sys.o parallel.o precision.o
m_getopts.o: f2kcli.o
m_getopts.o: f2kcli.o
moreParallelSubs.o: alloc.o local_sys.o m_io.o parallel.o precision.o
moreParallelSubs.o: alloc.o local_sys.o m_io.o parallel.o precision.o
units.o: precision.o
units.o: precision.o
io_hs.o: main_vars.o precision.o
io_hs.o: main_vars.o precision.o
main_vars.o: m_getopts.o precision.o subs.o units.o
main_vars.o: m_getopts.o precision.o subs.o units.o
orbital_set.o: main_vars.o subs.o
orbital_set.o: main_vars.o subs.o
read_curves.o: orbital_set.o precision.o
read_curves.o: orbital_set.o precision.o
subs.o: precision.o
subs.o: precision.o
sys.o: local_sys.o
