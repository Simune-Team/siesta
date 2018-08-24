# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---

# Makefile for pvtsp

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

default: what pvtsp

VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

what:
	@echo
	@echo "Compilation architecture to be used: $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo

# Set the default non-mpi routines
FC_DEFAULT := $(FC)
FC_SERIAL  ?= $(FC_DEFAULT)
FC := $(FC_SERIAL)

# We force the compilation of the graphviz library
FPPFLAGS += -DGRAPHVIZ

SIESTA_OBJS  = precision.o parallel.o alloc.o local_sys.o m_io.o
SIESTA_OBJS += pxf.o reclat.o io.o memory.o memoryinfo.o

# Add specific class objects
SIESTA_OBJS += m_uuid.o object_debug.o class_Sparsity.o \
	class_OrbitalDistribution.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o

# Some of NPA modules
SIESTA_OBJS += intrinsic_missing.o geom_helper.o m_io_s.o m_os.o
SIESTA_OBJS += create_Sparsity_SC.o m_ts_io.o
SIESTA_OBJS += m_pivot.o m_pivot_methods.o
SIESTA_OBJS += m_sparse.o m_sparsity_handling.o
SIESTA_OBJS += m_region.o m_interpolate.o

OBJS = $(SIESTA_OBJS) pvtsp.o

pvtsp: FPPFLAGS+=$(DEFS_PREFIX) -UMPI -UNCDF -UNCDF_4 -UNCDF_PARALLEL
pvtsp: COMP_LIBS := $(filter-out libncdf.a,$(COMP_LIBS))
pvtsp: COMP_LIBS := $(filter-out libfdict.a,$(COMP_LIBS))
pvtsp: LIBS := $(filter-out libncdf.a,$(LIBS))
pvtsp: LIBS := $(filter-out libfdict.a,$(LIBS))
pvtsp: $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o pvtsp \
              $(OBJS) $(COMP_LIBS) $(LIBS) || echo "***Please see top of Makefile"

clean: 
	@echo "==> Cleaning object, library, and executable files"
	-rm -f pvtsp *.o *.a *.mod

# Dependencies
dep:
	-sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
                $(VPATH)/*.T90 \
		*.f90 *.F90 *.f *.F
	@sed -i -e "/^m_ts_io.o:/s/m_ncdf_io.o//" Makefile
	@sed -i -e "/^siesta_cmlsubs.o:/s/version.o//" Makefile
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# DO NOT DELETE THIS LINE - used by make depend
