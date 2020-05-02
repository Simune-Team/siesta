# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone get_chem_labels

OBJDIR=Obj

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

VPATH:=$(shell pwd)/../../../Src


EXE = ts2ts
default: $(EXE)

dep:
	sfmakedepend --depend=obj  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 $(VPATH)/*.T90
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

override WITH_MPI=

ARCH_MAKE=../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)

# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)

SYSOBJ=$(SYS).o

# Note that machine-specific files are now in top Src directory.
DEP_OBJS = precision.o f2kcli.o parallel.o \
          alloc.o m_io.o pxf.o units.o $(SYSOBJ)

LOCAL_OBJS += ts2ts.o handlers.o

OBJS = $(DEP_OBJS) $(LOCAL_OBJS)

# Use the makefile in Src/fdf and all the sources there.
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

$(EXE): $(FDF) $(OBJS)
	$(FC) -o $(EXE) \
	       $(LDFLAGS) $(OBJS) $(FDF)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f $(EXE) *.o  *.a *.mod
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true

# DO NOT DELETE THIS LINE - used by make depend
units.o: precision.o
ts2ts.o: f2kcli.o precision.o units.o
