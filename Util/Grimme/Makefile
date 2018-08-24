# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---

OBJDIR=Obj

.SUFFIXES: .f .F .o .a  .f90 .F90

VPATH:=$(shell pwd)/../../Src


EXE = fdf2grimme
default: $(EXE)

dep:
	sfmakedepend --depend=obj  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 *.f90

ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:=$(DEFS_PREFIX) -UMPI $(DEFS) $(DEFS_PREFIX) -UMPI
FPPFLAGS:=$(DEFS_PREFIX) -UMPI $(FPPFLAGS) $(DEFS_PREFIX) -UMPI
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)

# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)

SYSOBJ=$(SYS).o

# Note that machine-specific files are now in top Src directory.
OBJS = precision.o f2kcli.o parallel.o \
          chemical.o io.o m_io.o alloc.o \
          memoryinfo.o memory.o periodic_table.o \
          pxf.o atom_options.o local_sys.o $(SYSOBJ)

OBJS += grimme.o


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
	-rm -rf ./fdf

# DO NOT DELETE THIS LINE - used by make depend
