# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
.SUFFIXES: .f .F .o .a .f90 .F90

VPATH:=$(shell pwd)/../../../Src

OBJDIR=Obj

ARCH_MAKE=../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

# Set the default non-mpi routines
FC_DEFAULT  := $(FC)
FC_SERIAL   ?= $(FC_DEFAULT)

EXE = tshs2tshs
default: $(EXE)

dep:
	-sfmakedepend --depend=obj  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 $(VPATH)/*.T90 *.f90
	@sed -i -e "/^siesta_cmlsubs.o:/s/version.o//" Makefile
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# Only because the NetCDF libraries could be parallel, in
# which case we need the mpi-compiler :(
# However, we still compile everything else as though it did not exist
#FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:=$(DEFS_PREFIX) $(DEFS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UNCDF_PARALLEL
FPPFLAGS:=$(DEFS_PREFIX) $(FPPFLAGS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UNCDF_PARALLEL
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)

# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)

SYSOBJ=$(SYS).o

# systematic siesta files
OBJS = precision.o parallel.o \
        alloc.o pxf.o units.o local_sys.o \
	memory.o memoryinfo.o m_io.o io.o \
	m_sparse.o intrinsic_missing.o geom_helper.o reclat.o \
	m_io_s.o m_ncdf_io.o m_ts_io.o 

# Add class objects
OBJS += m_uuid.o object_debug.o \
	class_Sparsity.o class_OrbitalDistribution.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o

# Add sys object
OBJS += $(SYSOBJ)

OBJS += m_os.o
OBJS += m_ts_io_version.o
OBJS += tshs2tshs.o

#######################
#  TS-TBT libs:
#
# Libraries used by transiesta and tbtrans
# They are currently supplied by this distribution
# All rights reserved for the author of these libraries.
# Author: Nick Papior Andersen

# This library can be found at:
#  https://github.com/zerothi/fdict
FDICT=libfdict.a
FDICT_MAKEFILE=$(VPATH)/fdict/Makefile
$(FDICT): 
	-mkdir -p fdict
	(cd fdict ; \
	echo "TOP_DIR=$(VPATH)/fdict" > Makefile ; echo "include $(VPATH)/fdict/Makefile" >> Makefile ; \
	echo "FC=$(FC_SERIAL)" > setup.make ; \
	echo "AR=$(AR)" >> setup.make ; \
	echo "RANLIB=$(RANLIB)" >> setup.make ; \
	echo "FFLAGS=$(FFLAGS)" >> setup.make ; \
	$(MAKE) -j 1 copy ; $(MAKE) -j 1 lib SHARED=0)
	@cp -p fdict/*.mod ./
	cp fdict/$(FDICT) ./

# This library can be found at:
#  https://github.com/zerothi/ncdf
NCDF=libncdf.a
NCDF_MAKEFILE=$(VPATH)/NCDF/src/Makefile
$(NCDF): COMP_LIBS += $(FDICT)
$(NCDF): $(FDICT)
	-mkdir -p ncdf
	(cd ncdf ; \
	echo "TOP_DIR=$(VPATH)/ncdf" > Makefile ; echo "include $(VPATH)/ncdf/Makefile" >> Makefile ; \
	echo "settings.bash: FORCE" >> Makefile ; \
	printf '%s\n' "	-cp ../fdict/settings.bash ." >> Makefile ; \
	echo "FC=$(FC)" > setup.make ; \
	echo "AR=$(AR)" >> setup.make ; \
	echo "RANLIB=$(RANLIB)" >> setup.make ; \
	echo "FFLAGS=$(FFLAGS)" >> setup.make ; \
	echo "FPPFLAGS = $(FPPFLAGS)" >> setup.make ; \
	echo "INCLUDES = $(INCFLAGS) -I../" >> setup.make ; \
	echo "LIBS = $(LIBS) -lfdict" >> setup.make ; \
	$(MAKE) -j 1 copy ; $(MAKE) -j 1 lib SHARED=0)
	@cp -p ncdf/*.mod ./
	cp ncdf/$(NCDF) ./

#######################

# Use the makefile in Src/fdf and all the sources there.
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)


$(OBJS): $(FDF) $(COMP_LIBS) $(TSHS_COMP_LIBS)
$(EXE): $(FDF) $(COMP_LIBS) $(TSHS_COMP_LIBS) $(OBJS) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXE) \
		$(OBJS) $(COMP_LIBS) $(TSHS_COMP_LIBS) $(FDF) $(LIBS) 

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f $(EXE) *.o *.a *.mod
	-rm -rf ./fdf
	-rm -rf ./fdict
	-rm -rf ./ncdf


# DO NOT DELETE THIS LINE - used by make depend
