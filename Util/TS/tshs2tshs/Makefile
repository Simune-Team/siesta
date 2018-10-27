# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
.SUFFIXES: .f .F .o .a .f90 .F90

EXE = tshs2tshs
default: $(EXE)

VPATH:=$(shell pwd)/../../../Src

OBJDIR=Obj

MAIN_OBJDIR:=$(shell cd -P -- ../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)

# Set the default non-mpi routines
FC_DEFAULT  := $(FC)
FC_SERIAL   ?= $(FC_DEFAULT)


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
DEP_OBJS = precision.o parallel.o \
        alloc.o pxf.o units.o  \
	memory_log.o m_io.o io.o \
	m_sparse.o intrinsic_missing.o geom_helper.o reclat.o \
	m_io_s.o m_ncdf_io.o m_ts_io.o 

# Add class objects
DEP_OBJS += m_uuid.o object_debug.o \
	class_Sparsity.o class_OrbitalDistribution.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o

# Add sys object
DEP_OBJS += $(SYSOBJ)

DEP_OBJS += m_os.o
DEP_OBJS += m_ts_io_version.o

LOCAL_OBJS = tshs2tshs.o handlers.o local_sys.o

OBJS=$(DEP_OBJS) $(LOCAL_OBJS)

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
	(cd fdict ; $(MAKE) -f $(FDICT_MAKEFILE) \
		"VPATH=$(VPATH)/fdict" \
		"FC=$(FC_SERIAL)" \
		"ARCH_MAKE=$(ARCH_MAKE)" \
		"FPPFLAGS=$(FPPFLAGS)" \
		"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" lib)
	@cp -p fdict/*.*d ./
	cp fdict/$(FDICT) ./

# This library can be found at:
#  https://github.com/zerothi/ncdf
NCDF=libncdf.a
NCDF_MAKEFILE=$(VPATH)/NCDF/src/Makefile
$(NCDF): COMP_LIBS += $(FDICT)
$(NCDF): $(FDICT)
	-mkdir -p NCDF/src
	(cd NCDF/src ; $(MAKE) -f $(NCDF_MAKEFILE) \
		"VPATH=$(VPATH)/NCDF/src" "LIBFDICT=$(FDICT)" \
		"DIR_FDICT=$(VPATH)/fdict" \
		"ARCH_MAKE=$(ARCH_MAKE)" \
		"FPPFLAGS=$(FPPFLAGS)" \
		"FFLAGS=$(FFLAGS:$(IPO_FLAG)=) -I../../" lib)
	@cp -p NCDF/src/*.*d ./
	cp NCDF/src/$(NCDF) ./

#######################

# Use the makefile in Src/fdf and all the sources there.
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a

$(OBJS): $(FDF) $(COMP_LIBS) $(TSHS_COMP_LIBS)
$(EXE): $(FDF) $(NCPS) $(COMP_LIBS) $(TSHS_COMP_LIBS) $(OBJS) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXE) \
		$(OBJS) $(COMP_LIBS) $(TSHS_COMP_LIBS) $(FDF) $(LIBS) 

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f $(EXE) *.o *.a *.mod
	-rm -rf ./fdict
	-rm -rf ./NCDF

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(VPATH)/class*.T90 \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# DO NOT DELETE THIS LINE - used by make depend
class_Data1D.o: alloc.o
class_Data2D.o: alloc.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o
class_SpData1D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o
class_SpData2D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o
class_Sparsity.o: alloc.o
class_Sparsity.o: alloc.o
class_TriMat.o: alloc.o intrinsic_missing.o
geom_helper.o: intrinsic_missing.o
geom_helper.o: intrinsic_missing.o
io.o: m_io.o
io.o: m_io.o
m_io_s.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_io_s.o: class_Sparsity.o parallel.o precision.o
m_io_s.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_io_s.o: class_Sparsity.o parallel.o precision.o
m_ncdf_io.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_ncdf_io.o: class_Sparsity.o m_io_s.o parallel.o precision.o
m_ncdf_io.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_ncdf_io.o: class_Sparsity.o m_io_s.o parallel.o precision.o
m_sparse.o: alloc.o class_OrbitalDistribution.o class_SpData2D.o
m_sparse.o: class_Sparsity.o geom_helper.o intrinsic_missing.o parallel.o
m_sparse.o: precision.o
m_sparse.o: alloc.o class_OrbitalDistribution.o class_SpData2D.o
m_sparse.o: class_Sparsity.o geom_helper.o intrinsic_missing.o parallel.o
m_sparse.o: precision.o
m_ts_io.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_io.o: class_SpData2D.o class_Sparsity.o geom_helper.o local_sys.o m_io_s.o
m_ts_io.o: m_ncdf_io.o m_os.o m_sparse.o memory_log.o parallel.o precision.o
m_ts_io.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_io.o: class_SpData2D.o class_Sparsity.o geom_helper.o local_sys.o m_io_s.o
m_ts_io.o: m_ncdf_io.o m_os.o m_sparse.o memory_log.o parallel.o precision.o
memory_log.o: m_io.o parallel.o precision.o
memory_log.o: m_io.o parallel.o precision.o
units.o: precision.o
units.o: precision.o
tshs2tshs.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
tshs2tshs.o: class_Sparsity.o geom_helper.o m_sparse.o m_ts_io.o parallel.o
tshs2tshs.o: precision.o units.o
tshs2tshs.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
tshs2tshs.o: class_Sparsity.o geom_helper.o m_sparse.o m_ts_io.o parallel.o
tshs2tshs.o: precision.o units.o
class_cdata1d.o: class_Data1D.o
class_cdata2d.o: class_Data2D.o
class_cspdata1d.o: class_SpData1D.o
class_cspdata2d.o: class_SpData2D.o
class_ddata1d.o: class_Data1D.o
class_ddata2d.o: class_Data2D.o
class_dspdata1d.o: class_SpData1D.o
class_dspdata2d.o: class_SpData2D.o
class_gdata1d.o: class_Data1D.o
class_gspdata1d.o: class_SpData1D.o
class_idata1d.o: class_Data1D.o
class_idata2d.o: class_Data2D.o
class_ispdata1d.o: class_SpData1D.o
class_ispdata2d.o: class_SpData2D.o
class_ldata1d.o: class_Data1D.o
class_ldata2d.o: class_Data2D.o
class_lspdata1d.o: class_SpData1D.o
class_lspdata2d.o: class_SpData2D.o
class_sdata1d.o: class_Data1D.o
class_sdata2d.o: class_Data2D.o
class_sspdata1d.o: class_SpData1D.o
class_sspdata2d.o: class_SpData2D.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
m_object_debug.o: object_debug.o
mtprng.o: m_uuid.o
sys.o: local_sys.o
