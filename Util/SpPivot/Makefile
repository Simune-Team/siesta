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

SIESTA_OBJS  = precision.o parallel.o alloc.o m_io.o
SIESTA_OBJS += pxf.o reclat.o io.o memory.o memory_log.o

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

LOCAL_OBJS= pvtsp.o handlers.o local_sys.o

OBJS = $(SIESTA_OBJS) $(LOCAL_OBJS)

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
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.f) $(SIESTA_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.F) $(SIESTA_OBJS:.o=.F90)) \
		$(VPATH)/class*.T90 \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true

	@sed -i -e "/^m_ts_io.o:/s/m_ncdf_io.o//" Makefile
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
create_Sparsity_SC.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
create_Sparsity_SC.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
geom_helper.o: intrinsic_missing.o
geom_helper.o: intrinsic_missing.o
io.o: local_sys.o
io.o: local_sys.o
m_io_s.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_io_s.o: class_Sparsity.o parallel.o precision.o
m_io_s.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_io_s.o: class_Sparsity.o parallel.o precision.o
m_pivot.o: class_Sparsity.o m_pivot_methods.o m_region.o
m_pivot.o: class_Sparsity.o m_pivot_methods.o m_region.o
m_pivot_methods.o: class_Sparsity.o m_region.o precision.o
m_pivot_methods.o: class_Sparsity.o m_region.o precision.o
m_region.o: class_OrbitalDistribution.o class_Sparsity.o geom_helper.o
m_region.o: intrinsic_missing.o
m_region.o: class_OrbitalDistribution.o class_Sparsity.o geom_helper.o
m_region.o: intrinsic_missing.o
m_sparse.o: alloc.o class_OrbitalDistribution.o class_SpData2D.o
m_sparse.o: class_Sparsity.o geom_helper.o intrinsic_missing.o parallel.o
m_sparse.o: precision.o
m_sparse.o: alloc.o class_OrbitalDistribution.o class_SpData2D.o
m_sparse.o: class_Sparsity.o geom_helper.o intrinsic_missing.o parallel.o
m_sparse.o: precision.o
m_sparsity_handling.o: class_Data1D.o class_Data2D.o
m_sparsity_handling.o: class_OrbitalDistribution.o class_SpData1D.o
m_sparsity_handling.o: class_SpData2D.o class_Sparsity.o geom_helper.o
m_sparsity_handling.o: intrinsic_missing.o m_interpolate.o m_region.o
m_sparsity_handling.o: precision.o
m_sparsity_handling.o: class_Data1D.o class_Data2D.o
m_sparsity_handling.o: class_OrbitalDistribution.o class_SpData1D.o
m_sparsity_handling.o: class_SpData2D.o class_Sparsity.o geom_helper.o
m_sparsity_handling.o: intrinsic_missing.o m_interpolate.o m_region.o
m_sparsity_handling.o: precision.o
m_ts_io.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_io.o: class_SpData2D.o class_Sparsity.o geom_helper.o local_sys.o m_io_s.o
m_ts_io.o: m_os.o m_sparse.o memory_log.o parallel.o precision.o
m_ts_io.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_io.o: class_SpData2D.o class_Sparsity.o geom_helper.o local_sys.o m_io_s.o
m_ts_io.o: m_os.o m_sparse.o memory_log.o parallel.o precision.o
memory.o: memory_log.o parallel.o
memory.o: memory_log.o parallel.o
memory_log.o: m_io.o parallel.o precision.o
memory_log.o: m_io.o parallel.o precision.o
pvtsp.o: class_OrbitalDistribution.o class_Sparsity.o create_Sparsity_SC.o
pvtsp.o: geom_helper.o m_io_s.o m_os.o m_pivot.o m_pivot_methods.o m_region.o
pvtsp.o: m_sparsity_handling.o m_ts_io.o precision.o
pvtsp.o: class_OrbitalDistribution.o class_Sparsity.o create_Sparsity_SC.o
pvtsp.o: geom_helper.o m_io_s.o m_os.o m_pivot.o m_pivot_methods.o m_region.o
pvtsp.o: m_sparsity_handling.o m_ts_io.o precision.o
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
