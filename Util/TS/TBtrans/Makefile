# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---

# Makefile for tbtrans

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

.PHONY: default what dep version clean-tbt clean

default: what tbtrans 

VPATH:=$(shell pwd)/../../../Src
OBJDIR=Obj
ARCH_MAKE=../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

# Set the default non-mpi routines
FC_DEFAULT := $(FC)
FC_SERIAL  ?= $(FC_DEFAULT)

# For removing any -ipo compilation in sub directories...
# Later compilers are also having this enabled, and it lets you
# specify what to not be used in libraries.
# This is so because inter-procedural compilation should NOT be performed in
# libraries (it removes names when a routine can be fully moved into others)
IPO_FLAG?= -ipo

SIESTA_SYS=$(SIESTA_ARCH)

what:
	@echo
	@echo "Compilation architecture to be used: $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo

SYSOBJ=$(SYS).o

SIESTA_OBJS= bye.o precision.o m_io.o alloc.o memory_log.o memory.o pxf.o \
             moreParallelSubs.o m_mpi_utils.o \
             reclat.o idiag.o minvec.o \
             parallel.o io.o files.o cellsubs.o sorting.o \
	     kpoint_convert.o timestamp.o m_wallclock.o m_spin.o \
             diag.o diag_option.o \
	     debugmpi.o find_kgrid.o units.o

# Siesta timings
SIESTA_OBJS += m_walltime.o timer_tree.o m_timer.o timer.o

# Add specific class objects
# class Objects
SIESTA_OBJS += class_OrbitalDistribution.o \
	class_Sparsity.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
	class_TriMat.o \
        m_uuid.o object_debug.o

# As the classes are based on inputs we force the dependency to
# alloc
class_Sparsity.o class_Data1D.o class_Data2D.o class_TriMat.o: alloc.o

# Some of NPA modules
SIESTA_OBJS += m_char.o cross.o m_os.o
SIESTA_OBJS += intrinsic_missing.o geom_helper.o m_io_s.o m_sparse.o
SIESTA_OBJS += m_handle_sparse.o m_iodm.o
SIESTA_OBJS += m_integrate.o m_interpolate.o m_mat_invert.o
SIESTA_OBJS += m_iterator.o fdf_extra.o
SIESTA_OBJS += create_Sparsity_SC.o create_Sparsity_Union.o
SIESTA_OBJS += m_region.o m_sparsity_handling.o m_pivot_array.o
SIESTA_OBJS += m_pivot.o m_pivot_methods.o
SIESTA_OBJS += m_ncdf_io.o

# Geometry modules
SIESTA_OBJS += m_geom_aux.o m_geom_objects.o
SIESTA_OBJS += m_geom_box.o m_geom_coord.o m_geom_square.o m_geom_plane.o

# The needed transiesta objects
TS_OBJS  = m_ts_electype.o m_ts_method.o m_ts_chem_pot.o
TS_OBJS += m_ts_electrode.o m_ts_elec_se.o m_ts_gf.o
TS_OBJS += m_ts_cctype.o m_ts_sparse_helper.o
TS_OBJS += m_ts_io.o m_ts_iodm.o m_ts_io_ctype.o m_ts_debug.o
TS_OBJS += m_ts_io_contour.o m_gauss_quad.o

# These could be left out, however, there are dependencies
TS_OBJS += m_ts_contour_eq.o m_ts_contour_neq.o
TS_OBJS += m_gauss_fermi_inf.o m_gauss_fermi_30.o m_gauss_fermi_28.o
TS_OBJS += m_gauss_fermi_26.o m_gauss_fermi_24.o m_gauss_fermi_22.o
TS_OBJS += m_gauss_fermi_20.o m_gauss_fermi_19.o m_gauss_fermi_18.o
TS_OBJS += m_gauss_fermi_17.o m_ts_aux.o

# The block-tri-diagonal matrix routines
TS_OBJS += m_ts_tri_scat.o m_trimat_invert.o m_ts_trimat_invert.o
TS_OBJS += m_ts_tri_common.o m_ts_rgn2trimat.o

# For creating the global-unit-cell sparsity pattern
TS_OBJS += m_ts_sparse.o m_ts_tdir.o m_ts_pivot.o

DEP_OBJS = $(SIESTA_OBJS) $(TS_OBJS) $(SYSOBJ) $(EXTRA)

# To separate tbtrans objects from siesta/transiesta
TBTOBJS  = version.o m_verbosity.o

# Add generic tbtrans routines
TBTOBJS += tbt_init.o tbt_end.o tbt_reinit.o m_tbt_hs.o m_tbt_regions.o
TBTOBJS += m_tbt_options.o m_tbt_tri_init.o
TBTOBJS += m_tbt_contour.o m_tbt_kpoint.o m_tbt_gf.o

# Different solution methods
TBTOBJS += m_tbtrans.o m_tbt_trik.o m_tbt_tri_scat.o
TBTOBJS += m_tbt_proj.o m_tbt_diag.o m_tbt_dH.o

# The saving methods
TBTOBJS += m_tbt_sigma_save.o m_tbt_save.o

# The k-region files
TBTOBJS += m_tbt_kregions.o m_tbt_sparse_helper.o

# The main program
TBTOBJS += tbtrans.o 

# Handlers and resolvers of last resort
TBTOBJS += handlers.o local_sys.o

#
# The version program
.PHONY: version.o
version.o: version
version: version.F90
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;\
		s'TBTRANS_VERSION'$$(cat ${<D}/../../../version.info)'g;\
		s'LIBS'$(LIBS)'g;\
		s'FPPFLAGS'$(FPPFLAGS)'g" $<  > tmp.F90
	@awk '{if (length>80) { cur=78; \
		printf "%s&\n",substr($$0,0,78); \
		while(length-cur>78){ cur=cur+76 ; \
		printf "&%s&\n",substr($$0,cur-76,76) \
		} printf "&%s\n",substr($$0,cur)} else { print $$0 }}' tmp.F90 > compinfo.F90
	@($(MAKE) "FPPFLAGS=$(FPPFLAGS)" compinfo.o)
	@mv compinfo.o version.o
	@echo


# Interfaces to libraries
MPI_MAKEFILE=$(VPATH)/MPI/Makefile
libmpi_f90.a: 
	-mkdir -p MPI
	(cd MPI ; \
	$(MAKE) -f $(MPI_MAKEFILE) -j 1 "VPATH=$(VPATH)/MPI" \
	"MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=../$(ARCH_MAKE)" \
	"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module_built)

FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf -I../MPI $(INCFLAGS)
$(FDF): $(MPI_INTERFACE)
	-mkdir -p fdf
	(cd fdf ; \
	$(MAKE) -f $(FDF_MAKEFILE) -j 1 "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
	"ARCH_MAKE=../$(ARCH_MAKE)" \
	"INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module)

#######################
#  TS-TBT libs:
#
# Libraries used by transiesta and tbtrans
# They are currently supplied by this distribution
# All rights reserved for the author of these libraries.
# Author: Nick R. Papior

# This library can be found at:
#  https://github.com/zerothi/fdict
FDICT=libfdict.a
$(FDICT): 
	-mkdir -p fdict
	(cd fdict ; \
	echo "TOP_DIR=$(VPATH)/fdict" > Makefile ; echo "include $(VPATH)/fdict/Makefile" >> Makefile ; \
	echo "FC=$(FC_SERIAL)" > setup.make ; \
	echo "FFLAGS=$(FFLAGS)" >> setup.make ; \
	$(MAKE) -j 1 copy ; $(MAKE) -j 1 lib SHARED=0)
	@cp -p fdict/*.mod ./
	cp fdict/$(FDICT) ./

# This library can be found at:
#  https://github.com/zerothi/ncdf
NCDF=libncdf.a
$(NCDF): $(FDICT)
	-mkdir -p ncdf
	(cd ncdf ; \
	echo "TOP_DIR=$(VPATH)/ncdf" > Makefile ; echo "include $(VPATH)/ncdf/Makefile" >> Makefile ; \
	echo "settings.bash: FORCE" >> Makefile ; \
	printf '%s\n' "	-cp ../fdict/settings.bash ." >> Makefile ; \
	echo "FC=$(FC)" > setup.make ; \
	echo "FFLAGS=$(FFLAGS)" >> setup.make ; \
	echo "FPPFLAGS = $(FPPFLAGS)" >> setup.make ; \
	echo "INCLUDES = $(INCFLAGS) -I../" >> setup.make ; \
	echo "LIBS = $(LIBS) -lfdict" >> setup.make ; \
	$(MAKE) -j 1 copy ; $(MAKE) -j 1 lib SHARED=0)
	@cp -p ncdf/*.mod ./
	cp ncdf/$(NCDF) ./

#######################

#
# Linear algebra routines
#
libsiestaBLAS.a: $(VPATH)/Libs/blas.F
	@echo "==> Compiling libsiestaBLAS.a in Libs..."
	-mkdir -p Libs
	@(cd Libs ; $(MAKE) -f $(VPATH)/Libs/Makefile -j 1 \
		"VPATH=$(VPATH)/Libs" \
		"ARCH_MAKE=$(VPATH)/../$(OBJDIR)/arch.make" \
		"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libsiestaBLAS.a)
	@cp Libs/libsiestaBLAS.a .
libsiestaLAPACK.a: $(VPATH)/Libs/lapack.F
	@echo "==> Compiling libsiestaLAPACK.a in Libs..."
	-mkdir -p Libs
	@(cd Libs ; $(MAKE) -f $(VPATH)/Libs/Makefile -j 1 \
		"VPATH=$(VPATH)/Libs" \
		"ARCH_MAKE=$(VPATH)/../$(OBJDIR)/arch.make" \
		"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libsiestaLAPACK.a)
	@cp Libs/libsiestaLAPACK.a .


# All objects require the FDF library
# In practice they do not but it greatly simplifies the logic
# here.
# This will setup the correct compilation sequence as
# FDF => MPI
$(DEP_OBJS): $(FDF) $(COMP_LIBS) $(FDICT)
#####$(TBTOBJS): $(DEP_OBJS)

.PHONY: check-tbtrans
check-tbtrans:
	@if [ -f .phtrans ]; then \
              echo "** You were building phtrans in this directory." ;\
              echo "** Please type 'make clean' first " ;\
              false ;\
	fi
	@touch .tbtrans
tbtrans: check-tbtrans
.PHONY: tb tbt
tb: tbtrans
tbt: tbtrans


.PHONY: check-phtrans
check-phtrans:
	@if [ -f .tbtrans ]; then \
              echo "** You were building tbtrans in this directory." ;\
              echo "** Please type 'make clean' first " ;\
              false ;\
	fi
	@touch .phtrans
phtrans: check-phtrans
.PHONY: ph pht
ph: phtrans
pht: phtrans



# We append the -DTBTRANS as some special routines
# which could take advantage of this preprocessing flag
tbtrans: FPPFLAGS += -DTBTRANS
tbtrans: version $(MPI_INTERFACE) $(FDF) $(COMP_LIBS) $(TBT_COMP_LIBS) $(FDICT) $(DEP_OBJS) $(TBTOBJS) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o tbtrans \
              $(DEP_OBJS) $(TBTOBJS) $(MPI_INTERFACE) \
              $(FDF) $(COMP_LIBS) $(TBT_COMP_LIBS) $(FDICT) $(LIBS)

# We append the -DTBTRANS as some special routines
# which could take advantage of this preprocessing flag
# Tbtrans has TBT_PHONON flags to compile special capabilities
phtrans: FPPFLAGS += -DTBTRANS -DTBT_PHONON
phtrans: version $(MPI_INTERFACE) $(FDF) $(COMP_LIBS) $(TBT_COMP_LIBS) $(FDICT) $(DEP_OBJS) $(TBTOBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o phtrans \
              $(DEP_OBJS) $(TBTOBJS) $(MPI_INTERFACE) \
              $(FDF) $(COMP_LIBS) $(TBT_COMP_LIBS) $(FDICT) $(LIBS)

.PHONY: clean-tbt
clean-tbt:
	@echo "==> Cleaning tb/ph-trans only object, and executable files"
	-rm -f m_ts_chem_pot.o m_ts_elec_se.o m_ts_electype.o m_ts_gf.o
	-rm -f tbtrans phtrans $(TBTOBJS)
	-rm -f .phtrans .tbtrans

.PHONY: clean
clean: 
	@echo "==> Cleaning object, library, and executable files"
	-rm -f tbtrans phtrans *.o *.a *.mod
	-rm -f .phtrans .tbtrans
	-rm -rf ./fdf
	-rm -rf ./MPI
	-rm -rf ./fdict
	-rm -rf ./ncdf
	-rm -rf ./Libs

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(TBTOBJS:.o=.f90) $(TBTOBJS:.o=.F90) \
		$(TBTOBJS:.o=.f) $(TBTOBJS:.o=.F) || true

# Assert that FDF is compiled on before hand
$(DEP_OBJS): $(FDF)

# DO NOT DELETE THIS LINE - used by make depend
cellsubs.o: precision.o
cellsubs.o: precision.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o
class_Sparsity.o: alloc.o
class_Sparsity.o: alloc.o
create_Sparsity_SC.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
create_Sparsity_SC.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
create_Sparsity_Union.o: class_OrbitalDistribution.o class_Sparsity.o
create_Sparsity_Union.o: m_region.o parallel.o precision.o
create_Sparsity_Union.o: class_OrbitalDistribution.o class_Sparsity.o
create_Sparsity_Union.o: m_region.o parallel.o precision.o
debugmpi.o: parallel.o
debugmpi.o: parallel.o
diag.o: alloc.o diag_option.o local_sys.o parallel.o precision.o
diag.o: alloc.o diag_option.o local_sys.o parallel.o precision.o
diag_option.o: parallel.o precision.o
diag_option.o: parallel.o precision.o
fdf_extra.o: m_region.o
fdf_extra.o: m_region.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
geom_helper.o: intrinsic_missing.o
geom_helper.o: intrinsic_missing.o
idiag.o: local_sys.o parallel.o
idiag.o: local_sys.o parallel.o
io.o: local_sys.o
io.o: local_sys.o
kpoint_convert.o: local_sys.o precision.o units.o
kpoint_convert.o: local_sys.o precision.o units.o
m_geom_box.o: intrinsic_missing.o m_geom_aux.o
m_geom_box.o: intrinsic_missing.o m_geom_aux.o
m_geom_coord.o: intrinsic_missing.o m_geom_aux.o parallel.o
m_geom_coord.o: intrinsic_missing.o m_geom_aux.o parallel.o
m_geom_objects.o: m_geom_aux.o m_geom_box.o m_geom_coord.o m_geom_plane.o
m_geom_objects.o: m_geom_square.o
m_geom_objects.o: m_geom_aux.o m_geom_box.o m_geom_coord.o m_geom_plane.o
m_geom_objects.o: m_geom_square.o
m_geom_plane.o: intrinsic_missing.o m_geom_aux.o
m_geom_plane.o: intrinsic_missing.o m_geom_aux.o
m_geom_square.o: intrinsic_missing.o m_geom_aux.o m_geom_plane.o
m_geom_square.o: intrinsic_missing.o m_geom_aux.o m_geom_plane.o
m_handle_sparse.o: class_OrbitalDistribution.o class_SpData1D.o
m_handle_sparse.o: class_SpData2D.o class_Sparsity.o geom_helper.o m_iodm.o
m_handle_sparse.o: m_os.o m_ts_io.o parallel.o precision.o
m_handle_sparse.o: class_OrbitalDistribution.o class_SpData1D.o
m_handle_sparse.o: class_SpData2D.o class_Sparsity.o geom_helper.o m_iodm.o
m_handle_sparse.o: m_os.o m_ts_io.o parallel.o precision.o
m_integrate.o: precision.o
m_integrate.o: precision.o
m_io_s.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_io_s.o: class_Sparsity.o parallel.o precision.o
m_io_s.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_io_s.o: class_Sparsity.o parallel.o precision.o
m_iodm.o: class_OrbitalDistribution.o class_SpData2D.o class_Sparsity.o
m_iodm.o: m_io_s.o m_os.o parallel.o
m_iodm.o: class_OrbitalDistribution.o class_SpData2D.o class_Sparsity.o
m_iodm.o: m_io_s.o m_os.o parallel.o
m_mat_invert.o: intrinsic_missing.o m_pivot_array.o precision.o
m_mat_invert.o: intrinsic_missing.o m_pivot_array.o precision.o
m_mpi_utils.o: precision.o
m_mpi_utils.o: precision.o
m_ncdf_io.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_ncdf_io.o: class_Sparsity.o m_io_s.o parallel.o precision.o
m_ncdf_io.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_ncdf_io.o: class_Sparsity.o m_io_s.o parallel.o precision.o
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
m_spin.o: alloc.o local_sys.o parallel.o precision.o units.o
m_spin.o: alloc.o local_sys.o parallel.o precision.o units.o
m_timer.o: local_sys.o m_io.o m_walltime.o moreParallelSubs.o parallel.o
m_timer.o: precision.o
m_timer.o: local_sys.o m_io.o m_walltime.o moreParallelSubs.o parallel.o
m_timer.o: precision.o
m_trimat_invert.o: class_TriMat.o intrinsic_missing.o m_pivot_array.o
m_trimat_invert.o: precision.o
m_trimat_invert.o: class_TriMat.o intrinsic_missing.o m_pivot_array.o
m_trimat_invert.o: precision.o
m_ts_aux.o: precision.o
m_ts_aux.o: precision.o
m_ts_cctype.o: m_gauss_fermi_inf.o m_ts_io_ctype.o precision.o
m_ts_cctype.o: m_gauss_fermi_inf.o m_ts_io_ctype.o precision.o
m_ts_chem_pot.o: m_ts_cctype.o m_ts_io_ctype.o parallel.o precision.o units.o
m_ts_chem_pot.o: m_ts_cctype.o m_ts_io_ctype.o parallel.o precision.o units.o
m_ts_contour_eq.o: m_gauss_fermi_17.o m_gauss_fermi_18.o m_gauss_fermi_19.o
m_ts_contour_eq.o: m_gauss_fermi_20.o m_gauss_fermi_22.o m_gauss_fermi_24.o
m_ts_contour_eq.o: m_gauss_fermi_26.o m_gauss_fermi_28.o m_gauss_fermi_30.o
m_ts_contour_eq.o: m_gauss_fermi_inf.o m_gauss_quad.o m_integrate.o m_ts_aux.o
m_ts_contour_eq.o: m_ts_cctype.o m_ts_chem_pot.o m_ts_electype.o
m_ts_contour_eq.o: m_ts_io_contour.o m_ts_io_ctype.o parallel.o precision.o
m_ts_contour_eq.o: units.o
m_ts_contour_eq.o: m_gauss_fermi_17.o m_gauss_fermi_18.o m_gauss_fermi_19.o
m_ts_contour_eq.o: m_gauss_fermi_20.o m_gauss_fermi_22.o m_gauss_fermi_24.o
m_ts_contour_eq.o: m_gauss_fermi_26.o m_gauss_fermi_28.o m_gauss_fermi_30.o
m_ts_contour_eq.o: m_gauss_fermi_inf.o m_gauss_quad.o m_integrate.o m_ts_aux.o
m_ts_contour_eq.o: m_ts_cctype.o m_ts_chem_pot.o m_ts_electype.o
m_ts_contour_eq.o: m_ts_io_contour.o m_ts_io_ctype.o parallel.o precision.o
m_ts_contour_eq.o: units.o
m_ts_contour_neq.o: m_gauss_quad.o m_integrate.o m_ts_aux.o m_ts_cctype.o
m_ts_contour_neq.o: m_ts_chem_pot.o m_ts_electype.o m_ts_io_contour.o
m_ts_contour_neq.o: m_ts_io_ctype.o parallel.o precision.o units.o
m_ts_contour_neq.o: m_gauss_quad.o m_integrate.o m_ts_aux.o m_ts_cctype.o
m_ts_contour_neq.o: m_ts_chem_pot.o m_ts_electype.o m_ts_io_contour.o
m_ts_contour_neq.o: m_ts_io_ctype.o parallel.o precision.o units.o
m_ts_debug.o: class_Sparsity.o class_TriMat.o geom_helper.o parallel.o
m_ts_debug.o: precision.o
m_ts_debug.o: class_Sparsity.o class_TriMat.o geom_helper.o parallel.o
m_ts_debug.o: precision.o
m_ts_elec_se.o: intrinsic_missing.o m_ts_cctype.o m_ts_electype.o precision.o
m_ts_elec_se.o: units.o
m_ts_elec_se.o: intrinsic_missing.o m_ts_cctype.o m_ts_electype.o precision.o
m_ts_elec_se.o: units.o
m_ts_electrode.o: alloc.o class_SpData1D.o class_SpData2D.o class_Sparsity.o
m_ts_electrode.o: geom_helper.o intrinsic_missing.o local_sys.o m_iterator.o
m_ts_electrode.o: m_mat_invert.o m_pivot_array.o m_ts_elec_se.o m_ts_electype.o
m_ts_electrode.o: parallel.o precision.o units.o
m_ts_electrode.o: alloc.o class_SpData1D.o class_SpData2D.o class_Sparsity.o
m_ts_electrode.o: geom_helper.o intrinsic_missing.o local_sys.o m_iterator.o
m_ts_electrode.o: m_mat_invert.o m_pivot_array.o m_ts_elec_se.o m_ts_electype.o
m_ts_electrode.o: parallel.o precision.o units.o
m_ts_electype.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_ts_electype.o: class_Sparsity.o create_Sparsity_SC.o geom_helper.o
m_ts_electype.o: intrinsic_missing.o m_geom_box.o m_geom_plane.o
m_ts_electype.o: m_handle_sparse.o m_iodm.o m_os.o m_region.o m_ts_chem_pot.o
m_ts_electype.o: m_ts_io.o m_ts_io_ctype.o m_ts_iodm.o parallel.o precision.o
m_ts_electype.o: units.o
m_ts_electype.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_ts_electype.o: class_Sparsity.o create_Sparsity_SC.o geom_helper.o
m_ts_electype.o: intrinsic_missing.o m_geom_box.o m_geom_plane.o
m_ts_electype.o: m_handle_sparse.o m_iodm.o m_os.o m_region.o m_ts_chem_pot.o
m_ts_electype.o: m_ts_io.o m_ts_io_ctype.o m_ts_iodm.o parallel.o precision.o
m_ts_electype.o: units.o
m_ts_gf.o: local_sys.o m_os.o m_ts_cctype.o m_ts_contour_eq.o
m_ts_gf.o: m_ts_contour_neq.o m_ts_electrode.o m_ts_electype.o parallel.o
m_ts_gf.o: precision.o units.o
m_ts_gf.o: local_sys.o m_os.o m_ts_cctype.o m_ts_contour_eq.o
m_ts_gf.o: m_ts_contour_neq.o m_ts_electrode.o m_ts_electype.o parallel.o
m_ts_gf.o: precision.o units.o
m_ts_io.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_io.o: class_SpData2D.o class_Sparsity.o geom_helper.o local_sys.o m_io_s.o
m_ts_io.o: m_ncdf_io.o m_os.o m_sparse.o memory_log.o parallel.o precision.o
m_ts_io.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_io.o: class_SpData2D.o class_Sparsity.o geom_helper.o local_sys.o m_io_s.o
m_ts_io.o: m_ncdf_io.o m_os.o m_sparse.o memory_log.o parallel.o precision.o
m_ts_io_contour.o: precision.o units.o
m_ts_io_contour.o: precision.o units.o
m_ts_io_ctype.o: parallel.o precision.o units.o
m_ts_io_ctype.o: parallel.o precision.o units.o
m_ts_iodm.o: class_OrbitalDistribution.o class_SpData2D.o class_Sparsity.o
m_ts_iodm.o: m_io_s.o m_os.o parallel.o precision.o
m_ts_iodm.o: class_OrbitalDistribution.o class_SpData2D.o class_Sparsity.o
m_ts_iodm.o: m_io_s.o m_os.o parallel.o precision.o
m_ts_method.o: alloc.o fdf_extra.o geom_helper.o m_region.o m_ts_electype.o
m_ts_method.o: alloc.o fdf_extra.o geom_helper.o m_region.o m_ts_electype.o
m_ts_pivot.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o
m_ts_pivot.o: create_Sparsity_Union.o geom_helper.o intrinsic_missing.o
m_ts_pivot.o: m_char.o m_pivot.o m_region.o m_sparsity_handling.o m_ts_debug.o
m_ts_pivot.o: m_ts_electype.o parallel.o precision.o
m_ts_pivot.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o
m_ts_pivot.o: create_Sparsity_Union.o geom_helper.o intrinsic_missing.o
m_ts_pivot.o: m_char.o m_pivot.o m_region.o m_sparsity_handling.o m_ts_debug.o
m_ts_pivot.o: m_ts_electype.o parallel.o precision.o
m_ts_rgn2trimat.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o
m_ts_rgn2trimat.o: create_Sparsity_Union.o geom_helper.o m_region.o
m_ts_rgn2trimat.o: m_ts_electype.o m_ts_method.o m_ts_tri_common.o parallel.o
m_ts_rgn2trimat.o: precision.o
m_ts_rgn2trimat.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o
m_ts_rgn2trimat.o: create_Sparsity_Union.o geom_helper.o m_region.o
m_ts_rgn2trimat.o: m_ts_electype.o m_ts_method.o m_ts_tri_common.o parallel.o
m_ts_rgn2trimat.o: precision.o
m_ts_sparse.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_sparse.o: class_Sparsity.o create_Sparsity_SC.o geom_helper.o
m_ts_sparse.o: intrinsic_missing.o m_region.o m_sparsity_handling.o
m_ts_sparse.o: m_ts_debug.o m_ts_electype.o m_ts_method.o parallel.o
m_ts_sparse.o: precision.o
m_ts_sparse.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_sparse.o: class_Sparsity.o create_Sparsity_SC.o geom_helper.o
m_ts_sparse.o: intrinsic_missing.o m_region.o m_sparsity_handling.o
m_ts_sparse.o: m_ts_debug.o m_ts_electype.o m_ts_method.o parallel.o
m_ts_sparse.o: precision.o
m_ts_sparse_helper.o: class_OrbitalDistribution.o class_SpData1D.o
m_ts_sparse_helper.o: class_SpData1D.o class_SpData2D.o class_SpData2D.o
m_ts_sparse_helper.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
m_ts_sparse_helper.o: m_region.o m_ts_electype.o m_ts_method.o precision.o
m_ts_sparse_helper.o: class_OrbitalDistribution.o class_SpData1D.o
m_ts_sparse_helper.o: class_SpData1D.o class_SpData2D.o class_SpData2D.o
m_ts_sparse_helper.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
m_ts_sparse_helper.o: m_region.o m_ts_electype.o m_ts_method.o precision.o
m_ts_tri_common.o: intrinsic_missing.o m_region.o m_ts_electype.o precision.o
m_ts_tri_common.o: intrinsic_missing.o m_region.o m_ts_electype.o precision.o
m_ts_tri_scat.o: alloc.o class_TriMat.o m_region.o m_trimat_invert.o
m_ts_tri_scat.o: m_ts_electype.o m_ts_method.o m_ts_trimat_invert.o precision.o
m_ts_tri_scat.o: alloc.o class_TriMat.o m_region.o m_trimat_invert.o
m_ts_tri_scat.o: m_ts_electype.o m_ts_method.o m_ts_trimat_invert.o precision.o
m_ts_trimat_invert.o: class_TriMat.o m_pivot_array.o m_region.o
m_ts_trimat_invert.o: m_trimat_invert.o m_ts_electype.o m_ts_method.o
m_ts_trimat_invert.o: precision.o
m_ts_trimat_invert.o: class_TriMat.o m_pivot_array.o m_region.o
m_ts_trimat_invert.o: m_trimat_invert.o m_ts_electype.o m_ts_method.o
m_ts_trimat_invert.o: precision.o
m_wallclock.o: m_walltime.o
m_wallclock.o: m_walltime.o
memory.o: memory_log.o parallel.o
memory.o: memory_log.o parallel.o
memory_log.o: m_io.o parallel.o precision.o
memory_log.o: m_io.o parallel.o precision.o
minvec.o: cellsubs.o local_sys.o precision.o sorting.o
minvec.o: cellsubs.o local_sys.o precision.o sorting.o
moreParallelSubs.o: alloc.o local_sys.o m_io.o parallel.o precision.o
moreParallelSubs.o: alloc.o local_sys.o m_io.o parallel.o precision.o
timer.o: local_sys.o m_timer.o parallel.o timer_tree.o
timer.o: local_sys.o m_timer.o parallel.o timer_tree.o
timer_tree.o: m_walltime.o
timer_tree.o: m_walltime.o
units.o: precision.o
units.o: precision.o
m_tbt_contour.o: m_gauss_quad.o m_integrate.o m_tbt_save.o m_ts_aux.o
m_tbt_contour.o: m_ts_cctype.o m_ts_chem_pot.o m_ts_electype.o
m_tbt_contour.o: m_ts_io_contour.o m_ts_io_ctype.o parallel.o precision.o
m_tbt_contour.o: units.o
m_tbt_contour.o: m_gauss_quad.o m_integrate.o m_tbt_save.o m_ts_aux.o
m_tbt_contour.o: m_ts_cctype.o m_ts_chem_pot.o m_ts_electype.o
m_tbt_contour.o: m_ts_io_contour.o m_ts_io_ctype.o parallel.o precision.o
m_tbt_contour.o: units.o
m_tbt_dH.o: class_OrbitalDistribution.o class_SpData1D.o class_Sparsity.o
m_tbt_dH.o: class_TriMat.o intrinsic_missing.o m_ncdf_io.o m_os.o m_region.o
m_tbt_dH.o: m_sparsity_handling.o m_tbt_save.o parallel.o precision.o
m_tbt_dH.o: class_OrbitalDistribution.o class_SpData1D.o class_Sparsity.o
m_tbt_dH.o: class_TriMat.o intrinsic_missing.o m_ncdf_io.o m_os.o m_region.o
m_tbt_dH.o: m_sparsity_handling.o m_tbt_save.o parallel.o precision.o
m_tbt_diag.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_tbt_diag.o: class_Sparsity.o intrinsic_missing.o m_region.o
m_tbt_diag.o: m_ts_sparse_helper.o parallel.o precision.o
m_tbt_diag.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_tbt_diag.o: class_Sparsity.o intrinsic_missing.o m_region.o
m_tbt_diag.o: m_ts_sparse_helper.o parallel.o precision.o
m_tbt_gf.o: local_sys.o m_os.o m_tbt_contour.o m_ts_cctype.o m_ts_electrode.o
m_tbt_gf.o: m_ts_electype.o m_ts_gf.o parallel.o precision.o
m_tbt_gf.o: local_sys.o m_os.o m_tbt_contour.o m_ts_cctype.o m_ts_electrode.o
m_tbt_gf.o: m_ts_electype.o m_ts_gf.o parallel.o precision.o
m_tbt_hs.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_tbt_hs.o: class_Sparsity.o files.o m_handle_sparse.o m_interpolate.o
m_tbt_hs.o: m_sparsity_handling.o m_ts_io.o m_ts_io_ctype.o parallel.o
m_tbt_hs.o: precision.o units.o
m_tbt_hs.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_tbt_hs.o: class_Sparsity.o files.o m_handle_sparse.o m_interpolate.o
m_tbt_hs.o: m_sparsity_handling.o m_ts_io.o m_ts_io_ctype.o parallel.o
m_tbt_hs.o: precision.o units.o
m_tbt_kpoint.o: fdf_extra.o files.o find_kgrid.o intrinsic_missing.o
m_tbt_kpoint.o: m_gauss_quad.o m_integrate.o m_os.o m_tbt_save.o m_ts_tdir.o
m_tbt_kpoint.o: m_verbosity.o parallel.o precision.o units.o
m_tbt_kpoint.o: fdf_extra.o files.o find_kgrid.o intrinsic_missing.o
m_tbt_kpoint.o: m_gauss_quad.o m_integrate.o m_os.o m_tbt_save.o m_ts_tdir.o
m_tbt_kpoint.o: m_verbosity.o parallel.o precision.o units.o
m_tbt_kregions.o: class_OrbitalDistribution.o class_Sparsity.o fdf_extra.o
m_tbt_kregions.o: m_region.o m_tbt_kpoint.o m_ts_cctype.o m_ts_electrode.o
m_tbt_kregions.o: m_ts_electype.o m_verbosity.o parallel.o precision.o
m_tbt_kregions.o: class_OrbitalDistribution.o class_Sparsity.o fdf_extra.o
m_tbt_kregions.o: m_region.o m_tbt_kpoint.o m_ts_cctype.o m_ts_electrode.o
m_tbt_kregions.o: m_ts_electype.o m_verbosity.o parallel.o precision.o
m_tbt_options.o: files.o intrinsic_missing.o m_os.o m_tbt_contour.o m_tbt_dH.o
m_tbt_options.o: m_tbt_diag.o m_tbt_hs.o m_tbt_save.o m_tbt_sigma_save.o
m_tbt_options.o: m_ts_chem_pot.o m_ts_electype.o m_ts_method.o m_ts_tdir.o
m_tbt_options.o: parallel.o precision.o units.o
m_tbt_options.o: files.o intrinsic_missing.o m_os.o m_tbt_contour.o m_tbt_dH.o
m_tbt_options.o: m_tbt_diag.o m_tbt_hs.o m_tbt_save.o m_tbt_sigma_save.o
m_tbt_options.o: m_ts_chem_pot.o m_ts_electype.o m_ts_method.o m_ts_tdir.o
m_tbt_options.o: parallel.o precision.o units.o
m_tbt_proj.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData1D.o
m_tbt_proj.o: class_SpData2D.o class_Sparsity.o fdf_extra.o intrinsic_missing.o
m_tbt_proj.o: m_ncdf_io.o m_os.o m_region.o m_tbt_diag.o m_tbt_hs.o
m_tbt_proj.o: m_tbt_save.o m_ts_electype.o parallel.o precision.o timestamp.o
m_tbt_proj.o: units.o
m_tbt_proj.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData1D.o
m_tbt_proj.o: class_SpData2D.o class_Sparsity.o fdf_extra.o intrinsic_missing.o
m_tbt_proj.o: m_ncdf_io.o m_os.o m_region.o m_tbt_diag.o m_tbt_hs.o
m_tbt_proj.o: m_tbt_save.o m_ts_electype.o parallel.o precision.o timestamp.o
m_tbt_proj.o: units.o
m_tbt_regions.o: class_OrbitalDistribution.o class_Sparsity.o
m_tbt_regions.o: create_Sparsity_SC.o fdf_extra.o geom_helper.o
m_tbt_regions.o: intrinsic_missing.o m_char.o m_pivot.o m_pivot_methods.o
m_tbt_regions.o: m_region.o m_sparsity_handling.o m_ts_debug.o m_ts_electype.o
m_tbt_regions.o: m_ts_method.o m_ts_pivot.o m_ts_sparse.o m_verbosity.o
m_tbt_regions.o: parallel.o precision.o
m_tbt_regions.o: class_OrbitalDistribution.o class_Sparsity.o
m_tbt_regions.o: create_Sparsity_SC.o fdf_extra.o geom_helper.o
m_tbt_regions.o: intrinsic_missing.o m_char.o m_pivot.o m_pivot_methods.o
m_tbt_regions.o: m_region.o m_sparsity_handling.o m_ts_debug.o m_ts_electype.o
m_tbt_regions.o: m_ts_method.o m_ts_pivot.o m_ts_sparse.o m_verbosity.o
m_tbt_regions.o: parallel.o precision.o
m_tbt_save.o: class_OrbitalDistribution.o class_SpData1D.o class_Sparsity.o
m_tbt_save.o: files.o m_interpolate.o m_ncdf_io.o m_os.o m_region.o m_tbt_hs.o
m_tbt_save.o: m_ts_electype.o m_verbosity.o parallel.o precision.o timestamp.o
m_tbt_save.o: units.o
m_tbt_save.o: class_OrbitalDistribution.o class_SpData1D.o class_Sparsity.o
m_tbt_save.o: files.o m_interpolate.o m_ncdf_io.o m_os.o m_region.o m_tbt_hs.o
m_tbt_save.o: m_ts_electype.o m_verbosity.o parallel.o precision.o timestamp.o
m_tbt_save.o: units.o
m_tbt_sigma_save.o: m_os.o m_region.o m_tbt_hs.o m_tbt_save.o m_ts_electype.o
m_tbt_sigma_save.o: parallel.o timestamp.o units.o
m_tbt_sigma_save.o: m_os.o m_region.o m_tbt_hs.o m_tbt_save.o m_ts_electype.o
m_tbt_sigma_save.o: parallel.o timestamp.o units.o
m_tbt_sparse_helper.o: class_OrbitalDistribution.o class_SpData1D.o
m_tbt_sparse_helper.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
m_tbt_sparse_helper.o: m_region.o m_tbt_kregions.o m_ts_electype.o
m_tbt_sparse_helper.o: m_ts_method.o m_ts_sparse_helper.o parallel.o
m_tbt_sparse_helper.o: precision.o
m_tbt_sparse_helper.o: class_OrbitalDistribution.o class_SpData1D.o
m_tbt_sparse_helper.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
m_tbt_sparse_helper.o: m_region.o m_tbt_kregions.o m_ts_electype.o
m_tbt_sparse_helper.o: m_ts_method.o m_ts_sparse_helper.o parallel.o
m_tbt_sparse_helper.o: precision.o
m_tbt_tri_init.o: class_OrbitalDistribution.o class_Sparsity.o
m_tbt_tri_init.o: create_Sparsity_Union.o m_region.o m_sparsity_handling.o
m_tbt_tri_init.o: m_tbt_options.o m_tbt_regions.o m_ts_debug.o m_ts_electype.o
m_tbt_tri_init.o: m_ts_method.o m_ts_rgn2trimat.o m_ts_tri_common.o
m_tbt_tri_init.o: m_verbosity.o parallel.o precision.o
m_tbt_tri_init.o: class_OrbitalDistribution.o class_Sparsity.o
m_tbt_tri_init.o: create_Sparsity_Union.o m_region.o m_sparsity_handling.o
m_tbt_tri_init.o: m_tbt_options.o m_tbt_regions.o m_ts_debug.o m_ts_electype.o
m_tbt_tri_init.o: m_ts_method.o m_ts_rgn2trimat.o m_ts_tri_common.o
m_tbt_tri_init.o: m_verbosity.o parallel.o precision.o
m_tbt_tri_scat.o: class_SpData1D.o class_SpData1D.o class_Sparsity.o
m_tbt_tri_scat.o: class_TriMat.o geom_helper.o intrinsic_missing.o m_region.o
m_tbt_tri_scat.o: m_tbt_proj.o m_ts_cctype.o m_ts_electype.o m_ts_tri_common.o
m_tbt_tri_scat.o: m_ts_tri_scat.o m_ts_trimat_invert.o precision.o units.o
m_tbt_tri_scat.o: class_SpData1D.o class_SpData1D.o class_Sparsity.o
m_tbt_tri_scat.o: class_TriMat.o geom_helper.o intrinsic_missing.o m_region.o
m_tbt_tri_scat.o: m_tbt_proj.o m_ts_cctype.o m_ts_electype.o m_ts_tri_common.o
m_tbt_tri_scat.o: m_ts_tri_scat.o m_ts_trimat_invert.o precision.o units.o
m_tbt_trik.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData1D.o
m_tbt_trik.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o class_TriMat.o
m_tbt_trik.o: intrinsic_missing.o m_iterator.o m_mat_invert.o m_pivot_array.o
m_tbt_trik.o: m_region.o m_tbt_contour.o m_tbt_dH.o m_tbt_gf.o m_tbt_hs.o
m_tbt_trik.o: m_tbt_kpoint.o m_tbt_kregions.o m_tbt_options.o m_tbt_proj.o
m_tbt_trik.o: m_tbt_regions.o m_tbt_save.o m_tbt_sigma_save.o
m_tbt_trik.o: m_tbt_sparse_helper.o m_tbt_tri_init.o m_tbt_tri_scat.o m_timer.o
m_tbt_trik.o: m_trimat_invert.o m_ts_cctype.o m_ts_elec_se.o m_ts_electype.o
m_tbt_trik.o: m_ts_method.o m_ts_sparse_helper.o m_ts_tri_common.o
m_tbt_trik.o: m_ts_trimat_invert.o m_verbosity.o parallel.o precision.o units.o
m_tbt_trik.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData1D.o
m_tbt_trik.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o class_TriMat.o
m_tbt_trik.o: intrinsic_missing.o m_iterator.o m_mat_invert.o m_pivot_array.o
m_tbt_trik.o: m_region.o m_tbt_contour.o m_tbt_dH.o m_tbt_gf.o m_tbt_hs.o
m_tbt_trik.o: m_tbt_kpoint.o m_tbt_kregions.o m_tbt_options.o m_tbt_proj.o
m_tbt_trik.o: m_tbt_regions.o m_tbt_save.o m_tbt_sigma_save.o
m_tbt_trik.o: m_tbt_sparse_helper.o m_tbt_tri_init.o m_tbt_tri_scat.o m_timer.o
m_tbt_trik.o: m_trimat_invert.o m_ts_cctype.o m_ts_elec_se.o m_ts_electype.o
m_tbt_trik.o: m_ts_method.o m_ts_sparse_helper.o m_ts_tri_common.o
m_tbt_trik.o: m_ts_trimat_invert.o m_verbosity.o parallel.o precision.o units.o
m_tbtrans.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_tbtrans.o: class_SpData2D.o class_Sparsity.o m_region.o m_sparsity_handling.o
m_tbtrans.o: m_tbt_contour.o m_tbt_dH.o m_tbt_hs.o m_tbt_kpoint.o
m_tbtrans.o: m_tbt_kregions.o m_tbt_options.o m_tbt_proj.o m_tbt_regions.o
m_tbtrans.o: m_tbt_save.o m_tbt_sigma_save.o m_tbt_tri_init.o m_tbt_trik.o
m_tbtrans.o: m_ts_electype.o m_ts_gf.o m_verbosity.o parallel.o precision.o
m_tbtrans.o: units.o
m_tbtrans.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_tbtrans.o: class_SpData2D.o class_Sparsity.o m_region.o m_sparsity_handling.o
m_tbtrans.o: m_tbt_contour.o m_tbt_dH.o m_tbt_hs.o m_tbt_kpoint.o
m_tbtrans.o: m_tbt_kregions.o m_tbt_options.o m_tbt_proj.o m_tbt_regions.o
m_tbtrans.o: m_tbt_save.o m_tbt_sigma_save.o m_tbt_tri_init.o m_tbt_trik.o
m_tbtrans.o: m_ts_electype.o m_ts_gf.o m_verbosity.o parallel.o precision.o
m_tbtrans.o: units.o
tbt_end.o: m_wallclock.o memory_log.o parallel.o timestamp.o
tbt_end.o: m_wallclock.o memory_log.o parallel.o timestamp.o
tbt_init.o: class_SpData1D.o class_SpData2D.o class_Sparsity.o files.o
tbt_init.o: local_sys.o m_sparsity_handling.o m_spin.o m_tbt_contour.o
tbt_init.o: m_tbt_dH.o m_tbt_gf.o m_tbt_hs.o m_tbt_kpoint.o m_tbt_kregions.o
tbt_init.o: m_tbt_options.o m_tbt_proj.o m_tbt_regions.o m_tbt_save.o m_timer.o
tbt_init.o: m_ts_electrode.o m_ts_electype.o m_wallclock.o memory_log.o
tbt_init.o: parallel.o precision.o timestamp.o
tbt_init.o: class_SpData1D.o class_SpData2D.o class_Sparsity.o files.o
tbt_init.o: local_sys.o m_sparsity_handling.o m_spin.o m_tbt_contour.o
tbt_init.o: m_tbt_dH.o m_tbt_gf.o m_tbt_hs.o m_tbt_kpoint.o m_tbt_kregions.o
tbt_init.o: m_tbt_options.o m_tbt_proj.o m_tbt_regions.o m_tbt_save.o m_timer.o
tbt_init.o: m_ts_electrode.o m_ts_electype.o m_wallclock.o memory_log.o
tbt_init.o: parallel.o precision.o timestamp.o
tbt_reinit.o: m_verbosity.o parallel.o
tbt_reinit.o: m_verbosity.o parallel.o
tbtrans.o: m_tbt_hs.o m_tbt_options.o m_tbtrans.o
tbtrans.o: m_tbt_hs.o m_tbt_options.o m_tbtrans.o
class_cdata1d.o: class_Data1D.o
class_cdata2d.o: class_Data2D.o
class_cspdata1d.o: class_SpData1D.o
class_cspdata2d.o: class_SpData2D.o
class_ctrimat.o: class_TriMat.o
class_ddata1d.o: class_Data1D.o
class_ddata2d.o: class_Data2D.o
class_dspdata1d.o: class_SpData1D.o
class_dspdata2d.o: class_SpData2D.o
class_dtrimat.o: class_TriMat.o
class_gdata1d.o: class_Data1D.o
class_gspdata1d.o: class_SpData1D.o
class_idata1d.o: class_Data1D.o
class_idata2d.o: class_Data2D.o
class_ispdata1d.o: class_SpData1D.o
class_ispdata2d.o: class_SpData2D.o
class_itrimat.o: class_TriMat.o
class_ldata1d.o: class_Data1D.o
class_ldata2d.o: class_Data2D.o
class_lspdata1d.o: class_SpData1D.o
class_lspdata2d.o: class_SpData2D.o
class_ltrimat.o: class_TriMat.o
class_sdata1d.o: class_Data1D.o
class_sdata2d.o: class_Data2D.o
class_sspdata1d.o: class_SpData1D.o
class_sspdata2d.o: class_SpData2D.o
class_strimat.o: class_TriMat.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
class_ztrimat.o: class_TriMat.o
m_diag.o: diag.o
m_diag_option.o: diag_option.o
m_find_kgrid.o: find_kgrid.o
m_minvec.o: minvec.o
m_object_debug.o: object_debug.o
m_timer_tree.o: timer_tree.o
m_timestamp.o: timestamp.o
mtprng.o: m_uuid.o
sys.o: local_sys.o
t_spin.o: m_spin.o
timer_options.o: timer.o
version_info.o: version.o
