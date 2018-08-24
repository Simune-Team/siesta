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

SIESTA_OBJS= precision.o local_sys.o m_io.o alloc.o memoryinfo.o memory.o pxf.o \
             moreParallelSubs.o m_mpi_utils.o \
             reclat.o idiag.o minvec.o \
             parallel.o io.o files.o cellsubs.o sorting.o \
	     kpoint_convert.o timestamp.o m_wallclock.o m_spin.o \
             diag.o diag_option.o \
	     debugmpi.o find_kgrid.o

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

# Output stuff
SIESTA_OBJS += init_output.o

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
TBTOBJS += sorted_search_t.o

# Different solution methods
TBTOBJS += m_tbtrans.o m_tbt_trik.o m_tbt_tri_scat.o
TBTOBJS += m_tbt_proj.o m_tbt_diag.o
TBTOBJS += m_tbt_delta.o m_tbt_dH.o m_tbt_dSE.o

# The saving methods
TBTOBJS += m_tbt_sigma_save.o m_tbt_save.o

# The k-region files
TBTOBJS += m_tbt_kregions.o m_tbt_sparse_helper.o

# The main program
TBTOBJS += tbtrans.o

# Extract information of the compiler version used.
# This is very basic but for intel, gcc, pgi and open64
# this at least contain the version string of the compiler.
# Possibly later the used compiler version may be extracted.
COMPILER_VERSION := $(shell $(FC) --version | sed '/^[ ]*$$/d;q')
ifeq ($(COMPILER_VERSION),)
  # Just try and capture the case were it can not figure out anything
  COMPILER_VERSION := Unknown version
endif

# The version program
.PHONY: version.o
version.o: version
version: version.F90
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;\
		s^COMPILER_VERSION^$(COMPILER_VERSION)^;\
		s'TBTRANS_VERSION'$$(cat ${<D}/../../../version.info)'g;\
		s'LIBS'$(COMP_LIBS) $(LIBS)'g;\
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
	echo "AR=$(AR)" >> setup.make ; \
	echo "RANLIB=$(RANLIB)" >> setup.make ; \
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
$(TBTOBJS): $(DEP_OBJS)

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
	-sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
		$(VPATH)/*.T90 \
		*.f90 *.F90 *.f *.F
	@sed -i -e "/^siesta_cmlsubs.o:/s/version.o//" Makefile
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# Assert that FDF is compiled on before hand
$(DEP_OBJS): $(FDF)

# DO NOT DELETE THIS LINE - used by make depend
