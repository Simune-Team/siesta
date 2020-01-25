# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---

# Makefile for unfold

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

.PHONY: default what dep version clean-tbt clean

override WITH_NCDF=
override WITH_FLOOK=

default: what unfold 

VPATH:=$(shell pwd)/../../../Src
OBJDIR=Obj
MAIN_OBJDIR:=$(shell cd -P -- ../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
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


NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src
INCFLAGS:=$(NETCDF_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSML_INCFLAGS) $(PSOP_INCFLAGS) \
          $(INCFLAGS)

XC = $(GRIDXC_LIBS)
NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
PSOP=$(MAIN_OBJDIR)/psoplib/src/libpsop.a

INCFLAGS += $(GRIDXC_INCFLAGS)

SYSOBJ=$(SYS).o

SIESTA_OBJS= bye.o precision.o m_io.o alloc.o memory_log.o memory.o pxf.o \
             moreParallelSubs.o m_mpi_utils.o atmfuncs.o atm_types.o \
             radial.o reclat.o idiag.o minvec.o \
             parallel.o io.o files.o cellsubs.o sorting.o \
	     kpoint_convert.o timestamp.o m_wallclock.o m_spin.o \
             diag.o diag_option.o bloch_unfold.o \
	     debugmpi.o find_kgrid.o units.o posix_calls.o xml.o \
             interpolation.o spher_harm.o basis_types.o atmparams.o basis_io.o \
             atom_options.o chemical.o \
             get_kpoints_scale.o siesta_geom.o \
             radfft.o bessph.o m_fft_gpfa.o \
             broadcast_basis.o ldau_specs.o basis_specs.o \
             atom.o hamann.o m_filter.o \
             m_cite.o m_spin_orbit_potentials.o \
             coor.o chkdim.o semicore_info_froyen.o \
             zmatrix.o m_cell.o  periodic_table.o volcel.o digcel.o \
             arw.o redcel.o
# Siesta timings
SIESTA_OBJS += m_walltime.o timer_tree.o m_timer.o timer.o

# Add specific class objects
# class Objects
SIESTA_OBJS += class_OrbitalDistribution.o \
	class_Sparsity.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
        m_uuid.o object_debug.o

# Output stuff
SIESTA_OBJS += init_output.o

# As the classes are based on inputs we force the dependency to
# alloc
class_Sparsity.o class_Data1D.o class_Data2D.o: alloc.o

DEP_OBJS = $(SIESTA_OBJS) $(TS_OBJS) $(SYSOBJ) $(EXTRA)

# To separate unfold objects from siesta/transiesta
UNFOLD_OBJS  = array.o hsx_m.o
UNFOLD_OBJS += unfold.o 

# Handlers and resolvers of last resort
UNFOLD_OBJS += handlers.o local_sys.o
#
# Extract information of the compiler version used.
# This is very basic but for intel, gcc, pgi and open64
# this at least contain the version string of the compiler.
# Possibly later the used compiler version may be extracted.
COMPILER_VERSION := $(shell $(FC) --version | sed '/^[ ]*$$/d;q')
ifeq ($(COMPILER_VERSION),)
  # Just try and capture the case were it can not figure out anything
  COMPILER_VERSION := Unknown version
endif

# Interfaces to libraries
MPI_MAKEFILE=$(VPATH)/MPI/Makefile
libmpi_f90.a: 
	-mkdir -p MPI
	(cd MPI ; \
	$(MAKE) -f $(MPI_MAKEFILE) -j 1 "VPATH=$(VPATH)/MPI" \
	"MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=$(ARCH_MAKE)" \
	"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module_built)

FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf -I../MPI $(INCFLAGS)
$(FDF): $(MPI_INTERFACE)
	-mkdir -p fdf
	(cd fdf ; \
	$(MAKE) -f $(FDF_MAKEFILE) -j 1 "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
	"ARCH_MAKE=$(ARCH_MAKE)" \
	"INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module)

#######################
#
# Linear algebra routines
#
libsiestaBLAS.a: $(VPATH)/Libs/blas.F
	@echo "==> Compiling libsiestaBLAS.a in Libs..."
	-mkdir -p Libs
	@(cd Libs ; $(MAKE) -f $(VPATH)/Libs/Makefile -j 1 \
		"VPATH=$(VPATH)/Libs" \
		"ARCH_MAKE=$(ARCH_MAKE)" \
		"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libsiestaBLAS.a)
	@cp Libs/libsiestaBLAS.a .
libsiestaLAPACK.a: $(VPATH)/Libs/lapack.F
	@echo "==> Compiling libsiestaLAPACK.a in Libs..."
	-mkdir -p Libs
	@(cd Libs ; $(MAKE) -f $(VPATH)/Libs/Makefile -j 1 \
		"VPATH=$(VPATH)/Libs" \
		"ARCH_MAKE=$(ARCH_MAKE)" \
		"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libsiestaLAPACK.a)
	@cp Libs/libsiestaLAPACK.a .


# All objects require the FDF library
# In practice they do not but it greatly simplifies the logic
# here.
# This will setup the correct compilation sequence as
# FDF => MPI
$(DEP_OBJS): $(FDF) $(COMP_LIBS) 

unfold: $(MPI_INTERFACE) $(FDF) $(COMP_LIBS) $(DEP_OBJS) $(UNFOLD_OBJS) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o unfold \
              $(DEP_OBJS) $(UNFOLD_OBJS) $(MPI_INTERFACE) \
              $(FDF) $(XC) $(NCPS) $(PSOP) $(PSML_LIBS) $(XMLF90_LIBS) $(COMP_LIBS) $(LIBS)

.PHONY: clean
clean: 
	@echo "==> Cleaning object, library, and executable files"
	-rm -f unfold *.o *.a *.mod
	-rm -rf ./fdf
	-rm -rf ./MPI
	-rm -rf ./Libs

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(VPATH)/class*.T90 \
		$(UNFOLD_OBJS:.o=.f90) $(UNFOLD_OBJS:.o=.F90) \
		$(UNFOLD_OBJS:.o=.f) $(UNFOLD_OBJS:.o=.F) || true

	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# Assert that FDF is compiled on before hand
$(DEP_OBJS): $(FDF)

# DO NOT DELETE THIS LINE - used by make depend
arw.o: alloc.o local_sys.o parallel.o precision.o
arw.o: alloc.o local_sys.o parallel.o precision.o
atm_types.o: precision.o radial.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o local_sys.o precision.o radial.o spher_harm.o
atmfuncs.o: atm_types.o local_sys.o precision.o radial.o spher_harm.o
atom.o: atm_types.o atmparams.o atom_options.o basis_specs.o basis_types.o
atom.o: hamann.o interpolation.o local_sys.o m_filter.o periodic_table.o
atom.o: precision.o radial.o
atom.o: atm_types.o atmparams.o atom_options.o basis_specs.o basis_types.o
atom.o: hamann.o interpolation.o local_sys.o m_filter.o periodic_table.o
atom.o: precision.o radial.o
atom_options.o: local_sys.o
atom_options.o: local_sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: local_sys.o precision.o radial.o xml.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: local_sys.o precision.o radial.o xml.o
basis_specs.o: basis_types.o chemical.o local_sys.o m_spin.o
basis_specs.o: m_spin_orbit_potentials.o periodic_table.o precision.o
basis_specs.o: semicore_info_froyen.o
basis_specs.o: basis_types.o chemical.o local_sys.o m_spin.o
basis_specs.o: m_spin_orbit_potentials.o periodic_table.o precision.o
basis_specs.o: semicore_info_froyen.o
basis_types.o: alloc.o atmparams.o local_sys.o precision.o
basis_types.o: alloc.o atmparams.o local_sys.o precision.o
bessph.o: local_sys.o precision.o
bessph.o: local_sys.o precision.o
bloch_unfold.o: units.o
bloch_unfold.o: units.o
broadcast_basis.o: atm_types.o chemical.o ldau_specs.o m_mpi_utils.o parallel.o
broadcast_basis.o: radial.o
broadcast_basis.o: atm_types.o chemical.o ldau_specs.o m_mpi_utils.o parallel.o
broadcast_basis.o: radial.o
cellsubs.o: precision.o
cellsubs.o: precision.o
chemical.o: local_sys.o parallel.o precision.o
chemical.o: local_sys.o parallel.o precision.o
chkdim.o: local_sys.o
chkdim.o: local_sys.o
class_Data1D.o: alloc.o
class_Data2D.o: alloc.o
class_Data3D.o: alloc.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o
class_SpData1D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData3D.o: class_OrbitalDistribution.o class_Sparsity.o
class_Sparsity.o: alloc.o
class_Sparsity.o: alloc.o
class_TriMat.o: alloc.o
coor.o: alloc.o local_sys.o parallel.o precision.o siesta_geom.o units.o
coor.o: zmatrix.o
coor.o: alloc.o local_sys.o parallel.o precision.o siesta_geom.o units.o
coor.o: zmatrix.o
debugmpi.o: parallel.o
debugmpi.o: parallel.o
diag.o: alloc.o diag_option.o local_sys.o parallel.o precision.o
diag.o: alloc.o diag_option.o local_sys.o parallel.o precision.o
diag_option.o: parallel.o precision.o
diag_option.o: parallel.o precision.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
get_kpoints_scale.o: parallel.o siesta_geom.o
get_kpoints_scale.o: parallel.o siesta_geom.o
idiag.o: local_sys.o parallel.o
idiag.o: local_sys.o parallel.o
init_output.o: files.o
init_output.o: files.o
io.o: m_io.o
io.o: m_io.o
kpoint_convert.o: local_sys.o precision.o units.o
kpoint_convert.o: local_sys.o precision.o units.o
ldau_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
ldau_specs.o: basis_specs.o basis_types.o interpolation.o local_sys.o m_cite.o
ldau_specs.o: parallel.o precision.o radial.o units.o
ldau_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
ldau_specs.o: basis_specs.o basis_types.o interpolation.o local_sys.o m_cite.o
ldau_specs.o: parallel.o precision.o radial.o units.o
m_cell.o: precision.o siesta_geom.o units.o
m_cell.o: precision.o siesta_geom.o units.o
m_filter.o: bessph.o local_sys.o precision.o radfft.o
m_filter.o: bessph.o local_sys.o precision.o radfft.o
m_mpi_utils.o: precision.o
m_mpi_utils.o: precision.o
m_spin.o: alloc.o files.o local_sys.o m_cite.o parallel.o precision.o units.o
m_spin.o: alloc.o files.o local_sys.o m_cite.o parallel.o precision.o units.o
m_timer.o: local_sys.o m_io.o m_walltime.o moreParallelSubs.o parallel.o
m_timer.o: precision.o
m_timer.o: local_sys.o m_io.o m_walltime.o moreParallelSubs.o parallel.o
m_timer.o: precision.o
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
periodic_table.o: local_sys.o precision.o
periodic_table.o: local_sys.o precision.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
radial.o: alloc.o interpolation.o precision.o xml.o
redcel.o: local_sys.o parallel.o precision.o
redcel.o: local_sys.o parallel.o precision.o
siesta_geom.o: precision.o
siesta_geom.o: precision.o
spher_harm.o: alloc.o local_sys.o precision.o
spher_harm.o: alloc.o local_sys.o precision.o
timer.o: local_sys.o m_timer.o parallel.o timer_tree.o
timer.o: local_sys.o m_timer.o parallel.o timer_tree.o
timer_tree.o: m_walltime.o
timer_tree.o: m_walltime.o
units.o: precision.o
units.o: precision.o
xml.o: precision.o
xml.o: precision.o
zmatrix.o: alloc.o local_sys.o m_cell.o parallel.o precision.o siesta_geom.o
zmatrix.o: units.o
zmatrix.o: alloc.o local_sys.o m_cell.o parallel.o precision.o siesta_geom.o
zmatrix.o: units.o
array.o: local_sys.o precision.o
array.o: local_sys.o precision.o
unfold.o: alloc.o array.o atmfuncs.o basis_io.o basis_types.o cellsubs.o
unfold.o: diag_option.o get_kpoints_scale.o hsx_m.o local_sys.o m_io.o
unfold.o: m_mpi_utils.o m_timer.o memory_log.o parallel.o precision.o radfft.o
unfold.o: siesta_geom.o spher_harm.o
unfold.o: alloc.o array.o atmfuncs.o basis_io.o basis_types.o cellsubs.o
unfold.o: diag_option.o get_kpoints_scale.o hsx_m.o local_sys.o m_io.o
unfold.o: m_mpi_utils.o m_timer.o memory_log.o parallel.o precision.o radfft.o
unfold.o: siesta_geom.o spher_harm.o
aux_proj.o: atom.o
bloch_unfold_m.o: bloch_unfold.o
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
class_lspdata1d.o: class_SpData1D.o
class_sdata1d.o: class_Data1D.o
class_sspdata1d.o: class_SpData1D.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
m_array.o: array.o
m_bessph.o: bessph.o
m_diag.o: diag.o
m_diag_option.o: diag_option.o
m_find_kgrid.o: find_kgrid.o
m_get_kpoints_scale.o: get_kpoints_scale.o
m_hamann.o: hamann.o
m_minvec.o: minvec.o
m_object_debug.o: object_debug.o
m_radfft.o: radfft.o
m_semicore_info_froyen.o: semicore_info_froyen.o
m_timer_tree.o: timer_tree.o
m_timestamp.o: timestamp.o
mtprng.o: m_uuid.o
sys.o: local_sys.o
t_spin.o: m_spin.o
timer_options.o: timer.o
