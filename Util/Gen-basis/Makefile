# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone Gen-basis and ioncat
#
# The idea is to use the code in the top Src directory as much as possible.
# This is achieved by the VPATH directive below.
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

OBJDIR=Obj

.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

VPATH:=$(shell pwd)/../../Src

default: gen-basis ioncat

MAIN_OBJDIR:=$(shell cd -P -- ../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

GRIDXC_SERIAL_ROOT ?= $(GRIDXC_ROOT)
#
DEFS:= $(DEFS) $(DEFS_PREFIX) -UMPI
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX) -UMPI
#
include $(PSML_ROOT)/psml.mk
include $(XMLF90_ROOT)/xmlf90.mk
#
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src
INCFLAGS:=$(NETCDF_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSOP_INCFLAGS) $(PSML_INCFLAGS) \
          $(INCFLAGS)
SYSOBJ=$(SYS).o

#
# Note that machine-specific files are now in top Src directory.
#
DEPS_GEN-BASIS=m_io.o io.o alloc.o basis_types.o precision.o  parallel.o \
          parsing.o basis_io.o chemical.o atm_transfer.o atm_types.o\
          atmparams.o old_atmfuncs.o radial.o interpolation.o \
	  bessph.o radfft.o  sorting.o m_filter.o m_spin.o \
          basis_specs.o ldau_specs.o atom.o periodic_table.o\
          pxf.o dot.o  atom_options.o arw.o  m_spin_orbit_potentials.o \
          sys.o xml.o m_walltime.o read_xc_info.o  $(SYSOBJ) \
          xc.o bsc_xcmod.o  m_fft_gpfa.o m_cite.o reclat.o units.o \
          semicore_info_froyen.o hamann.o
GEN-BASIS_LOCAL_OBJS=gen-basis.o handlers.o
OBJS_GEN-BASIS=$(DEPS_GEN-BASIS) $(GEN-BASIS_LOCAL_OBJS)

DEPS_IONCAT=f2kcli.o m_getopts.o basis_types.o precision.o  parallel.o \
          parsing.o m_io.o io.o alloc.o atom_options.o basis_io.o atm_transfer.o atm_types.o \
          atmparams.o radial.o  chkdim.o m_filter.o \
	  chemical.o m_spin.o m_spin_orbit_potentials.o reclat.o \
          basis_specs.o periodic_table.o pxf.o semicore_info_froyen.o \
          atom.o atmfuncs.o spher_harm.o interpolation.o old_atmfuncs.o arw.o \
          sys.o xc.o bsc_xcmod.o xml.o m_walltime.o ldau_specs.o\
          hamann.o m_cite.o  units.o $(SYSOBJ)
IONCAT_LOCAL_OBJS=ioncat.o handlers.o
OBJS_IONCAT=$(DEPS_IONCAT) $(IONCAT_LOCAL_OBJS)


#
include $(GRIDXC_SERIAL_ROOT)/gridxc.mk
XC = $(GRIDXC_LIBS)
INCFLAGS:= $(INCFLAGS) $(GRIDXC_INCFLAGS)

#
NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
PSOP=$(MAIN_OBJDIR)/psoplib/src/libpsop.a
#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

#
# Simple hack to link to any compiled lapack/blas code
# If COMP_LIBS is empty, we do a no-op, else link the library file here.
#
link_linalg:
	for llibs in $(COMP_LIBS); do \
   [ -f ../../$(OBJDIR)/$$llibs ] && ln -sf ../../$(OBJDIR)/$$llibs . ; \
	done

gen-basis: $(FDF) $(NCPS) $(PSOP) $(OBJS_GEN-BASIS) link_linalg
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(OBJS_GEN-BASIS) $(FDF) $(XC) \
              $(NCPS) $(PSOP) $(PSML_LIBS) $(XMLF90_LIBS) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)
ioncat:  $(FDF) $(NCPS) $(PSOP) $(OBJS_IONCAT)
	$(FC) -o ioncat \
	       $(LDFLAGS) $(OBJS_IONCAT) $(FDF) \
              $(NCPS) $(PSOP) $(XC) $(PSML_LIBS) $(XMLF90_LIBS) $(NETCDF_LIBS) $(LIBS)
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f gen-basis ioncat *.o  *.a
	rm -f *.mod
	rm -rf fdf
	rm -f _tmp_deps deps.list  protomake*

ALL_OBJS=$(OBJS_IONCAT) $(OBJS_GEN-BASIS)

DEP_OBJS= $(DEPS_GEN-BASIS) $(DEPS_IONCAT)
LOCAL_OBJS=$(IONCAT_LOCAL_OBJS) $(GEN-BASIS_LOCAL_OBJS)

dep:
	sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) \
                || true

# DO NOT DELETE THIS LINE - used by make depend
arw.o: alloc.o parallel.o precision.o sys.o
arw.o: alloc.o parallel.o precision.o sys.o
atm_transfer.o: atm_types.o atmparams.o ldau_specs.o old_atmfuncs.o
atm_transfer.o: periodic_table.o radial.o sys.o
atm_transfer.o: atm_types.o atmparams.o ldau_specs.o old_atmfuncs.o
atm_transfer.o: periodic_table.o radial.o sys.o
atm_types.o: precision.o radial.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o bsc_xcmod.o
atom.o: interpolation.o m_filter.o old_atmfuncs.o periodic_table.o precision.o
atom.o: sys.o hamann.o
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o bsc_xcmod.o
atom.o: interpolation.o m_filter.o old_atmfuncs.o periodic_table.o precision.o
atom.o: sys.o
atom_options.o: sys.o
atom_options.o: sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o radial.o sys.o xml.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o radial.o sys.o xml.o
basis_specs.o: basis_types.o chemical.o m_spin.o m_spin_orbit_potentials.o
basis_specs.o: periodic_table.o precision.o semicore_info_froyen.o sys.o
basis_specs.o: basis_types.o chemical.o m_spin.o m_spin_orbit_potentials.o
basis_specs.o: periodic_table.o precision.o semicore_info_froyen.o sys.o
basis_types.o: alloc.o atmparams.o precision.o sys.o
basis_types.o: alloc.o atmparams.o precision.o sys.o
bessph.o: precision.o sys.o
bessph.o: precision.o sys.o
bsc_xcmod.o: parallel.o precision.o sys.o
bsc_xcmod.o: parallel.o precision.o sys.o
chemical.o: parallel.o precision.o sys.o
chemical.o: parallel.o precision.o sys.o
chkdim.o: sys.o
chkdim.o: sys.o
io.o: sys.o
io.o: sys.o
ldau_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
ldau_specs.o: basis_specs.o basis_types.o interpolation.o m_cite.o parallel.o
ldau_specs.o: precision.o radial.o sys.o units.o
ldau_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
ldau_specs.o: basis_specs.o basis_types.o interpolation.o m_cite.o parallel.o
ldau_specs.o: precision.o radial.o sys.o units.o
m_filter.o: bessph.o precision.o radfft.o sorting.o sys.o
m_filter.o: bessph.o precision.o radfft.o sorting.o sys.o
m_getopts.o: f2kcli.o
m_getopts.o: f2kcli.o
m_io.o: sys.o
m_io.o: sys.o
m_spin.o: alloc.o parallel.o precision.o sys.o units.o
m_spin.o: alloc.o parallel.o precision.o sys.o units.o
old_atmfuncs.o: alloc.o atmparams.o precision.o sys.o
old_atmfuncs.o: alloc.o atmparams.o precision.o sys.o
periodic_table.o: precision.o sys.o
periodic_table.o: precision.o sys.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
radial.o: alloc.o interpolation.o precision.o xml.o
read_xc_info.o: parallel.o precision.o sys.o
read_xc_info.o: parallel.o precision.o sys.o
spher_harm.o: alloc.o precision.o sys.o
spher_harm.o: alloc.o precision.o sys.o
units.o: precision.o
units.o: precision.o
xc.o: alloc.o bsc_xcmod.o precision.o sys.o
xc.o: alloc.o bsc_xcmod.o precision.o sys.o
xml.o: precision.o
xml.o: precision.o
gen-basis.o: atom.o atom_options.o basis_io.o basis_specs.o basis_types.o
gen-basis.o: old_atmfuncs.o sys.o
gen-basis.o: atom.o atom_options.o basis_io.o basis_specs.o basis_types.o
gen-basis.o: old_atmfuncs.o sys.o
ioncat.o: atm_types.o atmfuncs.o basis_io.o m_getopts.o precision.o
ioncat.o: atm_types.o atmfuncs.o basis_io.o m_getopts.o precision.o
aux_proj.o: atom.o
m_bessph.o: bessph.o
m_radfft.o: radfft.o
m_semicore_info_froyen.o: semicore_info_froyen.o
t_spin.o: m_spin.o
