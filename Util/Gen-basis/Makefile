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

ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:= $(DEFS) $(DEFS_PREFIX) -UMPI
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX) -UMPI
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)

SYSOBJ=$(SYS).o

#
# Note that machine-specific files are now in top Src directory.
#
OBJS_GEN-BASIS=m_io.o io.o alloc.o basis_types.o precision.o  parallel.o \
          parsing.o basis_io.o  chemical.o atm_transfer.o atm_types.o\
          atmparams.o old_atmfuncs.o radial.o interpolation.o \
	      bessph.o radfft.o sorting.o m_filter.o \
          basis_specs.o atom.o memoryinfo.o memory.o periodic_table.o\
          pseudopotential.o pxf.o dot.o  atom_options.o ldau_specs.o arw.o  \
          timer.o  xml.o m_walltime.o read_xc_info.o gen-basis.o $(SYSOBJ) \
          xc.o bsc_xcmod.o m_cite.o local_sys.o

OBJS_IONCAT=f2kcli.o m_getopts.o basis_types.o precision.o  parallel.o \
          parsing.o m_io.o io.o alloc.o atom_options.o basis_io.o atm_transfer.o atm_types.o \
          atmparams.o radial.o  memory.o chkdim.o \
	      pseudopotential.o chemical.o m_filter.o \
          basis_specs.o atom.o memoryinfo.o  m_memory.o periodic_table.o pxf.o \
          atmfuncs.o spher_harm.o interpolation.o old_atmfuncs.o arw.o \
          local_sys.o xc.o bsc_xcmod.o xml.o m_walltime.o ioncat.o ldau_specs.o m_cite.o $(SYSOBJ)


#
# Compile XC in a private directory, without MPI
#
XC=libSiestaXC.a
XC_MAKEFILE=$(VPATH)/SiestaXC/makefile
XC_INCFLAGS:=-I $(VPATH)/SiestaXC ##$(INCFLAGS)
$(XC): 
	(mkdir -p SiestaXC ; cd SiestaXC ; $(MAKE) -f $(XC_MAKEFILE) \
                          "FC=$(FC)" "VPATH_ROOT=$(VPATH)" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(XC_INCFLAGS)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "FFLAGS=$(FFLAGS)"  module )

#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
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

gen-basis: $(FDF) $(XC) $(OBJS_GEN-BASIS) link_linalg
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(OBJS_GEN-BASIS) $(FDF) $(XC) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)

ioncat:  $(FDF) $(XC) $(OBJS_IONCAT)
	$(FC) -o ioncat \
	       $(LDFLAGS) $(OBJS_IONCAT) $(FDF) $(XC) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f gen-basis ioncat *.o  *.a
	rm -f *.mod
	rm -rf fdf
	rm -rf SiestaXC
	rm -f _tmp_deps deps.list  protomake*

ALL_OBJS=$(OBJS_IONCAT) $(OBJS_GEN-BASIS)

dep:
	-sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
		*.f *.f90 *.F *.F90
	@sed '/^siesta_cmlsubs.o:/s/version.o//'  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile

# DO NOT DELETE THIS LINE - used by make depend
