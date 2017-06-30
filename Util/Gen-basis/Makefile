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
#  
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
TOP_LEVEL=$(shell pwd)/../../Src
VPATH:=$(TOP_LEVEL)
#
default: gen-basis ioncat
#
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:= $(DEFS) $(DEFS_PREFIX)-UMPI
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS_GEN-BASIS=m_io.o io.o alloc.o basis_types.o precision.o  parallel.o \
          parsing.o basis_io.o  chemical.o atm_transfer.o atm_types.o\
          atmparams.o old_atmfuncs.o radial.o interpolation.o \
	  paste.o   bessph.o radfft.o  sorting.o m_filter.o \
          basis_specs.o atom.o memoryinfo.o memory.o periodic_table.o\
          pseudopotential.o pxf.o dot.o  atom_options.o arw.o  \
          timer.o  xml.o m_walltime.o read_xc_info.o gen-basis.o $(SYSOBJ) \
          xc.o bsc_xcmod.o local_sys.o
#
IONCAT_OBJS=f2kcli.o m_getopts.o basis_types.o precision.o  parallel.o \
          parsing.o m_io.o io.o alloc.o atom_options.o basis_io.o  atm_types.o\
          atmparams.o radial.o  memory.o chkdim.o\
	  paste.o pseudopotential.o chemical.o \
          basis_specs.o memoryinfo.o  m_memory.o periodic_table.o pxf.o \
           atmfuncs.o spher_harm.o interpolation.o\
          xml.o m_walltime.o ioncat.o  local_sys.o $(SYSOBJ)
#
# Compile XC in a private directory, without MPI
#
XC=libSiestaXC.a
XC_MAKEFILE=$(TOP_LEVEL)/SiestaXC/makefile
XC_INCFLAGS:=-I $(TOP_LEVEL)/SiestaXC ##$(INCFLAGS)
$(XC): 
	(mkdir -p SiestaXC ; cd SiestaXC ; $(MAKE) -f $(XC_MAKEFILE) \
                          "FC=$(FC)" "VPATH_ROOT=$(TOP_LEVEL)" \
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
FDF_MAKEFILE=$(TOP_LEVEL)/fdf/makefile
FDF_INCFLAGS:=-I $(TOP_LEVEL)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(TOP_LEVEL)/fdf" \
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
#
gen-basis: $(FDF) $(XC) $(OBJS_GEN-BASIS) link_linalg
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(OBJS_GEN-BASIS) $(FDF) $(XC) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)
#
ioncat:  $(FDF) $(IONCAT_OBJS)
	$(FC) -o ioncat \
	       $(LDFLAGS) $(IONCAT_OBJS) $(FDF) $(NETCDF_LIBS)
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f gen-basis ioncat *.o  *.a
	rm -f *.mod
	rm -rf fdf
	rm -rf SiestaXC
	rm -f _tmp_deps deps.list  protomake*
#
ALL_OBJS=$(IONCAT_OBJS) $(OBJS_GEN-BASIS)
dep:
	rm -f _tmp_deps deps.list  protomake*
	sfmakedepend --depend=obj --file _tmp_defs  --modext=o \
          $(TOP_LEVEL)/*.f $(TOP_LEVEL)/*.f90 $(TOP_LEVEL)/*.F $(TOP_LEVEL)/*.F90 \
          *.f *.F *.f90 *.F90
	@echo "Ignore errors in the following command"
	@echo "They appear if the last grep does not match anything"
	@-for i in $(ALL_OBJS) ; do grep "^$$i: " _tmp_defs  ; done > deps.list
	@echo "ok"
	@csplit -k -f "protomake" Makefile '/^# -- specially prepared/+1'
	@cat protomake00 deps.list >| Makefile
#
# -- specially prepared dependency list -- do not remove this line!
m_getopts.o: f2kcli.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o local_sys.o
m_io.o: local_sys.o
io.o: local_sys.o
alloc.o: debugmpi.o m_io.o parallel.o precision.o local_sys.o
atom_options.o: local_sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o local_sys.o xml.o
atm_types.o: precision.o radial.o
radial.o: alloc.o interpolation.o precision.o xml.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
chkdim.o: local_sys.o
pseudopotential.o: alloc.o atom_options.o interpolation.o precision.o local_sys.o
chemical.o: parallel.o precision.o local_sys.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o local_sys.o
memoryinfo.o: precision.o
m_memory.o: memoryinfo.o parallel.o local_sys.o
periodic_table.o: precision.o local_sys.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o local_sys.o
spher_harm.o: alloc.o precision.o local_sys.o
xml.o: precision.o
ioncat.o: atm_types.o atmfuncs.o basis_io.o m_getopts.o precision.o
m_io.o: local_sys.o
io.o: local_sys.o
alloc.o: debugmpi.o m_io.o parallel.o precision.o local_sys.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o local_sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o local_sys.o xml.o
chemical.o: parallel.o precision.o local_sys.o
atm_transfer.o: atm_types.o atmparams.o old_atmfuncs.o periodic_table.o
atm_transfer.o: radial.o local_sys.o
atm_types.o: precision.o radial.o
old_atmfuncs.o: alloc.o atmparams.o precision.o local_sys.o
radial.o: alloc.o interpolation.o precision.o xml.o
bessph.o: precision.o local_sys.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
m_filter.o: bessph.o precision.o radfft.o sorting.o local_sys.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o local_sys.o
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o bsc_xcmod.o
atom.o: interpolation.o m_filter.o old_atmfuncs.o periodic_table.o precision.o
atom.o: pseudopotential.o local_sys.o
memoryinfo.o: precision.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
periodic_table.o: precision.o local_sys.o
pseudopotential.o: alloc.o atom_options.o interpolation.o precision.o local_sys.o
atom_options.o: local_sys.o
arw.o: alloc.o parallel.o precision.o local_sys.o
timer.o: m_timer.o parallel.o local_sys.o timer_tree.o
xml.o: precision.o
read_xc_info.o: parallel.o precision.o local_sys.o
gen-basis.o: atom.o atom_options.o basis_io.o basis_specs.o basis_types.o
gen-basis.o: old_atmfuncs.o parallel.o precision.o local_sys.o
xc.o: alloc.o bsc_xcmod.o precision.o local_sys.o
bsc_xcmod.o: parallel.o precision.o local_sys.o
