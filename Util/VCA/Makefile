# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone mixps and fractional
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
#  The dependency list at the end is overly large, but harmless
#
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
VPATH:=$(shell pwd)/../../Src
#
default: what mixps fractional
#
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS_PREFIX)-UMPI $(DEFS) $(DEFS_PREFIX)-UMPI
FPPFLAGS:=$(DEFS_PREFIX)-UMPI $(FPPFLAGS) $(DEFS_PREFIX)-UMPI
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2
#
SYSOBJ=$(SYS).o
#
#------------------------------------------------------------------------
#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF):
	(mkdir -p fdf; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
MIXPS_OBJS= f2kcli.o interpolation.o atom_options.o pseudopotential.o \
            local_sys.o  m_io.o io.o alloc.o pxf.o parallel.o periodic_table.o \
            memory.o memoryinfo.o mixps.o

mixps:  $(FDF) $(MIXPS_OBJS)
	$(FC) -o mixps \
	       $(LDFLAGS) $(MIXPS_OBJS)   $(FDF) 
#----------------------------------------------------------------------------
FRACTIONAL_OBJS= f2kcli.o m_io.o alloc.o interpolation.o fractional.o \
                 atom_options.o  memory.o memoryinfo.o pseudopotential.o  \
                 io.o pxf.o parallel.o periodic_table.o local_sys.o

fractional: $(FDF) $(FRACTIONAL_OBJS)
	$(FC) -o fractional \
	       $(LDFLAGS) $(FRACTIONAL_OBJS)  $(FDF) 
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f mixps fractional *.o  *.a
	rm -f *.mod
	rm -rf ./fdf
	rm -f _tmp_deps deps.list  protomake*
#
ALL_OBJS= $(MIXPS_OBJS) $(FRACTIONAL_OBJS)
#
dep:
	rm -f _tmp_deps deps.list  protomake*
	sfmakedepend --depend=obj --file _tmp_defs  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
          *.f *.F *.f90 *.F90
	@echo "Ignore errors in the following command"
	@echo "They appear if the last grep does not match anything"
	@-for i in $(ALL_OBJS) ; do grep "^$$i: " _tmp_defs  ; done > deps.list
	@echo "ok"
	@csplit -k -f "protomake" Makefile '/^# -- specially prepared/+1'
	@cat protomake00 deps.list >| Makefile
#
# -- specially prepared dependency list -- do not remove this line!
m_io.o: local_sys.o
io.o: local_sys.o
alloc.o: m_io.o parallel.o precision.o local_sys.o
periodic_table.o: precision.o local_sys.o
mixps.o: f2kcli.o interpolation.o periodic_table.o precision.o
mixps.o: pseudopotential.o
fractional.o: f2kcli.o periodic_table.o precision.o pseudopotential.o
atom_options.o: local_sys.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
memoryinfo.o: precision.o
pseudopotential.o: alloc.o atom_options.o interpolation.o precision.o local_sys.o

