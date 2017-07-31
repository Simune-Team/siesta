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

VPATH:=$(shell pwd)/../../Src

default: what mixps fractional

MAIN_OBJDIR:=$(shell cd -P -- ../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:=$(DEFS_PREFIX)-UMPI $(DEFS) $(DEFS_PREFIX)-UMPI
#
include $(PSML_ROOT)/psml.mk
include $(XMLF90_ROOT)/xmlf90.mk
#
FPPFLAGS:=$(DEFS_PREFIX)-UMPI $(FPPFLAGS) $(DEFS_PREFIX)-UMPI
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src

INCFLAGS:=$(NETCDF_INCFLAGS)\
          $(NCPS_INCFLAGS) \
          $(PSML_INCFLAGS) \
          $(INCFLAGS)
what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2

SYSOBJ=$(SYS).o

#------------------------------------------------------------------------
#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF):
	(mkdir -p fdf; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

MIXPS_DEPS= f2kcli.o interpolation.o atom_options.o \
            m_io.o io.o alloc.o periodic_table.o  precision.o
MIXPS_LOCAL_OBJS= mixps.o handlers.o local_sys.o
MIXPS_OBJS= $(MIXPS_DEPS) $(MIXPS_LOCAL_OBJS)
#
NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
#
mixps:  $(FDF) $(MIXPS_OBJS) $(NCPS) 
	$(FC) -o mixps \
	       $(LDFLAGS) $(MIXPS_OBJS) $(NCPS) $(FDF) $(PSML_LIBS) $(XMLF90_LIBS)

#----------------------------------------------------------------------------
FRACTIONAL_DEPS= f2kcli.o m_io.o alloc.o interpolation.o \
                 atom_options.o io.o periodic_table.o precision.o
FRACTIONAL_LOCAL_OBJS= fractional.o handlers.o local_sys.o
FRACTIONAL_OBJS= $(FRACTIONAL_DEPS) $(FRACTIONAL_LOCAL_OBJS)

fractional: $(FDF) $(FRACTIONAL_OBJS) $(NCPS) 
	$(FC) -o fractional \
	       $(LDFLAGS) $(FRACTIONAL_OBJS) $(NCPS) $(PSML) $(FDF)\
               $(PSML_LIBS) $(XMLF90_LIBS)
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f mixps fractional *.o  *.a
	rm -f *.mod
	rm -rf ./fdf

DEP_OBJS= $(MIXPS_DEPS) $(FRACTIONAL_DEPS)
LOCAL_OBJS=$(FRACTIONAL_LOCAL_OBJS) $(MIXPS_LOCAL_OBJS)
dep:
	sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) \
                || true

# DO NOT DELETE THIS LINE - used by make depend
atom_options.o: local_sys.o
atom_options.o: local_sys.o
io.o: local_sys.o
io.o: local_sys.o
m_io.o: local_sys.o
m_io.o: local_sys.o
periodic_table.o: local_sys.o precision.o
periodic_table.o: local_sys.o precision.o
fractional.o: f2kcli.o periodic_table.o precision.o
fractional.o: f2kcli.o periodic_table.o precision.o
mixps.o: f2kcli.o interpolation.o periodic_table.o precision.o
mixps.o: f2kcli.o interpolation.o periodic_table.o precision.o
sys.o: local_sys.o
