# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone get_chem_labels
#
OBJDIR=Obj
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
VPATH:=$(shell pwd)/../../Src
#
default: get_chem_labels
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
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
SIESTA_OBJS=precision.o f2kcli.o parallel.o \
          chemical.o io.o m_io.o alloc.o \
          memory_log.o memory.o periodic_table.o\
          pxf.o atom_options.o  \
          $(SYSOBJ)

LOCAL_OBJS=get_chem_labels.o handlers.o local_sys.o
OBJS_GET_CHEM_LABELS=$(SIESTA_OBJS) $(LOCAL_OBJS)
#
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
get_chem_labels: $(FDF)  $(OBJS_GET_CHEM_LABELS)
	$(FC) -o get_chem_labels \
	       $(LDFLAGS) $(OBJS_GET_CHEM_LABELS) $(FDF) $(LIBS)
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f get_chem_labels ioncat *.o  *.a
	rm -f *.mod .tmp_fdflog
	rm -rf fdf

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.f) $(SIESTA_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.F) $(SIESTA_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true

# DO NOT DELETE THIS LINE - used by make depend
atom_options.o: local_sys.o
atom_options.o: local_sys.o
chemical.o: local_sys.o parallel.o precision.o
chemical.o: local_sys.o parallel.o precision.o
io.o: local_sys.o
io.o: local_sys.o
memory.o: memory_log.o parallel.o
memory.o: memory_log.o parallel.o
memory_log.o: m_io.o parallel.o precision.o
memory_log.o: m_io.o parallel.o precision.o
periodic_table.o: local_sys.o precision.o
periodic_table.o: local_sys.o precision.o
get_chem_labels.o: chemical.o f2kcli.o
get_chem_labels.o: chemical.o f2kcli.o
sys.o: local_sys.o
