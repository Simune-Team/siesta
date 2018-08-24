# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone Denchar
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

OBJDIR=Obj

.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

VPATH:=$(shell pwd)/../../../Src

default: what denchar

dep:
	-sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
		$(VPATH)/*.T90 \
		*.f *.f90 *.F *.F90
	@sed -i -e "/^siesta_cmlsubs.o:/s/version.o//" Makefile
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile


ARCH_MAKE=../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC := $(FC_SERIAL)         # Make it non-recursive

DENCHAR_UNDEF  = $(DEFS_PREFIX)-UMPI
DENCHAR_UNDEF += $(DEFS_PREFIX)-UCDF
DENCHAR_UNDEF += $(DEFS_PREFIX)-USIESTA__METIS
DENCHAR_UNDEF += $(DEFS_PREFIX)-UON_DOMAIN_DECOMP

DEFS := $(DEFS) $(DENCHAR_UNDEF)
FPPFLAGS := $(FPPFLAGS) $(DENCHAR_UNDEF)

# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)

what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2

SYSOBJ=$(SYS).o

# Note that machine-specific files are now in top Src directory.
OBJS_DENCHAR = atompla.o planed.o\
        bessph.o chkdim.o dismin.o \
	dot.o readpla.o readsts.o \
        io.o iodm.o itochar.o length.o \
        matvect.o memory.o colinear.o\
	neighb.o pxf.o \
	radfft.o m_fft_gpfa.o ranger.o \
        readwavesx.o rhoofr.o wavofr.o \
        redata_denchar.o local_reinit.o denchar.o \
	volcel.o wrout.o stsofr.o local_sys.o

# Use the makefile in Src/fdf and all the sources there.
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	-mkdir -p fdf
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
		"ARCH_MAKE=../$(ARCH_MAKE)" \
		"DEFS=$(DEFS)" \
		"FPPFLAGS=$(FPPFLAGS)" \
		"MPI_INTERFACE= " \
		"INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

# This is crude but will have to do for now.
# Note : precision must be the first module

MOD_OBJS=precision.o parallel.o debugmpi.o m_io.o alloc.o memoryinfo.o listsc.o  \
         atmparams.o atom_options.o atmfuncs.o atm_types.o \
         radial.o spher_harm.o basis_io.o basis_types.o parallelsubs.o \
         schecomm.o sparse_matrices.o interpolation.o\
         files.o pseudopotential.o chemical.o xml.o \
         domain_decom.o printmatrix.o qsort.o mmio.o pspltm1.o lenstr.o \
         spatial.o timer_local.o 

# Additional class objects
CLS_OBJS = m_uuid.o object_debug.o class_Sparsity.o \
	class_OrbitalDistribution.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
        class_Geometry.o \
	class_Pair_Data1D.o \
	class_Fstack_Pair_Data1D.o \
	class_Pair_Geometry_SpData2D.o \
	class_Fstack_Pair_Geometry_SpData2D.o

COM_OBJS_DENCHAR=$(OBJS_DENCHAR) $(SYSOBJ)
ALL_OBJS_DENCHAR= $(CLS_OBJS) $(MOD_OBJS) $(COM_OBJS_DENCHAR)


denchar: $(FDF) $(ALL_OBJS_DENCHAR)
	$(FC) -o denchar \
	       $(LDFLAGS) $(ALL_OBJS_DENCHAR) $(FDF) $(METIS_LIB)


clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f denchar *.o  *.a
	rm -f *.mod
	rm -f _tmp_deps deps.list _deps.list* protomake*
	rm -rf fdf

# DO NOT DELETE THIS LINE - used by make depend
