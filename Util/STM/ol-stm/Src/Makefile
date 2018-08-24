# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for Ordejon-Lorente STM program
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
#------------- NOTE
# Special needs for the use of FFTW
# You must set these symbols properly, either here or 
# (better) in your top arch.make
# 
FFTW_ROOT ?=/usr/local
FFTW_INCFLAGS ?= -I$(FFTW_ROOT)/include
FFTW_LIBS ?= -L$(FFTW_ROOT)/lib -lfftw3f -lfftw3
#-----------------------------------------

OBJDIR=Obj

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

VPATH:=$(shell pwd)/../../../../Src

default: what stm

dep:
	-sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
                $(VPATH)/*.T90 \
		*.f *.f90 *.F *.F90
	@sed -i -e "/^siesta_cmlsubs.o:/s/version.o//" Makefile
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile


ARCH_MAKE=../../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

STM_UNDEF  = $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX) -UCDF
STM_UNDEF += $(DEFS_PREFIX)-UCDF
STM_UNDEF += $(DEFS_PREFIX)-USIESTA__METIS
STM_UNDEF += $(DEFS_PREFIX)-UON_DOMAIN_DECOMP

DEFS := $(DEFS) $(STM_UNDEF)
FPPFLAGS := $(FPPFLAGS) $(STM_UNDEF)

#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)

what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@if [ -z "$(FFTW_LIBS)" ] ; then \
           echo "You must set the FFTW symbols in your arch.make" ; \
           false ; \
        fi
	@sleep 2

SYSOBJ=$(SYS).o

# Note that machine-specific files are now in top Src directory.

OBJS_STM = atompla.o  \
        bessph.o chkdim.o dismin.o \
	dot.o fftw3_mymod.o extrapolate.o \
        io.o iodm.o length.o \
        mainstm.o matvect.o memory.o \
	neighb.o pxf.o \
	radfft.o ranger.o \
        readstm.o m_wfsx.o \
        redata.o  reinitstm.o stm.o \
	vacpot.o volcel.o handlers.o local_sys.o

INCFLAGS:= $(FFTW_INCFLAGS) $(NETCDF_INCFLAGS) $(INCFLAGS)

# Use the makefile in Src/fdf and all the sources there.

FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)"\
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

# This is crude but will have to do for now.
# Note : precision must be the first module

MOD_OBJS=precision.o parallel.o m_io.o alloc.o memoryinfo.o listsc.o  \
         atmparams.o atom_options.o atmfuncs.o atm_types.o \
         radial.o spher_harm.o basis_io.o basis_types.o \
         domain_decom.o schecomm.o printmatrix.o qsort.o mmio.o pspltm1.o lenstr.o \
         sparse_matrices.o interpolation.o m_fft_gpfa.o \
         files.o pseudopotential.o chemical.o xml.o \
         parallelsubs.o spatial.o

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

COM_OBJS_STM=$(OBJS_STM) $(SYSOBJ)
ALL_OBJS_STM= $(CLS_OBJS) $(MOD_OBJS) $(COM_OBJS_STM)

stm: $(FDF) $(ALL_OBJS_STM)
	$(FC) -o stm \
	       $(LDFLAGS) $(ALL_OBJS_STM) $(FDF) $(FFTW_LIBS) 

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f stm *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)

# DO NOT DELETE THIS LINE - used by make depend
