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

OBJDIR=Obj

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

VPATH:=$(shell pwd)/../../../../Src
VPATH_GRID:=$(shell pwd)/../../../../Util/Grid

default: what stm

MAIN_OBJDIR:=$(shell cd -P -- ../../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make

include $(ARCH_MAKE)

#------------- NOTE
# Special needs for the use of FFTW
# You must set these symbols properly, either here or 
# (better) in your top arch.make
# 
#FFTW_ROOT ?=/usr/local
FFTW_INCFLAGS ?= -I$(FFTW_ROOT)/include
FFTW_LIBS ?= -L$(FFTW_ROOT)/lib -lfftw3f -lfftw3
#-----------------------------------------

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

STM_UNDEF  = $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX) -UCDF
STM_UNDEF += $(DEFS_PREFIX)-UCDF
STM_UNDEF += $(DEFS_PREFIX)-USIESTA__METIS
STM_UNDEF += $(DEFS_PREFIX)-UON_DOMAIN_DECOMP

#
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


NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
#PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src
INCFLAGS:=$(NETCDF_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSOP_INCFLAGS) \
          $(PSML_INCFLAGS) \
          $(INCFLAGS)
INCFLAGS:= $(FFTW_INCFLAGS) $(INCFLAGS)
#
NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
# Use the makefile in Src/fdf and all the sources there.

FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)"\
                          "ARCH_MAKE=$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

LOCAL_OBJS = fftw3_mymod.o extrapolate.o \
        mainstm.o neighb.o ranger.o \
        readstm.o redata.o  reinitstm.o stm.o sts.o \
	vacpot.o handlers.o local_sys.o

SIESTA_OBJS= bessph.o chkdim.o dismin.o dot.o io.o iodm.o memory.o \
             pxf.o radfft.o volcel.o

SIESTA_OBJS +=precision.o parallel.o m_io.o alloc.o listsc.o  \
         atmparams.o atom_options.o atmfuncs.o atm_types.o \
         radial.o spher_harm.o basis_io.o basis_types.o \
         domain_decom.o schecomm.o printmatrix.o qsort.o mmio.o pspltm1.o \
         sparse_matrices.o interpolation.o m_fft_gpfa.o \
         files.o   chemical.o xml.o \
         parallelsubs.o spatial.o memory_log.o m_gridfunc.o units.o

# Additional class objects
# These are just needed because sparse_matrices has "delete" operations for
# some pieces of data. (And sparse_matrices is needed because of the index arrays)
SIESTA_OBJS += m_uuid.o object_debug.o class_Sparsity.o \
	class_OrbitalDistribution.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
        class_Geometry.o \
	class_Pair_Data1D.o \
	class_Fstack_Pair_Data1D.o \
	class_Pair_Geometry_SpData2D.o \
	class_Fstack_Pair_Geometry_SpData2D.o

ALL_OBJS_STM= $(SIESTA_OBJS) $(LOCAL_OBJS)

m_gridfunc.o: $(VPATH_GRID)/m_gridfunc.F90
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
vacpot.o: m_gridfunc.o

stm: $(FDF)  $(NCPS) $(ALL_OBJS_STM)
	$(FC) -o stm \
	       $(LDFLAGS) $(ALL_OBJS_STM) $(FDF) $(NCPS) \
               $(PSML_LIBS) $(XMLF90_LIBS) $(FFTW_LIBS) 
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f stm *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "ARCH_MAKE=$(ARCH_MAKE)" clean)

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.f) $(SIESTA_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.F) $(SIESTA_OBJS:.o=.F90)) \
		$(VPATH)/class*.T90 \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# DO NOT DELETE THIS LINE - used by make depend
atm_types.o: precision.o radial.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o local_sys.o precision.o radial.o spher_harm.o
atmfuncs.o: atm_types.o local_sys.o precision.o radial.o spher_harm.o
atom_options.o: local_sys.o
atom_options.o: local_sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: local_sys.o precision.o radial.o xml.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: local_sys.o precision.o radial.o xml.o
basis_types.o: alloc.o atmparams.o local_sys.o precision.o
basis_types.o: alloc.o atmparams.o local_sys.o precision.o
bessph.o: local_sys.o precision.o
bessph.o: local_sys.o precision.o
chemical.o: local_sys.o parallel.o precision.o
chemical.o: local_sys.o parallel.o precision.o
chkdim.o: local_sys.o
chkdim.o: local_sys.o
class_Data1D.o: alloc.o
class_Data2D.o: alloc.o
class_Data3D.o: alloc.o
class_Fstack_Pair_Data1D.o: class_Pair_Data1D.o class_Pair_Data1D.o
class_Fstack_Pair_Data1D.o: class_Pair_Data1D.o class_Pair_Data1D.o
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
class_Geometry.o: alloc.o
class_Geometry.o: alloc.o
class_Pair_Data1D.o: class_Data1D.o class_Data1D.o
class_Pair_Data1D.o: class_Data1D.o class_Data1D.o
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
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
domain_decom.o: alloc.o local_sys.o parallel.o precision.o printmatrix.o
domain_decom.o: schecomm.o sparse_matrices.o
domain_decom.o: alloc.o local_sys.o parallel.o precision.o printmatrix.o
domain_decom.o: schecomm.o sparse_matrices.o
io.o: m_io.o
io.o: m_io.o
iodm.o: alloc.o files.o local_sys.o parallel.o parallelsubs.o precision.o
iodm.o: alloc.o files.o local_sys.o parallel.o parallelsubs.o precision.o
listsc.o: alloc.o
listsc.o: alloc.o
memory.o: memory_log.o parallel.o
memory.o: memory_log.o parallel.o
memory_log.o: m_io.o parallel.o precision.o
memory_log.o: m_io.o parallel.o precision.o
parallelsubs.o: domain_decom.o local_sys.o parallel.o spatial.o
parallelsubs.o: domain_decom.o local_sys.o parallel.o spatial.o
printmatrix.o: alloc.o
printmatrix.o: alloc.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
radial.o: alloc.o interpolation.o precision.o xml.o
schecomm.o: alloc.o
schecomm.o: alloc.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpData2D.o
sparse_matrices.o: class_OrbitalDistribution.o class_SpData1D.o
sparse_matrices.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o
sparse_matrices.o: precision.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpData2D.o
sparse_matrices.o: class_OrbitalDistribution.o class_SpData1D.o
sparse_matrices.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o
sparse_matrices.o: precision.o
spatial.o: precision.o
spatial.o: precision.o
spher_harm.o: precision.o
units.o: precision.o
xml.o: precision.o
xml.o: precision.o
extrapolate.o: fftw3_mymod.o precision.o
extrapolate.o: fftw3_mymod.o precision.o
handlers.o: memory_log.o
handlers.o: memory_log.o
mainstm.o: basis_io.o listsc.o precision.o sts.o
mainstm.o: basis_io.o listsc.o precision.o sts.o
ranger.o: precision.o
ranger.o: precision.o
stm.o: atmfuncs.o chemical.o precision.o
stm.o: atmfuncs.o chemical.o precision.o
sts.o: atmfuncs.o chemical.o precision.o
sts.o: atmfuncs.o chemical.o precision.o
vacpot.o: precision.o units.o
vacpot.o: precision.o units.o
broadening.o: sts.o
class_ddata1d.o: class_Data1D.o
class_ddata2d.o: class_Data2D.o
class_dspdata1d.o: class_SpData1D.o
class_dspdata2d.o: class_SpData2D.o
class_fstack_pair_ddata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_gdata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_geometry_dspdata2d.o: class_Fstack_Pair_Geometry_SpData2D.o
class_fstack_pair_sdata1d.o: class_Fstack_Pair_Data1D.o
class_gdata1d.o: class_Data1D.o
class_gspdata1d.o: class_SpData1D.o
class_idata1d.o: class_Data1D.o
class_idata2d.o: class_Data2D.o
class_ispdata1d.o: class_SpData1D.o
class_ispdata2d.o: class_SpData2D.o
class_ldata1d.o: class_Data1D.o
class_lspdata1d.o: class_SpData1D.o
class_pair_ddata1d.o: class_Pair_Data1D.o
class_pair_gdata1d.o: class_Pair_Data1D.o
class_pair_geometry_dspdata2d.o: class_Pair_Geometry_SpData2D.o
class_pair_sdata1d.o: class_Pair_Data1D.o
class_sdata1d.o: class_Data1D.o
class_sspdata1d.o: class_SpData1D.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
listsc_module.o: listsc.o
m_bessph.o: bessph.o
m_object_debug.o: object_debug.o
m_radfft.o: radfft.o
mtprng.o: m_uuid.o
printmat.o: printmatrix.o
sys.o: local_sys.o
