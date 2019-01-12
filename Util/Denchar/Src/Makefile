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

MAIN_OBJDIR:=$(shell cd -P -- ../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make

override WITH_MPI=
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC := $(FC_SERIAL)         # Make it non-recursive
#

NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
INCFLAGS:=$(NCPS_INCFLAGS) $(PSML_INCFLAGS) $(INCFLAGS)

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

# Use the makefile in Src/fdf and all the sources there.
FDF=libfdf.a
FDF_MAKEFILE:=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	-mkdir -p fdf
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
		"ARCH_MAKE=$(ARCH_MAKE)" \
		"DEFS=$(DEFS)" \
		"FPPFLAGS=$(FPPFLAGS)" \
		"MPI_INTERFACE= " \
		"INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

# This is crude but will have to do for now.
# Note : precision must be the first module

# Note that machine-specific files are now in top Src directory.
SIESTA_OBJS =  bessph.o chkdim.o dismin.o \
               dot.o  io.o iodm.o memory.o memory_log.o\
	       pxf.o radfft.o m_fft_gpfa.o volcel.o

       
SIESTA_OBJS += precision.o parallel.o debugmpi.o m_io.o alloc.o listsc.o  \
         atmparams.o atom_options.o atmfuncs.o atm_types.o \
         radial.o spher_harm.o basis_io.o basis_types.o parallelsubs.o \
         schecomm.o sparse_matrices.o\
         files.o interpolation.o chemical.o xml.o \
         domain_decom.o printmatrix.o qsort.o mmio.o pspltm1.o \
         spatial.o 

# Additional class objects
SIESTA_OBJS += m_uuid.o object_debug.o class_Sparsity.o \
	class_OrbitalDistribution.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
        class_Geometry.o \
	class_Pair_Data1D.o \
	class_Fstack_Pair_Data1D.o \
	class_Pair_Geometry_SpData2D.o \
	class_Fstack_Pair_Geometry_SpData2D.o

LOCAL_OBJS= atompla.o colinear.o denchar.o itochar.o length.o \
            local_die.o local_reinit.o local_sys.o matvect.o neighb.o \
            planed.o ranger.o readpla.o readsts.o readwaves.o \
            readwavesx.o redata_denchar.o rhoofr.o stsofr.o \
	    timer_local.o wavofr.o wrout.o alloc_handlers.o

ALL_OBJS_DENCHAR= $(SIESTA_OBJS) $(LOCAL_OBJS)

NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
#
denchar: $(FDF)  $(ALL_OBJS_DENCHAR) $(NCPS) 
	$(FC) -o denchar \
	       $(LDFLAGS) $(ALL_OBJS_DENCHAR) $(FDF) $(NCPS) $(PSML_LIBS) \
               $(XMLF90_LIBS) $(METIS_LIB)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f denchar *.o  *.a
	rm -f *.mod
	rm -f _tmp_deps deps.list _deps.list* protomake*
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
debugmpi.o: parallel.o
debugmpi.o: parallel.o
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
spher_harm.o: alloc.o local_sys.o precision.o
spher_harm.o: alloc.o local_sys.o precision.o
xml.o: precision.o
xml.o: precision.o
atompla.o: local_sys.o precision.o
atompla.o: local_sys.o precision.o
colinear.o: precision.o
colinear.o: precision.o
denchar.o: basis_io.o listsc.o parallel.o precision.o
denchar.o: basis_io.o listsc.o parallel.o precision.o
length.o: precision.o
length.o: precision.o
local_reinit.o: files.o
local_reinit.o: files.o
matvect.o: precision.o
matvect.o: precision.o
neighb.o: precision.o
neighb.o: precision.o
planed.o: precision.o
planed.o: precision.o
ranger.o: precision.o
ranger.o: precision.o
readpla.o: precision.o
readpla.o: precision.o
readsts.o: precision.o
readsts.o: precision.o
readwaves.o: precision.o
readwaves.o: precision.o
readwavesx.o: precision.o
readwavesx.o: precision.o
redata_denchar.o: precision.o
redata_denchar.o: precision.o
rhoofr.o: atmfuncs.o chemical.o listsc.o planed.o precision.o
rhoofr.o: atmfuncs.o chemical.o listsc.o planed.o precision.o
stsofr.o: atmfuncs.o chemical.o listsc.o precision.o
stsofr.o: atmfuncs.o chemical.o listsc.o precision.o
wavofr.o: atmfuncs.o chemical.o listsc.o planed.o precision.o
wavofr.o: atmfuncs.o chemical.o listsc.o planed.o precision.o
wrout.o: precision.o
wrout.o: precision.o
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
listsc_module.o: listsc.o
m_bessph.o: bessph.o
m_object_debug.o: object_debug.o
m_radfft.o: radfft.o
mtprng.o: m_uuid.o
printmat.o: printmatrix.o
sys.o: local_sys.o
