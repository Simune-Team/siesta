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
FFTW_ROOT=/usr/local
FFTW_INCFLAGS= -I$(FFTW_ROOT)/include
FFTW_LIBS = -L$(FFTW_ROOT)/lib -lfftw3f -lfftw3
#-----------------------------------------
#
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH:=$(shell pwd)/../../../../Src
#
default: what stm
#
ARCH_MAKE=../../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
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
	@if [ -z "$(FFTW_LIBS)" ] ; then \
           echo "You must set the FFTW symbols in your arch.make" ; \
           false ; \
        fi
	@sleep 2
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS_STM = atompla.o  \
        bessph.o chkdim.o dismin.o \
	dot.o fftw3_mymod.o extrapolate.o \
        io.o iodm.o length.o \
        mainstm.o matvect.o memory.o \
	neighb.o paste.o pxf.o \
	 radfft.o ranger.o \
        readstm.o m_wfsx.o \
        redata.o  reinitstm.o stm.o \
	vacpot.o volcel.o handlers.o local_sys.o
#
INCFLAGS:= $(FFTW_INCFLAGS) $(NETCDF_INCFLAGS) $(INCFLAGS)
#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)"\
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=precision.o parallel.o m_io.o alloc.o memoryinfo.o listsc.o  \
         atmparams.o atom_options.o atmfuncs.o   atm_types.o \
         radial.o spher_harm.o  basis_io.o basis_types.o \
         domain_decom.o schecomm.o printmatrix.o qsort.o mmio.o pspltm1.o lenstr.o \
         sparse_matrices.o interpolation.o m_fft_gpfa.o\
         files.o   pseudopotential.o chemical.o xml.o \
         parallelsubs.o spatial.o
#
COM_OBJS_STM=$(OBJS_STM) $(SYSOBJ)
ALL_OBJS_STM=$(MOD_OBJS) $(COM_OBJS_STM)
#
stm: $(FDF)  $(ALL_OBJS_STM)
	$(FC) -o stm \
	       $(LDFLAGS) $(ALL_OBJS_STM) $(FDF) $(FFTW_LIBS) 
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f stm *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)
	rm -f _tmp_deps deps.list  protomake*
#
dep:
	rm -f _tmp_deps deps.list  protomake*
	sfmakedepend --depend=obj --file _tmp_defs  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
          *.f *.F *.f90 *.F90
	@echo "Ignore errors in the following command"
	@echo "They appear if the last grep does not match anything"
	@-for i in $(ALL_OBJS_STM) ; do grep "^$$i: " _tmp_defs  ; done > deps.list
	@echo "ok"
	@csplit -k -f "protomake" Makefile '/^# -- specially prepared/+1'
	@cat protomake00 deps.list >| Makefile
#
# -- specially prepared dependency list -- do not remove this line!
m_io.o: local_sys.o
alloc.o: debugmpi.o m_io.o parallel.o precision.o local_sys.o
memoryinfo.o: precision.o
listsc.o: alloc.o
atom_options.o: local_sys.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o local_sys.o
atm_types.o: precision.o radial.o
radial.o: alloc.o interpolation.o precision.o xml.o
spher_harm.o: alloc.o precision.o local_sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o local_sys.o xml.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o local_sys.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o local_sys.o
schecomm.o: alloc.o
printmatrix.o: alloc.o
sparse_matrices.o: alloc.o precision.o
pseudopotential.o: alloc.o atom_options.o interpolation.o precision.o local_sys.o
chemical.o: parallel.o precision.o local_sys.o
xml.o: precision.o
parallelsubs.o: domain_decom.o parallel.o spatial.o local_sys.o
spatial.o: precision.o
atompla.o: precision.o local_sys.o
bessph.o: precision.o local_sys.o
chkdim.o: local_sys.o
extrapolate.o: fftw3_mymod.o precision.o
io.o: local_sys.o
iodm.o: alloc.o files.o parallel.o parallelsubs.o precision.o local_sys.o
mainstm.o: basis_io.o listsc.o parallel.o precision.o m_wfsx.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
ranger.o: precision.o
stm.o: atmfuncs.o chemical.o precision.o
vacpot.o: precision.o
