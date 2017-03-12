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
#
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
TOP_LEVEL:=$(shell pwd)/../../../Src
VPATH:=.:$(TOP_LEVEL)
#
default: what denchar
#
#
ARCH_MAKE=../../../$(OBJDIR)/arch.make
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
	@sleep 2
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS_DENCHAR = atompla.o planed.o\
        bessph.o chkdim.o dismin.o \
	dot.o readpla.o readsts.o \
        io.o iodm.o itochar.o length.o \
        matvect.o memory.o colinear.o\
	neighb.o paste.o pxf.o \
	radfft.o m_fft_gpfa.o ranger.o \
        readwavesx.o rhoofr.o wavofr.o \
        redata_denchar.o  local_reinit.o denchar.o \
	volcel.o wrout.o stsofr.o local_sys.o
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
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=precision.o parallel.o debugmpi.o m_io.o alloc.o memoryinfo.o listsc.o  \
         atmparams.o atom_options.o atmfuncs.o   atm_types.o \
         radial.o spher_harm.o  basis_io.o basis_types.o parallelsubs.o \
         schecomm.o sparse_matrices.o interpolation.o\
         files.o pseudopotential.o chemical.o xml.o \
         domain_decom.o printmatrix.o qsort.o mmio.o pspltm1.o lenstr.o \
         spatial.o timer_local.o 
#
COM_OBJS_DENCHAR=$(OBJS_DENCHAR) $(SYSOBJ)
ALL_OBJS_DENCHAR=$(MOD_OBJS) $(COM_OBJS_DENCHAR)
#

denchar: $(FDF)  $(ALL_OBJS_DENCHAR)
	$(FC) -o denchar \
	       $(LDFLAGS) $(ALL_OBJS_DENCHAR) $(FDF) $(METIS_LIB)
#
#
dep:
	rm -f _tmp_deps deps.list _deps.list* protomake*
	sfmakedepend --depend=obj --file _tmp_defs  --modext=o \
          $(TOP_LEVEL)/*.f $(TOP_LEVEL)/*.f90 $(TOP_LEVEL)/*.F $(TOP_LEVEL)/*.F90 \
          *.f *.F *.f90 *.F90
	@echo "Ignore errors in the following command"
	@echo "They appear if the last grep does not match anything"
	@-for i in $(ALL_OBJS_DENCHAR) ; do grep "^$$i: " _tmp_defs  ; done > _deps.list.1
	@echo "ok"
	@sed '/^siesta_cmlsubs.o:/s/version.o//'  _deps.list.1 > deps.list
	@csplit -k -f "protomake" Makefile '/^# -- specially prepared/+1'
	@cat protomake00 deps.list >| Makefile
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f denchar *.o  *.a
	rm -f *.mod
	rm -f _tmp_deps deps.list _deps.list* protomake*
	rm -rf fdf

# -- specially prepared dependency list -- do not remove this line!
debugmpi.o: parallel.o
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
parallelsubs.o: domain_decom.o parallel.o spatial.o local_sys.o
schecomm.o: alloc.o
sparse_matrices.o: alloc.o precision.o
pseudopotential.o: alloc.o atom_options.o interpolation.o precision.o local_sys.o
chemical.o: parallel.o precision.o local_sys.o
xml.o: precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o local_sys.o
printmatrix.o: alloc.o
spatial.o: precision.o
die.o: parallel.o
atompla.o: precision.o local_sys.o
planed.o: precision.o
bessph.o: precision.o local_sys.o
chkdim.o: local_sys.o
readpla.o: precision.o
readsts.o: precision.o
io.o: local_sys.o
iodm.o: alloc.o files.o parallel.o parallelsubs.o precision.o local_sys.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
colinear.o: precision.o
neighb.o: precision.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
ranger.o: precision.o
readwavesx.o: precision.o
rhoofr.o: atmfuncs.o chemical.o listsc.o planed.o precision.o
wavofr.o: atmfuncs.o chemical.o listsc.o planed.o precision.o
redata_denchar.o: precision.o
local_reinit.o: files.o
denchar.o: basis_io.o listsc.o parallel.o precision.o
wrout.o: precision.o
stsofr.o: atmfuncs.o chemical.o listsc.o precision.o
