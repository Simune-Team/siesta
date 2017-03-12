# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for tbtrans
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
# Please see note about Linear-Algebra libraries below
#

default: what tbtrans 

VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
# For removing any -ipo compilation in sub directories...
# Later compilers are also having this enabled, and it lets you
# specify what to not be used in libraries.
# This is so because inter-procedural compilation should NOT be performed in
# libraries (it removes names when a routine can be fully moved into others)
IPO_FLAG?= -ipo

# If compiling without system lapack, you might need something like this
# 
#TBTRANS_LINALG=linalg.a # Or maybe only: zgeev_lapack.a zgesv_lapack.a
LINALG:= $(TBTRANS_LINALG) $(LINALG)

SIESTA_SYS=$(SIESTA_ARCH)

#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG) -g

what:
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo

SYSOBJ=$(SYS).o

# Tbtrans:
MOD_TOBJS=  m_tbt_options.o m_tbt_kpoints.o m_tbt_out.o m_tbt_read_tshs.o m_lowdin.o

SIESTA_OBJS= alloc.o cdiag.o cellsubs.o debugmpi.o domain_decom.o \
	files.o find_kgrid.o idiag.o io.o kpoint_convert.o m_fdf_global.o		\
	m_hs_matrix.o m_io.o m_mpi_utils.o m_spin.o m_timer.o  \
	m_ts_aux_rout.o m_ts_cctype.o m_ts_contour.o m_ts_electrode.o \
        m_ts_gf.o m_ts_io.o m_ts_scattering.o m_wallclock.o m_walltime.o memory.o\
        memoryinfo.o minvec.o moreParallelSubs.o parallel.o parallelsubs.o \
        paste.o precision.o printmatrix.o pxf.o reclat.o schecomm.o sorting.o sparse_matrices.o\
        spatial.o timer.o timer_tree.o timestamp.o units.o zheevds.o\
        siesta_options.o mmio.o pspltm1.o qsort.o lenstr.o

TBTOBJS= $(SIESTA_OBJS) $(MOD_TOBJS) $(SYSOBJ) $(EXTRA) \
         tbtrans.o transmission.o tbt_init.o tbt_reinit.o \
	 coop.o tt_eig.o atompdos.o local_sys.o


# Interfaces to libraries
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf -I../MPI $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module)

INCFLAGS:=$(INCFLAGS)

MPI_MAKEFILE=$(VPATH)/MPI/Makefile
$(MPI_INTERFACE): 
	(mkdir -p MPI; cd MPI ; $(MAKE) -f $(MPI_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/MPI" \
                          "MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module_built)

LINALG_MAKEFILE=$(VPATH)/Libs/makefile
LINALG_INCFLAGS:=-I$(VPATH)/Libs $(INCFLAGS)
$(LINALG): 
	(mkdir -p Libs; cd Libs ; $(MAKE) -f $(LINALG_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/Libs" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(LINALG_INCFLAGS)" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" $(LINALG))

tbtrans: FPPFLAGS+=$(DEFS_PREFIX) -DTBTRANS
tbtrans: DEFS+=$(DEFS_PREFIX) -DTBTRANS
tbtrans: $(MPI_INTERFACE) $(FDF) $(TBTOBJS) $(SYSOBJ) $(LINALG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o tbtrans \
              $(TBTOBJS) $(MPI_INTERFACE)  \
              $(FDF) $(LIBS) $(LINALG)  || echo "***Please see top of Makefile"

clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f tbtrans *.o *.a *.mod
	rm -rf ./fdf
	rm -rf ./MPI
	rm -rf ./Libs
	rm -f _tmp_deps deps.list _deps.list* protomake*
#
# Dependencies
dep:	
	rm -f _tmp_deps deps.list _deps.list* protomake*
	sfmakedepend --depend=obj --srcdir=$(VPATH)  --modext=o\
                   --file _tmp_defs *.f90 *.F90 *.f *.F \
                 $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90
	@echo "Ignore errors in the following command"
	@echo "They appear if the last grep does not match anything"
	@-for i in $(TBTOBJS) ; do grep "^$$i: " _tmp_defs ; done >| _deps.list.1
	@echo "ok"
	@sed "/^m_ts_io.o:/s/m_ts_kpoints.o//" _deps.list.1 >| _deps.list.2
	@sed "/^m_ts_electrode.o:/s/m_ts_kpoints.o//" _deps.list.2 >| deps.list
	@csplit -k -f "protomake" Makefile '/^# -- specially prepared/+1'
	@cat protomake00 deps.list >| Makefile
#
$(TBTOBJS): $(FDF)
#
# -- specially prepared dependency list -- do not remove this line!
alloc.o: m_io.o parallel.o precision.o local_sys.o
cdiag.o: alloc.o m_spin.o memoryinfo.o parallel.o precision.o local_sys.o
cellsubs.o: precision.o
debugmpi.o: parallel.o
die.o: parallel.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o local_sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
idiag.o: parallel.o local_sys.o
io.o: local_sys.o
kpoint_convert.o: precision.o local_sys.o units.o
m_fdf_global.o: m_mpi_utils.o parallel.o precision.o
m_hs_matrix.o: alloc.o precision.o local_sys.o
m_io.o: local_sys.o
m_mpi_utils.o: precision.o local_sys.o
m_spin.o: alloc.o parallel.o precision.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o local_sys.o
m_ts_aux_rout.o: precision.o local_sys.o
m_ts_cctype.o: precision.o
m_ts_contour.o: m_ts_aux_rout.o m_ts_cctype.o parallel.o precision.o local_sys.o
m_ts_contour.o: units.o
m_ts_electrode.o: files.o m_hs_matrix.o m_tbt_kpoints.o m_ts_aux_rout.o
m_ts_electrode.o: m_ts_cctype.o m_ts_io.o  parallel.o precision.o
m_ts_electrode.o: local_sys.o units.o
m_ts_gf.o: m_ts_cctype.o m_ts_electrode.o parallel.o precision.o local_sys.o units.o
m_ts_io.o: files.o m_spin.o m_tbt_kpoints.o  parallel.o
m_ts_io.o: parallelsubs.o precision.o local_sys.o
m_ts_scattering.o: m_ts_aux_rout.o parallel.o precision.o local_sys.o units.o
m_wallclock.o: m_walltime.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
memoryinfo.o: precision.o
minvec.o: cellsubs.o precision.o sorting.o local_sys.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o local_sys.o
parallelsubs.o: domain_decom.o parallel.o spatial.o local_sys.o
printmatrix.o: alloc.o
schecomm.o: alloc.o
sparse_matrices.o: alloc.o precision.o
spatial.o: precision.o
local_sys.o: parallel.o
timer.o: m_timer.o parallel.o local_sys.o timer_tree.o
timer_tree.o: m_walltime.o
units.o: precision.o
m_tbt_options.o: files.o m_fdf_global.o m_ts_io.o memoryinfo.o parallel.o
m_tbt_options.o: precision.o local_sys.o units.o
m_tbt_kpoints.o: alloc.o minvec.o parallel.o precision.o local_sys.o units.o
m_tbt_out.o: parallel.o precision.o units.o
m_tbt_read_tshs.o: m_ts_io.o parallel.o precision.o
m_lowdin.o: precision.o
tbtrans.o: alloc.o files.o m_hs_matrix.o m_tbt_kpoints.o m_tbt_options.o
tbtrans.o: m_tbt_out.o m_tbt_read_tshs.o m_ts_contour.o m_ts_electrode.o
tbtrans.o: m_ts_gf.o m_ts_io.o m_ts_scattering.o m_wallclock.o parallel.o
tbtrans.o: precision.o local_sys.o timestamp.o units.o
transmission.o: m_ts_aux_rout.o precision.o local_sys.o
tbt_init.o: alloc.o files.o m_spin.o m_tbt_kpoints.o m_tbt_options.o m_timer.o
tbt_init.o: m_ts_contour.o m_wallclock.o parallel.o precision.o local_sys.o
tbt_init.o: timestamp.o
tbt_reinit.o: files.o parallel.o
coop.o: m_tbt_out.o parallel.o precision.o local_sys.o units.o
tt_eig.o: precision.o local_sys.o
atompdos.o: m_lowdin.o m_tbt_out.o parallel.o precision.o local_sys.o  units.o
