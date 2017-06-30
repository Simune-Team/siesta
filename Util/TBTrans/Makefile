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
default: what  tbtrans 
#
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

# If compiling without system lapack, you might need something like this
# 
#TBTRANS_LINALG=linalg.a # Or maybe only: zgeev_lapack.a zgesv_lapack.a
LINALG:= $(TBTRANS_LINALG) $(LINALG)

SIESTA_SYS=${SIESTA_ARCH}

#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG) -g
#
what:
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo

SYSOBJ=$(SYS).o
# Tbtrans:
MOD_TOBJS=  m_tbt_options.o m_tbt_gf.o m_tbt_kpts.o hamil2.o tsread2.o

SIESTA_OBJS= precision.o m_io.o alloc.o  memoryinfo.o memory.o pxf.o \
             moreParallelSubs.o m_walltime.o m_timer.o timer.o m_mpi_utils.o \
             m_fdf_global.o paste.o timer_tree.o \
             parallel.o  io.o files.o  zheevds.o debugmpi.o 

TBTOBJS= tbtrans.o transmission.o   \
         calc_green.o cdiag_local.o getsfe.o \
         mkRealContour.o reinit_tb.o setupkham2.o TSiohs.o \
         csolve.o teigchana.o helpers.o local_sys.o \
         $(SIESTA_OBJS) $(MOD_TOBJS) $(SYSOBJ) $(EXTRA)


#
# Interfaces to libraries
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf -I../MPI $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
INCFLAGS:=$(INCFLAGS)
#
MPI_MAKEFILE=$(VPATH)/MPI/Makefile
$(MPI_INTERFACE): 
	(mkdir -p MPI; cd MPI ; $(MAKE) -f $(MPI_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/MPI" \
                          "MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "FFLAGS=$(FFLAGS)" module_built)
#
LINALG_MAKEFILE=$(VPATH)/Libs/makefile
LINALG_INCFLAGS:=-I $(VPATH)/Libs $(INCFLAGS)
$(LINALG): 
	(mkdir -p Libs; cd Libs ; $(MAKE) -f $(LINALG_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/Libs" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(LINALG_INCFLAGS)" "FFLAGS=$(FFLAGS)" $(LINALG))
#
tbtrans: $(MPI_INTERFACE) $(TBTOBJS) $(FDF) $(SYSOBJ)  $(LINALG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o tbtrans \
              $(TBTOBJS) $(MPI_INTERFACE)  \
              $(FDF) $(LIBS) $(LINALG)  || echo "***Please see top of Makefile"
#	@(cd ../../bin ; rm -f tsScan ; ln -s tbtrans tsScan)

#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f tbtrans *.o *.a *.mod
	rm -rf ./fdf
	rm -rf ./MPI
	rm -rf ./Libs
	rm -f _tmp_deps deps.list  protomake*
#
# Dependencies
dep:
	rm -f _tmp_deps deps.list  protomake*
	sfmakedepend --silent --depend=obj -I $(VPATH)  --modext=o\
                   --file _tmp_defs *.f90 *.F90 *.f *.F \
         $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 #$(VPATH)/*.T90 $(VPATH)/*.inc
	@echo "Ignore errors in the following command"
	@echo "They appear if the last grep does not match anything"
	@-for i in $(TBTOBJS) ; do grep "^$$i: " _tmp_defs ; done > deps.list
	@echo "ok"
	@csplit -k -f "protomake" Makefile '/^# -- specially prepared/+1'
	@cat protomake00 deps.list >| Makefile
#
$(TBTOBJS)      : $(FDF)

reinit.o        : $(MPI_INTERFACE)
local_sys.o           : $(MPI_INTERFACE) parallel.o
#
# -- specially prepared dependency list -- do not remove this line!
tbtrans.o: files.o hamil2.o m_mpi_utils.o m_tbt_gf.o m_tbt_kpts.o
tbtrans.o: m_tbt_options.o parallel.o local_sys.o
calc_green.o: m_tbt_gf.o
getsfe.o: m_tbt_gf.o
mkRealContour.o: parallel.o
reinit_tb.o: files.o parallel.o
setupkham2.o: m_tbt_kpts.o tsread2.o
TSiohs.o: m_mpi_utils.o parallel.o local_sys.o
local_sys.o: parallel.o
m_io.o: local_sys.o
alloc.o: m_io.o parallel.o precision.o local_sys.o
memoryinfo.o: precision.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o local_sys.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o local_sys.o
timer.o: m_timer.o parallel.o local_sys.o timer_tree.o
m_mpi_utils.o: precision.o local_sys.o
m_fdf_global.o: m_mpi_utils.o parallel.o precision.o
timer_tree.o: m_walltime.o
io.o: local_sys.o
debugmpi.o: parallel.o
die.o: parallel.o
m_tbt_options.o: files.o parallel.o precision.o
m_tbt_gf.o: m_tbt_kpts.o m_tbt_options.o parallel.o local_sys.o tsread2.o
m_tbt_kpts.o: m_mpi_utils.o m_tbt_options.o parallel.o local_sys.o
