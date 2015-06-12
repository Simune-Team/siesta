#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
ARCH_MAKE_DEFAULT=lib-arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
VPATH=.
#
default: module
#
module:  libGridXC.a
#
dep:
	sfmakedepend --depend=obj  --modext=o *.f *.f90 *.F *.F90
#
GRID_SRCS= alloc.F90 bessph.f cellsubs.f \
             chkgmx.f  interpolation.f90  \
             minvec.f m_io.f moreParallelSubs.F90 \
             m_fft_gpfa.F  m_walltime.f90 \
	     precision.F radfft.f sorting.f sys.F
GRID_OBJS:= $(GRID_SRCS:.f=.o)
GRID_OBJS:= $(GRID_OBJS:.F=.o)
GRID_OBJS:= $(GRID_OBJS:.f90=.o)
GRID_OBJS:= $(GRID_OBJS:.F90=.o)
#
LOCAL_SRCS= atomxc.F90 cellxc.F90 debugxc.f90 ggaxc.F ldaxc.F \
            gridxc.F90 vdwxc.F90 xcmod.f90  fft3d.F90 fftr.F90 \
            mesh1d.f90 mesh3d.F90 am05.f90 vv_vdwxc.F90 \
            gridxc_config.F90
LOCAL_OBJS:= $(LOCAL_SRCS:.f=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.F=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.f90=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.F90=.o)
#
ALL_OBJS=$(GRID_OBJS) $(LOCAL_OBJS)
#
#
INCFLAGS:= -I./MPI_instr $(INCFLAGS)
#
libmpi_f90.a: 
	@(cd MPI_instr ; $(MAKE) \
                        "ARCH_MAKE=../$(ARCH_MAKE)"\
                        "VPATH=$(VPATH)/MPI_instr" )
#
libGridXC.a: $(MPI_INTERFACE) $(ALL_OBJS)
	ar $(ARFLAGS_EXTRA) cru libGridXC.a $(ALL_OBJS) 
	@if [ ! -z "$(MPI_INTERFACE)" ] ; then \
            ar $(ARFLAGS_EXTRA) cru libGridXC.a MPI_instr/*.o ; fi
	-$(RANLIB) libGridXC.a
#
clean:
	rm -f *.o *.*d *.a
	@(cd MPI_instr ; $(MAKE) \
                        "ARCH_MAKE=../$(ARCH_MAKE)"\
                        "VPATH=$(VPATH)/MPI_instr" clean)
	@(cd Testers ; $(MAKE) \
                   "ARCH_MAKE=../$(ARCH_MAKE)"\
                        "VPATH=$(VPATH)/Testers" clean)
#
# DO NOT DELETE THIS LINE - used by make depend
am05.o: precision.o sys.o
atomxc.o: alloc.o debugxc.o ggaxc.o gridxc_config.o ldaxc.o mesh1d.o
atomxc.o: precision.o radfft.o sys.o vdwxc.o xcmod.o
bessph.o: precision.o sys.o
cellsubs.o: precision.o
cellxc.o: alloc.o cellsubs.o chkgmx.o debugxc.o fftr.o ggaxc.o gridxc_config.o
cellxc.o: ldaxc.o  mesh3d.o moreParallelSubs.o precision.o sys.o
cellxc.o: vdwxc.o xcmod.o
chkgmx.o: cellsubs.o minvec.o precision.o sys.o
debugxc.o: m_io.o moreParallelSubs.o
fft3d.o: alloc.o debugxc.o m_fft_gpfa.o  mesh3d.o precision.o sys.o
fftr.o: alloc.o fft3d.o mesh3d.o precision.o sys.o
ggaxc.o: am05.o ldaxc.o precision.o sys.o
ldaxc.o: precision.o sys.o
m_io.o: sys.o
mesh1d.o: interpolation.o precision.o
mesh3d.o: alloc.o debugxc.o gridxc_config.o precision.o sorting.o sys.o
minvec.o: cellsubs.o precision.o sorting.o sys.o
moreParallelSubs.o: alloc.o m_io.o precision.o sys.o gridxc_config.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
gridxc.o: atomxc.o cellxc.o debugxc.o ggaxc.o gridxc_config.o ldaxc.o
gridxc.o: m_fft_gpfa.o  mesh3d.o precision.o xcmod.o
vdwxc.o: alloc.o debugxc.o ggaxc.o interpolation.o ldaxc.o mesh1d.o precision.o
vdwxc.o: radfft.o sys.o vv_vdwxc.o
vv_vdwxc.o: alloc.o debugxc.o interpolation.o mesh1d.o precision.o radfft.o
vv_vdwxc.o: sys.o
xcmod.o: precision.o sys.o vdwxc.o
m_atomxc.o: atomxc.o
m_bessph.o: bessph.o
m_cellxc.o: cellxc.o
m_chkgmx.o: chkgmx.o
m_fft3d.o: fft3d.o
m_ggaxc.o: ggaxc.o
m_ldaxc.o: ldaxc.o
m_minvec.o: minvec.o
m_radfft.o: radfft.o
m_vdwxc.o: vdwxc.o
m_vv_vdwxc.o: vv_vdwxc.o
