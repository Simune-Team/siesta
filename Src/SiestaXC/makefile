OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
include ../../$(OBJDIR)/arch.make
#
VPATH:=$(VPATH_ROOT)/SiestaXC:$(VPATH_ROOT)
#
default: libSiestaXC.a
#
SIESTA_OBJS= alloc.o bessph.o cellsubs.o  \
             chkgmx.o fft1d.o  flib_spline.o  \
             m_filter.o minvec.o m_io.o moreParallelSubs.o \
             m_recipes.o m_timer.o m_walltime.o parallel.o \
	     precision.o radfft.o sorting.o sys.o 

LOCAL_OBJS= atomxc.o cellxc.o debugxc.o ggaxc.o ldaxc.o \
            siestaxc.o vdwxc.o xcmod.o  fftr.o \
            mesh1d.o mesh3d.o 

module:  libSiestaXC.a
	cp *.mod *.a ..
#
ALL_FILES=$(SIESTA_OBJS) $(LOCAL_OBJS)
#
libSiestaXC.a: $(ALL_FILES)
	ar $(ARFLAGS_EXTRA) cru libSiestaXC.a $(ALL_FILES)
	-$(RANLIB) libSiestaXC.a
#
clean:
	rm -f *.o *.*d *.a
# DO NOT DELETE THIS LINE - used by make depend
alloc.o: m_io.o parallel.o precision.o sys.o
atomxc.o: alloc.o ggaxc.o ldaxc.o m_filter.o mesh1d.o precision.o radfft.o
atomxc.o: sys.o vdwxc.o xcmod.o
bessph.o: precision.o sys.o
cellsubs.o: precision.o
cellxc.o: alloc.o cellsubs.o chkgmx.o debugxc.o fftr.o ggaxc.o ldaxc.o
cellxc.o: m_timer.o mesh3d.o moreParallelSubs.o parallel.o precision.o sys.o
cellxc.o: vdwxc.o xcmod.o
chkgmx.o: cellsubs.o minvec.o precision.o sys.o
debugxc.o: m_io.o moreParallelSubs.o parallel.o
fft1d.o: parallel.o precision.o sys.o
fft3d.o: alloc.o debugxc.o fft1d.o m_timer.o mesh3d.o precision.o sys.o
fftr.o: alloc.o fft3d.o mesh3d.o precision.o sys.o
flib_spline.o: precision.o
ggaxc.o: ldaxc.o precision.o sys.o
ldaxc.o: precision.o sys.o
m_io.o: sys.o
m_recipes.o: sys.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o sys.o
mesh1d.o: m_recipes.o precision.o
mesh3d.o: alloc.o debugxc.o parallel.o precision.o sorting.o sys.o
minvec.o: cellsubs.o precision.o sorting.o sys.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o sys.o
radfft.o: alloc.o bessph.o m_recipes.o precision.o
siestaxc.o: atomxc.o cellxc.o xcmod.o
sys.o: parallel.o
vdwxc.o: debugxc.o flib_spline.o ldaxc.o m_recipes.o mesh1d.o precision.o
vdwxc.o: radfft.o
xcmod.o: precision.o sys.o

