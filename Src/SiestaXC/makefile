.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
module:  libSiestaXC.a
include ../arch.make
#
OBJFILES= alloc.o atomxc.o bessph.o cellsubs.o cellxc.o chkgmx.o \
          debugxc.o fft1d.o fft3d.o fftr.o flib_spline.o \
          ggaxc.o ldaxc.o mesh1d.o mesh3d.o minvec.o moreParallelSubs.o \
          m_filter.o m_io.o m_recipes.o m_timer.o m_walltime.o \
          parallel.o precision.o radfft.o siestaxc.o sorting.o sys.o \
          vdwxc.o xcmod.o
#
libSiestaXC.a: $(OBJFILES)
	ar $(ARFLAGS_EXTRA) cru libSiestaXC.a $(OBJFILES)
	-$(RANLIB) libSiestaXC.a
	cp libSiestaXC.a ..
	@cp -p *.*d ..
	@echo "siestaXC/work.pc" >> ../work.pcl
#
#         (Note trick for Intel Fortran Compiler)
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
m_atomxc.mod: atomxc.o
m_bessph.mod: bessph.o
m_cellxc.mod: cellxc.o
m_chkgmx.mod: chkgmx.o
m_fft3d.mod: fft3d.o
m_ggaxc.mod: ggaxc.o
m_ldaxc.mod: ldaxc.o
m_minvec.mod: minvec.o
m_radfft.mod: radfft.o
m_vdwxc.mod: vdwxc.o
