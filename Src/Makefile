# Makefile for siesta
#
.SUFFIXES: .f .F .o .a  .f90
MOD_DIR=Modules
#
#mpi:
#	$(MAKE) "DEFS=-DMPI" siesta
#
default: what siesta
#
# System dependencies are dealt with by an include file. See Sys directory
# for examples and instructions. The variable SIESTA_SYS should be in the 
# shell environment. 
#
SIESTA_SYS?=unknown
#
include Sys/${SIESTA_SYS}.make
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
what:
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = arw.o atomlwf.o bands.o bessph.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o conjgr.o constr.o coxmol.o cross.o \
	denmat.o detover.o dfscf.o dhscf.o diagon.o digcel.o fft3d.o \
        diagg.o diagk.o diagkp.o diag2g.o diag2k.o diagpol.o \
	dipole.o dismin.o dnaefs.o dot.o dynamics.o \
	efield.o egandd.o ener3.o extrapol.o extrapolon.o \
	fermid.o fermispin.o fixed.o gauleg.o gradient.o grdsam.o \
	hsparse.o idiag.o  initdm.o inver.o  \
	iodm.o iohs.o iolwf.o iorder.o iorho.o ioxv.o  \
	kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o lofilm.o \
	madelung.o meshmatrix.o memory.o meshsubs.o minvec.o mulliken.o \
        naefs.o neighb.o nlefsm.o \
	on_subs.o order.o ordern.o ordvec.o outcell.o outcoor.o overfsm.o \
	parse.o paste.o phirphi.o pixmol.o timestamp.o prmem.o propor.o pulayx.o  \
	ranger.o ran3.o recipes.o reclat.o recoor.o redata.o redcel.o\
	reinit.o reord.o resize.o rhoofd.o rhooda.o rlylm.o \
	savepsi.o shaper.o superc.o superx.o timer.o \
	vmb.o vmat.o volcel.o xc.o xijorb.o ylmexp.o ylmylm.o \
	cellxc.o cdiag.o rdiag.o \
        cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o repol.o typecell.o \
        ofc.o poison.o radfft.o die.o symbol.o\
	siesta.o atom.o redbasis.o plcharge.o dump_atom.o version.o io.o
#
FDF= fdf/libfdf.a
$(FDF): fdf/fdf.f fdf/fdf_mod.f
	(cd fdf ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" module)
#
# Routines using fdf calls.
#
dhscf.o initdm.o iodm.o iohs.o iolwf.o iorho.o grdsam.o : $(FDF)
plcharge.o recoor.o outcoor.o ioxv.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o: $(FDF)
cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o repol.o ofc.o : $(FDF)
#
atom.o redbasis.o: $(FDF)
#
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=precision.o atmfuncs.o atminit.o listsc.o matel.o memoryinfo.o numbvect.o \
        parallel.o
COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
atminit.o atom.o atomlwf.o dhscf.o dnaefs.o       : atmfuncs.o
dump_atom.o efield.o matel.o mulliken.o outcoor.o : atmfuncs.o
overfsm.o plcharge.o redbasis.o shaper.o          : atmfuncs.o
kinefsm.o naefs.o nlefsm.o overfsm.o              : atmfuncs.o matel.o
siesta.o          : atminit.o atmfuncs.o matel.o
dfscf.o hsparse.o : atmfuncs.o listsc.o
phirphi.o         : matel.o
rhoofd.o vmat.o   : listsc.o
chempot.o         : numbvect.o
memory.o          : memoryinfo.o
dhscf.o dfscf.o rhooda.o rhoofd.o vmat.o : meshsubs.o
dfscf.o rhoofd.o vmat.o : meshmatrix.o
#
#
$(COM_OBJS): precision.o parallel.o $(MPILIB)
#
# Compile "AS IS" means basically not to use "-dp" on the Cray
#
precision.o: precision.f
	$(FC) $(FFLAGS) -c precision.f
#
libmpi_f90.a: precision.o
	@(cd MPI ; $(MAKE) )
#
# Libraries that need to be compiled
#
linalg.a: Libs/blas.f Libs/lapack.f
	@(cd Libs ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)")
#
version:
	sed  "s/SYS/${SIESTA_SYS}/g" version.f > temp0.f
	sed  "s/FFLAGS/${FC} ${FFLAGS}/g" temp0.f > temp1.f
	sed  "s/MPILIB/${MPILIB}/g" temp1.f > temp.f
	$(FC) -c $(FFLAGS) temp.f
	@mv temp.o version.o
	@rm -f temp.f temp0.f temp1.f
	@echo

siesta: version $(MPILIB) $(FDF) $(COMP_LIBS) $(ALL_OBJS)
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(MPILIB) $(COMP_LIBS) $(LIBS) 
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)
	(cd MPI ; $(MAKE) clean )
	(cd Libs ; $(MAKE) clean )
pristine:  clean
	cp -fp Include/constr.f .
#

