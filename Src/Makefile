# Makefile for siesta
#
.SUFFIXES: .f .F .o .a  .f90 .F90
MOD_DIR=Modules
#
default: what siesta
#
include arch.make
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
	@echo "Hit return to continue, or ^C to abort..."
	@read ans
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = arw.o  atomlwf.o bands.o bessph.o cgwf.o chkdim.o chkgmx.o \
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
	paste.o phirphi.o pixmol.o timestamp.o \
        propor.o pulayx.o  \
	ranger.o ran3.o recipes.o reclat.o redata.o redcel.o\
	reinit.o reord.o rhoofd.o rhooda.o rlylm.o \
	savepsi.o shaper.o superx.o timer.o \
	vmb.o vmat.o volcel.o xc.o xijorb.o ylmexp.o ylmylm.o \
	cellxc.o cdiag.o rdiag.o \
        cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o repol.o typecell.o \
        ofc.o poison.o radfft.o die.o symbol.o\
	siesta.o atom.o redbasis.o plcharge.o version.o io.o\
        spin_init.o coor.o dump2cdf.o transfer.o\
        read_user_basis.o broadcast_basis.o dump_ascii.o
#
FDF=libfdf.a
$(FDF): 
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
coor.o read_user_basis.o: $(FDF)
#
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=precision.o atmparams.o atmfuncs.o atminit.o listsc.o \
         matel.o memoryinfo.o numbvect.o  parallel.o \
	 atomlist.o ionew.o types.o old_atmfuncs.o radial.o parsing.o\
         interfaces.o  alloc.o
COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
atom.o redbasis.o plcharge.o transfer.o     : old_atmfuncs.o
atminit.o atom.o redbasis.o old_atmfuncs.o atmfuncs.o: atmparams.o
atomlwf.o dhscf.o dnaefs.o                 : atmfuncs.o
efield.o matel.o mulliken.o outcoor.o : atmfuncs.o
overfsm.o shaper.o                                : atmfuncs.o
atomlist.o                                        : atmfuncs.o 
kinefsm.o naefs.o nlefsm.o overfsm.o              : atmfuncs.o matel.o
siesta.o          : atminit.o atmfuncs.o matel.o
dfscf.o hsparse.o : atmfuncs.o listsc.o
phirphi.o         : matel.o
rhoofd.o vmat.o   : listsc.o
chempot.o         : numbvect.o
memory.o          : memoryinfo.o
dhscf.o dfscf.o rhooda.o rhoofd.o vmat.o : meshsubs.o
dfscf.o rhoofd.o vmat.o : meshmatrix.o
siesta.o: atomlist.o
old_atmfuncs.o atomlist.o siesta.o atom.o: ionew.o
coor.o: atomlist.o
types.o: atmparams.o
transfer.o: types.o atmfuncs.o
atmfuncs.o: types.o atmparams.o
#
#
dump2cdf.o: types.o atmfuncs.o
atmfuncs.o: radial.o
redbasis.o siesta.o bands.o efield.o fixed.o grdsam.o: parsing.o
initdm.o repol.o ksv.o: parsing.o
#
siesta.o: interfaces.o
siesta.o atomlist.o coor.o: alloc.o
#
read_user_basis.o: types.o atmfuncs.o radial.o parsing.o
broadcast_basis.o dump_ascii.o : types.o atmfuncs.o

read_user_basis.o dump2cdf.o: $(NETCDF_INTERFACE)
#
$(COM_OBJS): precision.o parallel.o $(MPI_INTERFACE)
#
# Interfaces to libraries
#
libmpi_f90.a: 
	@(cd MPI ; $(MAKE) )
libnetcdf_f90.a: 
	@(cd NetCDF ; $(MAKE) )
#
# Libraries that need to be compiled
#
linalg.a: Libs/blas.f Libs/lapack.f
	@(cd Libs ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)")
#
version:
	sed  "s/SIESTA_ARCH/${SIESTA_ARCH}/g" version.F > temp0.F
	sed  "s/FFLAGS/${FC} ${FFLAGS}/g" temp0.F > temp.F
	$(FC) -c $(FFLAGS) $(DEFS) temp.F
	@mv temp.o version.o
	@rm -f temp.F temp0.F 
	@echo

siesta: version $(MPI_INTERFACE) $(NETCDF_INTERFACE) $(FDF) \
                $(COMP_LIBS) $(ALL_OBJS)
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(MPI_INTERFACE)\
               $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)
	(cd MPI ; $(MAKE) clean )
	(cd NetCDF ; $(MAKE) clean )
	(cd Libs ; $(MAKE) clean )
pristine:  clean
	cp -fp Include/constr.f .
#

