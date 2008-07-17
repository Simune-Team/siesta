# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#
# Makefile for siesta
#
# Remove implicit rules (including for Modula-2, which just confuses matters)
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
VPATH=.
#
default: what siesta # gen-basis denchar fsiesta.o
#
configure: configure.ac aclocal.m4
	autoconf ${<D}/configure.ac
#
M4_FILES = m4/ac_config_aux_dirs.m4  m4/TW_CHECK_FC_FPP.m4      m4/TW_FC_ID_FLAGS.m4  \
 m4/ACX_MPI.m4             m4/TW_CHECK_FC_TR15580.m4  m4/TW_FC_ID.m4        \
 m4/fortran.m4             m4/TW_CHECK_FC_TR15581.m4  m4/TW_FIND_FC_BLAS.m4 \
 m4/TW_CHECK_BLACS.m4      m4/TW_CHECK_SCALAPACK.m4   m4/TW_FIND_LAPACK.m4  \
 m4/TW_CHECK_FC_90.m4      m4/TW_FC_CHECK_ABORT.m4    m4/TW_PATH_NETCDF.m4  \
 m4/TW_CHECK_FC_95.m4      m4/TW_FC_CHECK_DCFUNS.m4   m4/TW_TRY_DC_LAPACK.m4\
 m4/TW_CHECK_FC_FPP_90.m4  m4/TW_FC_CHECK_FLUSH.m4    m4/TW_FC_KINDS.m4
#
aclocal.m4: $(M4_FILES)
	cat ${<D}/$(M4_FILES)  > aclocal.m4
dep:
	@if [ $(VPATH) == "." ] ; then \
	sfmakedepend --depend=obj  --modext=o \
                              *.f *.f90 *.F *.F90 ; \
        else \
	sfmakedepend --depend=obj  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 ; \
        fi
	@sed '/^siesta_cmlsubs.o:/s/version.o//'  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile
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
	@echo "Hit ^C to abort..."
	@sleep 2
#
SYSOBJ=$(SYS).o pxf.o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS =  domain_decom.o automatic_cell.o \
        arw.o  atomlwf.o bands.o basis_enthalpy.o bessph.o bonds.o \
        born_charge.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o coxmol.o cross.o \
	denmat.o denmatlomem.o detover.o dfscf.o m_diagon.o diagon.o digcel.o \
        fft3d.o dhscf.o m_dhscf.o constr.o \
	diagg.o diagk.o diagkp.o diag2g.o diag2k.o diagpol.o \
	diagsprl.o dipole.o dismin.o dnaefs.o dot.o  \
	efield.o egandd.o ener3.o ener3lomem.o extrapol.o  extrapolon.o \
	fixed.o flib_spline.o gradient.o gradientlomem.o grdsam.o \
	hsparse.o idiag.o  initatom.o initdm.o inver.o \
	iodm.o iohs.o iolwf.o iozm.o ipack.o iopipes.o \
	kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o local_DOS.o\
	madelung.o matel.o meshmatrix.o memory.o meshsubs.o metaforce.o \
        m_check_supercell.o mulliken.o minvec.o m_pulay.o naefs.o \
	new_dm.o normalize_dm.o ordern.o outcell.o outcoor.o \
	paste.o pdos.o pdosg.o pdosk.o pdoskp.o phirphi.o pixmol.o plcharge.o \
	projected_DOS.o propor.o \
	ran3.o reclat.o redcel.o \
	reinit.o reord.o rhoofd.o rhoofdsp.o rhooda.o \
	savepsi.o shaper.o timer.o \
	vmb.o vmat.o vmatsp.o volcel.o xc.o \
	cellxc.o cdiag.o rdiag.o \
        cgvc.o cgvc_zmatrix.o  \
        iocg.o ioeig.o iofa.o iokp.o iomd.o typecell.o \
        ofc.o poison.o readsp.o radfft.o \
	write_md_record.o kpoint_grid.o find_kgrid.o proximity_check.o\
	state_init.o siesta_move.o setup_hamiltonian.o compute_dm.o mixer.o\
	scfconvergence_test.o post_scf_work.o state_analysis.o write_subs.o \
	siesta_init.o struct_init.o siesta_options.o siesta_geom.o \
        siesta_analysis.o siesta_end.o siesta_forces.o \
        siesta.o io.o sparse_matrices.o \
        coor.o atm_transfer.o \
        broadcast_basis.o eggbox.o dsyevds.o zheevds.o \
        optical.o phirphi_opt.o reoptical.o transition_rate.o \
        initparallel.o show_distribution.o setspatial.o setatomnodes.o uncell.o \
        cart2frac.o obc.o \
        precision.o sys.o m_cell.o recipes.o files.o \
        spatial.o parallel.o parallelsubs.o  parsing.o chemical.o xcmod.o atom.o atmparams.o \
        m_mpi_utils.o m_fdf_global.o m_history.o m_iorho.o  \
        atmfuncs.o listsc.o memoryinfo.o m_memory.o sorting.o \
        atomlist.o atm_types.o old_atmfuncs.o radial.o \
        alloc.o phonon.o spher_harm.o periodic_table.o version.o timestamp.o \
        basis_types.o xml.o pseudopotential.o basis_specs.o basis_io.o\
        onmod.o densematrix.o writewave.o on_subs.o fermid.o m_broyddj.o \
        electrostatic.o mneighb.o globalise.o siesta_cmlsubs.o siesta_cml.o \
        units.o zmatrix.o m_broyden_mixing.o forhar.o m_walltime.o m_wallclock.o m_iostruct.o \
        nlefsm.o overfsm.o overlap.o conjgr.o conjgr_old.o m_energies.o m_steps.o\
        m_broyddj_nocomm.o broyden_optim.o ioxv.o dynamics.o md_out.o \
        molecularmechanics.o filter.o filter_subs.o zm_broyden_optim.o \
        cell_broyden_optim.o remove_intramol_pressure.o \
        m_ntm.o m_dipol.o  m_kinetic.o  m_rmaxh.o \
        m_forces.o m_stress.o m_eo.o m_spin.o m_gamma.o setup_H0.o \
        final_H_f_stress.o phionmesh.o initmesh.o \
        debugmpi.o qsort.o distriphionmesh.o \
        dscfcomm.o schecomm.o moremeshsubs.o \
        cellxc_mod.o pspltm1.o lenstr.o printmatrix.o jacobi_davidson.o \
        matrixop.o precond.o gmres.o
#--------------------------------------------------------------
#
WXML=libwxml.a
$(WXML): 
	(cd wxml ; $(MAKE) "VPATH=$(VPATH)/wxml" module)
#
# Routines using wxml calls.
#
pdos.o : $(WXML)
#--------------------------------------------------------------
XMLPARSER=libxmlparser.a
$(XMLPARSER): 
	(cd xmlparser ; $(MAKE) "VPATH=$(VPATH)/xmlparser" module)
#
# Routines using xmlparser calls.
#
######  somefile.o : $(WXML)
#--------------------------------------------------------------
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "VPATH=$(VPATH)/fdf" module)
#
# Routines using fdf calls.
#
dhscf.o initdm.o iodm.o iohs.o iolwf.o  grdsam.o : $(FDF)
recoor.o outcoor.o ioxv.o m_iostruct.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o pdos.o plcharge.o: $(FDF)
coor.o cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o  ofc.o : $(FDF)
readsp.o eggbox.o iopipes.o : $(FDF)
#
atom.o basis_io.o basis_specs.o: $(FDF)
gen-basis.o: $(FDF)
#
#
ALL_OBJS=$(OBJS) $(SYSOBJ)
#
basis_io.o md_out.o: $(NETCDF_INTERFACE)
#
# Interfaces to libraries
#
libmpi_f90.a: 
	@(cd MPI ; $(MAKE) "VPATH=$(VPATH)/MPI")
libnetcdf_f90.a: 
	@(cd NetCDF ; $(MAKE) "VPATH=$(VPATH)/NetCDF" )
#
# Libraries that might need to be compiled
#
libblas.a: Libs/blas.f
	@echo "==> Compiling libblas.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" libblas.a)
liblapack.a: Libs/lapack.f
	@echo "==> Compiling liblapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" liblapack.a)
linalg.a: Libs/blas.f Libs/lapack.f
	@echo "==> Compiling linalg.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" linalg.a)
dc_lapack.a:  Libs/dc_lapack.f
	@echo "==> Compiling dc_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" dc_lapack.a)
#

version: version.F90
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;s'SIESTA_VERSION'\
                                       $$(cat ${<D}/../version.info)'g" $<  > compinfo.F90
	@($(MAKE) compinfo.o)
	@mv compinfo.o version.o
	#@rm -f compinfo.F90
	@echo

siesta: version $(MPI_INTERFACE) $(NETCDF_INTERFACE) $(FDF) $(WXML) $(XMLPARSER) \
                $(COMP_LIBS) $(ALL_OBJS) 
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(WXML) $(XMLPARSER) $(MPI_INTERFACE)\
               $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#
#------------------------------------------------------------------------
GEN-BASIS_OBJS=alloc.o gen-basis.o basis_types.o precision.o recipes.o parallel.o \
          parsing.o xcmod.o basis_io.o  chemical.o atm_transfer.o atm_types.o\
          atmparams.o old_atmfuncs.o radial.o io.o \
	  paste.o  flib_spline.o bessph.o filter.o filter_subs.o \
          basis_specs.o atom.o memoryinfo.o memory.o periodic_table.o\
          pseudopotential.o dot.o xc.o  arw.o  \
          sys.o timer.o  xml.o m_walltime.o $(SYSOBJ)
# precision must be the first one here...
gen-basis.o: $(FDF)
gen-basis: version $(NETCDF_INTERFACE) $(MPI_INTERFACE) $(FDF) $(WXML) \
                 $(GEN-BASIS_OBJS) 
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(GEN-BASIS_OBJS) $(MPI_INTERFACE) $(FDF) $(WXML) \
               $(NETCDF_INTERFACE) $(LIBS) 
#
#------------------------------------------------------------------------
IONCAT_OBJS=f2kcli.o m_getopts.o ioncat.o basis_types.o precision.o  parallel.o \
          parsing.o xcmod.o basis_io.o  atm_types.o\
          atmparams.o radial.o io.o memory.o chkdim.o\
	  paste.o pseudopotential.o chemical.o  \
          basis_specs.o memoryinfo.o  m_memory.o periodic_table.o\
          xc.o flib_spline.o atmfuncs.o spher_harm.o\
          sys.o timer.o  xml.o m_walltime.o $(SYSOBJ)
# precision must be the first one here...
ioncat.o: $(FDF)
ioncat: version $(NETCDF_INTERFACE) $(MPI_INTERFACE) $(FDF) $(WXML) \
                 $(IONCAT_OBJS) $(COMP_LIBS)
	$(FC) -o ioncat \
	       $(LDFLAGS) $(IONCAT_OBJS) $(MPI_INTERFACE) $(FDF) $(WXML) \
               $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#
#------------------------------------------------------------------------
# Denchar is now at the main source level...
#
DENCHAR_OBJS = m_denchar_init.o m_denchar_geom.o m_denchar_io.o \
               m_denchar_neighb.o m_denchar_work.o \
               denchar.o flib_spline.o
#
m_denchar_io.o m_denchar_init.o denchar.o : $(FDF)
#
OTHER_DENCHAR_OBJS=  precision.o recipes.o \
        f2kcli.o bessph.o chkdim.o dismin.o  dot.o iodm.o memory.o \
        paste.o  radfft.o io.o  spatial.o flib_spline.o\
        volcel.o parallel.o parallelsubs.o memoryinfo.o sys.o listsc.o  \
        atmparams.o atmfuncs.o atm_types.o m_memory.o \
        radial.o spher_harm.o  basis_io.o basis_types.o \
        pseudopotential.o chemical.o xml.o files.o $(SYSOBJ)

denchar: $(NETCDF_INTERFACE) $(MPI_INTERFACE) $(FDF) \
         $(OTHER_DENCHAR_OBJS) $(DENCHAR_OBJS) $(COMP_LIBS)
	$(FC) -o denchar \
	   $(LDFLAGS) $(DENCHAR_OBJS) $(OTHER_DENCHAR_OBJS) \
           $(MPI_INTERFACE) $(FDF) $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#
# ---------------------------------------------------------------------------
# Other auxiliary programs
#----------------------------------------------------------------------------
MT_OBJS=precision.o m_mpi_utils.o sys.o parallel.o 
mpi_test: $(MPI_INTERFACE) $(MT_OBJS) mpi_test.o
	$(FC) -o mpi_test \
	       $(LDFLAGS) mpi_test.o $(MT_OBJS) $(MPI_INTERFACE) $(LIBS) 
#----------------------------------------------------------------------------
MIXPS_OBJS= f2kcli.o flib_spline.o mixps.o  pseudopotential.o sys.o \
            io.o pxf.o parallel.o periodic_table.o

mixps: $(MPI_INTERFACE)  $(MIXPS_OBJS)
	$(FC) -o mixps \
	       $(LDFLAGS) $(MIXPS_OBJS) $(MPI_INTERFACE) $(LIBS) 
#----------------------------------------------------------------------------
FRACTIONAL_OBJS= f2kcli.o flib_spline.o fractional.o  pseudopotential.o sys.o \
            io.o pxf.o parallel.o periodic_table.o

fractional: $(MPI_INTERFACE)  $(FRACTIONAL_OBJS)
	$(FC) -o fractional \
	       $(LDFLAGS) $(FRACTIONAL_OBJS) $(MPI_INTERFACE) $(LIBS) 
#
#------------------------------------------------------------------------
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta gen-basis ioncat denchar *.o  *.a *.pcl *.pc
	rm -f mixps mpi_test
	rm -f *.mod
	rm -f *.bck *~
	rm -f aux_*.f aux_*.f90
	(cd fdf ; $(MAKE) clean)
	(cd wxml ; $(MAKE) clean)
	(cd xmlparser ; $(MAKE) clean)
	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi
	(cd NetCDF ; $(MAKE) clean )
	(cd Libs ; $(MAKE) clean )
#
##pristine:  clean
##	cp -fp Include/constr.f .
#

# DO NOT DELETE THIS LINE - used by make depend
alloc.o: debugmpi.o parallel.o precision.o sys.o
arw.o: alloc.o parallel.o precision.o sys.o
atm_transfer.o: atm_types.o atmparams.o old_atmfuncs.o periodic_table.o
atm_transfer.o: radial.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atmparams.o basis_specs.o basis_types.o old_atmfuncs.o periodic_table.o
atom.o: precision.o pseudopotential.o recipes.o sys.o xcmod.o
atomlist.o: alloc.o atm_types.o atmfuncs.o parallel.o precision.o siesta_geom.o
atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o sys.o
automatic_cell.o: atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o densematrix.o files.o parallel.o parallelsubs.o precision.o
bands.o: sys.o
basis_enthalpy.o: atmfuncs.o atomlist.o precision.o siesta_geom.o units.o
basis_io.o: atm_types.o atmparams.o basis_types.o chemical.o precision.o
basis_io.o: pseudopotential.o radial.o sys.o xml.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o sys.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o sys.o
bessph.o: parallel.o precision.o sys.o
bonds.o: alloc.o atmfuncs.o mneighb.o precision.o recipes.o sorting.o units.o
born_charge.o: atomlist.o ksv.o ksvinit.o m_mpi_utils.o m_spin.o parallel.o
born_charge.o: siesta_geom.o siesta_options.o sparse_matrices.o
broadcast_basis.o: atm_types.o parallel.o radial.o
broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
broyden_optim.o: units.o
cart2frac.o: sys.o
cdiag.o: alloc.o memoryinfo.o parallel.o precision.o sys.o
cell_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
cell_broyden_optim.o: units.o zmatrix.o
cellxc.o: alloc.o cellxc_mod.o meshsubs.o moremeshsubs.o parallel.o
cellxc.o: parallelsubs.o precision.o sys.o xcmod.o
cellxc_mod.o: sys.o xcmod.o
cgvc.o: alloc.o conjgr_old.o m_mpi_utils.o parallel.o precision.o sys.o units.o
cgvc_zmatrix.o: alloc.o conjgr.o m_mpi_utils.o parallel.o precision.o sys.o
cgvc_zmatrix.o: units.o zmatrix.o
cgwf.o: onmod.o precision.o sys.o
chemical.o: parallel.o precision.o sys.o
chempot.o: alloc.o mneighb.o parallel.o parallelsubs.o precision.o sys.o
chkdim.o: sys.o
chkgmx.o: minvec.o parallel.o precision.o sys.o
coceri.o: files.o periodic_table.o precision.o
compute_dm.o: atomlist.o domain_decom.o jacobi_davidson.o kpoint_grid.o
compute_dm.o: m_energies.o m_eo.o m_gamma.o m_mpi_utils.o m_rmaxh.o m_spin.o
compute_dm.o: m_steps.o ordern.o precision.o siesta_geom.o siesta_options.o
compute_dm.o: sparse_matrices.o sys.o
conjgr.o: precision.o
conjgr_old.o: precision.o
coor.o: alloc.o parallel.o precision.o siesta_geom.o sys.o zmatrix.o
coxmol.o: files.o periodic_table.o precision.o
debugmpi.o: parallel.o
denchar.o: basis_io.o f2kcli.o listsc.o m_denchar_geom.o m_denchar_init.o
denchar.o: m_denchar_io.o m_denchar_work.o parallel.o precision.o sys.o
denmat.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
denmatlomem.o: alloc.o globalise.o onmod.o precision.o
densematrix.o: alloc.o precision.o
detover.o: alloc.o parallel.o parallelsubs.o precision.o
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
dfscf.o: meshsubs.o meshsubs.o parallel.o precision.o sys.o
dhscf.o: alloc.o atmfuncs.o cellxc_mod.o forhar.o m_dhscf.o m_iorho.o
dhscf.o: meshsubs.o moremeshsubs.o parallel.o parsing.o precision.o
dhscf.o: siesta_options.o sys.o xcmod.o
diag2g.o: fermid.o parallel.o parallelsubs.o precision.o sys.o
diag2k.o: fermid.o parallel.o parallelsubs.o precision.o sys.o
diagg.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagg.o: writewave.o
diagk.o: fermid.o parallel.o parallelsubs.o precision.o sys.o writewave.o
diagkp.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagon.o: alloc.o densematrix.o m_diagon.o m_memory.o parallel.o parallelsubs.o
diagon.o: precision.o sys.o
diagpol.o: precision.o sys.o
diagsprl.o: atomlist.o fermid.o parallel.o parallelsubs.o precision.o
diagsprl.o: siesta_geom.o sys.o
dipole.o: meshsubs.o parallel.o precision.o sys.o
distriphionmesh.o: alloc.o atmfuncs.o cellxc_mod.o meshsubs.o moremeshsubs.o
distriphionmesh.o: parallel.o precision.o
dnaefs.o: atmfuncs.o mneighb.o precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o sys.o
dscfcomm.o: alloc.o parallel.o parallelsubs.o precision.o schecomm.o
dynamics.o: alloc.o atomlist.o files.o ioxv.o m_mpi_utils.o parallel.o
dynamics.o: precision.o sys.o units.o
efield.o: atmfuncs.o meshsubs.o parallel.o parallelsubs.o parsing.o precision.o
efield.o: siesta_cml.o sys.o
egandd.o: alloc.o globalise.o onmod.o onmod.o precision.o sys.o
eggbox.o: parallel.o precision.o
electrostatic.o: atm_types.o atmfuncs.o atmparams.o precision.o radfft.o
electrostatic.o: radial.o recipes.o sys.o
ener3.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
ener3lomem.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
extrapol.o: alloc.o parallel.o precision.o sys.o
extrapolon.o: parallel.o precision.o sys.o
fermid.o: parallel.o precision.o recipes.o sys.o
fft3d.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
final_H_f_stress.o: alloc.o atmfuncs.o atomlist.o dnaefs.o kinefsm.o m_dipol.o
final_H_f_stress.o: m_energies.o m_forces.o m_gamma.o m_mpi_utils.o m_ntm.o
final_H_f_stress.o: m_spin.o m_steps.o m_stress.o metaforce.o
final_H_f_stress.o: molecularmechanics.o naefs.o nlefsm.o normalize_dm.o
final_H_f_stress.o: overfsm.o siesta_cml.o siesta_geom.o siesta_options.o
final_H_f_stress.o: sparse_matrices.o sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
fixed.o: alloc.o parallel.o precision.o
flib_spline.o: precision.o
forhar.o: alloc.o meshsubs.o moremeshsubs.o precision.o
fractional.o: f2kcli.o periodic_table.o precision.o pseudopotential.o sys.o
gen-basis.o: atom.o basis_io.o basis_specs.o basis_types.o old_atmfuncs.o
gen-basis.o: parallel.o precision.o sys.o xcmod.o
globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o
gmres.o: alloc.o matrixop.o precision.o precond.o sys.o
gradient.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
gradientlomem.o: alloc.o globalise.o onmod.o precision.o
grdsam.o: alloc.o m_dhscf.o parallel.o precision.o sys.o
hsparse.o: alloc.o atmfuncs.o atomlist.o domain_decom.o listsc.o mneighb.o
hsparse.o: parallel.o parallelsubs.o precision.o sorting.o sparse_matrices.o
hsparse.o: sys.o
idiag.o: parallel.o sys.o
initatom.o: atmparams.o atom.o basis_io.o basis_specs.o basis_types.o
initatom.o: electrostatic.o old_atmfuncs.o precision.o
initdm.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o
initmesh.o: alloc.o atmfuncs.o meshsubs.o moremeshsubs.o parallel.o
initmesh.o: parallelsubs.o precision.o siesta_cml.o
initparallel.o: alloc.o domain_decom.o parallel.o parallelsubs.o precision.o
initparallel.o: spatial.o sys.o
io.o: sys.o
iocg.o: files.o precision.o sys.o
iodm.o: alloc.o files.o parallel.o parallelsubs.o precision.o sys.o
ioeig.o: files.o precision.o siesta_cml.o units.o
iofa.o: files.o precision.o
iohs.o: alloc.o files.o parallel.o parallelsubs.o precision.o sys.o
iokp.o: files.o precision.o
iolwf.o: alloc.o files.o onmod.o parallel.o parallelsubs.o precision.o sys.o
iomd.o: files.o precision.o
ioncat.o: atm_types.o atmfuncs.o basis_io.o m_getopts.o precision.o
iopipes.o: parallel.o precision.o sys.o
ioxv.o: files.o parallel.o precision.o
iozm.o: files.o parallel.o precision.o siesta_geom.o zmatrix.o
ipack.o: sys.o
jacobi_davidson.o: alloc.o domain_decom.o gmres.o m_mpi_utils.o matrixop.o
jacobi_davidson.o: parallel.o precision.o precond.o printmatrix.o sys.o
kgrid.o: parallel.o precision.o units.o
kgridinit.o: minvec.o parallel.o precision.o sys.o units.o
kinefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
kpoint_grid.o: find_kgrid.o minvec.o parallel.o precision.o siesta_cml.o
kpoint_grid.o: siesta_options.o sys.o units.o
ksv.o: alloc.o atmfuncs.o densematrix.o ksvinit.o parallel.o precision.o sys.o
ksvinit.o: alloc.o parallel.o precision.o sys.o
listsc.o: alloc.o
local_DOS.o: atomlist.o files.o kpoint_grid.o m_dhscf.o m_energies.o m_eo.o
local_DOS.o: m_forces.o m_gamma.o m_ntm.o m_spin.o parallel.o siesta_geom.o
local_DOS.o: siesta_options.o sparse_matrices.o sys.o
m_broyddj.o: alloc.o m_mpi_utils.o parallel.o precision.o sys.o
m_broyddj_nocomm.o: alloc.o parallel.o precision.o sys.o
m_broyden_mixing.o: alloc.o m_broyddj.o m_mpi_utils.o parallel.o precision.o
m_cell.o: precision.o siesta_geom.o units.o
m_check_supercell.o: minvec.o parallel.o precision.o
m_denchar_geom.o: precision.o sys.o
m_denchar_init.o: files.o
m_denchar_io.o: m_denchar_geom.o
m_denchar_neighb.o: precision.o sys.o
m_denchar_work.o: atmfuncs.o chemical.o listsc.o m_denchar_geom.o
m_denchar_work.o: m_denchar_io.o m_denchar_neighb.o parallel.o
m_dhscf.o: precision.o
m_diagon.o: alloc.o precision.o
m_dipol.o: precision.o
m_energies.o: precision.o
m_eo.o: precision.o
m_fdf_global.o: m_mpi_utils.o parallel.o precision.o
m_forces.o: precision.o
m_getopts.o: f2kcli.o
m_history.o: alloc.o
m_iorho.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
m_iostruct.o: alloc.o files.o m_mpi_utils.o parallel.o precision.o
m_iostruct.o: siesta_geom.o sys.o units.o
m_kinetic.o: precision.o
m_memory.o: memoryinfo.o parallel.o sys.o
m_mpi_utils.o: precision.o sys.o
m_options.o: precision.o
m_pulay.o: alloc.o atomlist.o files.o m_mpi_utils.o m_spin.o parallel.o
m_pulay.o: precision.o siesta_options.o sparse_matrices.o sys.o
m_rmaxh.o: precision.o
m_spin.o: alloc.o parallel.o precision.o
m_stress.o: precision.o
m_wallclock.o: m_walltime.o
madelung.o: parallel.o precision.o sys.o
matel.o: alloc.o atmfuncs.o precision.o radfft.o recipes.o spher_harm.o
matrixop.o: alloc.o domain_decom.o m_mpi_utils.o parallel.o precision.o
md_out.o: files.o precision.o sys.o units.o
md_utils.o: precision.o
memory.o: memoryinfo.o parallel.o precision.o
memoryinfo.o: precision.o
meshcomm.o: meshsubs.o precision.o
meshmatrix.o: alloc.o atomlist.o dscfcomm.o meshsubs.o parallel.o
meshmatrix.o: parallelsubs.o precision.o
meshsubs.o: alloc.o atm_types.o atmfuncs.o precision.o radial.o
metaforce.o: alloc.o parallel.o precision.o sys.o
minvec.o: precision.o sorting.o sys.o
mixer.o: atomlist.o m_broyden_mixing.o m_mpi_utils.o m_pulay.o m_spin.o
mixer.o: m_steps.o siesta_geom.o siesta_options.o sparse_matrices.o
mixps.o: f2kcli.o flib_spline.o periodic_table.o precision.o pseudopotential.o
mixps.o: sys.o
mneighb.o: alloc.o precision.o sys.o
molecularmechanics.o: alloc.o parallel.o precision.o sys.o units.o
moremeshsubs.o: alloc.o m_dhscf.o meshsubs.o parallel.o precision.o schecomm.o
moremeshsubs.o: sys.o
mpi_test.o: m_mpi_utils.o parallel.o precision.o
mulliken.o: alloc.o atmfuncs.o parallel.o parallelsubs.o precision.o
mulliken.o: siesta_cml.o
naefs.o: atmfuncs.o mneighb.o precision.o
new_dm.o: atomlist.o extrapol.o m_spin.o m_steps.o siesta_geom.o
new_dm.o: siesta_options.o sparse_matrices.o
nlefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
normalize_dm.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
normalize_dm.o: sparse_matrices.o
obc.o: alloc.o precision.o
ofc.o: alloc.o files.o precision.o units.o
old_atmfuncs.o: alloc.o atmparams.o precision.o sys.o
on_subs.o: alloc.o onmod.o onmod.o spatial.o
onmod.o: precision.o
optical.o: alloc.o atmfuncs.o densematrix.o fermid.o files.o parallel.o
optical.o: parallelsubs.o precision.o sys.o units.o
ordern.o: alloc.o chempot.o globalise.o m_mpi_utils.o onmod.o parallel.o
ordern.o: precision.o siesta_cml.o sys.o units.o
outcell.o: precision.o siesta_cml.o units.o
outcoor.o: alloc.o atmfuncs.o precision.o siesta_geom.o sys.o units.o
overfsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
overlap.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
parallelsubs.o: domain_decom.o parallel.o spatial.o sys.o
pdos.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o parallel.o
pdos.o: parallelsubs.o precision.o siesta_geom.o units.o xml.o
pdosg.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o
pdosk.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o
pdoskp.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o
periodic_table.o: precision.o sys.o
phionmesh.o: alloc.o atmfuncs.o meshsubs.o meshsubs.o parallel.o precision.o
phirphi.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
phirphi_opt.o: alloc.o atm_types.o atmfuncs.o atmparams.o mneighb.o parallel.o
phirphi_opt.o: parallelsubs.o precision.o spher_harm.o sys.o
phonon.o: alloc.o files.o parallel.o precision.o sys.o
pixmol.o: files.o periodic_table.o precision.o
plcharge.o: files.o old_atmfuncs.o parallel.o precision.o
poison.o: alloc.o parallel.o precision.o sys.o
post_scf_work.o: atomlist.o final_H_f_stress.o m_energies.o m_spin.o mneighb.o
post_scf_work.o: siesta_cml.o siesta_geom.o siesta_options.o sparse_matrices.o
precond.o: alloc.o matrixop.o precision.o
printmatrix.o: alloc.o
projected_DOS.o: atomlist.o kpoint_grid.o m_eo.o m_gamma.o m_spin.o parallel.o
projected_DOS.o: precision.o siesta_options.o sparse_matrices.o sys.o
propor.o: precision.o sys.o
proximity_check.o: mneighb.o parallel.o precision.o siesta_geom.o
proximity_check.o: siesta_options.o
pseudopotential.o: alloc.o flib_spline.o precision.o sys.o
radfft.o: alloc.o precision.o recipes.o
radial.o: precision.o xml.o
rdiag.o: alloc.o m_diagon.o memoryinfo.o parallel.o precision.o sys.o
readsp.o: parallel.o precision.o sys.o
recipes.o: sys.o
redcel.o: parallel.o precision.o sys.o
reinit.o: files.o parallel.o precision.o
remove_intramol_pressure.o: precision.o zmatrix.o
reoptical.o: parallel.o precision.o sys.o
reord.o: alloc.o precision.o
rhooda.o: atmfuncs.o atomlist.o meshsubs.o meshsubs.o precision.o
rhoofd.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
rhoofd.o: meshsubs.o meshsubs.o parallel.o precision.o sys.o
rhoofdsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
rhoofdsp.o: meshsubs.o meshsubs.o precision.o sys.o
savepsi.o: alloc.o parallel.o parallelsubs.o precision.o
scfconvergence_test.o: m_energies.o m_wallclock.o parallel.o siesta_cml.o
scfconvergence_test.o: siesta_options.o write_subs.o
schecomm.o: alloc.o
setatomnodes.o: alloc.o parallel.o precision.o spatial.o sys.o
setspatial.o: alloc.o parallel.o precision.o spatial.o
setup_H0.o: alloc.o atmfuncs.o atomlist.o dnaefs.o kinefsm.o m_dhscf.o
setup_H0.o: m_dipol.o m_energies.o m_forces.o m_gamma.o m_mpi_utils.o m_ntm.o
setup_H0.o: m_spin.o m_steps.o m_stress.o metaforce.o molecularmechanics.o
setup_H0.o: naefs.o nlefsm.o normalize_dm.o parallel.o siesta_geom.o
setup_H0.o: siesta_options.o sparse_matrices.o sys.o
setup_hamiltonian.o: alloc.o atmfuncs.o atomlist.o m_dhscf.o m_dipol.o
setup_hamiltonian.o: m_energies.o m_forces.o m_gamma.o m_mpi_utils.o m_ntm.o
setup_hamiltonian.o: m_spin.o m_steps.o m_stress.o metaforce.o
setup_hamiltonian.o: molecularmechanics.o nlefsm.o normalize_dm.o parallel.o
setup_hamiltonian.o: siesta_cml.o siesta_geom.o siesta_options.o
setup_hamiltonian.o: sparse_matrices.o sys.o
setup_kscell.o: minvec.o parallel.o precision.o units.o
shaper.o: atmfuncs.o mneighb.o precision.o
show_distribution.o: atomlist.o parallel.o parallelsubs.o siesta_geom.o sys.o
siesta.o: m_steps.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o
siesta_analysis.o: alloc.o atomlist.o bands.o basis_enthalpy.o files.o
siesta_analysis.o: kpoint_grid.o ksv.o ksvinit.o local_DOS.o m_dhscf.o
siesta_analysis.o: m_dipol.o m_energies.o m_eo.o m_forces.o m_gamma.o
siesta_analysis.o: m_mpi_utils.o m_ntm.o m_spin.o m_steps.o parallel.o
siesta_analysis.o: projected_DOS.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_analysis.o: sparse_matrices.o write_subs.o writewave.o zmatrix.o
siesta_cml.o: siesta_cmlsubs.o
siesta_cmlsubs.o: files.o parallel.o timestamp.o 
siesta_end.o: alloc.o densematrix.o fixed.o m_diagon.o m_pulay.o m_wallclock.o
siesta_end.o: meshmatrix.o moremeshsubs.o parallel.o siesta_cml.o
siesta_end.o: sparse_matrices.o timestamp.o
siesta_forces.o: alloc.o atomlist.o compute_dm.o m_spin.o m_steps.o mixer.o
siesta_forces.o: ordern.o parallel.o post_scf_work.o scfconvergence_test.o
siesta_forces.o: setup_H0.o setup_hamiltonian.o siesta_options.o
siesta_forces.o: sparse_matrices.o state_analysis.o state_init.o
siesta_geom.o: precision.o
siesta_init.o: alloc.o atomlist.o bands.o densematrix.o files.o ioxv.o
siesta_init.o: kpoint_grid.o ksvinit.o m_check_supercell.o m_energies.o m_eo.o
siesta_init.o: m_forces.o m_gamma.o m_iostruct.o m_mpi_utils.o m_rmaxh.o
siesta_init.o: m_spin.o m_steps.o m_wallclock.o metaforce.o
siesta_init.o: molecularmechanics.o parallel.o parallelsubs.o phonon.o
siesta_init.o: siesta_cmlsubs.o siesta_geom.o siesta_options.o
siesta_init.o: sparse_matrices.o struct_init.o sys.o timestamp.o writewave.o
siesta_init.o: xcmod.o zmatrix.o
siesta_move.o: atomlist.o broyden_optim.o cell_broyden_optim.o dynamics.o
siesta_move.o: iopipes.o ioxv.o m_energies.o m_forces.o m_kinetic.o m_steps.o
siesta_move.o: m_stress.o parallel.o phonon.o siesta_cml.o siesta_geom.o
siesta_move.o: siesta_options.o write_subs.o zm_broyden_optim.o zmatrix.o
siesta_options.o: files.o memoryinfo.o parallel.o precision.o siesta_cml.o
siesta_options.o: sys.o units.o
sparse_matrices.o: alloc.o precision.o
spatial.o: precision.o
spher_harm.o: alloc.o precision.o sys.o
state_analysis.o: atomlist.o born_charge.o fixed.o m_energies.o m_forces.o
state_analysis.o: m_iostruct.o m_ntm.o m_spin.o m_stress.o m_wallclock.o
state_analysis.o: parallel.o siesta_cml.o siesta_geom.o siesta_options.o
state_analysis.o: sparse_matrices.o write_subs.o zmatrix.o
state_init.o: alloc.o atomlist.o hsparse.o kpoint_grid.o m_check_supercell.o
state_init.o: m_energies.o m_eo.o m_gamma.o m_mpi_utils.o m_pulay.o m_rmaxh.o
state_init.o: m_spin.o m_steps.o new_dm.o overlap.o parallel.o
state_init.o: proximity_check.o siesta_cml.o siesta_cmlsubs.o siesta_geom.o
state_init.o: siesta_options.o sparse_matrices.o write_subs.o zmatrix.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o iopipes.o ioxv.o
struct_init.o: m_iostruct.o parallel.o periodic_table.o siesta_cml.o
struct_init.o: siesta_geom.o siesta_options.o zmatrix.o
sys.o: parallel.o
timer.o: m_walltime.o parallel.o precision.o
transition_rate.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
typecell.o: precision.o
uncell.o: precision.o
units.o: precision.o
vmat.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
vmat.o: meshsubs.o meshsubs.o parallel.o precision.o
vmatsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
vmatsp.o: meshsubs.o meshsubs.o precision.o
vmb.o: fixed.o parallel.o precision.o sys.o
write_md_record.o: atomlist.o m_energies.o m_kinetic.o m_steps.o m_stress.o
write_md_record.o: md_out.o parallel.o siesta_geom.o siesta_options.o
write_subs.o: alloc.o atomlist.o m_energies.o m_forces.o m_iostruct.o m_spin.o
write_subs.o: m_steps.o m_stress.o parallel.o precision.o siesta_cml.o
write_subs.o: siesta_geom.o siesta_options.o units.o zmatrix.o
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o parallel.o
writewave.o: parallelsubs.o precision.o siesta_geom.o sys.o
xc.o: alloc.o precision.o sys.o xcmod.o
xcmod.o: parallel.o precision.o sys.o
xml.o: precision.o
zm_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
zm_broyden_optim.o: units.o zmatrix.o
zmatrix.o: alloc.o m_cell.o parallel.o precision.o siesta_geom.o sys.o units.o
band.o: bands.o
diagmemory.o: memoryinfo.o
listsc_module.o: listsc.o
m_born_charge.o: born_charge.o
m_broyden_optim.o: broyden_optim.o
m_cell_broyden_optim.o: cell_broyden_optim.o
m_chempot.o: chempot.o
m_compute_dm.o: compute_dm.o
m_conjgr.o: conjgr.o
m_conjgr_old.o: conjgr_old.o
m_dnaefs.o: dnaefs.o
m_dscfcomm.o: dscfcomm.o
m_dynamics.o: dynamics.o
m_extrapol.o: extrapol.o
m_fermid.o: fermid.o
m_final_h_f_stress.o: final_H_f_stress.o
m_find_kgrid.o: find_kgrid.o
m_fixed.o: fixed.o
m_forhar.o: forhar.o
m_hsparse.o: hsparse.o
m_ioxv.o: ioxv.o
m_kinefsm.o: kinefsm.o
m_ksv.o: ksv.o
m_ksvinit.o: ksvinit.o
m_local_dos.o: local_DOS.o
m_minvec.o: minvec.o
m_mixer.o: mixer.o
m_naefs.o: naefs.o
m_new_dm.o: new_dm.o
m_nlefsm.o: nlefsm.o
m_normalize_dm.o: normalize_dm.o
m_ordern.o: ordern.o
m_overfsm.o: overfsm.o
m_overlap.o: overlap.o
m_post_scf_work.o: post_scf_work.o
m_projected_dos.o: projected_DOS.o
m_proximity_check.o: proximity_check.o
m_radfft.o: radfft.o
m_recipes.o: recipes.o
m_scf_options.o: m_options.o
m_scfconvergence_test.o: scfconvergence_test.o
m_setup_h0.o: setup_H0.o
m_setup_hamiltonian.o: setup_hamiltonian.o
m_siesta_analysis.o: siesta_analysis.o
m_siesta_end.o: siesta_end.o
m_siesta_forces.o: siesta_forces.o
m_siesta_init.o: siesta_init.o
m_siesta_move.o: siesta_move.o
m_state_analysis.o: state_analysis.o
m_state_init.o: state_init.o
m_struct_init.o: struct_init.o
m_timestamp.o: timestamp.o
m_xyz_phiatm.o: matel.o
m_zm_broyden_optim.o: zm_broyden_optim.o
mesh.o: meshsubs.o
meshdscf.o: meshmatrix.o
meshphi.o: meshsubs.o
neighbour.o: mneighb.o
on_core.o: onmod.o
on_main.o: onmod.o
printmat.o: printmatrix.o
version_info.o: version.o
