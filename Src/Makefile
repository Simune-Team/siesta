# 
# Siesta/Src/Makefile
# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#
# Makefile for siesta
#
# Remove implicit rules (including for Modula-2, which just confuses matters)
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90 .c
#
VPATH=.
#
default: what siesta #  fsiesta_pipes.o libSiestaForces.a
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
	@if [ $(VPATH) = "." ] ; then \
	sfmakedepend --depend=obj  --modext=o \
                              *.f *.f90 *.F *.F90 ; \
        else \
	sfmakedepend --depend=obj  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 ; \
        fi
	@sed "/^siesta_cmlsubs.o:/s/version.o//;\
		/^m_ts_io.o:/s/m_tbt_kpoints.o//"  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile
#
include arch.make
INCFLAGS:= -I. $(NETCDF_INCFLAGS) $(INCFLAGS)
# For removing any -ipo compilation in sub directories...
# Later compilers are also having this enabled, and it lets you
# specify what to not be used in libraries.
# This is so because inter-procedural compilation should NOT be performed in
# libraries (it removes names when a routine can be fully moved into others)
IPO_FLAG?= -ipo
#
# Uncomment the following line for debugging support
#
##FFLAGS=$(FFLAGS_DEBUG)
#
.PHONY: .siesta .transiesta what version
#
what:
	@echo
	@if [ "$(VPATH)" = . ]; then \
              echo "** You can no longer build SIESTA in Src." ;\
              echo "** Go to the Obj directory and see the README file." ;\
              false ;\
         fi
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Src/Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2
#
SYSOBJ=$(SYS).o pxf.o
#
# Note that machine-specific files are now in top Src directory.
#
MAIN_OBJ = siesta.o

OBJS =  automatic_cell.o atom_options.o \
        arw.o  atomlwf.o bands.o basis_enthalpy.o bessph.o bonds.o \
        born_charge.o cellxc_mod.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o coxmol.o cross.o compute_norm.o\
	denmat.o denmatlomem.o detover.o dfscf.o m_diagon.o diagon.o digcel.o \
        fft.o dhscf.o constr.o diagk_file.o \
	diagg.o diagk.o diagkp.o diag2g.o diag2k.o diagmemory.o diagpol.o \
	diagsprl.o dipole.o dismin.o dnaefs.o doping_uniform.o dot.o  \
	m_efield.o egandd.o ener3.o ener3lomem.o errorf.o extrapolon.o \
	fixed.o interpolation.o gradient.o gradientlomem.o grdsam.o \
	hsparse.o idiag.o  initatom.o inver.o iodm_netcdf.o \
	iodmhs_netcdf.o iogrid_netcdf.o iolwf.o iozm.o ipack.o iopipes.o \
	iowfs_netcdf.o kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o local_DOS.o\
	madelung.o matel.o mesh.o meshphi.o meshdscf.o memory.o meshsubs.o metaforce.o \
        m_check_supercell.o m_fft_gpfa.o mulliken.o minvec.o m_pulay.o naefs.o \
	new_dm.o normalize_dm.o ordern.o outcell.o outcoor.o \
	paste.o pdos.o pdosg.o pdosk.o pdoskp.o phirphi.o pixmol.o plcharge.o \
	projected_DOS.o propor.o m_partial_charges.o \
	randomg.o reclat.o redcel.o \
	reinit.o reord.o rhoofd.o rhoofdsp.o rhooda.o \
	savepsi.o shaper.o timer_tree.o timer.o \
	vmb.o vmat.o vmatsp.o volcel.o \
	cdiag.o rdiag.o \
        cgvc.o cgvc_zmatrix.o m_convergence.o \
        iocg.o ioeig.o iofa.o iokp.o iomd.o kpoint_pdos.o typecell.o \
        ofc.o poison.o readsp.o radfft.o \
	write_md_record.o kpoint_grid.o find_kgrid.o proximity_check.o\
	state_init.o siesta_move.o setup_hamiltonian.o compute_dm.o mixer.o\
	scfconvergence_test.o post_scf_work.o state_analysis.o write_subs.o \
	siesta_init.o struct_init.o siesta_options.o siesta_geom.o \
        siesta_analysis.o siesta_end.o siesta_forces.o \
        io.o sparse_matrices.o \
        coor.o atm_transfer.o \
        broadcast_basis.o eggbox.o dsyevds.o zheevds.o \
        optical.o phirphi_opt.o reoptical.o transition_rate.o \
        initparallel.o show_distribution.o setspatial.o setatomnodes.o uncell.o \
        cart2frac.o obc.o \
        precision.o sys.o m_cell.o files.o \
        spatial.o parallel.o parallelsubs.o  parsing.o chemical.o atom.o atmparams.o \
        m_mpi_utils.o m_fdf_global.o m_history.o m_iorho.o  \
        atmfuncs.o listsc.o memory_log.o sorting.o \
        atomlist.o atm_types.o old_atmfuncs.o radial.o \
        alloc.o phonon.o spher_harm.o periodic_table.o version.o timestamp.o \
        basis_types.o xml.o pseudopotential.o basis_specs.o basis_io.o\
        onmod.o densematrix.o writewave.o on_subs.o fermid.o m_broyddj.o \
        electrostatic.o mneighb.o globalise.o siesta_cmlsubs.o siesta_cml.o \
        units.o zmatrix.o m_broyden_mixing.o forhar.o m_walltime.o m_wallclock.o m_iostruct.o \
        nlefsm.o overfsm.o overlap.o conjgr.o conjgr_old.o m_energies.o m_steps.o\
        m_broyddj_nocomm.o broyden_optim.o ioxv.o dynamics.o md_out.o \
        molecularmechanics.o zm_broyden_optim.o \
        cell_broyden_optim.o remove_intramol_pressure.o \
        m_ntm.o m_dipol.o  m_kinetic.o  m_rmaxh.o \
        m_forces.o m_stress.o m_eo.o m_spin.o m_gamma.o m_hsx.o m_fire.o setup_H0.o \
        get_target_stress.o fire_optim.o zm_fire_optim.o \
        cell_fire_optim.o m_fire_para.o m_fire_mixing.o write_raw_efs.o \
	fsiesta_mpi.o \
        final_H_f_stress.o \
        debugmpi.o qsort.o \
        m_dscfcomm.o schecomm.o moremeshsubs.o \
        domain_decom.o printmatrix.o  mmio.o pspltm1.o lenstr.o \
        setup_ordern_indexes.o\
        m_iodm.o m_sparse.o \
	m_filter.o \
        m_io.o \
	m_timer.o \
	moreParallelSubs.o \
	read_xc_info.o \
	siesta_master.o \
        bsc_xcmod.o bsc_cellxc.o xc.o \
	vacuum_level.o \
	write_orb_indx.o \
	m_hs_matrix.o kpoint_convert.o m_glob_sparse.o\
        die.o m_dminim.o m_zminim.o m_svd.o compute_max_diff.o \
        save_density_matrix.o alloc_error_report.o\
        alloc_memory_event.o\
        m_rhog.o rhofft.o m_diis.o compute_energies.o \
        class_Vector.o class_Pair_Vectors.o class_Fstack_Pair_Vectors.o
#--------------------------------------------------------------
# FoX whatnot
#
# FoX build targets:
FoX_configured=FoX/.config
FoX_built=FoX/.FoX
#
# This is how we pick up modules and libraries for FoX:
FoX_FCFLAGS=`FoX/FoX-config --fcflags`
FoX_LIBS=`FoX/FoX-config --libs --wcml`
#
# And add them to global compiler flags:
INCFLAGS:=$(INCFLAGS) $(FoX_FCFLAGS) 
#
# First, it needs to be configured. This may have been done
# by the main SIESTA configure, but in case not:
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
#
$(FoX_configured):
	(cd FoX; touch arch.make ; \
         CONFIGURE="$(VPATH)/FoX/configure"; \
         $$CONFIGURE VPATH="$(VPATH)/FoX" \
         FC="$(FC_SERIAL)" FCFLAGS="$(FFLAGS:$(IPO_FLAG)=)" \
         --enable-wcml $(DUMMY_FOX) || false )
#
# Note we have to make sure to use the same compiler as SIESTA,
# and pick up all the same FFLAGS, and also remember to maybe ask for a dummy library.
# Note also that the "false" clause will stop the 'make' process if the configuration fails.
#
# Then we want to make FoX itself. Like so:
$(FoX_built): $(FoX_configured)
	(cd FoX; $(MAKE) )
#
# Now ensure that building FoX is inserted into the dependency chain.
# Routines using FoX:
# version (compinfo.o) does not actually need FoX, but making it build after
# FoX is done avoids a warning as we use the FoX-config script everywhere.
#
version $(SYSOBJ) $(OBJS) \
basis_io.o efield.o initatom.o ioeig.o kpoint_grid.o meshsubs.o \
mulliken.o ordern.o outcell.o pdos.o post_scf_work.o radial.o   \
scfconvergence_test.o setup_hamiltonian.o siesta_analysis.o     \
siesta_cml.o siesta_cmlsubs.o siesta_end.o siesta_forces.o      \
siesta_init.o siesta_move.o siesta_options.o state_analysis.o   \
state_init.o struct_init.o write_subs.o xml.o: $(FoX_built)
#
#--------------------------------------------------------------
XMLPARSER=libxmlparser.a
$(XMLPARSER): 
	(cd xmlparser ; $(MAKE) "VPATH=$(VPATH)/xmlparser" \
	FFFLAGS="$(FFLAGS:$(IPO_FLAG)=)" module)
#
#--------------------------------------------------------------
XC=gridxc/libGridXC.a
$(XC): 
	(cd gridxc ; $(MAKE) \
                    "ARCH_MAKE=../arch.make"  \
                    "VPATH=$(VPATH)/gridxc"  \
		    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module )
INCFLAGS:=  $(INCFLAGS) -I./gridxc
#--------------------------------------------------------------
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "VPATH=$(VPATH)/fdf" \
		    "FPPFLAGS=$(FPPFLAGS)" \
		    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module )
#
# Routines using fdf calls.
#
dhscf.o iolwf.o  grdsam.o : $(FDF)
recoor.o outcoor.o ioxv.o m_iostruct.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o pdos.o plcharge.o: $(FDF)
coor.o cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o  ofc.o : $(FDF)
readsp.o eggbox.o iopipes.o : $(FDF)
#
atom.o basis_io.o basis_specs.o: $(FDF)
#
#
ALL_OBJS=$(OBJS) $(SYSOBJ) $(MAIN_OBJ)
#
# Interfaces to libraries
#
libmpi_f90.a: 
	@(cd MPI ; $(MAKE) \
                    "VPATH=$(VPATH)/MPI" \
                    "FPP_FLAGS=$(FPP_FLAGS)" \
                    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" )
#
# Libraries that might need to be compiled
#
libblas.a: Libs/blas.f
	@echo "==> Compiling libblas.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libblas.a)
liblapack.a: Libs/lapack.f
	@echo "==> Compiling liblapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" liblapack.a)
linalg.a: Libs/blas.f Libs/lapack.f
	@echo "==> Compiling linalg.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" linalg.a)
dc_lapack.a:  Libs/dc_lapack.f
	@echo "==> Compiling dc_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" dc_lapack.a)
zgesv_lapack.a:  Libs/zgesv_lapack.f
	@echo "==> Compiling zgesv_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" zgesv_lapack.a)
#
essl_lapack.a: Libs/essl_lapack.f
	@echo "==> Compiling essl_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" essl_lapack.a)
#

version: version.F90
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;\
		s'SIESTA_VERSION'$$(cat ${<D}/../version.info)'g;\
		s'FPPFLAGS'$(FPPFLAGS) $(DEFS)'g" $<  > tmp.F90
	@awk '{if (length>80) { cur=78; \
		printf "%s&\n",substr($$0,0,78); \
		while(length-cur>78){ cur=cur+76 ; \
		printf "&%s&\n",substr($$0,cur-76,76) \
		} printf "&%s\n",substr($$0,cur)} else { print $$0 }}' tmp.F90 > compinfo.F90
	@($(MAKE) "FPPFLAGS=$(FPPFLAGS)" "DEFS=$(DEFS)" compinfo.o)
	@mv compinfo.o version.o
	@echo
#
check-siesta:
	@if [ -f .transiesta ]; then \
              echo "** You were building transiesta in this directory." ;\
              echo "** Please type 'make clean' first " ;\
              false ;\
	fi
	@touch .siesta
#
siesta: check-siesta what version \
                $(MPI_INTERFACE) $(FDF) $(WXML) $(XMLPARSER) $(XC) \
                $(COMP_LIBS) $(ALL_OBJS) 
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(WXML) $(XMLPARSER) \
               $(XC) $(MPI_INTERFACE) \
               $(COMP_LIBS) $(FoX_LIBS) $(LIBS) 
#
libSiestaForces.a: $(OBJS) $(SYSOBJ)
	ar $(ARFLAGS_EXTRA) cru libSiestaForces.a $(OBJS) $(SYSOBJ)
	@if [ ! -z "$(MPI_INTERFACE)" ] ; then \
            ar $(ARFLAGS_EXTRA) cru libSiestaForces.a MPI/*.o ; fi
	-$(RANLIB) libSiestaForces.a
#
#----------------------------------------------------------------------------
# TranSIESTA...
#
TS_OBJS= m_ts_global_vars.o m_ts_options.o m_ts_aux_rout.o m_ts_kpoints.o \
         m_ts_iodm.o m_ts_io.o m_ts_scattering.o m_ts_contour.o \
         m_ts_electrode.o m_ts_in_siesta.o m_ts_voltage.o m_ts_gf.o \
	 m_ts_cctype.o ts_init.o ts_show_regions.o
#
transiesta: FPPFLAGS+=$(DEFS_PREFIX)-DTRANSIESTA 
transiesta: DEFS+=$(DEFS_PREFIX)-DTRANSIESTA
check-ts:
	@if [ -f .siesta ]; then \
              echo "** You were building plain siesta in this directory." ;\
              echo "** Please type 'make clean' first" ;\
              false ;\
	fi
	@touch .transiesta
#
transiesta: check-ts what version \
            $(MPI_INTERFACE) $(NETCDF_INTERFACE) $(FDF) $(XC) \
            $(COMP_LIBS) $(ALL_OBJS) $(TS_OBJS)
	$(FC) -o transiesta \
               $(LDFLAGS) $(ALL_OBJS) $(TS_OBJS) $(FDF) $(XC) \
               $(MPI_INTERFACE) $(NETCDF_INTERFACE) \
               $(COMP_LIBS) $(FoX_LIBS) $(LIBS)
#
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta transiesta *.o  *.a *.pcl *.pc
	rm -f .siesta .transiesta
	rm -f *.mod
	rm -f *.bck *~
	rm -f aux_*.f aux_*.f90
	(cd fdf ; $(MAKE) clean)
	(cd FoX ; touch arch.make ; $(MAKE) distclean)
	(cd xmlparser ; $(MAKE) clean)
	(cd wxml ; $(MAKE) clean)
	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi
	(cd Libs ; $(MAKE) clean )
	(cd gridxc; $(MAKE) \
                    "ARCH_MAKE=../arch.make"  \
                    "VPATH=$(VPATH)/gridxc"  \
		    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" clean )

# DO NOT DELETE THIS LINE - used by make depend
alloc_error_report.o: parallel.o
alloc_memory_event.o: memory_log.o
arw.o: alloc.o parallel.o precision.o sys.o
atm_transfer.o: atm_types.o atmparams.o old_atmfuncs.o periodic_table.o
atm_transfer.o: radial.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o bsc_xcmod.o
atom.o: interpolation.o m_filter.o old_atmfuncs.o periodic_table.o precision.o
atom.o: pseudopotential.o sys.o
atom_options.o: sys.o
atomlist.o: alloc.o atm_types.o atmfuncs.o parallel.o precision.o siesta_geom.o
atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o sys.o
automatic_cell.o: atm_types.o atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o m_gamma.o
bands.o: parallel.o parallelsubs.o precision.o siesta_geom.o sys.o writewave.o
basis_enthalpy.o: atmfuncs.o atomlist.o precision.o siesta_geom.o units.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o sys.o xml.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o sys.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o sys.o
bessph.o: precision.o sys.o
bonds.o: alloc.o atmfuncs.o mneighb.o precision.o sorting.o units.o
born_charge.o: atomlist.o ksv.o ksvinit.o m_mpi_utils.o m_spin.o parallel.o
born_charge.o: siesta_geom.o siesta_options.o sparse_matrices.o
broadcast_basis.o: atm_types.o parallel.o radial.o
broyden_optim.o: m_broyddj_nocomm.o memory_log.o parallel.o precision.o
broyden_optim.o: siesta_options.o sys.o units.o
bsc_cellxc.o: alloc.o bsc_xcmod.o cellxc_mod.o mesh.o moremeshsubs.o parallel.o
bsc_cellxc.o: parallelsubs.o precision.o sys.o
bsc_xcmod.o: parallel.o precision.o sys.o
cart2frac.o: sys.o
cdiag.o: alloc.o diagmemory.o m_spin.o parallel.o precision.o sys.o
cell_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
cell_broyden_optim.o: units.o zmatrix.o
cell_fire_optim.o: alloc.o m_fdf_global.o m_fire.o memory_log.o parallel.o
cell_fire_optim.o: precision.o siesta_options.o sys.o zmatrix.o
cellsubs.o: precision.o
cellxc_mod.o: bsc_xcmod.o sys.o
cgvc.o: alloc.o conjgr_old.o m_mpi_utils.o parallel.o precision.o sys.o units.o
cgvc_zmatrix.o: alloc.o conjgr.o m_mpi_utils.o parallel.o precision.o sys.o
cgvc_zmatrix.o: units.o zmatrix.o
cgwf.o: egandd.o onmod.o parallel.o precision.o sys.o
chemical.o: parallel.o precision.o sys.o
chempot.o: alloc.o mneighb.o parallel.o parallelsubs.o precision.o sys.o
chkdim.o: sys.o
chkgmx.o: cellsubs.o minvec.o precision.o sys.o
class_Fstack_Pair_Vectors.o: Fstack.T90 basic_type.inc
class_Fstack_Pair_Vectors.o: class_Pair_Vectors.o
class_Pair_Vectors.o: Pair.T90 basic_type.inc
class_Pair_Vectors.o: class_Vector.o
class_Vector.o: basic_type.inc
class_Vector.o: alloc.o
coceri.o: files.o periodic_table.o precision.o
compute_dm.o: atomlist.o iodmhs_netcdf.o kpoint_grid.o m_dminim.o m_energies.o
compute_dm.o: m_eo.o m_gamma.o m_rmaxh.o m_spin.o m_steps.o m_ts_global_vars.o
compute_dm.o: m_ts_in_siesta.o m_ts_options.o m_zminim.o normalize_dm.o
compute_dm.o: ordern.o parallel.o precision.o siesta_geom.o siesta_options.o
compute_dm.o: sparse_matrices.o sys.o units.o
compute_energies.o: atomlist.o dhscf.o files.o m_dipol.o m_energies.o
compute_energies.o: m_mpi_utils.o m_ntm.o m_rhog.o m_spin.o precision.o
compute_energies.o: siesta_geom.o siesta_options.o sparse_matrices.o
compute_max_diff.o: atomlist.o m_mpi_utils.o m_spin.o precision.o
compute_max_diff.o: sparse_matrices.o
compute_norm.o: m_mpi_utils.o precision.o sparse_matrices.o
compute_rhog.o: alloc.o atomlist.o dhscf.o files.o m_dipol.o m_energies.o
compute_rhog.o: m_gamma.o m_ntm.o m_spin.o m_stress.o siesta_geom.o
compute_rhog.o: siesta_options.o sparse_matrices.o
conjgr.o: precision.o
conjgr_old.o: precision.o
coor.o: alloc.o parallel.o precision.o siesta_geom.o sys.o zmatrix.o
coxmol.o: files.o periodic_table.o precision.o
debugmpi.o: parallel.o
denmat.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
denmatlomem.o: alloc.o globalise.o onmod.o precision.o
densematrix.o: alloc.o precision.o
detover.o: alloc.o parallel.o parallelsubs.o precision.o
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
dfscf.o: meshphi.o parallel.o parallelsubs.o precision.o sys.o
dhscf.o: alloc.o atmfuncs.o bsc_xcmod.o cellxc_mod.o doping_uniform.o files.o
dhscf.o: forhar.o iogrid_netcdf.o m_efield.o m_iorho.o m_partial_charges.o
dhscf.o: m_rhog.o m_ts_global_vars.o m_ts_in_siesta.o m_ts_options.o
dhscf.o: m_ts_voltage.o mesh.o meshdscf.o meshsubs.o moremeshsubs.o parallel.o
dhscf.o: parsing.o precision.o rhofft.o siesta_options.o sys.o units.o vmat.o
diag2g.o: fermid.o m_spin.o parallel.o parallelsubs.o precision.o sys.o
diag2k.o: fermid.o m_spin.o parallel.o parallelsubs.o precision.o sys.o
diagg.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagg.o: writewave.o
diagk.o: compute_norm.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagk.o: writewave.o
diagk_file.o: fermid.o iowfs_netcdf.o parallel.o parallelsubs.o precision.o
diagk_file.o: sys.o writewave.o
diagkp.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagmemory.o: precision.o
diagon.o: alloc.o densematrix.o m_diagon.o memory_log.o parallel.o
diagon.o: parallelsubs.o precision.o siesta_options.o sys.o
diagpol.o: precision.o sys.o
diagsprl.o: atomlist.o fermid.o parallel.o parallelsubs.o precision.o
diagsprl.o: siesta_geom.o sys.o
die.o: parallel.o
dipole.o: mesh.o parallel.o precision.o sys.o
dnaefs.o: atmfuncs.o mneighb.o precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o sys.o
doping_uniform.o: alloc.o m_ntm.o mesh.o parallel.o precision.o sys.o
dynamics.o: alloc.o atomlist.o files.o ioxv.o m_mpi_utils.o parallel.o
dynamics.o: precision.o sys.o units.o
egandd.o: alloc.o denmat.o ener3.o globalise.o gradient.o on_subs.o onmod.o
egandd.o: onmod.o precision.o sys.o
eggbox.o: parallel.o precision.o
electrostatic.o: atm_types.o atmfuncs.o atmparams.o bessph.o interpolation.o
electrostatic.o: precision.o radfft.o radial.o sys.o
ener3.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
ener3lomem.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
extrapolon.o: parallel.o precision.o sys.o
fermid.o: errorf.o parallel.o precision.o sys.o
fft.o: alloc.o m_fft_gpfa.o m_timer.o mesh.o parallel.o parallelsubs.o
fft.o: precision.o sys.o
fft1d.o: parallel.o precision.o sys.o
final_H_f_stress.o: alloc.o atmfuncs.o atomlist.o dnaefs.o files.o grdsam.o
final_H_f_stress.o: kinefsm.o m_dipol.o m_energies.o m_forces.o m_gamma.o
final_H_f_stress.o: m_hsx.o m_mpi_utils.o m_ntm.o m_spin.o m_steps.o m_stress.o
final_H_f_stress.o: m_ts_global_vars.o m_ts_io.o m_ts_kpoints.o m_ts_options.o
final_H_f_stress.o: metaforce.o molecularmechanics.o naefs.o nlefsm.o overfsm.o
final_H_f_stress.o: parallel.o siesta_geom.o siesta_options.o sparse_matrices.o
final_H_f_stress.o: sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_mpi_utils.o parallel.o
fire_optim.o: precision.o siesta_options.o units.o
fixed.o: alloc.o parallel.o precision.o
forhar.o: alloc.o mesh.o moremeshsubs.o parallel.o precision.o
fsiesta_mpi.o: precision.o siesta_analysis.o siesta_end.o siesta_forces.o
fsiesta_mpi.o: siesta_init.o siesta_master.o siesta_move.o sys.o
get_target_stress.o: parallel.o precision.o sys.o units.o
globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o
gradient.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
gradientlomem.o: alloc.o globalise.o onmod.o precision.o
grdsam.o: alloc.o dhscf.o files.o m_partial_charges.o parallel.o precision.o
grdsam.o: siesta_geom.o siesta_options.o sys.o
hsparse.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mneighb.o
hsparse.o: parallel.o parallelsubs.o precision.o siesta_options.o sorting.o
hsparse.o: sparse_matrices.o sys.o
idiag.o: parallel.o sys.o
initatom.o: atmparams.o atom.o atom_options.o basis_io.o basis_specs.o
initatom.o: basis_types.o electrostatic.o old_atmfuncs.o precision.o
initparallel.o: alloc.o domain_decom.o m_fdf_global.o parallel.o parallelsubs.o
initparallel.o: precision.o siesta_options.o spatial.o sys.o
io.o: sys.o
iocg.o: files.o precision.o sys.o
iodm.o: alloc.o files.o parallel.o parallelsubs.o precision.o sys.o
iodm_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
iodmhs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
ioeig.o: files.o precision.o siesta_cml.o units.o
iofa.o: files.o precision.o
iogrid_netcdf.o: alloc.o parallel.o precision.o sys.o
iokp.o: files.o precision.o
iolwf.o: alloc.o files.o onmod.o parallel.o parallelsubs.o precision.o sys.o
iomd.o: files.o precision.o
iopipes.o: parallel.o precision.o sys.o
iowfs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
ioxv.o: files.o parallel.o precision.o
iozm.o: files.o parallel.o precision.o siesta_geom.o zmatrix.o
ipack.o: sys.o
kgrid.o: parallel.o precision.o units.o
kgridinit.o: minvec.o parallel.o precision.o sys.o units.o
kinefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
kpoint_convert.o: precision.o sys.o units.o
kpoint_grid.o: find_kgrid.o minvec.o parallel.o precision.o siesta_cml.o
kpoint_grid.o: siesta_options.o sys.o units.o
kpoint_pdos.o: find_kgrid.o minvec.o parallel.o precision.o siesta_options.o
kpoint_pdos.o: sys.o units.o
ksv.o: alloc.o atmfuncs.o densematrix.o ksvinit.o parallel.o precision.o sys.o
ksvinit.o: alloc.o parallel.o precision.o sys.o
listsc.o: alloc.o
local_DOS.o: atomlist.o dhscf.o files.o kpoint_grid.o m_energies.o m_eo.o
local_DOS.o: m_forces.o m_gamma.o m_ntm.o m_spin.o parallel.o siesta_geom.o
local_DOS.o: siesta_options.o sparse_matrices.o sys.o
m_broyddj.o: alloc.o m_mpi_utils.o parallel.o precision.o sys.o
m_broyddj_nocomm.o: alloc.o m_fdf_global.o parallel.o precision.o sys.o
m_broyden_mixing.o: alloc.o m_broyddj.o m_mpi_utils.o parallel.o precision.o
m_cell.o: precision.o siesta_geom.o units.o
m_check_supercell.o: minvec.o parallel.o precision.o
m_diagon.o: precision.o
m_diis.o: class_Fstack_Pair_Vectors.o class_Pair_Vectors.o class_Vector.o
m_diis.o: m_svd.o parallel.o precision.o
m_dipol.o: precision.o
m_dminim.o: densematrix.o files.o parallel.o parallelsubs.o precision.o
m_dminim.o: siesta_options.o sys.o
m_dscfcomm.o: alloc.o parallel.o parallelsubs.o precision.o schecomm.o
m_efield.o: atmfuncs.o mesh.o parallel.o precision.o siesta_cml.o siesta_geom.o
m_efield.o: sys.o units.o
m_energies.o: precision.o
m_eo.o: precision.o
m_fdf_global.o: m_mpi_utils.o parallel.o precision.o
m_filter.o: bessph.o precision.o radfft.o sorting.o sys.o
m_fire.o: parallel.o precision.o
m_fire_mixing.o: alloc.o m_fdf_global.o m_fire_para.o m_mpi_utils.o parallel.o
m_fire_mixing.o: precision.o
m_fire_para.o: m_mpi_utils.o parallel.o precision.o
m_forces.o: precision.o
m_getopts.o: f2kcli.o
m_glob_sparse.o: parallel.o parallelsubs.o precision.o
m_history.o: alloc.o
m_hs_matrix.o: alloc.o precision.o sys.o
m_hsx.o: atm_types.o atmfuncs.o atomlist.o files.o parallel.o parallelsubs.o
m_hsx.o: precision.o siesta_geom.o sys.o
m_io.o: sys.o
m_iodm.o: alloc.o files.o m_mpi_utils.o parallel.o parallelsubs.o precision.o
m_iodm.o: sys.o
m_iorho.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
m_iostruct.o: alloc.o files.o m_mpi_utils.o parallel.o precision.o
m_iostruct.o: siesta_geom.o sys.o units.o
m_kinetic.o: precision.o
m_mpi_utils.o: precision.o sys.o
m_options.o: precision.o
m_partial_charges.o: atmfuncs.o atomlist.o mesh.o meshphi.o parallel.o
m_partial_charges.o: precision.o siesta_geom.o siesta_options.o sys.o
m_pulay.o: alloc.o atomlist.o m_mpi_utils.o m_spin.o m_svd.o parallel.o
m_pulay.o: precision.o siesta_options.o sparse_matrices.o sys.o
m_rhog.o: alloc.o atomlist.o class_Fstack_Pair_Vectors.o class_Pair_Vectors.o
m_rhog.o: class_Vector.o m_diis.o m_mpi_utils.o m_spin.o parallel.o precision.o
m_rhog.o: siesta_options.o sorting.o
m_rmaxh.o: precision.o
m_sparse.o: alloc.o precision.o
m_spin.o: alloc.o parallel.o precision.o
m_stress.o: precision.o
m_svd.o: precision.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o sys.o
m_ts_aux_rout.o: precision.o sys.o
m_ts_cctype.o: precision.o
m_ts_contour.o: m_ts_aux_rout.o m_ts_cctype.o parallel.o precision.o sys.o
m_ts_contour.o: units.o
m_ts_electrode.o: files.o m_hs_matrix.o m_ts_aux_rout.o m_ts_cctype.o m_ts_io.o
m_ts_electrode.o: m_ts_kpoints.o parallel.o precision.o sys.o units.o
m_ts_gf.o: m_ts_cctype.o m_ts_electrode.o parallel.o precision.o sys.o units.o
m_ts_global_vars.o: precision.o
m_ts_in_siesta.o: m_glob_sparse.o m_hs_matrix.o m_ts_cctype.o m_ts_contour.o
m_ts_in_siesta.o: m_ts_gf.o m_ts_kpoints.o m_ts_options.o m_ts_scattering.o
m_ts_in_siesta.o: parallel.o parallelsubs.o precision.o sys.o units.o
m_ts_io.o: files.o m_spin.o m_ts_kpoints.o parallel.o parallelsubs.o
m_ts_io.o: precision.o sys.o
m_ts_iodm.o: alloc.o files.o m_mpi_utils.o m_ts_global_vars.o parallel.o
m_ts_iodm.o: parallelsubs.o precision.o sys.o
m_ts_kpoints.o: files.o find_kgrid.o minvec.o parallel.o precision.o
m_ts_kpoints.o: siesta_cml.o siesta_options.o sys.o
m_ts_options.o: files.o m_fdf_global.o m_ts_contour.o m_ts_global_vars.o
m_ts_options.o: m_ts_io.o parallel.o precision.o siesta_options.o sys.o units.o
m_ts_scattering.o: m_ts_aux_rout.o parallel.o precision.o sys.o units.o
m_ts_voltage.o: m_ts_options.o mesh.o parallel.o precision.o units.o
m_wallclock.o: m_walltime.o
m_zminim.o: files.o parallel.o parallelsubs.o precision.o siesta_options.o
m_zminim.o: sys.o
madelung.o: parallel.o precision.o sys.o
matel.o: alloc.o atmfuncs.o errorf.o interpolation.o parallel.o precision.o
matel.o: radfft.o spher_harm.o
md_out.o: files.o m_energies.o precision.o sys.o units.o
md_utils.o: precision.o
memory.o: memory_log.o parallel.o
memory_log.o: m_io.o parallel.o precision.o
mesh.o: precision.o
meshcomm.o: mesh.o precision.o
meshdscf.o: alloc.o atomlist.o m_dscfcomm.o meshphi.o parallel.o parallelsubs.o
meshdscf.o: precision.o
meshphi.o: alloc.o precision.o
meshsubs.o: alloc.o atm_types.o atmfuncs.o cellsubs.o cellxc_mod.o chkgmx.o
meshsubs.o: m_fft_gpfa.o mesh.o meshphi.o moremeshsubs.o parallel.o
meshsubs.o: parallelsubs.o precision.o radial.o siesta_cml.o sys.o
metaforce.o: alloc.o parallel.o precision.o sys.o
minvec.o: cellsubs.o precision.o sorting.o sys.o
mixer.o: atomlist.o m_broyden_mixing.o m_fire_mixing.o m_pulay.o m_spin.o
mixer.o: m_ts_global_vars.o parallel.o siesta_geom.o siesta_options.o
mixer.o: sparse_matrices.o
mneighb.o: alloc.o precision.o sys.o
molecularmechanics.o: alloc.o parallel.o precision.o sys.o units.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o sys.o
moremeshsubs.o: alloc.o mesh.o parallel.o precision.o schecomm.o sys.o
mulliken.o: alloc.o atmfuncs.o parallel.o parallelsubs.o precision.o
mulliken.o: siesta_cml.o
naefs.o: atmfuncs.o mneighb.o precision.o
new_dm.o: alloc.o atomlist.o m_energies.o m_fdf_global.o m_iodm.o m_mpi_utils.o
new_dm.o: m_sparse.o m_spin.o m_steps.o m_ts_global_vars.o m_ts_iodm.o
new_dm.o: m_ts_options.o normalize_dm.o parallel.o parallelsubs.o parsing.o
new_dm.o: precision.o siesta_geom.o siesta_options.o sparse_matrices.o sys.o
new_dm.o: units.o
nlefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
normalize_dm.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
normalize_dm.o: siesta_options.o sparse_matrices.o sys.o
obc.o: alloc.o precision.o
ofc.o: alloc.o files.o precision.o units.o
old_atmfuncs.o: alloc.o atmparams.o precision.o sys.o
on_subs.o: alloc.o onmod.o onmod.o
onmod.o: precision.o
optical.o: alloc.o atomlist.o densematrix.o fermid.o files.o parallel.o
optical.o: parallelsubs.o precision.o sys.o units.o
ordern.o: alloc.o cgwf.o chempot.o globalise.o m_mpi_utils.o onmod.o parallel.o
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
phirphi.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
phirphi_opt.o: alloc.o atm_types.o atmfuncs.o atmparams.o mneighb.o parallel.o
phirphi_opt.o: parallelsubs.o precision.o spher_harm.o sys.o
phonon.o: alloc.o files.o parallel.o precision.o sys.o
pixmol.o: files.o periodic_table.o precision.o
plcharge.o: files.o old_atmfuncs.o parallel.o precision.o
poison.o: alloc.o cellsubs.o chkgmx.o fft.o parallel.o precision.o sys.o
post_scf_work.o: atomlist.o compute_dm.o final_H_f_stress.o kpoint_grid.o
post_scf_work.o: m_dminim.o m_energies.o m_eo.o m_gamma.o m_spin.o m_steps.o
post_scf_work.o: m_zminim.o mneighb.o parallel.o siesta_geom.o siesta_options.o
post_scf_work.o: sparse_matrices.o
printmatrix.o: alloc.o
projected_DOS.o: alloc.o atomlist.o kpoint_grid.o kpoint_pdos.o m_eo.o
projected_DOS.o: m_gamma.o m_spin.o parallel.o precision.o siesta_geom.o
projected_DOS.o: siesta_options.o sparse_matrices.o sys.o
propor.o: precision.o sys.o
proximity_check.o: m_ts_options.o mneighb.o parallel.o precision.o
proximity_check.o: siesta_geom.o siesta_options.o
pseudopotential.o: alloc.o atom_options.o interpolation.o precision.o sys.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
rdiag.o: alloc.o diagmemory.o m_diagon.o parallel.o precision.o sys.o
read_xc_info.o: parallel.o precision.o sys.o
readsp.o: parallel.o precision.o sys.o
redcel.o: parallel.o precision.o sys.o
reinit.o: files.o parallel.o siesta_master.o
remove_intramol_pressure.o: precision.o zmatrix.o
reoptical.o: parallel.o precision.o sys.o
reord.o: alloc.o precision.o
rhofft.o: cellsubs.o fft.o precision.o sys.o
rhooda.o: atmfuncs.o atomlist.o mesh.o meshphi.o precision.o
rhoofd.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
rhoofd.o: meshphi.o parallel.o parallelsubs.o precision.o sys.o
rhoofdsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o
rhoofdsp.o: meshdscf.o meshphi.o precision.o sys.o
save_density_matrix.o: atomlist.o iodm_netcdf.o m_energies.o m_iodm.o m_spin.o
save_density_matrix.o: m_steps.o m_ts_global_vars.o m_ts_iodm.o precision.o
save_density_matrix.o: siesta_options.o sparse_matrices.o
savepsi.o: alloc.o parallel.o parallelsubs.o precision.o
scfconvergence_test.o: m_convergence.o m_energies.o m_wallclock.o parallel.o
scfconvergence_test.o: siesta_cml.o siesta_options.o units.o write_subs.o
schecomm.o: alloc.o
setatomnodes.o: alloc.o parallel.o precision.o spatial.o sys.o
setspatial.o: alloc.o parallel.o precision.o spatial.o
setup_H0.o: alloc.o atmfuncs.o atomlist.o dhscf.o dnaefs.o kinefsm.o
setup_H0.o: m_energies.o m_ntm.o m_spin.o metaforce.o molecularmechanics.o
setup_H0.o: naefs.o nlefsm.o siesta_geom.o siesta_options.o sparse_matrices.o
setup_hamiltonian.o: alloc.o atmfuncs.o atomlist.o dhscf.o files.o m_dipol.o
setup_hamiltonian.o: m_energies.o m_gamma.o m_hsx.o m_mpi_utils.o m_ntm.o
setup_hamiltonian.o: m_partial_charges.o m_rhog.o m_spin.o m_steps.o m_stress.o
setup_hamiltonian.o: m_ts_global_vars.o metaforce.o molecularmechanics.o
setup_hamiltonian.o: parallel.o siesta_geom.o siesta_options.o
setup_hamiltonian.o: sparse_matrices.o sys.o
setup_kscell.o: minvec.o parallel.o precision.o units.o
setup_ordern_indexes.o: alloc.o domain_decom.o parallel.o spatial.o
shaper.o: atmfuncs.o mneighb.o precision.o
show_distribution.o: atomlist.o parallel.o parallelsubs.o siesta_geom.o sys.o
siesta.o: alloc.o m_steps.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o
siesta_analysis.o: alloc.o atomlist.o bands.o basis_enthalpy.o dhscf.o files.o
siesta_analysis.o: kpoint_grid.o ksv.o ksvinit.o local_DOS.o m_dipol.o
siesta_analysis.o: m_energies.o m_eo.o m_forces.o m_gamma.o m_mpi_utils.o
siesta_analysis.o: m_ntm.o m_partial_charges.o m_spin.o m_steps.o parallel.o
siesta_analysis.o: projected_DOS.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_analysis.o: sparse_matrices.o write_subs.o writewave.o zmatrix.o
siesta_cmlsubs.o: files.o parallel.o siesta_cml.o timestamp.o 
siesta_end.o: bands.o densematrix.o fixed.o m_dscfcomm.o m_pulay.o m_rhog.o
siesta_end.o: m_wallclock.o memory_log.o meshdscf.o meshphi.o moremeshsubs.o
siesta_end.o: parallel.o siesta_cml.o siesta_cmlsubs.o siesta_master.o
siesta_end.o: sparse_matrices.o timestamp.o
siesta_forces.o: compute_dm.o compute_energies.o compute_max_diff.o
siesta_forces.o: final_H_f_stress.o m_convergence.o m_energies.o m_forces.o
siesta_forces.o: m_pulay.o m_rhog.o m_steps.o m_stress.o m_ts_global_vars.o
siesta_forces.o: m_ts_options.o mixer.o parallel.o post_scf_work.o
siesta_forces.o: save_density_matrix.o scfconvergence_test.o setup_H0.o
siesta_forces.o: setup_hamiltonian.o siesta_cml.o siesta_geom.o siesta_master.o
siesta_forces.o: siesta_options.o sparse_matrices.o state_analysis.o
siesta_forces.o: state_init.o sys.o units.o write_subs.o
siesta_geom.o: precision.o
siesta_init.o: alloc.o atomlist.o bands.o bsc_xcmod.o densematrix.o files.o
siesta_init.o: ioxv.o kpoint_grid.o kpoint_pdos.o ksvinit.o m_check_supercell.o
siesta_init.o: m_energies.o m_eo.o m_forces.o m_gamma.o m_iostruct.o
siesta_init.o: m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o m_timer.o
siesta_init.o: m_wallclock.o memory_log.o metaforce.o molecularmechanics.o
siesta_init.o: parallel.o parallelsubs.o phonon.o projected_DOS.o
siesta_init.o: siesta_cmlsubs.o siesta_geom.o siesta_options.o
siesta_init.o: sparse_matrices.o struct_init.o sys.o timer.o timestamp.o
siesta_init.o: ts_init.o writewave.o zmatrix.o
siesta_master.o: iopipes.o precision.o sys.o
siesta_move.o: atomlist.o broyden_optim.o cell_broyden_optim.o
siesta_move.o: cell_fire_optim.o dynamics.o fire_optim.o ioxv.o m_energies.o
siesta_move.o: m_forces.o m_kinetic.o m_steps.o m_stress.o parallel.o phonon.o
siesta_move.o: siesta_cml.o siesta_geom.o siesta_master.o siesta_options.o
siesta_move.o: sys.o write_subs.o zm_broyden_optim.o zm_fire_optim.o zmatrix.o
siesta_options.o: diagmemory.o files.o parallel.o precision.o siesta_cml.o
siesta_options.o: sys.o units.o
sparse_matrices.o: alloc.o precision.o
spatial.o: precision.o
spher_harm.o: alloc.o precision.o sys.o
state_analysis.o: atomlist.o born_charge.o fixed.o m_energies.o m_forces.o
state_analysis.o: m_iostruct.o m_ntm.o m_spin.o m_stress.o m_wallclock.o
state_analysis.o: parallel.o remove_intramol_pressure.o siesta_cml.o
state_analysis.o: siesta_geom.o siesta_options.o sparse_matrices.o write_subs.o
state_analysis.o: zmatrix.o
state_init.o: alloc.o atomlist.o domain_decom.o files.o hsparse.o iodm_netcdf.o
state_init.o: iodmhs_netcdf.o kpoint_grid.o m_check_supercell.o m_energies.o
state_init.o: m_eo.o m_gamma.o m_mpi_utils.o m_pulay.o m_rmaxh.o m_spin.o
state_init.o: m_steps.o m_ts_io.o m_ts_options.o new_dm.o overlap.o parallel.o
state_init.o: proximity_check.o siesta_cml.o siesta_geom.o siesta_options.o
state_init.o: sparse_matrices.o sys.o write_subs.o zmatrix.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o ioxv.o m_iostruct.o
struct_init.o: m_mpi_utils.o parallel.o periodic_table.o siesta_cml.o
struct_init.o: siesta_geom.o siesta_master.o siesta_options.o zmatrix.o
sys.o: parallel.o
timer.o: m_timer.o parallel.o sys.o timer_tree.o
timer_tree.o: m_walltime.o
transition_rate.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
ts_init.o: files.o m_ts_contour.o m_ts_gf.o m_ts_kpoints.o m_ts_options.o
ts_init.o: parallel.o
ts_show_regions.o: parallel.o precision.o units.o
typecell.o: precision.o
uncell.o: precision.o
units.o: precision.o
vacuum_level.o: moreParallelSubs.o precision.o
vmat.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
vmat.o: meshphi.o parallel.o parallelsubs.o precision.o
vmatsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
vmatsp.o: meshphi.o precision.o
vmb.o: fixed.o parallel.o precision.o sys.o
write_md_record.o: atomlist.o m_energies.o m_kinetic.o m_steps.o m_stress.o
write_md_record.o: md_out.o parallel.o siesta_geom.o siesta_options.o
write_orb_indx.o: atmfuncs.o cellsubs.o files.o precision.o
write_raw_efs.o: atmfuncs.o atomlist.o precision.o siesta_geom.o
write_subs.o: alloc.o atomlist.o m_energies.o m_forces.o m_iostruct.o m_spin.o
write_subs.o: m_steps.o m_stress.o m_ts_global_vars.o parallel.o precision.o
write_subs.o: siesta_cml.o siesta_geom.o siesta_options.o units.o zmatrix.o
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o kpoint_grid.o
writewave.o: parallel.o parallelsubs.o precision.o siesta_geom.o sys.o units.o
xc.o: alloc.o bsc_xcmod.o precision.o sys.o
xml.o: precision.o
zm_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
zm_broyden_optim.o: units.o zmatrix.o
zm_fire_optim.o: m_fdf_global.o m_fire.o parallel.o precision.o
zm_fire_optim.o: siesta_options.o units.o zmatrix.o
zmatrix.o: alloc.o m_cell.o parallel.o precision.o siesta_geom.o sys.o units.o
band.o: bands.o
fsiesta.o: fsiesta_pipes.o
listsc_module.o: listsc.o
m_bessph.o: bessph.o
m_born_charge.o: born_charge.o
m_broyden_optim.o: broyden_optim.o
m_cell_broyden_optim.o: cell_broyden_optim.o
m_cell_fire_optim.o: cell_fire_optim.o
m_cgwf.o: cgwf.o
m_chempot.o: chempot.o
m_chkgmx.o: chkgmx.o
m_compute_dm.o: compute_dm.o
m_compute_energies.o: compute_energies.o
m_compute_max_diff.o: compute_max_diff.o
m_compute_rhog.o: compute_rhog.o
m_conjgr.o: conjgr.o
m_conjgr_old.o: conjgr_old.o
m_denmat.o: denmat.o
m_dhscf.o: dhscf.o
m_dnaefs.o: dnaefs.o
m_doping_uniform.o: doping_uniform.o
m_dynamics.o: dynamics.o
m_eandg.o: egandd.o
m_ener3.o: ener3.o
m_errorf.o: errorf.o
m_fermid.o: fermid.o
m_fft.o: fft.o
m_final_h_f_stress.o: final_H_f_stress.o
m_find_kgrid.o: find_kgrid.o
m_fire_optim.o: fire_optim.o
m_fixed.o: fixed.o
m_forhar.o: forhar.o
m_gradient.o: gradient.o
m_grdsam.o: grdsam.o
m_hsparse.o: hsparse.o
m_intramol_pressure.o: remove_intramol_pressure.o
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
m_norm.o: compute_norm.o
m_normalize_dm.o: normalize_dm.o
m_on_subs.o: on_subs.o
m_ordern.o: ordern.o
m_overfsm.o: overfsm.o
m_overlap.o: overlap.o
m_post_scf_work.o: post_scf_work.o
m_projected_dos.o: projected_DOS.o
m_proximity_check.o: proximity_check.o
m_radfft.o: radfft.o
m_rhofft.o: rhofft.o
m_save_density_matrix.o: save_density_matrix.o
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
m_timer_tree.o: timer_tree.o
m_timestamp.o: timestamp.o
m_ts_init.o: ts_init.o
m_vmat.o: vmat.o
m_xyz_phiatm.o: matel.o
m_zm_broyden_optim.o: zm_broyden_optim.o
m_zm_fire_optim.o: zm_fire_optim.o
neighbour.o: mneighb.o
on_core.o: onmod.o
on_main.o: onmod.o
printmat.o: printmatrix.o
timer_options.o: timer.o
version_info.o: version.o
write_subs_energies.o: write_subs.o
write_subs_positions.o: write_subs.o
write_subs_pressure.o: write_subs.o
