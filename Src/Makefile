# 
# Siesta/Src/Makefile
# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
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
default: what siesta  # fsiesta_pipes.o fsiesta_sockets.o libSiestaForces.a
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
INCFLAGS:= $(NETCDF_INCFLAGS) $(INCFLAGS)
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
SYSOBJ=$(SYS).o pxf.o sockets.o fsockets.o # includes C socket library and F90 interface
#
# Note that machine-specific files are now in top Src directory.
#
MAIN_OBJ = siesta.o

OBJS =  automatic_cell.o atom_options.o \
        arw.o  atomlwf.o bands.o basis_enthalpy.o bessph.o bonds.o \
        born_charge.o cellxc_mod.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o coxmol.o cross.o compute_norm.o\
	denmat.o denmatlomem.o detover.o dfscf.o m_diagon_opt.o diagon.o digcel.o \
        fft.o dhscf.o constr.o diagk_file.o \
	diagg.o diagk.o diagkp.o diag2g.o diag2k.o diag3g.o diag3k.o diagpol.o \
	diagsprl.o dipole.o dismin.o dnaefs.o doping_uniform.o dot.o  \
	m_efield.o egandd.o ener3.o ener3lomem.o errorf.o extrapolon.o \
	fixed.o interpolation.o gradient.o gradientlomem.o grdsam.o \
	hsparse.o idiag.o  initatom.o inver.o iodm_netcdf.o \
	iodmhs_netcdf.o iogrid_netcdf.o iolwf.o iozm.o ipack.o iopipes.o iosockets.o \
	iowfs_netcdf.o kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o local_DOS.o\
	madelung.o mesh.o meshphi.o meshdscf.o memory.o meshsubs.o metaforce.o \
        m_check_supercell.o m_fft_gpfa.o mulliken.o minvec.o m_pulay.o naefs.o \
	new_dm.o normalize_dm.o ordern.o outcell.o outcoor.o \
	paste.o pdos.o pdosg.o pdosk.o pdoskp.o phirphi.o pixmol.o plcharge.o \
	projected_DOS.o propor.o m_partial_charges.o \
	randomg.o reclat.o redcel.o \
	reinit.o reord.o reordpsi.o rhoofd.o rhoofdsp.o rhooda.o \
	savepsi.o shaper.o timer_tree.o timer.o \
	vmb.o vmat.o vmatsp.o volcel.o \
	cdiag.o rdiag.o \
        cgvc.o cgvc_zmatrix.o m_convergence.o \
        iocg.o ioeig.o iofa.o iokp.o iomd.o kpoint_pdos.o typecell.o \
        ofc.o poison.o readsp.o radfft.o \
	write_md_record.o kpoint_grid.o find_kgrid.o proximity_check.o\
	state_init.o siesta_move.o setup_hamiltonian.o compute_dm.o mixer.o\
	scfconvergence_test.o post_scf_work.o state_analysis.o write_subs.o \
	siesta_init.o struct_init.o siesta_options.o read_options.o siesta_geom.o \
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
        atmfuncs.o listsc.o memoryinfo.o m_memory.o sorting.o \
        atomlist.o atm_types.o old_atmfuncs.o radial.o \
        alloc.o spher_harm.o periodic_table.o version.o timestamp.o \
        basis_types.o xml.o pseudopotential.o basis_specs.o basis_io.o\
	ldau.o ldau_specs.o \
        onmod.o densematrix.o writewave.o on_subs.o fermid.o m_broyddj.o \
        electrostatic.o mneighb.o globalise.o siesta_cmlsubs.o siesta_cml.o \
        units.o zmatrix.o m_broyden_mixing.o forhar.o m_walltime.o m_wallclock.o m_iostruct.o \
        nlefsm.o overfsm.o overlap.o conjgr.o conjgr_old.o m_energies.o m_steps.o\
        m_broyddj_nocomm.o broyden_optim.o ioxv.o dynamics.o md_out.o \
        molecularmechanics.o zm_broyden_optim.o \
        cell_broyden_optim.o remove_intramol_pressure.o \
        m_ntm.o m_dipol.o  m_kinetic.o  m_rmaxh.o \
        m_forces.o m_stress.o m_eo.o m_spin.o m_gamma.o m_hsx.o m_fire.o setup_H0.o \
        get_kpoints_scale.o get_target_stress.o fire_optim.o zm_fire_optim.o \
        cell_fire_optim.o m_fire_para.o m_fire_mixing.o write_raw_efs.o \
        pdos2g.o pdos2k.o pdos3g.o pdos3k.o spinorbit.o moments.o \
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
        save_density_matrix.o \
        m_rhog.o rhofft.o m_diis.o compute_energies.o \
        class_Vector.o class_Pair_Vectors.o class_Fstack_Pair_Vectors.o\
        m_target_stress.o\
	matel_registry.o register_rfs.o new_matel.o\
	m_trialorbitalclass.o siesta2wannier90.o m_planewavematrixvar.o \
	delk.o m_planewavematrix.o m_digest_nnkp.o broadcast_projections.o \
	compute_pw_matrix.o m_writedelk.o mmn.o m_noccbands.o amn.o \
	m_overkkneig.o write_inp_wannier.o diagonalizeHk.o m_orderbands.o

# Add citation package
OBJS += m_cite.o

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
INCFLAGS:=$(INCFLAGS) $(FoX_FCFLAGS) $(NETCDF_INCFLAGS)
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
XC=libSiestaXC.a
$(XC): 
	(cd SiestaXC ; $(MAKE) \
                    "VPATH_ROOT=$(VPATH)"  \
	            "INCFLAGS=-I../MPI" \
		    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module )
#--------------------------------------------------------------
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "VPATH=$(VPATH)/fdf" \
		    "FPPFLAGS=$(FPPFLAGS)" \
		    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module )
#
# Routines using fdf calls. Woefully incomplete, but the $(FDF) target is made first.
#
dhscf.o iolwf.o  grdsam.o : $(FDF)
recoor.o outcoor.o ioxv.o m_iostruct.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o pdos.o plcharge.o: $(FDF)
coor.o cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o  ofc.o : $(FDF)
readsp.o eggbox.o iopipes.o : $(FDF)
#
atom.o basis_io.o basis_specs.o: $(FDF)
ldau_specs.o: $(FDF)
#
# C/Fortran sockets library 
#
# As it stands, it will use the CC in the default rules
# Define CC (and maybe CFLAGS) in arch.make for full control
#
sockets.o:  sockets.c   # Will be picked-up in the VPATH
	$(CC) -c $(CFLAGS) $<
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
	(cd SiestaXC ; $(MAKE) clean )

# DO NOT DELETE THIS LINE - used by make depend
alloc.o: debugmpi.o m_io.o parallel.o precision.o sys.o
amn.o: alloc.o atmfuncs.o atomlist.o m_orderbands.o m_trialorbitalclass.o
amn.o: matel_registry.o mneighb.o new_matel.o parallel.o parallelsubs.o
amn.o: precision.o siesta2wannier90.o siesta_geom.o sys.o
arw.o: alloc.o parallel.o precision.o sys.o
atm_transfer.o: atm_types.o atmparams.o ldau_specs.o old_atmfuncs.o
atm_transfer.o: periodic_table.o radial.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o bsc_xcmod.o
atom.o: interpolation.o m_filter.o old_atmfuncs.o periodic_table.o precision.o
atom.o: pseudopotential.o sys.o
atom_options.o: sys.o
atomlist.o: alloc.o atm_types.o atmfuncs.o ldau_specs.o parallel.o precision.o
atomlist.o: radial.o siesta_geom.o
