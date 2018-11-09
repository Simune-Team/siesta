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
# Makefile for SIESTA
#
# Remove implicit rules
.SUFFIXES:

VPATH=.

# Default compiler variables
AR ?= ar
RANLIB ?= ranlib

default: siesta
include arch.make

# Ensure the correct implicit rules for the compilation
# SIESTA relies on these extensions and thus we
# force them to be here.
# Importantly these should be after the include statement.
# Note that this will not help the utilities.
.SUFFIXES: .f .f90 .F .F90 .c .o .a


lib: libSiestaXC.a libfdf.a libSiestaForces.a

dist:
	make -f makefile.lit dist
dep:
	-if [ $(VPATH) = "." ] ; then \
	sfmakedepend --depend=obj --modext=o \
                              *.f *.f90 *.F *.F90 *.T90 ; \
        else \
	sfmakedepend --depend=obj --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 ; \
        fi
	@sed -i -e "/^siesta_cmlsubs.o:/s/version.o//" Makefile
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

INCFLAGS:= $(NETCDF_INCFLAGS) $(PEXSI_INCFLAGS) $(INCFLAGS)


# For removing any -ipo compilation in sub directories...
# Later compilers are also having this enabled, and it lets you
# specify what to not be used in libraries.
# This is so because inter-procedural compilation should NOT be performed in
# libraries (it removes names when a routine can be fully moved into others)
IPO_FLAG?= -ipo


.PHONY: what version

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


SYSOBJ=$(SYS).o pxf.o sockets.o fsockets.o # includes C socket library and F90 interface

# Note that machine-specific files are now in top Src directory.
MAIN_OBJ = siesta.o


OBJS =  automatic_cell.o atom_options.o \
        arw.o  atomlwf.o bands.o basis_enthalpy.o bessph.o bonds.o \
        born_charge.o cellxc_mod.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o coxmol.o cross.o compute_norm.o\
	denmat.o denmatlomem.o detover.o dfscf.o diagon.o digcel.o \
        fft.o dhscf.o constr.o diagk_file.o \
	diagg.o diagk.o diagkp.o diag2g.o diag2k.o diag3g.o diag3k.o diagpol.o \
	diagsprl.o dipole.o dismin.o dnaefs.o doping_uniform.o dot.o  \
	m_efield.o egandd.o ener3.o ener3lomem.o errorf.o extrapolon.o \
	m_fixed.o interpolation.o gradient.o gradientlomem.o grdsam.o \
	hsparse.o idiag.o  initatom.o inver.o iodm_netcdf.o \
	iodmhs_netcdf.o iogrid_netcdf.o iolwf.o iozm.o ipack.o iopipes.o iosockets.o \
	iowfs_netcdf.o kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o local_DOS.o\
	madelung.o mesh.o meshphi.o meshdscf.o memory.o meshsubs.o metaforce.o \
        m_supercell.o mulliken.o minvec.o naefs.o \
	m_new_dm.o normalize_dm.o ordern.o outcell.o outcoor.o \
	m_fft_gpfa.o \
	pdos.o pdosg.o pdosk.o pdoskp.o phirphi.o pixmol.o plcharge.o \
	projected_DOS.o propor.o m_partial_charges.o \
	randomg.o reclat.o redcel.o \
	reinit.o reord.o reordpsi.o rhoofd.o rhoofdsp.o rhooda.o \
	savepsi.o shaper.o timer_tree.o timer.o \
	vmb.o vmat.o vmatsp.o volcel.o \
        cgvc.o cgvc_zmatrix.o m_convergence.o \
        iocg.o ioeig.o iofa.o iokp.o iomd.o typecell.o \
        ofc.o poison.o readsp.o radfft.o \
	write_md_record.o find_kgrid.o proximity_check.o\
	state_init.o siesta_move.o setup_hamiltonian.o compute_dm.o mixer.o\
	scfconvergence_test.o post_scf_work.o state_analysis.o write_subs.o \
	siesta_init.o struct_init.o siesta_options.o read_options.o siesta_geom.o \
        siesta_analysis.o siesta_end.o siesta_forces.o \
	siesta_tddft.o \
        io.o sparse_matrices.o \
        coor.o \
        broadcast_basis.o eggbox.o \
        optical.o phirphi_opt.o reoptical.o transition_rate.o \
        initparallel.o show_distribution.o setspatial.o setatomnodes.o uncell.o \
        cart2frac.o obc.o \
        precision.o sys.o m_cell.o files.o \
        spatial.o parallel.o parallelsubs.o  parsing.o chemical.o atom.o atmparams.o \
        m_mpi_utils.o fdf_extra.o m_iorho.o  \
        atmfuncs.o listsc.o memoryinfo.o m_memory.o sorting.o \
        atomlist.o atm_types.o radial.o \
        alloc.o spher_harm.o periodic_table.o version.o timestamp.o \
        basis_types.o xml.o pseudopotential.o basis_specs.o basis_io.o\
	ldau.o ldau_specs.o \
        onmod.o densematrix.o writewave.o on_subs.o fermid.o m_broyddj.o \
        electrostatic.o mneighb.o globalise.o siesta_cmlsubs.o siesta_cml.o \
        units.o zmatrix.o m_broyden_mixing.o forhar.o m_walltime.o m_wallclock.o m_iostruct.o \
        nlefsm.o overfsm.o overlap.o overlap_gradient.o conjgr.o conjgr_old.o m_energies.o m_steps.o\
        m_broyddj_nocomm.o broyden_optim.o ioxv.o dynamics.o md_out.o \
        molecularmechanics.o zm_broyden_optim.o \
        cell_broyden_optim.o remove_intramol_pressure.o \
        m_ntm.o m_dipol.o m_kinetic.o m_rmaxh.o \
        m_forces.o m_stress.o m_eo.o m_spin.o m_hsx.o m_fire.o setup_H0.o \
        get_kpoints_scale.o get_target_stress.o fire_optim.o zm_fire_optim.o \
        cell_fire_optim.o m_fire_para.o m_fire_mixing.o write_raw_efs.o \
        pdos2g.o pdos2k.o pdos3g.o pdos3k.o spinorbit.o moments.o \
	fsiesta_mpi.o \
        final_H_f_stress.o \
        debugmpi.o qsort.o save_density_matrix.o \
        m_dscfcomm.o schecomm.o moremeshsubs.o \
        domain_decom.o printmatrix.o mmio.o pspltm1.o lenstr.o \
        setup_ordern_indexes.o\
        m_iodm.o m_iodm_old.o \
	m_filter.o \
        m_io.o m_io_yaml.o \
	m_timer.o \
	extrae_module.o extrae_eventllist.o\
	moreParallelSubs.o \
	read_xc_info.o \
	siesta_master.o \
        bsc_xcmod.o bsc_cellxc.o xc.o \
	vacuum_level.o \
	write_orb_indx.o \
        die.o m_pexsi.o m_pexsi_driver.o m_pexsi_dos.o m_pexsi_local_dos.o\
		m_chess.o \
        m_redist_spmatrix.o \
        class_Distribution.o\
        m_dminim.o m_zminim.o \
        m_getopts.o f2kcli.o \
        m_svd.o m_matio.o\
        rusage.o memory_snapshot.o compute_ebs_shift.o memory_all.o \
	matel_registry.o register_rfs.o new_matel.o \
	kpoint_convert.o \
	m_target_stress.o compute_max_diff.o \
	m_initwf.o wavefunctions.o \
	m_evolve.o cranknic_evolk.o cranknic_evolg.o sankey_change_basis.o m_iotddft.o \
	m_matdiag.o m_matswinvers.o iotdxv.o  m_inversemm.o

OBJS += kpoint_t.o kpoint_scf.o kpoint_dos.o kpoint_pdos.o kpoint_ldos.o

OBJS += init_output.o

# Diagonalization algorithm
OBJS += diag_option.o diag.o

# Mixing algorithms
OBJS += m_mixing.o m_mixing_scf.o

# Wannier90 interface
OBJS += m_trialorbitalclass.o siesta2wannier90.o m_planewavematrixvar.o \
	delk.o m_planewavematrix.o m_digest_nnkp.o broadcast_projections.o \
	compute_pw_matrix.o m_writedelk.o mmn.o m_noccbands.o amn.o \
	m_overkkneig.o write_inp_wannier.o diagonalizeHk.o m_orderbands.o\
        m_check_walltime.o print_spin.o

# class Objects
OBJS += class_OrbitalDistribution.o \
	class_Sparsity.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
        restructSpData2D.o extrapolateSpData2D.o \
        class_Geometry.o \
	class_Fstack_Data1D.o \
	class_Pair_Data1D.o \
	class_Fstack_Pair_Data1D.o \
	class_Pair_Geometry_SpData2D.o \
	class_Fstack_Pair_Geometry_SpData2D.o \
	class_Pair_SpData1D.o \
	class_Fstack_Pair_SpData1D.o \
	class_TriMat.o m_trimat_invert.o \
        m_uuid.o object_debug.o

# As the classes are based on inputs we force the dependency to alloc
class_Sparsity.o class_Data1D.o class_Data2D.o class_TriMat.o: alloc.o

OBJS += m_rhog.o rhofft.o m_diis.o compute_energies.o

# Generic modules
OBJS += m_char.o m_os.o
OBJS += intrinsic_missing.o geom_helper.o m_sparse.o
OBJS += m_handle_sparse.o m_mesh_node.o
OBJS += create_Sparsity_SC.o create_Sparsity_Union.o
OBJS += m_gauss_quad.o
OBJS += m_gauss_fermi_inf.o m_gauss_fermi_30.o m_gauss_fermi_28.o
OBJS += m_gauss_fermi_26.o m_gauss_fermi_24.o m_gauss_fermi_22.o
OBJS += m_gauss_fermi_20.o m_gauss_fermi_19.o m_gauss_fermi_18.o
OBJS += m_gauss_fermi_17.o
OBJS += m_integrate.o m_interpolate.o m_mat_invert.o m_monitor.o
OBJS += m_iterator.o m_pivot_array.o

OBJS += dm_charge.o

# Pivot (of sparsity pattern) routines
OBJS += m_pivot.o m_pivot_methods.o

# Super-cell calculations/information modules
OBJS += atom_graph.o

# Geometry modules
OBJS += m_geom_aux.o m_geom_objects.o
OBJS += m_geom_box.o m_geom_coord.o m_geom_square.o m_geom_plane.o

# Charge and Hartree Geometries
OBJS += m_charge_add.o m_hartree_add.o

# New io-routines
OBJS += m_io_s.o
ifneq (,$(findstring -DTEST_IO,$(FPPFLAGS)))
OBJS += m_test_io.o m_matio.o
endif

# Region objects
OBJS += m_region.o m_sparsity_handling.o

# IO-modules
OBJS += m_ncdf_siesta.o m_ncdf_io.o m_exp_coord.o

# Lua interface
OBJS += flook_siesta.o siesta_dicts.o

# Add citation package
OBJS += m_cite.o

# Special dependency that sfmakedepend cannot grasp for now
m_pexsi.o m_pexsi_driver.o m_pexsi_dos.o m_pexsi_local_dos.o: f_ppexsi_interface.o

#--------------------------------------------------------------
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
#
# These subordinate make processes take options from the top level.
# In some cases they are superfluous, since they can be read
# from the arch.make file, which is 'included' by all subordinate
# makefiles. But in some other cases they are needed: VPATH, for
# example, needs to be passed always.
#
#--------------------------------------------------------------
XMLPARSER=libxmlparser.a
$(XMLPARSER): 
	(cd xmlparser ; $(MAKE) -j 1 "VPATH=$(VPATH)/xmlparser" \
	FFFLAGS="$(FFLAGS:$(IPO_FLAG)=)" module)

WXML=libwxml.a
$(WXML):
	(cd wxml ; $(MAKE) -j 1 "VPATH=$(VPATH)/wxml" \
	FFFLAGS="$(FFLAGS:$(IPO_FLAG)=)" module)

#--------------------------------------------------------------
libmpi_f90.a: 
	@(cd MPI ; $(MAKE) -j 1 \
                    "VPATH=$(VPATH)/MPI" \
                    "FPPFLAGS=$(FPPFLAGS)" \
                    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" )

#--------------------------------------------------------------
FDF=libfdf.a
$(FDF): $(MPI_INTERFACE)
	(cd fdf ; $(MAKE) -j 1 "VPATH=$(VPATH)/fdf" \
		    "FPPFLAGS=$(FPPFLAGS)" \
		    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module )

#--------------------------------------------------------------
MS=MatrixSwitch.a
$(MS): $(MPI_INTERFACE)
	(cd MatrixSwitch/src ; \
	 $(MAKE) -j 1 \
	            "VPATH=$(VPATH)/MatrixSwitch/src" \
		    "FPPFLAGS=$(FPPFLAGS)" \
		    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module)
#
m_initwf.o sankey_change_basis.o m_evolve.o cranknic_evolg.o cranknic_evolk.o: $(MS)
wavefunctions.o m_matdiag.o m_matswinvers.o m_inversemm.o: $(MS)
#--------------------------------------------------------------
XC=libSiestaXC.a
$(XC): $(FDF)
	(cd SiestaXC ; $(MAKE) -j 1 \
                    "VPATH_ROOT=$(VPATH)"  \
	            "INCFLAGS=-I../MPI" \
		    "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module )

#######################
#  TS-TBT libs:
#
# Libraries used by transiesta and tbtrans
# They are currently supplied by this distribution
# All rights reserved for the author of these libraries.
# Author: Nick R. Papior

# This library can be found at:
#  https://github.com/zerothi/fdict
FDICT=libfdict.a
$(FDICT): 
	(cd fdict ; rm -rf obj ; mkdir obj ; cd obj ; \
	echo "TOP_DIR=$(VPATH)/fdict" > Makefile ; echo "include $(VPATH)/fdict/Makefile" >> Makefile ; \
	echo "FC=$(FC_SERIAL)" > setup.make ; \
	echo "AR=$(AR)" >> setup.make ; \
	echo "RANLIB=$(RANLIB)" >> setup.make ; \
	echo "FFLAGS=$(FFLAGS)" >> setup.make ; \
	echo "LINK_FLAGS=$(FFLAGS)" >> setup.make ; \
	$(MAKE) -j 1 copy ; $(MAKE) -j 1 lib SHARED=0)
	@cp -p fdict/obj/*.mod ./
	cp fdict/obj/$(FDICT) ./

# This library can be found at:
#  https://github.com/zerothi/ncdf
NCDF=libncdf.a
$(NCDF): $(FDICT)
	(cd ncdf ; rm -rf obj ; mkdir obj ; cd obj ; \
	echo "TOP_DIR=$(VPATH)/ncdf" > Makefile ; echo "include $(VPATH)/ncdf/Makefile" >> Makefile ; \
	echo "settings.bash: FORCE" >> Makefile ; \
	printf '%s\n' "	-cp ../../fdict/obj/settings.bash ." >> Makefile ; \
	echo "FC=$(FC)" > setup.make ; \
	echo "AR=$(AR)" >> setup.make ; \
	echo "RANLIB=$(RANLIB)" >> setup.make ; \
	echo "FFLAGS=$(FFLAGS)" >> setup.make ; \
	echo "LINK_FLAGS=$(FFLAGS)" >> setup.make ; \
	echo "FPPFLAGS = $(FPPFLAGS)" >> setup.make ; \
	echo "INCLUDES = $(INCFLAGS) -I../../" >> setup.make ; \
	echo "LIBS = $(LIBS) -lfdict" >> setup.make ; \
	$(MAKE) -j 1 copy ; $(MAKE) -j 1 lib SHARED=0)
	@cp -p ncdf/obj/*.mod ./
	cp ncdf/obj/$(NCDF) ./

#######################


# Routines using fdf calls.
# Woefully incomplete, but the $(FDF) target is made first.
dhscf.o iolwf.o grdsam.o : $(FDF)
recoor.o outcoor.o ioxv.o m_iostruct.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o pdos.o plcharge.o: $(FDF)
coor.o cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o  ofc.o : $(FDF)
readsp.o eggbox.o iopipes.o : $(FDF)
atom.o basis_io.o basis_specs.o: $(FDF)
ldau_specs.o: $(FDF)

cranknic_evolk.o: $(FDF)
cranknic_evolg.o: $(FDF)
sankey_change_basis.o iotdxv.o: $(FDF)

#
# Linear algebra routines
#
libsiestaBLAS.a: Libs/blas.F
	@echo "==> Compiling libsiestaBLAS.a in Libs..."
	@(cd Libs ; $(MAKE) -j 1 "VPATH=$(VPATH)/Libs" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libsiestaBLAS.a)
	@cp Libs/libsiestaBLAS.a .

ifneq (,$(findstring libsiestaLAPACK.a,$(COMP_LIBS)))
 # Automatically add the DIAG_2STAGE solvers.
 FPPFLAGS += -DSIESTA__DIAG_2STAGE
 ifeq (,$(findstring MPI,$(FPPFLAGS)))
  FPPFLAGS += -DSIESTA__MRRR
 endif
endif
libsiestaLAPACK.a: Libs/lapack.F
	@echo "==> Compiling libsiestaLAPACK.a in Libs..."
	@(cd Libs ; $(MAKE) -j 1 "VPATH=$(VPATH)/Libs" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libsiestaLAPACK.a)
	@cp Libs/libsiestaLAPACK.a .

libblas.a:
	@echo ""
	@echo "SIESTA: $@ has been deprecated in favor of:"
	@echo "  COMP_LIBS += libsiestaBLAS.a"
	@echo "Please update your arch.make file accordingly..."
	@echo ""
	@false
liblapack.a dc_lapack.a zgesv_lapack.a essl_lapack.a:
	@echo ""
	@echo "SIESTA: $@ has been deprecated in favor of:"
	@echo "  COMP_LIBS += libsiestaLAPACK.a"
	@echo "Please update your arch.make file accordingly..."
	@echo ""
	@false
linalg.a:
	@echo ""
	@echo "SIESTA: $@ has been deprecated in favor of:"
	@echo "  COMP_LIBS += libsiestaLAPACK.a libsiestaBLAS.a"
	@echo "Please update your arch.make file accordingly..."
	@echo ""
	@false


# Extract information of the compiler version used.
# This is very basic but for intel, gcc, pgi and open64
# this at least contain the version string of the compiler.
# Possibly later the used compiler version may be extracted.
COMPILER_VERSION := $(shell $(FC) --version | sed '/^[ ]*$$/d;q')
ifeq ($(COMPILER_VERSION),)
  # Just try and capture the case were it can not figure out anything
  COMPILER_VERSION := Unknown version
endif

.PHONY: version version.o
version.o: version
version: version.F90
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;\
		s^COMPILER_VERSION^$(COMPILER_VERSION)^;\
		s'SIESTA_VERSION'$$(cat ${<D}/../version.info)'g;\
		s'LIBS'$(COMP_LIBS) $(LIBS)'g;\
		s'FPPFLAGS'$(FPPFLAGS)'g" $<  > tmp.F90
	@awk '{if (length>80) { cur=78; \
		printf "%s&\n",substr($$0,0,78); \
		while(length-cur>78){ cur=cur+76 ; \
		printf "&%s&\n",substr($$0,cur-76,76) \
		} printf "&%s\n",substr($$0,cur)} else { print $$0 }}' tmp.F90 > compinfo.F90
	$(MAKE) "FPPFLAGS=$(FPPFLAGS)" compinfo.o
	@mv compinfo.o version.o
	@echo

#----------------------------------------------------------------------------
# TranSIESTA objects needed for siesta
OBJS += m_ts_io.o

# TranSIESTA...
TS_OBJS= m_ts_global_vars.o m_ts_options.o m_ts_aux.o ts_kpoint_scf.o \
         m_ts_iodm.o \
         m_ts_electrode.o m_ts_gf.o \
	 m_ts_cctype.o ts_init.o ts_show_regions.o m_ts_electype.o \
	 m_ts_method.o m_ts_elec_se.o \
	 m_ts_debug.o m_ts_tdir.o \
	 m_ts_chem_pot.o m_transiesta.o m_ts_hartree.o m_ts_voltage.o

# contour objects
TS_OBJS += m_ts_contour_eq.o m_ts_contour_neq.o m_ts_io_ctype.o \
	   m_ts_contour.o m_ts_io_contour.o

# Transiesta objects (new memory reduced objects)
TS_OBJS += m_ts_rgn2trimat.o \
	m_ts_dm_update.o m_ts_weight.o \
	m_ts_sparse.o m_ts_sparse_helper.o \
	m_ts_charge.o m_ts_pivot.o

# Sparse full-matrix method
TS_OBJS += m_ts_fullg.o m_ts_fullk.o m_ts_full_scat.o

# Sparse tri-diagonal matrix method
TS_OBJS += m_ts_tri_init.o m_ts_tri_common.o \
	m_ts_trimat_invert.o \
	m_ts_tri_scat.o m_ts_trig.o m_ts_trik.o

# MUMPS method
TS_OBJS += m_ts_mumps_init.o m_ts_mumpsg.o m_ts_mumpsk.o m_ts_mumps_scat.o

OBJS += $(TS_OBJS)

# Create all objects
ALL_OBJS=$(OBJS) $(SYSOBJ) $(MAIN_OBJ)

# All objects require the XC library (in theory)
# In practice they do not but it greatly simplifies the logic
# here.
# This will setup the correct compilation sequence as
#  WXML => XC => FDF => MPI
$(ALL_OBJS): $(WXML) $(XC) $(MS) $(COMP_LIBS)

# Build options

libSiestaForces.a: $(OBJS) $(SYSOBJ)
	$(AR) $(ARFLAGS_EXTRA) cru libSiestaForces.a $(OBJS) $(SYSOBJ)
	@if [ ! -z "$(MPI_INTERFACE)" ] ; then \
            $(AR) $(ARFLAGS_EXTRA) cru libSiestaForces.a MPI/*.o ; fi
	-$(RANLIB) libSiestaForces.a

siesta: what version \
                $(MPI_INTERFACE) $(FDF) $(WXML) $(XMLPARSER) \
		$(XC) $(MS) $(COMP_LIBS) $(ALL_OBJS) 
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(WXML) $(XMLPARSER) \
               $(XC) $(MS) $(MPI_INTERFACE) $(COMP_LIBS) $(LIBS) 

# Notify users of erroneous transiesta compilation
transiesta:
	@echo ""
	@echo " TranSiesta executable is deprecated!"
	@echo ""
	@echo " TranSiesta method is enabled in the siesta executable"
	@echo ""

clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta *.o  *.a *.pcl *.pc
	rm -f *.mod
	rm -f compinfo.F90 tmp.F90
	rm -f *.bck *~
	rm -f aux_*.f aux_*.f90
	(cd fdf ; $(MAKE) clean)
	(cd xmlparser ; $(MAKE) clean)
	(cd wxml ; $(MAKE) clean)
	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi
	@if [ -d fdict/obj ] ; then $(MAKE) -C fdict/obj clean ; fi
	@if [ -d ncdf/obj ] ; then $(MAKE) -C ncdf/obj clean ; fi
	(cd Libs ; $(MAKE) clean )
	(cd SiestaXC ; $(MAKE) clean )
	(cd MatrixSwitch/src ; $(MAKE) clean )


# DO NOT DELETE THIS LINE - used by make depend
Fstack.o: basic_type.inc
Pair.o: basic_type.inc
alloc.o: m_io.o parallel.o precision.o sys.o
amn.o: alloc.o atmfuncs.o atomlist.o m_orderbands.o m_trialorbitalclass.o
amn.o: matel_registry.o mneighb.o new_matel.o parallel.o parallelsubs.o
amn.o: precision.o siesta2wannier90.o siesta_geom.o sys.o
arw.o: alloc.o parallel.o precision.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atm_types.o atmparams.o atom_options.o basis_specs.o basis_types.o
atom.o: bsc_xcmod.o interpolation.o m_filter.o periodic_table.o precision.o
atom.o: pseudopotential.o radial.o sys.o
atom_graph.o: alloc.o atm_types.o atmfuncs.o class_OrbitalDistribution.o
atom_graph.o: class_SpData1D.o class_SpData2D.o class_SpData2D.o
atom_graph.o: class_Sparsity.o intrinsic_missing.o ldau_specs.o mneighb.o
atom_graph.o: precision.o radial.o sorting.o sys.o
atom_options.o: sys.o
atomlist.o: alloc.o atm_types.o atmfuncs.o ldau_specs.o parallel.o precision.o
atomlist.o: radial.o siesta_geom.o
atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o sys.o
automatic_cell.o: atm_types.o atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o atmfuncs.o atomlist.o densematrix.o diag.o diag_option.o
bands.o: files.o get_kpoints_scale.o m_spin.o parallel.o parallelsubs.o
bands.o: precision.o siesta_geom.o sys.o writewave.o
basis_enthalpy.o: atmfuncs.o atomlist.o precision.o siesta_geom.o units.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o sys.o xml.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o sys.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o sys.o
bessph.o: precision.o sys.o
bonds.o: alloc.o atmfuncs.o mneighb.o precision.o sorting.o units.o
born_charge.o: atomlist.o ksv.o ksvinit.o m_mpi_utils.o m_spin.o parallel.o
born_charge.o: precision.o siesta_geom.o siesta_options.o sparse_matrices.o
broadcast_basis.o: atm_types.o chemical.o ldau_specs.o m_mpi_utils.o parallel.o
broadcast_basis.o: radial.o
broadcast_projections.o: m_trialorbitalclass.o parallel.o siesta2wannier90.o
broyden_optim.o: m_broyddj_nocomm.o m_memory.o parallel.o precision.o
broyden_optim.o: siesta_options.o sys.o units.o
bsc_cellxc.o: alloc.o bsc_xcmod.o cellxc_mod.o mesh.o moremeshsubs.o parallel.o
bsc_cellxc.o: parallelsubs.o precision.o sys.o
bsc_xcmod.o: parallel.o precision.o sys.o
cart2frac.o: precision.o sys.o
cell_broyden_optim.o: m_broyddj_nocomm.o parallel.o precision.o sys.o units.o
cell_broyden_optim.o: zmatrix.o
cell_fire_optim.o: alloc.o m_fire.o m_memory.o parallel.o precision.o
cell_fire_optim.o: siesta_options.o sys.o zmatrix.o
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
class_Data1D.o: class_Data1D.T90 basic_type.inc
class_Data1D.o: alloc.o
class_Data1D.o: basic_type.inc
class_Data1D.o: alloc.o
class_Data2D.o: class_Data2D.T90 basic_type.inc
class_Data2D.o: alloc.o
class_Data2D.o: basic_type.inc
class_Data2D.o: alloc.o
class_Distribution.o: basic_type.inc
class_Fstack_Data1D.o: Fstack.T90 basic_type.inc
class_Fstack_Data1D.o: class_Data1D.o class_Data1D.o
class_Fstack_Pair_Data1D.o: Fstack.T90 basic_type.inc
class_Fstack_Pair_Data1D.o: class_Pair_Data1D.o class_Pair_Data1D.o
class_Fstack_Pair_Geometry_SpData2D.o: Fstack.T90 basic_type.inc
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
class_Fstack_Pair_SpData1D.o: Fstack.T90 basic_type.inc
class_Fstack_Pair_SpData1D.o: class_Pair_SpData1D.o class_Pair_SpData1D.o
class_Geometry.o: basic_type.inc
class_Geometry.o: alloc.o
class_OrbitalDistribution.o: basic_type.inc
class_Pair_Data1D.o: Pair.T90 basic_type.inc
class_Pair_Data1D.o: class_Data1D.o class_Data1D.o
class_Pair_Geometry_SpData2D.o: Pair.T90 basic_type.inc
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
class_Pair_SpData1D.o: Pair.T90 basic_type.inc
class_Pair_SpData1D.o: class_SpData1D.o class_SpData1D.o
class_SpData1D.o: class_SpData1D.T90 basic_type.inc
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_OrbitalDistribution.o
class_SpData1D.o: class_Sparsity.o
class_SpData1D.o: basic_type.inc
class_SpData1D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData2D.o: class_SpData2D.T90 basic_type.inc
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_OrbitalDistribution.o
class_SpData2D.o: class_Sparsity.o
class_SpData2D.o: basic_type.inc
class_SpData2D.o: class_OrbitalDistribution.o class_Sparsity.o
class_Sparsity.o: basic_type.inc
class_Sparsity.o: alloc.o
class_TriMat.o: class_TriMat.T90 basic_type.inc
class_TriMat.o: alloc.o intrinsic_missing.o
class_TriMat.o: basic_type.inc
class_TriMat.o: alloc.o intrinsic_missing.o
coceri.o: files.o periodic_table.o precision.o units.o
compute_dm.o: atomlist.o class_SpData1D.o compute_ebs_shift.o diagon.o
compute_dm.o: iodmhs_netcdf.o kpoint_scf.o m_chess.o m_dminim.o m_energies.o
compute_dm.o: m_eo.o m_hsx.o m_pexsi_driver.o m_rmaxh.o m_spin.o m_steps.o
compute_dm.o: m_transiesta.o m_ts_global_vars.o m_zminim.o normalize_dm.o
compute_dm.o: ordern.o parallel.o precision.o siesta_geom.o siesta_options.o
compute_dm.o: sparse_matrices.o sys.o units.o
compute_ebs_shift.o: m_mpi_utils.o parallel.o precision.o
compute_energies.o: atomlist.o class_SpData1D.o class_SpData2D.o
compute_energies.o: class_SpData2D.o dhscf.o files.o m_dipol.o m_energies.o
compute_energies.o: m_mpi_utils.o m_ntm.o m_rhog.o m_spin.o parallel.o
compute_energies.o: precision.o siesta_geom.o siesta_options.o
compute_energies.o: sparse_matrices.o
compute_max_diff.o: m_mpi_utils.o precision.o
compute_norm.o: m_mpi_utils.o m_spin.o precision.o sparse_matrices.o
compute_pw_matrix.o: alloc.o m_planewavematrix.o m_planewavematrixvar.o
compute_pw_matrix.o: parallel.o precision.o siesta2wannier90.o
compute_pw_matrix.o: sparse_matrices.o
compute_rhog.o: alloc.o atomlist.o dhscf.o files.o m_dipol.o m_energies.o
compute_rhog.o: m_ntm.o m_spin.o m_stress.o siesta_geom.o siesta_options.o
compute_rhog.o: sparse_matrices.o
conjgr.o: precision.o
conjgr_old.o: precision.o
coor.o: alloc.o parallel.o precision.o siesta_geom.o sys.o units.o zmatrix.o
coxmol.o: files.o periodic_table.o precision.o units.o
cranknic_evolg.o: alloc.o atomlist.o m_energies.o m_eo.o m_inversemm.o
cranknic_evolg.o: m_matswinvers.o m_spin.o m_steps.o parallel.o parallelsubs.o
cranknic_evolg.o: precision.o siesta_options.o sparse_matrices.o sys.o units.o
cranknic_evolg.o: wavefunctions.o
cranknic_evolk.o: atomlist.o cranknic_evolg.o kpoint_scf.o m_energies.o m_eo.o
cranknic_evolk.o: m_spin.o parallel.o parallelsubs.o precision.o
cranknic_evolk.o: siesta_options.o sparse_matrices.o units.o wavefunctions.o
create_Sparsity_SC.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
create_Sparsity_Union.o: class_OrbitalDistribution.o class_Sparsity.o
create_Sparsity_Union.o: m_region.o parallel.o precision.o
debugmpi.o: parallel.o
delk.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o
delk.o: m_planewavematrixvar.o mesh.o meshdscf.o meshphi.o parallel.o
delk.o: parallelsubs.o precision.o siesta_geom.o
denmat.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
denmatlomem.o: alloc.o globalise.o onmod.o precision.o
densematrix.o: alloc.o precision.o
detover.o: alloc.o linpack.o parallel.o parallelsubs.o precision.o
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
dfscf.o: meshphi.o parallel.o parallelsubs.o precision.o sys.o
dhscf.o: alloc.o atmfuncs.o bsc_xcmod.o cellxc_mod.o delk.o dfscf.o
dhscf.o: doping_uniform.o files.o forhar.o iogrid_netcdf.o m_charge_add.o
dhscf.o: m_efield.o m_hartree_add.o m_iorho.o m_iotddft.o m_mesh_node.o
dhscf.o: m_ncdf_io.o m_ncdf_siesta.o m_partial_charges.o m_rhog.o m_spin.o
dhscf.o: m_ts_global_vars.o m_ts_hartree.o m_ts_options.o m_ts_voltage.o mesh.o
dhscf.o: meshdscf.o meshsubs.o moremeshsubs.o parallel.o parsing.o precision.o
dhscf.o: rhofft.o rhoofd.o siesta_options.o sys.o units.o vmat.o
diag.o: alloc.o diag_option.o parallel.o precision.o sys.o
diag2g.o: fermid.o intrinsic_missing.o parallel.o parallelsubs.o precision.o
diag2g.o: sys.o
diag2k.o: fermid.o parallel.o parallelsubs.o precision.o sys.o
diag3g.o: fermid.o intrinsic_missing.o parallel.o parallelsubs.o precision.o
diag3g.o: sys.o
diag3k.o: fermid.o parallel.o parallelsubs.o precision.o sys.o
diag_option.o: parallel.o precision.o
diagg.o: alloc.o fermid.o intrinsic_missing.o m_spin.o parallel.o
diagg.o: parallelsubs.o precision.o sys.o writewave.o
diagk.o: compute_norm.o fermid.o m_spin.o parallel.o parallelsubs.o precision.o
diagk.o: sys.o writewave.o
diagk_file.o: fermid.o iowfs_netcdf.o parallel.o parallelsubs.o precision.o
diagk_file.o: sys.o writewave.o
diagkp.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagon.o: alloc.o densematrix.o diag_option.o m_memory.o m_spin.o parallel.o
diagon.o: parallelsubs.o precision.o siesta_options.o sys.o
diagonalizeHk.o: alloc.o atomlist.o densematrix.o m_orderbands.o m_spin.o
diagonalizeHk.o: parallel.o parallelsubs.o precision.o siesta2wannier90.o
diagonalizeHk.o: sparse_matrices.o units.o
diagpol.o: parallel.o precision.o sys.o
diagsprl.o: atomlist.o fermid.o parallel.o parallelsubs.o precision.o
diagsprl.o: siesta_geom.o sys.o
die.o: parallel.o siesta_cml.o
die.o: parallel.o
dipole.o: mesh.o parallel.o precision.o sys.o
dm_charge.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
dm_charge.o: class_Sparsity.o m_mpi_utils.o m_spin.o precision.o
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
extrapolateSpData2D.o: class_Data2D.o class_OrbitalDistribution.o
extrapolateSpData2D.o: class_SpData2D.o class_Sparsity.o restructSpData2D.o
extrapolon.o: parallel.o precision.o sys.o
fdf_extra.o: m_region.o
fermid.o: errorf.o parallel.o precision.o sys.o
fft.o: alloc.o fft1d.o m_timer.o mesh.o parallel.o parallelsubs.o precision.o
fft.o: sys.o
fft1d.o: parallel.o precision.o sys.o
