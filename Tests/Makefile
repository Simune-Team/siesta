#
# Makefile for SIESTA tests
#
# Usage:
#
#          make                     --- Run all tests
#          make check               --- Run all tests and check with reference file
#          make show                --- Show the list of tests
#
#          make tests=fe_broyden { check | }     --- Run a single test
#          (Not implemented) make batch                --- Dispatch to batch system
#                                        (site-dependent)
#                                        (needs extra configuration)
#
#          By default, the version of SIESTA used is the one in the
#          parent directory (Obj or custom OBJDIR), and is executed
#          in serial mode.
#
#          To use another one, or a parallel version
#
#          make SIESTA=/some/path/to/siesta-program ....
#
#           or
#
#          make SIESTA="mpiexec -np 4 /path/to/parallel-siesta-program" ....
#
#          It is also possible to have separate working directories,
#          by using the a "label". For example:
#
#          make label=finer 
#     
#          The default label is "work".
#
# See also test.mk for details, and the scripts in Scripts for 
# examples of other options for batch processing.
#
label = work

tests = 32_h2o ag anneal-cont ar2_vdw batio3 benzene bessel born \
	born_spin carbon_nanoscroll ch4 chargeconf-h2o constant_volume \
	dipole_correction fe fe_broyden fe_clust_noncollinear fe_cohp \
	fen fe_noncol_kp fire_benzene floating force_2 force_constants \
	ge111 graphite_c6 graphite_c6_full graphite_vdw_df h2_bessel \
	h2o h2o_2 h2o_am05 h2o_basis h2o_bands h2o_bands_polarized h2o_bands_nc \
	h2o_coop h2o_dipole h2o_dipole2 \
	h2o_dos h2o_filteret_basis h2o_findp_bug h2o_netcdf h2o_op_broyden \
	h2o_orderN h2o-pop h2o_radialgrid h2o_reparam h2o-rhog-mixing h2oZ \
	h3po4_1 h3po4_2 h3po4_3 h_chain h_chain2 \
	md_anneal md-anneal-cont md-cont md-cont-no-restart-file \
	md_nose md_npr md_pr md-single-run md_verlet mgc-force mgco3 \
	nanotube-c-5-0 n_chain OMM_h2o OMM_si oxyn partial pb_bulk pb_bulk_1000 \
	pb_filter_cutoff pb_filter_tol pd_ice ptcda ptcda-au si001 si001-ldos \
	si111-spinpol si2x1h si2x1h-dipole si2x1h-dipole-gcs si2x1h-quench \
	si64 si64_coop si_bandpoints si_coop sic-slab si_fatbands sih sih_fire \
	sih_op_broyden sinw sinw_2 si-optical si_pdos_gamma si_pdos_kgrid \
	var_cell wannier zmatrix

# Currently removed tests:
#  h2o_4 (ON, does not work in parallel)
#  ipi-test (requires ipi)
#  TranSiesta-TBTrans (requires new compilation)

default: all

all:
	@for dir in $(tests); do\
	 $(MAKE) -C $$dir completed; \
	done

check:
	@for dir in $(tests); do\
	 $(MAKE) -C $$dir check; \
	done

show:
	@echo $(tests)

clean:
	@for i in $(tests) ; do \
         (cd $$i ; make clean) ; \
	done


