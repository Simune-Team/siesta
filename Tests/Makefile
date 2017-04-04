#
# Makefile for SIESTA tests
#
# Usage:
#
#          make                      --- Run all tests
#          make  check               --- Run all tests and check with reference file
#          make  check-only          --- Run all checks with reference file (run make before)
#          make  list                --- Show the list of tests
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

# Default target
default: all

label=work


tests = ag anneal-cont ar2_vdw \
	batio3 benzene bessel \
	born born_spin \
	carbon_nanoscroll \
	ch4 chargeconf-h2o constant_volume \
	dipole_correction \
	fe fe_broyden fe_clust_noncollinear fe_clust_noncollinear-gga \
	fe_cohp fen fe_noncol_kp \
	fire_benzene floating force_2 force_constants \
    gate_G_charge gate_G_hartree \
	ge111 graphite_c6 graphite_vdw_df \
	h2_bessel \
	h2o h2o_2 h2o_basis h2o_dipole h2o_dipole2 h2o_dos \
	h2o_findp_bug h2o_op_broyden h2o-pop h2o_radialgrid \
	h2o_reparam h2o-rhog-mixing h2oZ \
	h_chain \
	md_anneal md_nose md_npr md_pr md_verlet \
	mgc-force mgco3 \
	n_chain \
	OMM_h2o OMM_si \
	oxyn \
	pb_bulk pb_filter_cutoff pb_filter_tol \
	scfmix \
	si001 si111-spinpol si2x1h si2x1h-dipole \
	si64 si64_coop \
	si_bandpoints si_coop \
	sic-slab \
	sih sih_fire sih_op_broyden \
	si-optical si_pdos_kgrid \
	sinw sinw_2 \
	offsite_SOC_Pt2_xx offsite_SOC_Pt2_xz offsite_SOC_Pt2_zy offsite_SOC_Pt2_zz \
	onsite_SOC_Pt2_xx onsite_SOC_Pt2_xz onsite_SOC_Pt2_zy onsite_SOC_Pt2_zz \
	var_cell \
	wannier \
	zmatrix TDDFT_h2o

# These tests are extremely time consuming and
# should only be runned sometimes
# Currently they may be runned individually.
#	offsite_SOC_FePt_xx offsite_SOC_FePt_xz offsite_SOC_FePt_zy offsite_SOC_FePt_zz
#	onsite_SOC_FePt_xx onsite_SOC_FePt_xz onsite_SOC_FePt_zy onsite_SOC_FePt_zz
tests_soc_offsite  = offsite_SOC_Pt2_xx offsite_SOC_Pt2_xz offsite_SOC_Pt2_zy offsite_SOC_Pt2_zz
tests_soc_offsite += offsite_SOC_FePt_xx offsite_SOC_FePt_xz offsite_SOC_FePt_zy offsite_SOC_FePt_zz
tests_soc_onsite  = onsite_SOC_Pt2_xx onsite_SOC_Pt2_xz onsite_SOC_Pt2_zy onsite_SOC_Pt2_zz
tests_soc_onsite += onsite_SOC_FePt_xx onsite_SOC_FePt_xz onsite_SOC_FePt_zy onsite_SOC_FePt_zz

# Tests only applicable for FLOS
tests_flos  = flos_h2o_lbfgs flos_h2o_fire flos_h2o_meshcutoff flos_h2o_neb

# Accumulated variable for ALL tests
tests_all = $(tests) $(tests_soc_offsite) $(tests_soc_onsite) $(tests_flos)

all: $(tests)

# Define test specific targets:
tests: $(tests)
tests-soc: $(tests_soc_offsite) $(tests_soc_onsite)
tests-flos: $(tests_flos)

# Ensure that tests may be called individually
# (and that they are run, PHONY is required)
# This ensures one can call:
#   make h2o
# and the h2o test will be runned.
.NOTPARALLEL: $(tests_all)
.PHONY: $(tests_all)
$(tests_all):
	$(MAKE) -C $@ completed

.PHONY: check
check:
	@for dir in $(tests) ; do\
	 $(MAKE) -C $$dir check; \
	done

check-only:
	@for dir in $(tests) ; do\
	 $(MAKE) -C $$dir check-only ; \
	done

.PHONY: show
show:
	@echo $(tests)

clean:
	@for i in $(tests_all) ; do \
		(cd $$i ; make clean) ; \
	done
