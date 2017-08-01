#
# Single-test makefile template (TDDFT version)
#
# You can edit the SIESTA macro here, or pass it on the command line

MPI=mpirun -np 2
SIESTA=../../../siesta

# Example for BSC runs
#
#MPI=mpirun -np 2
#SIESTA= ../../../siesta

# Make compatibility layer for old test-runs
ifeq ($(strip $(firstword $(SIESTA))),mpirun)
MPI=
endif
ifeq ($(strip $(firstword $(SIESTA))),mpiexec)
MPI=
endif

#----------------------------------------------------------------------------
REFERENCE_DIR?=../../../Tests/Reference
REFERENCE_CHECKER?=../cmp_digest.sh

label=work
#
.PHONY: completed
completed: completed_$(label)
#
completed_$(label):
	@echo ">>>> Running $(name)1 test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) $(label) ; fi
	@for i in `cat $(name).pseudos` ; do \
	    if [[ $$i =~ .*\.(psf|psml) ]] ; then fps=$$i ; \
                    else fps=$${i}.psf ; fi; \
            echo "    ==> Copying pseudopotential file for $$fps ..." ;\
          ln ../Pseudos/$$fps $(label)/$$fps ;\
         done
	@echo "    ==> Running SIESTA as ${SIESTA}"
	@(cd $(label) ; ${SIESTA} 2>&1 > $(name)1.out < ../$(name)1.fdf) \
          && touch completed_$(label)
	@if [ -f completed_$(label) ] ; then cp $(label)/$(name)1.out $(label)/$(name).xml .;\
           echo "    ===> SIESTA finished successfully";\
         else \
           echo " **** Test $(name)1 did not complete successfully";\
         fi
	@echo ">>>> Running $(name)2 test..."
	@echo "    ==> Running SIESTA as ${SIESTA}"
	@(rm -f completed* $(name).xml)
	@(cd $(label) ; ${SIESTA} 2>&1 > $(name)2.out < ../$(name)2.fdf) \
          && touch completed_$(label)
	@if [ -f completed_$(label) ] ; then cp $(label)/$(name)2.out $(label)/$(name).xml .;\
           echo "    ===> SIESTA finished successfully";\
         else \
           echo " **** Test $(name)2 did not complete successfully";\
         fi
#
check: completed check-only

check-only:
	@echo "    ==> Running check for system $(name)"
	@REFERENCE_DIR=$(REFERENCE_DIR) sh $(REFERENCE_CHECKER) $(name).out
#
clean:
	@echo ">>>> Cleaning $(name)1 test..."
	rm -rf $(label) completed* $(name)1.out $(name)2.out $(name).xml *.dat 