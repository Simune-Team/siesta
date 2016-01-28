#
# Single-test makefile template
#
# You can edit the SIESTA macro here, or pass it on the command line
#
SIESTA=../../../siesta
#
# Example for BSC runs
#
#SIESTA= mpirun -np 4 ../../../siesta
#
#----------------------------------------------------------------------------
REFERENCE_DIR?=../../../Tests/Reference
REFERENCE_CHECKER?=../cmp_digest.sh
#
label=work
#
completed: completed_$(label)
#
completed_$(label):
	@echo ">>>> Running $(name) test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) $(label) ; fi
	@for i in `cat $(name).pseudos` ; do \
          echo "    ==> Copying pseudopotential file for $$i..." ;\
          ln ../Pseudos/$$i.psf $(label)/$$i.psf ;\
         done
	@echo "    ==> Running SIESTA as ${SIESTA}"
	@(cd $(label) ; ${SIESTA} 2>&1 > $(name).out < ../$(name).fdf) \
          && touch completed_$(label)
	@if [ -f completed_$(label) ] ; then cp $(label)/$(name).out $(label)/$(name).xml .;\
           echo "    ===> SIESTA finished successfully";\
         else \
           echo " **** Test $(name) did not complete successfully";\
         fi
#
check: completed
	@echo "    ==> Running check for system $(name)"
	@REFERENCE_DIR=$(REFERENCE_DIR) sh $(REFERENCE_CHECKER) $(name).out
#
clean:
	@echo ">>>> Cleaning $(name) test..."
	rm -rf $(label) completed_$(label) $(name).out $(name).xml
