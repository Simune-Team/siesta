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
XML-TESTER=../../Util/test-xml/test-xml
XML-REFERENCE=../../Tests/Reference-xml
#
label=work
#
completed: completed_$(label)
#
completed_$(label):
	@echo ">>>> Running $(name)1 test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) $(label) ; fi
	@for i in `cat $(name).pseudos` ; do \
          echo "    ==> Copying pseudopotential file for $$i..." ;\
          ln ../Pseudos/$$i.psf $(label)/$$i.psf ;\
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
xmlcheck: completed
	@echo "    ==> Running xmlcheck for system $(name)2"
	@ln -sf ../tolerances.dat ./tolerances.dat
	$(XML-TESTER) $(XML-REFERENCE)/$(name)2.xml $(label)/$(name).xml | tee $(label).diff-xml
        # The following line erases the file if it is empty
	@if [ ! -s $(label).diff-xml ] ; then rm -f $(label).diff-xml ; fi
#
clean:
	@echo ">>>> Cleaning $(name)1 test..."
	rm -rf $(label) completed* $(name)1.out $(name)2.out $(name).xml *.dat *diff-xml
