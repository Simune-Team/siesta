#
# Single-test makefile template
#
SIESTA= mpirun -np 1 -machinefile /home/bsc21/bsc21246/siesta-2.1.29_BSC/Tests/hfile /home/bsc21/bsc21246/siesta-2.1.29_BSC/Src/siesta
#SIESTA=../../../Src/siesta
#
completed:
	@echo ">>>> Running $(name) test..."
	@if [ -d work ] ; then rm -rf work ; fi; mkdir work
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) work ; fi
	@for i in `cat $(name).pseudos` ; do \
          echo "    ==> Copying pseudopotential file for $$i..." ;\
          ln ../Pseudos/$$i.psf work/$$i.psf ;\
         done
	@echo "    ==> Running SIESTA as ${SIESTA}"
	@(cd work ; ${SIESTA} 2>&1 > $(name).out < ../$(name).fdf) \
          && touch completed
	@if [ -f completed ] ; then cp work/$(name).out work/$(name).xml .;\
           echo "    ===> SIESTA finished successfully";\
         else \
           echo " **** Test $(name) did not complete successfully";\
         fi
#
xmlcheck: completed
	@echo "---- xmllint check $(name).xml ..."
	xmllint $(name).xml > /dev/null
#
clean:
	@echo ">>>> Cleaning $(name) test..."
	rm -rf work completed $(name).out $(name).xml
