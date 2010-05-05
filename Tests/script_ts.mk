#
# Single-test makefile template for script usage
#
OBJDIR=Obj_TS
#
completed:
	@echo ">>>> Running $(name) test..."
	@if [ -d work ] ; then rm -rf work ; fi; mkdir work
	@if [ -f script.sh ] ; then cp -f script.sh work ; fi
	@echo "    ==> Running script with TranSIESTA as ${TS}"
	@(cd work ; sh script.sh "${OBJDIR}" )
	@if [ -f completed ] ; then \
           echo "    ===> Script finished successfully";\
         else \
           echo " **** Test $(name) did not complete successfully";\
         fi
#
clean:
	@echo ">>>> Cleaning $(name) test..."
	rm -rf work completed *.out 
