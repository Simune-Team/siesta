#
# Makefile for test_filter
#
FC=gfortran -DGFORTRAN
# 
# FC = ifort
#
#
test_filter: f2kcli.F90  m_getopts.f90 m_filter.f90 test_filter.f90  
	$(FC) -o test_filter f2kcli.F90  m_getopts.f90 m_filter.f90 test_filter.f90 

#
#
run_tests: test_filter
	test_filter -l 0 -f core.dat -o core.log
        #
	test_filter -l 1 -k 0.7 -n 2 -f ti-3p-orb.dat -o ti-3p-orb.log
        #
	# Expect trouble here...
        #
	test_filter -l 1 -k 0.7 -n 2 -f ti-4p-orb.dat -o ti-4p-orb.log
        #
	test_filter -l 0 -k 0.7 -n 2 -f ba-6s-orb.dat -o ba-6s-orb.log
	# A given file could be used with several values of L...
	test_filter -l 1 -f core.dat -o core-with-l1.log
	# Ba 6s also works with L=1 
	test_filter -l 1 -k 0.7 -n 2 -f ba-6s-orb.dat -o ba-6s-orb-with-l1.log

clean:
	rm -f test_filter *.o *.mod *.log fort.33

