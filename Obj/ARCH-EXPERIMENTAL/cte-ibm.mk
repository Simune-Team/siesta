#
# Defs for CTE-Power at BSC. IBM xl-16 compilers
#
SCALAPACK_ROOT=/apps/SCALAPACK/2.0.2/IBM/SPECTRUM
LAPACK_ROOT=/apps/LAPACK/3.8.0/IBM
ESSL_ROOT=/usr
SCALAPACK_LIBS=-L${SCALAPACK_ROOT}/lib -lscalapack
LAPACK_LIBS=-L${LAPACK_ROOT}/lib64 -llapack  -L${ESSL_ROOT}/lib64 -lessl 
#
# Uses headers -- you still need the LD_LIBRARY_PATH at execution time
FFTW_ROOT="/apps/FFTW/3.3.8/GCC/OPENMPI"
#
# Needed for PEXSI (ELSI) support
LIBS_CPLUS=-lstdc++ 
#--------------------------------------------------------
#
# Define compiler names and flags
#
FCFLAGS_fixed_f=-qfixed -qsuffix=cpp=f
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=-qfixed -qsuffix=cpp=F
FPPFLAGS_free_F90=
#
FPP_PREFIX=-WF,
DEFS_PREFIX=-WF,
OMPI_FC=xlf2008_r   # To use this in the mpixlf wrapper
FC_PARALLEL=mpixlf
FC_SERIAL=xlf2008_r
#
FPP = $(FC_SERIAL) -E -P -x c
FFLAGS= -O0 -g -qarch=pwr9 -qstrict
FFLAGS_DEBUG= -g -O0
RANLIB=echo
# ----------------------------------------------------------
