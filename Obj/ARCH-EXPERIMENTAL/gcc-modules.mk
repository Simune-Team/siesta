#
# Defs for a GCC system, with all library symbols loaded from modules
#
# ml flook netcdf openmpi
# ml scalapack      # defines SCALAPACK_LIBS
# ml lapack   (Mac with homebrew: ml veclibfort)  # defines LAPACK_LIBS
# ml fftw
# ml elpa
# ------------ for later versions
# ml gridxc-multi libpsml
# ml elsi-ext-elpa
#--------------------------------------------------------
# Use these symbols to request particular features
# To turn on, set '=1'.
#--------------
#
WITH_EXTERNAL_ELPA=1
WITH_FLOOK=1
WITH_MPI=1
WITH_NETCDF=1
WITH_SEPARATE_NETCDF_FORTRAN=
WITH_NCDF=1
WITH_LEGACY_GRIDXC_INSTALL=
WITH_GRID_SP=
WITH_ELSI=
# These are mandatory for PSML and MaX Versions,
# but they should be turned off for 4.1
WITH_PSML=
WITH_GRIDXC=
#-------------
#
# Needed for PEXSI (ELSI) support
LIBS_CPLUS=-lstdc++ 
#--------------------------------------------------------
#
# Define compiler names and flags
#
FC_PARALLEL=mpif90
FC_SERIAL=gfortran
#
FPP = $(FC_SERIAL) -E -P -x c
FFLAGS= -O2 -g 
FFLAGS_DEBUG= -g -O0 -fcheck=all
RANLIB=echo
# ----------------------------------------------------------
