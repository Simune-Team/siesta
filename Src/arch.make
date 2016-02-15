#                                                                   
# This file is part of the SIESTA package.                          
#                                                                   
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:   
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .                                              
#                                                                     
# Use of this software constitutes agreement with the full conditions 
# given in the SIESTA license, as signed by all legitimate users.     
#                                                                     
# ------------                                                        
# Description:                                                        
#              Intel compiler/mkl  V12 on atto                        
#              OpenMPI with support for Intel compiler V12            
#              MKL 10.3 libraries, including a version of BLACS provided
#                             by Intel for the openmpi framework,       
#                             and Intel's own Scalapack, Lapack, and BLAS.
#                                                                         
# Execution:                                                              
#                                                                         
#           $(OPENMPI_ROOT)/bin in PATH                                   
#           $(OPENMPI_ROOT)/lib in LD_LIBRARY_PATH                        
#                                                                         
#           mpirun -np NPROCS siesta ....                                 
#                                                                         
SIESTA_ARCH=femto-intel-openmpi
OPENMPI_ROOT=/share/apps/openmpi-1.6.5-intel/
#                                                                         
#                                                                         
#FC=$(OPENMPI_ROOT)/bin/mpif90
FC=/share/apps/openmpi-1.6.5-intel/bin/mpif90
#
#                                                                         
#  You can play with other optimization options                           
#  I am not sure whether the compiler attempts to multithread the code    
#                                                                         
FFLAGS= -w  -O2 -mp                                                       
#FFLAGS=-g -O0 -debug full -traceback -C                                  
FFLAGS_CHECKS=-g -O0 -debug full -traceback -C                            
FFLAGS_DEBUG= -g                                                          
LDFLAGS=                                                                  
COMP_LIBS=                                                                
RANLIB=echo                                                               
#                                                                         
# You might want to turn off FoX                                          
#                                                                         
#DUMMY_FOX=--enable-dummy                                               
#                                                                         
#                                                                         
#NETCDF_ROOT=/share/apps/netcdf-3.6.2-ifort                                
#NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include                                  
#FPPFLAGS_CDF=-DCDF                                                        
#                                                                         
MPI_INTERFACE=libmpi_f90.a                                                
MPI_INCLUDE=.      # Note . for no-op                                     
FPPFLAGS_MPI=-DMPI                                                        
#                                                                         
#NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
#
# From the "Intel advisor"
#
MKLPATH=/share/apps/intel/composerxe-2011.0.084/mkl
#MKLPATH=/share/apps/intel/composer_xe_2013.3.163/mkl/
SUGGESTED_LIBS=-L$(MKLPATH)/lib/intel64 -lmkl_scalapack_lp64 \
               -Wl,--start-group \
                   -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
                   -lmkl_blacs_openmpi_lp64 \
               -Wl,--end-group \
               -lpthread
#
LIBS=$(SUGGESTED_LIBS) $(NETCDF_LIBS) $(METIS_LIB)
#
SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#

