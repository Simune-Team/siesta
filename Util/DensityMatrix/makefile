# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#
# Makefile for DM to DM-netCDF converter
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: dm2cdf cdf2dm
include ../../Src/arch.make
#
# We assume that the top arch.make defines the relevant NetCDF symbols,
# including the INCFLAGS line. If not, uncomment the following lines, 
# using appropriate information
#
# NETCDF_ROOT=/some/path
# INCFLAGS=-I$(NETCDF_ROOT)/include
#
LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf
#
     
dm2cdf: dm2cdf.o
	$(FC) $(LDFLAGS) -o $@ dm2cdf.o  $(LIBS)
#
cdf2dm: cdf2dm.o
	$(FC) $(LDFLAGS) -o $@ cdf2dm.o  $(LIBS)
#
clean: 
	rm -f *.o dm2cdf cdf2dm *.o *.*d
#










