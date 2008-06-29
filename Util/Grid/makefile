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
# Makefile for Grid utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: grid2cdf 
all: grid2cdf grid2val grdi2cube
#
include ../../Src/arch.make
FFLAGS=$(FFLAGS_DEBUG)
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
#------------------
grid2cdf: grid2cdf.o
	$(FC) $(LDFLAGS) -o $@  grid2cdf.o  $(LIBS)
#
cdf2grid: cdf2grid.o
	$(FC) $(LDFLAGS) -o $@  cdf2grid.o  $(LIBS)
#
cdf2xsf: cdf2xsf.o
	$(FC) $(LDFLAGS) -o $@  cdf2xsf.o  $(LIBS)
#
cdfget: cdfget.o
	$(FC) $(LDFLAGS) -o $@  cdfget.o  $(LIBS)
cdf_diff: m_grid.o cdf_diff.o 
	$(FC) $(LDFLAGS) -o $@  m_grid.o cdf_diff.o  $(LIBS)
cdf_diff.o: m_grid.o
#------------------
cdf_fft: cdf_fft.o libfft.o
	$(FC) $(LDFLAGS) -o $@  cdf_fft.o  libfft.o $(LIBS)
#
#------------------
get_f2kcli:
	ln -sf ../../Src/f2kcli.F90 ./f2kcli.F90
#
grid2val.o: f2kcli.o
grid2val: get_f2kcli grid2val.o
	$(FC) $(LDFLAGS) -o $@ f2kcli.o grid2val.o 
#
#------------------
grid2cube: grid2cube.o
	$(FC) $(LDFLAGS) -o $@ grid2cube.o 
#
#------------------
clean: 
	rm -f *.o grid2cdf grid2val grid2cube *.o *.*d
#










