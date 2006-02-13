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
default: libnetcdf_f90.a
#
include ../arch.make
#FFLAGS=$(FFLAGS) $(FFLAGS_NETCDF)
SRCS=netcdf.f90 netcdf_sizes.f90 netcdf_sizesOK.f90 \
     netcdf_constants.f90 netcdf_externals.f90\
     netcdf_overloads.f90 netcdf_visibility.f90 \
     netcdf_file.f90 netcdf_dims.f90 netcdf_attributes.f90 \
     netcdf_variables.f90 netcdf_text_variables.f90 netcdf_expanded.f90

#
libnetcdf_f90.a: netcdf.o
	ar $(ARFLAGS_EXTRA) cru libnetcdf_f90.a netcdf.o
	-$(RANLIB) libnetcdf_f90.a
	@mv libnetcdf_f90.a ..
	@cp *.mod ..
#
#netcdf_expanded.f90: gen.m4 aux.m4 nves.m4 nvea.m4
#	m4 gen.m4 > netcdf_expanded.f90
#
netcdf.o: netcdf.f90 $(SRCS)

clean:
	/bin/rm -f *.o *.mod *.MOD module_built
