# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=pathf90-2.0-serial
#
# For an opteron cluster with Pathscale compiler at UniCan
#
FC=pathf90-2.0 -ff2c-abi /localhome/jjunquer/acml2.5.1/pathscale64 -m64 -ipa
FC_ASIS=$(FC)
#
FFLAGS= -Ofast -Wl,-R/opt/PathScale/lib/2.0/
FFLAGS_DEBUG= -g 
LDFLAGS=-static 
RANLIB=echo
COMP_LIBS=
#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
FPPFLAGS_CDF=            #  -DCDF
#
MPI_INTERFACE=
MPI_INCLUDE=
FPPFLAGS_MPI=
#
LIBS=/opt/PathScale/lib/2.0/libpathfortran.a \
	/localhome/jjunquer/acml2.5.1/pathscale64/lib/libacml.a 
SYS=cpu_time
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)  # -DGRID_DP
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<

