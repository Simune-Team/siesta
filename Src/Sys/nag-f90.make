# $Id: nag-f90.make,v 1.1.2.1 1999/06/04 14:47:33 wdpgaara Exp $
#
# For NAG  V.2.x compiler. Your mileage may vary.
#
# -C=array will enable array bounds checking.
# -dusty and -x77 are useful for legacy codes
#
# It should be possible to construct timing routines using
# the Fortran90 intrinsics.
#
FC=f90
FFLAGS= -O4 -C=none -dusty -x77
FFLAGS_DEBUG= -g -O0 -C=none -dusty -x77
LIBS=
SYS=general
RANLIB=@echo No ranlib required for
#
