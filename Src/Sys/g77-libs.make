# $Id: g77-libs.make,v 1.1.2.1 1999/06/04 14:47:32 wdpgaara Exp $
#
#  This file will actually work for any system running modern
#  versions of g77, and with the BLAS and LAPACK libraries.
#
FC=g77
FFLAGS=  -O6 -malign-double -funroll-loops -fomit-frame-pointer -Wall
FFLAGS_DEBUG= -g -O0
LIBS= -llapack -lblas
SYS=bsd-libs
RANLIB=ranlib
#
