# $Id: g77.make,v 1.1.2.1 1999/06/04 14:47:11 wdpgaara Exp $
#
#  This file will actually work for any system running modern
#  versions of g77.
#
FC=g77
FFLAGS=  -O6 -malign-double -funroll-loops -fomit-frame-pointer -Wall
FFLAGS_DEBUG= -g -O0
LIBS=
SYS=bsd
RANLIB=ranlib
#
