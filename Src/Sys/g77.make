# $Id: g77.make,v 1.1 1999/02/26 10:22:37 wdpgaara Exp $
#
#  This file will actually work for any system running modern
#  versions of g77.
#
FC=g77
FFLAGS=  -O6 -malign-double -funroll-loops -fomit-frame-pointer -Wall
FFLAGS_DEBUG= -g -O0
LIBS=
SYS=bsd
RANLIB=@echo No ranlib required for
#
