# $Id: linux.make,v 1.6 1999/04/26 10:24:12 wdpgaara Exp $
#
# Use this file if you have an old version of g77 which does not
# come with the standard timing routines. If you have a modern
# version, use g77.make instead.
#
#
FC=g77
FFLAGS=  -O6 -malign-double -funroll-loops -fomit-frame-pointer -Wall
FFLAGS_DEBUG= -g -O0
LIBS=
SYS=general
RANLIB=@echo No ranlib required for
#
