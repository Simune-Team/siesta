# $Id: nag-f95.make,v 1.1.2.1 1999/06/04 14:47:33 wdpgaara Exp $
#
# (Experimental) For NAG (f95) V.4. compiler.
#
# dcfuns allows the standard f77 behavior of the complex*16 intrinsics
# -C=all will make the program stop and complain about inconsistent
#        interfaces ('legal' in f77 and used in some instances for
#        efficiency in siesta)
# -C=array will enable array bounds checking.
#
#
FC=f95
FFLAGS= -O4 -C=none -dcfuns
FFLAGS_DEBUG= -g -O0 -C=none -dcfuns
LIBS=
SYS=general
RANLIB=@echo No ranlib required for
#
