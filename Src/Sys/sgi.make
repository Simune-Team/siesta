# $Id: sgi.make,v 1.5 1999/04/26 13:12:16 wdpgaara Exp $
#
FC=f77 -mips4 -static
# Careful with -mips4! If turning single into double precision, switch to -32
#
FFLAGS= -O
# -O3 gives problems with -c in some versions
#
FFLAGS_DEBUG= -g -O0
LIBS=
SYS=sgi
# To use the sgimath library, use sgimath.make
#
RANLIB=@echo No ranlib required for 
COMMENTS=" ** Make sure you are using the right ABI and ISA (see man abi(5))"
#
