# $Id: sgimath.make,v 1.1 1999/04/23 10:10:51 emilio Exp $
#
FC=f77
FFLAGS= -mips4 -O3 
# -mips4 can turn single precision to double, change to -32 if it does
# if -O3 gives problems with -c, change to -O
#
FFLAGS_DEBUG= -mips4 -g -O0
LIBS=-lcomplib.sgimath 
# CAREFUL!! In some versions, it is quite unreliable for diagonalization
#
SYS=sgi
RANLIB=@echo No ranlib required for 
COMMENTS=" ** Make sure you are using the right ABI and ISA (see man abi(5))"
#
