#
FC=f77
FFLAGS= -mips4 -O3 
FFLAGS_DEBUG= -mips4 -g -O0
LIBS=-lcomplib.sgimath 
SYS=sgi
RANLIB=@echo No ranlib required for 
COMMENTS=" ** Make sure you are using the right ABI and ISA (see man abi(5))"
#
