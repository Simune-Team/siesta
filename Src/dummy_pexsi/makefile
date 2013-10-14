# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#
#  Makefile for Libraries
#  Some systems only need the stuff in dc_lapack.f (Divide and Conquer)
#  TranSiesta needs zgesv_lapack.f
#
default: libdummy_pexsi.a
#
ARCH_MAKE?=../arch.make
include $(ARCH_MAKE)
#
ALL_OBJS=raw_inertiacount.o inertiacount.o solver.o local_dos.o
RAW_ONLY_OBJS=raw_inertiacount.o
#
libdummy_pexsi.a: $(ALL_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru libdummy_pexsi.a $(ALL_OBJS)
	$(RANLIB) libdummy_pexsi.a
	mv libdummy_pexsi.a ..
#
# This for those without an updated interface that includes the
# raw inertia counts routine
#
libdummy_raw_pexsi.a: $(RAW_ONLY_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru libdummy_raw_pexsi.a $(RAW_ONLY_OBJS)
	$(RANLIB) libdummy_raw_pexsi.a
	mv libdummy_raw_pexsi.a ..
#
clean:
	rm -f $(ALL_OBJS) *.a
#
.c.o:
	$(CC) -c $(CFLAGS) $<




