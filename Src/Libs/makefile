# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.


default: libsiestaLAPACK.a libsiestaBLAS.a

ARCH_MAKE?=../arch.make
include $(ARCH_MAKE)

###########
# These are now the preferred complete BLAS and LAPACK
# sources.
# They are extracted using the linalg2file.py script.
#
# See create.sh for how the lapack.f and blas.f files are created.
#
# Note that the python script is extremely basic, but functions well for
# F77 code.

libsiestaBLAS.a: blas.o
	$(AR) $(ARFLAGS_EXTRA) cru $@ $<
	$(RANLIB) $@

libsiestaLAPACK.a: lapack.o
	$(AR) $(ARFLAGS_EXTRA) cru $@ $<
	$(RANLIB) $@


########################
#
# Below are old (pre 4.1) temporary libraries for
# linear algebra.
# They should probably not be used as it is better
# to use a consistent library with all functionality.
#
# We advice on using the newer
#   libsiestaLAPACK.a
# and/or
#   libsiestaBLAS.a
#
########################


ALL_OBJS = dc_lapack.o zgesv_lapack.o zgeev_lapack.o \
           svd_omm_lapack.o old_lapack.o old_blas.o machine.o

# Add machine.o here only if needed...
DC_OBJS =  dc_lapack.o # machine.o
# Used in TranSIESTA
ZGESV_OBJS =  zgesv_lapack.o # machine.o
# Used in TBTrans
ZGEEV_OBJS = zgeev_lapack.o

ESSL_OBJS = essl_lapack.o

libblas.a: old_blas.o
	$(AR) $(ARFLAGS_EXTRA) cru libblas.a old_blas.o
	$(RANLIB) libblas.a
	mv libblas.a ..

liblapack.a: old_lapack.o machine.o
	$(AR) $(ARFLAGS_EXTRA) cru liblapack.a old_lapack.o machine.o
	$(RANLIB) liblapack.a
	mv liblapack.a ..

linalg.a: $(ALL_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru linalg.a $(ALL_OBJS)
	$(RANLIB) linalg.a
	mv linalg.a ..

dc_lapack.a: $(DC_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru dc_lapack.a $(DC_OBJS)
	$(RANLIB) dc_lapack.a
	mv dc_lapack.a ..

zgesv_lapack.a: $(ZGESV_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru zgesv_lapack.a $(ZGESV_OBJS)
	$(RANLIB) zgesv_lapack.a
	mv zgesv_lapack.a ..

zgeev_lapack.a: $(ZGEEV_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru zgeev_lapack.a $(ZGEEV_OBJS)
	$(RANLIB) zgeev_lapack.a
	mv zgeev_lapack.a ..

essl_lapack.a: $(ESSL_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru essl_lapack.a $(ESSL_OBJS)
	$(RANLIB) essl_lapack.a
	mv essl_lapack.a ..

clean:
	rm -f $(ALL_OBJS) *.a




