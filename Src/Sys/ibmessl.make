# $Id: ibmessl.make,v 1.5 1999/01/31 12:05:53 emilio Exp $
#
FC=f77
# To trap floating exceptions:(normally it runs quicker but with possible NaNs)
# FFLAGS= -O -qmaxmem=6144 -qtkq_size=25000 -qflttrap=ov:und:zero:en
FFLAGS= -O -qmaxmem=6144 -qtkq_size=25000
FFLAGS_DEBUG= -g -qmaxmem=6144 -qtkq_size=25000
LIBS= -lessl
SYS=ibmessl
RANLIB=ranlib
#
