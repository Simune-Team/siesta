# $Id: ibmesslp3.make,v 1.1 1999/02/21 16:13:41 emilio Exp $
#
FC=xlf
#
# To trap floating exceptions:(normally it runs quicker but with possible NaNs)
# FFLAGS= -O -qmaxmem=6144 -qtkq_size=25000 -qflttrap=ov:und:zero:en
#
FFLAGS= -O4 -qmaxmem=16000 -qtkq_size=25000 -qnolm
FFLAGS_DEBUG= -g -qmaxmem=6144 -qtkq_size=25000
#
# -qipa gives speed optimization of a few percent, but takes long to link
#
#LIBS= -Wl,-bD:912000000 -qipa -qnolm -lessl
LIBS= -Wl,-bD:2000000000  -lessl -qnolm
SYS=ibmessl
RANLIB=ranlib
#
