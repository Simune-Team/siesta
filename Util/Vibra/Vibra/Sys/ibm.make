#
FC=f77
# To trap floating exceptions:(normally it runs quicker but with possible NaNs)
# FFLAGS= -O -qmaxmem=6144 -qtkq_size=25000 -qflttrap=ov:und:zero:en
FFLAGS= -O -qmaxmem=6144 -qtkq_size=25000
FFLAGS_DEBUG= -g -qmaxmem=6144 -qtkq_size=25000
LIBS= 
SYS=ibm
RANLIB=ranlib
#
