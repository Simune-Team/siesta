#
# Makefile for ATM
# ---------------------------------------------------------------
#
FC=f77
LDR=f77
#
# EXTRA: Options that must be consistent for the compiler AND linker
# eg: -g for debugging,  -p or -pg for profiling...
#
EXTRA=  
FFLAGS=  $(EXTRA)
LDFLAGS= $(EXTRA) 
#
# System differences are hopefully accounted for in the "aux???.f file
#
AUX=osf
#
# scilib.f emulates Cray's SCILIB. In particular, the LINPACK routines
# have single precision names but are actually double precision 
# implementations.
#
AUX_OBJS= aux$(AUX).o
#
SRCS=   atm.f  nucl_z.f   difnrl.f   difrel.f   dmixp.f   dsolv1.f   \
        pseudo.f  wf.f  hsc.f  excorr.f \
        tm2.f  ker.f  wrapup.f   splift.f   spliq.f  \
        velect.f   vionic.f   wtrans.f   v0pp.f  chg_mism.f  \
        dsolv2.f   etotal.f   ext.f  input.f   orban.f   \
        potran.f   potrv.f   potrw.f   prdiff.f  string.f  root.f  \
        genrhs.f  logder.f ode.f denplot.f scilib.f  recipes.f \
	atomxc.f  ldaxc.f ggaxc.f pbexc.f pw92xc.f pzxc.f \
	pw92c.f exchng.f compat_params.f
#
OBJS=$(SRCS:.f=.o)
#
atm:  $(OBJS) $(AUX_OBJS)   
	$(LDR) $(LDFLAGS) -o atm $(OBJS) $(AUX_OBJS) $(LIBS)
#
clean:	
	rm -f *.o atm
#
atm.o: banner.h radial.h orbital.h param.h charge.h elecpot.h energy.h
#
chg_mism.o: radial.h nonlinear.h linear.h coeffs.h
compat_params.o: compat.h
denplot.o:  radial.h charge.h
difnrl.o: radial.h param.h ion.h elecpot.h
difrel.o: radial.h param.h ion.h elecpot.h
dsolv1.o: radial.h orbital.h ion.h charge.h elecpot.h energy.h
dsolv2.o: radial.h orbital.h param.h ion.h charge.h elecpot.h energy.h
etotal.o: orbital.h param.h energy.h
excorr.o: radial.h param.h
genrhs.o: nonlinear.h
input.o:  radial.h orbital.h param.h compat.h
logder.o: radial.h orbital.h param.h ion.h elecpot.h energy.h ode_blk.h \
          ode_path.h
ode.o	: radial.h ode_blk.h
orban.o : radial.h param.h ion.h elecpot.h
pseudo.o: radial.h orbital.h param.h elecpot.h
tm2.o   : radial.h orbital.h ion.h param.h elecpot.h energy.h tm2_blk.h
v0pp.o  : radial.h tm2_blk.h
velect.o: radial.h orbital.h param.h charge.h elecpot.h energy.h compat.h
vionic.o: radial.h orbital.h param.h ion.h charge.h 
wf.o hsc.o ker.o : radial.h orbital.h ion.h param.h elecpot.h energy.h
wrapup.o: radial.h orbital.h param.h ion.h charge.h elecpot.h energy.h \
          compat.h
wtrans.o potran.o: plot.h






