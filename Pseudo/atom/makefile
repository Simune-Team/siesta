#
# Makefile for ATM
# ---------------------------------------------------------------
#
.SUFFIXES: .f .o .a .prj .smb .log
#
#----------------------------------------------$$$$$
# EDIT the following lines to suit your system
#
FC= f77
FFLAGS= -O
AUX_OBJS= aux_general.o
#----------------------------------------------$$$$$
#
SRCS=   atm.f  nucl_z.f   difnrl.f   difrel.f   dmixp.f   dsolv1.f   \
        pseudo.f  wf.f  hsc.f  excorr.f \
        tm2.f  ker.f  wrapup.f   splift.f   spliq.f  \
        velect.f   vionic.f   wtrans.f   v0pp.f  chg_mism.f  \
        dsolv2.f   etotal.f   ext.f  input.f   orban.f   \
        potran.f   potrv.f   potrw.f   prdiff.f  string.f  root.f  \
        genrhs.f  logder.f ode.f denplot.f scilib.f  recipes.f \
	atomxc.f  ldaxc.f ggaxc.f pbexc.f pw92xc.f pzxc.f \
	pw92c.f exchng.f compat_params.f blypxc.f
#
OBJS=$(SRCS:.f=.o)
#
atm:  $(OBJS) $(AUX_OBJS)   
	$(FC) $(LDFLAGS) -o atm $(OBJS) $(AUX_OBJS) $(LIBS)
#
clean:	
	rm -f *.o atm
#
atm.o: radial.h orbital.h param.h charge.h elecpot.h energy.h
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
#
# FTNCHEK and FORCHECK section
#
FILES=$(OBJS) $(AUX_OBJS)
#
# set up macro PRJS containing project filenames
#
PRJS= $(FILES:.o=.prj)
SMBS= $(FILES:.o=.smb)
#
.f.smb:
	-spag log=$*.log $<
#
chk:    $(SMBS)
	-gxchk $(SMBS)
#

# tell make how to create a .prj file from a .f file
#
.f.prj:
	ftnchek -project -noextern -library -array=0 \
		-quiet -sixchar  -output=$*.ana $<

# "make check" will check everything that has been changed.
check: $(PRJS)
	ftnchek -output=Global.chk -quiet $(PRJS)
#
clean_check:
	rm -f *.prj *.ana Global.chk
clean_chk:
	rm -f *.log *.smb *.spg Global.chk gxchk.out gxchk.sav
#
veryclean:  clean clean_check clean_chk







