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
LDFLAGS=
AUX_OBJS= auxf95.o   # Change if you do not have
#                    Â#a fortran-90 compiler
#----------------------------------------------$$$$$
#
#
XC_OBJS= xc.o
#
# Replace definition above with line below to use the 
# old Soler-Balbas XC routines.
#
#XC_OBJS= atomxc.o  ldaxc.o ggaxc.o pbexc.o pw92xc.o pzxc.o \
#        pw92c.o exchng.o ## Not really old: rpbexc.o revpbexc.o
#
MAIN_OBJS= atm.o  nucl_z.o   difnrl.o   difrel.o   dmixp.o   dsolv1.o   \
        pseudo.o  wf.o  hsc.o  excorr.o \
        tm2.o  ker.o  wrapup.o   splift.o   spliq.o  \
        velect.o   vionic.o   wtrans.o   v0pp.o  chg_mism.o  \
        dsolv2.o   etotal.o   ext.o  input.o  header.o orban.o   \
        potran.o   potrv.o potrvs.o  potrw.o   prdiff.o  string.o  root.o  \
        genrhs.o  logder.o ode.o denplot.o scilib.o  recipes.o \
	compat_params.o  pcc_exp.o  change_valence.o \
	coreq.o get_unit.o gauleg.o divdif.o symbols.o prversion.o
#
OBJS= $(MAIN_OBJS) $(XC_OBJS)
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
prversion.o: version.h
pseudo.o: radial.h orbital.h param.h elecpot.h
symbols.o: set2.h
tm2.o   : radial.h orbital.h ion.h param.h elecpot.h energy.h tm2_blk.h
v0pp.o  : radial.h tm2_blk.h
velect.o: radial.h orbital.h param.h charge.h elecpot.h energy.h compat.h
vionic.o: radial.h orbital.h param.h ion.h charge.h 
change_valence.o: radial.h param.h  charge.h ion.h
wf.o hsc.o ker.o : radial.h orbital.h ion.h param.h elecpot.h energy.h
wrapup.o: radial.h orbital.h param.h ion.h charge.h elecpot.h energy.h \
          compat.h
wtrans.o potran.o: plot.h
ninput.o header.o: input.h
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







