#
# Makefile for ATM
#
# Please note that ATM now depends on the SiestaXC and FoX libraries
# for its correct compilation. It is not currently possible to have
# a "standalone" version independent of the rest of the Siesta package.
# ---------------------------------------------------------------
#
.SUFFIXES:
.SUFFIXES: .F90 .f .F .f90 .o .a 
#
default: atm
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Src, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
OBJDIR=Obj
#
#
MAIN_OBJDIR:=$(shell cd -P -- ../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
# We support both DEFS and FPPFLAGS (crazy -- deprecate DEFS soon)
DEFS:=$(DEFS_PREFIX)-UMPI $(DEFS) $(DEFS_PREFIX)-UMPI 
FPPFLAGS:=$(DEFS_PREFIX)-UMPI $(FPPFLAGS) $(DEFS_PREFIX)-UMPI 
#
#----------------------------------------------$$$$$
# EDIT the following lines to suit your system
#
AUX_OBJS= auxf95.o
LIBS= $(EXTRA_LIBS)
INCFLAGS=
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
#
TOP_LEVEL=$(shell pwd)/../../Src
#
pseudoXML.o: $(WXML)
write_psml.o: $(WXML)
#
#--------------------------------------------------------------
dep:
	sfmakedepend --depend=obj  --modext=o *.f *.f90 *.F *.F90
#----------------------------------------------$$$$$
#
WXML=libwxml.a
#
$(WXML): 
	(cd wxml ; $(MAKE) \
               FFLAGS="$(FFLAGS:$(IPO_FLAG)=)" \
               "ARCH_MAKE=$(ARCH_MAKE)"  module)
#
# Compile XC in a private directory, without MPI
##
XC=libSiestaXC.a
XC_MAKEFILE=$(TOP_LEVEL)/SiestaXC/makefile
XC_INCFLAGS:=  #   -I $(TOP_LEVEL)/SiestaXC ##$(INCFLAGS)
$(XC): 
	(mkdir -p SiestaXC; cd SiestaXC ; $(MAKE) -f $(XC_MAKEFILE) \
                          "FC=$(FC)" "VPATH_ROOT=$(TOP_LEVEL)" \
                          "ARCH_MAKE=$(ARCH_MAKE)" \
                          "INCFLAGS=$(XC_INCFLAGS)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "FFLAGS=$(FFLAGS)"  module )
#
#
MAIN_OBJS= atm.o  nucl_z.o   difnrl.o   difrel.o   dmixp.o   dsolv1.o   \
        pseudo.o  wf.o  hsc.o  excorr.o \
        tm2.o  ker.o  wrapup.o   splift.o   spliq.o  \
        velect.o   vionic.o   wtrans.o   v0pp.o  chg_mism.o  \
        dsolv2.o  etotal.o   ext.o  input.o  header.o orban.o   \
        potran.o   potrv.o potrvs.o  potrw.o   prdiff.o  string.o  root.o  \
        genrhs.o  logder.o ode.o denplot.o scilib.o  recipes.o \
	compat_params.o  pcc_exp.o  change_valence.o pswf_store.o \
	coreq.o get_unit.o gauleg.o divdif.o symbols.o \
        prversion.o pseudoXML.o write_psml.o pxf.o force_underflow.o dsolv2_single.o dsolv2_save_wf.o
#
OBJS= $(MAIN_OBJS) 
#
pseudoXML.o: $(WXML)
write_psml.o: $(WXML)
#
velect.o: $(XC)
atm:	$(WXML) $(XC) $(OBJS) $(AUX_OBJS) 
	$(FC) $(LDFLAGS) -o atm $(OBJS) $(AUX_OBJS) $(XC) $(LIBS) $(WXML)

fhi2psml: $(WXML) fhi2psml.o
	$(FC) $(LDFLAGS) -o fhi2psml fhi2psml.o $(WXML)

#
clean:	
	rm -f *.o *.mod atm *.a
	rm -rf SiestaXC
	(cd wxml ; $(MAKE) clean)
#
# DO NOT DELETE THIS LINE - used by make depend
atm.o: charge.h param.h orbital.h radial.h energy.h elecpot.h
change_valence.o: charge.h param.h radial.h ion.h
chg_mism.o: nonlinear.h coeffs.h radial.h linear.h
compat_params.o: compat.h
coreq.o: charge.h radial.h
denplot.o: param.h charge.h radial.h
difnrl.o: param.h radial.h ion.h elecpot.h
difrel.o: param.h radial.h ion.h elecpot.h
dsolv1.o: charge.h orbital.h radial.h energy.h ion.h elecpot.h
dsolv2.o: charge.h param.h orbital.h radial.h energy.h ion.h elecpot.h
dsolv2_save_wf.o: charge.h param.h orbital.h radial.h energy.h ion.h elecpot.h
dsolv2_single.o: charge.h param.h orbital.h radial.h energy.h ion.h elecpot.h
dummy_pseudoXML.o: charge.h param.h orbital.h radial.h pseudowave.h ion.h
etotal.o: param.h orbital.h energy.h
excorr.o: param.h radial.h
genrhs.o: nonlinear.h
header.o: param.h orbital.h input.h compat.h radial.h version.h
hsc.o: param.h orbital.h radial.h energy.h ion.h elecpot.h
input.o: param.h orbital.h input.h compat.h radial.h version.h
ker.o: param.h orbital.h radial.h energy.h ion.h elecpot.h
logder.o: param.h orbital.h radial.h energy.h elecpot.h ode_blk.h ode_path.h
logder.o: ion.h
nv0pp.o: coeffs.h
ode.o: ode_blk.h radial.h
orban.o: param.h radial.h ion.h elecpot.h
potran.o: plot.h
prversion.o: version.h
pseudo.o: param.h orbital.h radial.h elecpot.h
pseudoXML.o: charge.h param.h orbital.h radial.h pseudowave.h ion.h
pswf_store.o: radial.h pseudowave.h ion.h
recipes.o: ode_path.h
symbols.o: set2.h
tm2.o: param.h orbital.h radial.h energy.h ion.h tm2_blk.h elecpot.h
v0pp.o: radial.h tm2_blk.h
velect.o: charge.h param.h orbital.h compat.h radial.h energy.h elecpot.h
vionic.o: charge.h param.h orbital.h radial.h ion.h
wf.o: param.h orbital.h radial.h energy.h ion.h elecpot.h
wrapup.o: param.h corecorr.h orbital.h radial.h energy.h version.h elecpot.h
wrapup.o: charge.h compat.h pseudowave.h ion.h
write_psml.o: charge.h param.h corecorr.h orbital.h radial.h pseudowave.h ion.h
wtrans.o: plot.h
