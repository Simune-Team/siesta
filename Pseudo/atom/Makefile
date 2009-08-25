#
# Makefile for ATM
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
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#
OBJDIR=Obj
#
include ../../$(OBJDIR)/arch.make
#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
#----------------------------------------------$$$$$
# EDIT the following lines to suit your system
#
AUX_OBJS= auxf95.o
LIBS= $(EXTRA_LIBS)
INCFLAGS=
#--------------------------------------------------------------
# FoX whatnot
#
# FoX build targets:
FoX_configured=../../$(OBJDIR)/FoX/.config
FoX_built=../../$(OBJDIR)/FoX/.FoX
#
# This is how we pick up modules and libraries for FoX:
FoX_FCFLAGS=`../../$(OBJDIR)/FoX/FoX-config --fcflags`
FoX_LIBS=`../../$(OBJDIR)/FoX/FoX-config --libs --wxml`
#
# And add them to global compiler flags:
INCFLAGS=$(FoX_FCFLAGS)
#
# First, it needs to be configured. This may have been done
# by the main SIESTA configure, but in case not:
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
#
# Note that we hardwire the need to use the configure script in the main Src directory.
#
TOP_LEVEL=$(shell pwd)/../../Src
$(FoX_configured):
	(cd ../../$(OBJDIR)/FoX; CONFIGURE="$(TOP_LEVEL)/FoX/configure"; $$CONFIGURE VPATH="$(TOP_LEVEL)/FoX" FC="$(FC_SERIAL)" FCFLAGS="$(FFLAGS)" --enable-wcml $(DUMMY_FOX) || false )
#
# Note we have to make sure to use the same compiler as SIESTA,
# and pick up all the same FFLAGS, and also remember to maybe ask for a dummy library.
# Note also that the "false" clause will stop the 'make' process if the configuration fails.
#
# Then we want to make FoX itself. Like so:
$(FoX_built): $(FoX_configured)
	(cd ../../$(OBJDIR)/FoX; $(MAKE))
#
# Now ensure that building FoX is inserted into the dependency chain.
# Routines using FoX:
#
pseudoXML.o: $(FoX_built)
#
#--------------------------------------------------------------
dep:
	sfmakedepend --depend=obj  --modext=o *.f *.f90 *.F *.F90
#----------------------------------------------$$$$$
#
INCFLAGS:=$(INCFLAGS)  ##### -I../../$(OBJDIR)
#
#
XC_OBJS= precision.o alloc.o radfft.o m_debug.o flib_spline.o m_filter.o \
         xcmod.o  atomxc.o ldaxc.o ggaxc.o vdwxc.o mesh1d.o \
         m_recipes.o memory.o memoryinfo.o sys.o timer.o parallel.o bessph.o io.o
#
#
MAIN_OBJS= atm.o  nucl_z.o   difnrl.o   difrel.o   dmixp.o   dsolv1.o   \
        pseudo.o  wf.o  hsc.o  excorr.o \
        tm2.o  ker.o  wrapup.o   splift.o   spliq.o  \
        velect.o   vionic.o   wtrans.o   v0pp.o  chg_mism.o  \
        dsolv2.o   etotal.o   ext.o  input.o  header.o orban.o   \
        potran.o   potrv.o potrvs.o  potrw.o   prdiff.o  string.o  root.o  \
        genrhs.o  logder.o ode.o denplot.o scilib.o  recipes.o \
	compat_params.o  pcc_exp.o  change_valence.o pswf_store.o \
	coreq.o get_unit.o gauleg.o divdif.o symbols.o \
        prversion.o pseudoXML.o pxf.o force_underflow.o
#
OBJS= $(MAIN_OBJS) $(XC_OBJS) 
#
provide_xml:
	@echo "Make sure that the wxml library is already compiled."
	@echo "Otherwise the following command might fail"
	@echo
#
pseudoXML.o: provide_xml
#
atm:	$(FoX_built) $(OBJS) $(AUX_OBJS)
	$(FC) $(LDFLAGS) -o atm $(OBJS) $(AUX_OBJS) $(LIBS) $(FoX_LIBS)

#
clean:	
	rm -f *.o *.mod atm
#
# DO NOT DELETE THIS LINE - used by make depend
alloc.o: parallel.o precision.o sys.o
atm.o: charge.h param.h orbital.h radial.h energy.h elecpot.h
atomxc.o: alloc.o m_filter.o mesh1d.o precision.o radfft.o sys.o vdwxc.o
atomxc.o: xcmod.o
bessph.o: parallel.o precision.o sys.o
change_valence.o: charge.h param.h radial.h ion.h
chg_mism.o: nonlinear.h coeffs.h radial.h linear.h
compat_params.o: compat.h
coreq.o: charge.h radial.h
denplot.o: param.h charge.h radial.h
difnrl.o: param.h radial.h ion.h elecpot.h
difrel.o: param.h radial.h ion.h elecpot.h
dsolv1.o: charge.h orbital.h radial.h energy.h ion.h elecpot.h
dsolv2.o: charge.h param.h orbital.h radial.h energy.h ion.h elecpot.h
dummy_pseudoXML.o: charge.h param.h orbital.h radial.h pseudowave.h ion.h
etotal.o: param.h orbital.h energy.h
excorr.o: param.h radial.h
flib_spline.o: precision.o
genrhs.o: nonlinear.h
ggaxc.o: precision.o sys.o
header.o: param.h orbital.h input.h compat.h radial.h version.h
hsc.o: param.h orbital.h radial.h energy.h ion.h elecpot.h
input.o: param.h orbital.h input.h compat.h radial.h version.h
io.o: sys.o
ker.o: param.h orbital.h radial.h energy.h ion.h elecpot.h
ldaxc.o: precision.o sys.o
logder.o: ode_blk.h param.h ode_path.h orbital.h radial.h energy.h ion.h
logder.o: elecpot.h
m_debug.o: parallel.o
m_recipes.o: sys.o
memory.o: memoryinfo.o parallel.o precision.o
memoryinfo.o: precision.o
mesh1d.o: m_recipes.o precision.o
nv0pp.o: coeffs.h
ode.o: ode_blk.h radial.h
orban.o: param.h radial.h ion.h elecpot.h
potran.o: plot.h
prversion.o: version.h
pseudo.o: param.h orbital.h radial.h elecpot.h
pseudoXML.o: charge.h param.h orbital.h radial.h pseudowave.h ion.h
pswf_store.o: radial.h pseudowave.h ion.h
radfft.o: alloc.o m_recipes.o precision.o
recipes.o: ode_path.h
symbols.o: set2.h
sys.o: parallel.o
tm2.o: param.h orbital.h radial.h energy.h ion.h tm2_blk.h elecpot.h
v0pp.o: radial.h tm2_blk.h
vdwxc.o: flib_spline.o m_debug.o m_recipes.o mesh1d.o precision.o radfft.o
velect.o: charge.h param.h orbital.h compat.h radial.h energy.h elecpot.h
velect.o: precision.o xcmod.o
vionic.o: charge.h param.h orbital.h radial.h ion.h
wf.o: param.h orbital.h radial.h energy.h ion.h elecpot.h
wrapup.o: param.h orbital.h radial.h energy.h version.h elecpot.h charge.h
wrapup.o: compat.h pseudowave.h ion.h
wtrans.o: plot.h
xcmod.o: precision.o
diagmemory.o: memoryinfo.o
m_radfft.o: radfft.o
m_vdwxc.o: vdwxc.o
