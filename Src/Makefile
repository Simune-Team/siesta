# $Id: Makefile,v 1.47 1999/05/18 17:17:33 ordejon Exp $
#
# Makefile for siesta
#
.SUFFIXES: .f .o .a .prj .smb .log
#
default: what siesta
#
# System dependencies are dealt with by an include file. See Sys directory
# for examples and instructions. The variable SIESTA_SYS should be in the 
# shell environment. 
#
SIESTA_SYS?=unknown
#
include Sys/${SIESTA_SYS}.make
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
what:
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo
#
SYSOBJ=Sys/$(SYS).o
$(SYSOBJ): Sys/$(SYS).f
	(cd Sys ; $(FC) -c $(FFLAGS) $(SYS).f)
#
#  Some routines do not depend on program dimensions. They are put
#  in a library. Note that only changed .o files need to be re-inserted
#  in it.
#
SUBS =  arw.f atmass.f atomlwf.f bessph.f cgwf.f chkdim.f chkgmx.f \
	cgvc.f coceri.f conjgr.f constr.f coxmol.f cross.f \
	dfscf.f diagg.f diagk.f diag2g.f diag2k.f digcel.f  \
	dipole.f dismin.f dnaefs.f dot.f dynamics.f \
	efield.f ener3.f extrapol.f \
	fermid.f fermispin.f fixed.f gauleg.f grdsam.f hsparse.f \
	idiag.f initdm.f inver.f iocg.f iodm.f ioeig.f iofa.f iohs.f \
	iokp.f iolwf.f iomd.f iorder.f iorho.f ioxv.f \
	kgrid.f kinefsm.f listsc.f lofilm.f madelung.f\
	minvec.f mulliken.f naefs.f neighb.f nlefsm.f \
	ofc.f on_subs.f order.f ordvec.f outcell.f outcoor.f overfsm.f \
	parse.f paste.f prdate.f prmem.f propor.f prversion.f pulayx.f  \
	pixmol.f radfft.f ran3.f recipes.f reclat.f recoor.f redata.f \
	redcel.f remass.f reord.f rhoofd.f rhooda.f rlylm.f \
	shaper.f superc.f superx.f timer.f transp.f typecell.f \
	vmb.f vlist.f vmat.f volcel.f xc.f xijorb.f ylmexp.f ylmylm.f
#
LIBOBJS=$(SUBS:.f=.o)
libsiesta.a: $(LIBOBJS)
	ar rv $@ $?
	$(RANLIB) $@
#
#  These modules are dimensioned.
#
DIM_SRCS= atom.f atom_subs.f atom_functions.f bands.f \
          dhscf.f diagon.f ordern.f chempot.f denmat.f \
          egandd.f gradient.f initatom.f matel.f \
          plcharge.f ranger.f redbasis.f redbasis_subs.f siesta.f
DIM_OBJS=$(DIM_SRCS:.f=.o)
#
atom.o           : atom.f      atom.h
plcharge.o       : plcharge.f  atom.h
atom_subs.o      : atom_subs.f atom.h
atom_functions.o : atom_functions.f atom.h 
initatom.o       : initatom.f atom.h
redbasis.o       : redbasis.f atom.h
redbasis_subs.o  : redbasis_subs.f atom.h
dhscf.o          : dhscf.f     dhscf.h
diagon.o         : diagon.f    diagon.h
bands.o          : bands.f     diagon.h
ordern.o         : ordern.f    ordern.h
chempot.o        : chempot.f   ordern.h
denmat.o         : denmat.f    ordern.h
egandd.o         : egandd.f    ordern.h
gradient.o       : gradient.f  ordern.h
ranger.o         : ranger.f    ranger.h
matel.o          : matel.f     matel.h
siesta.o         : siesta.f    siesta.h
#
# Routines using fdf calls.
#
atom.o dhscf.o initdm.o iodm.o iohs.o iolwf.o iorho.o grdsam.o : fdf/fdfdefs.h
plcharge.o recoor.o redata.o remass.o outcoor.o ioeig.o kgrid.o : fdf/fdfdefs.h
siesta.o redbasis.o redbasis_subs.o iomd.o ioxv.o iokp.o iofa.o : fdf/fdfdefs.h
cgvc.o iocg.o ofc.o : fdf/fdfdefs.h
#
FDF= fdf/fdf.o fdf/io.o
$(FDF): fdf/fdf.f fdf/fdf.h fdf/io.f
	(cd fdf ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" fdf.o io.o)
#
prversion.o: version.h
#
siesta: $(DIM_OBJS) $(SYSOBJ) $(FDF) libsiesta.a
	$(FC) -o siesta \
	       $(LDFLAGS) $(DIM_OBJS) $(SYSOBJ) $(FDF) libsiesta.a $(LIBS)
#
# Basic recipe to compile
#
.f.o:
	$(FC) -c $(FFLAGS) $<
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta *.o libsiesta.a
	rm -f Sys/*.o
	(cd fdf ; make clean)
pristine:  clean
	@echo "==> Copying original .h files from Include/"
	cp -fp Include/*.h .
	cp -fp Include/constr.f .
setup:
	@echo "==> Copying standard .h files from Include/"
	cp -fp Include/*.h .
#
# FTNCHEK and plusFORT section
#
FILES=$(SUBS) $(DIM_SRCS)
#
# set up macro PRJS containing project filenames
#
PRJS= $(FILES:.f=.prj)
SMBS= $(FILES:.f=.smb)
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
	rm -f *.log *.smb *.spg gxchk.sav gxchk.out
#
veryclean:  clean clean_check





