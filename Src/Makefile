# Makefile for siesta
#
# J.M.Soler. November 1996, June 1997
# Alberto Garcia, December 1996
# Pablo Ordejon, December 1996
#
# Restructured by Alberto Garcia, July 1997
#
.SUFFIXES: .f .o .a .prj
#
default: what siesta
#
# System dependencies are dealt with by an include file. See Sys directory
# for examples and instructions. The variable ARCH should be in the 
# shell environment. 
#
ARCH?=unknown
#
include Sys/${ARCH}.make
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
what:
	@echo
	@echo "Compilation architecture to be used:  $(ARCH)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable ARCH in your shell environment."
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
SUBS =  atmass.f atomlwf.f bessph.f cgwf.f chkdim.f chkgmx.f \
	coceri.f conjgr.f constr.f coxmol.f cross.f \
	dfscf.f dismin.f dnaefs.f dot.f dynamics.f efield.f ener3.f \
	extrapol.f fermid.f fixed.f four1.f gauleg.f hsparse.f \
	idiag.f initatom.f initdm.f inver.f iodm.f iolwf.f iorder.f \
	kinefsm.f lofilm.f \
	minvec.f mixing.f mulliken.f naefs.f neighb.f \
	on_subs.f order.f ordix.f ordvec.f outcell.f outcoor.f overfsm.f \
        parse.f prmem.f propor.f prversion.f recoor.f radfft.f ran3.f \
	reclat.f redata.f remass.f reord.f rhoofd.f rhooda.f rlylm.f \
	splin.f timer.f transp.f vmb.f vlist.f vmat.f volcel.f \
	xc.f ylmexp.f ylmylm.f \
	cdiag.f dipole.f iohs.f ioxv.f iorho.f kgrid.f nlefsm.f \
	prdate.f shaper.f xijorb.f pulayx.f
#
LIBOBJS=$(SUBS:.f=.o)
libsiesta.a: $(LIBOBJS)
	ar rv $@ $?
	$(RANLIB) $@
#
#  These modules are dimensioned.
#
DIM_SRCS= atom.f atom_subs.f bands.f chcore.f epskb.f \
          phiatm.f psover.f rcore.f rcut.f uion.f \
          dhscf.f diagon.f ordern.f denmat.f \
          egandd.f gradient.f ranger.f matel.f siesta.f
DIM_OBJS=$(DIM_SRCS:.f=.o)
#
atom.o      : atom.f      atom.h
atom_subs.o : atom_subs.f atom.h
chcore.o    : chcore.f    atom.h
epskb.o     : epskb.f     atom.h
phiatm.o    : phiatm.f    atom.h
psover.o    : psover.f    atom.h
rcore.o     : rcore.f     atom.h
rcut.o      : rcut.f      atom.h
uion.o      : uion.f      atom.h
dhscf.o     : dhscf.f     dhscf.h
diagon.o    : diagon.f    diagon.h
bands.o     : bands.f     diagon.h
ordern.o    : ordern.f    ordern.h
denmat.o    : denmat.f    ordern.h
egandd.o    : egandd.f    ordern.h
gradient.o  : gradient.f  ordern.h
ranger.o    : ranger.f    ranger.h
matel.o     : matel.f     matel.h
siesta.o    : siesta.f    siesta.h
#
# Routines using fdf calls.
#
atom.o dhscf.o initdm.o iodm.o iohs.o iolwf.o iorho.o : fdf/fdfdefs.h
recoor.o redata.o remass.o outcoor.o ioxv.o kgrid.o siesta.o : fdf/fdfdefs.h
#
FDF= fdf/fdf.o fdf/io.o
$(FDF):
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
#
# FTNCHEK section
#
FILES=$(SUBS) $(DIM_SRCS)
#
# set up macro PRJS containing project filenames
#
PRJS= $(FILES:.f=.prj)

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
veryclean:  clean clean_check



