#
# Makefile for MatrixSwitch in Siesta style
#
.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

LIBNAME=MatrixSwitch.a
all: $(LIBNAME)

ARCH_MAKE_DEFAULT=../../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)



# Add -DSLAP for MPI compilations only
#
EXTFLAGS=$(FPPFLAGS)
#
# Use the override command in case we still want
# to pass FPPFLAGS in the command line. (Now it
# is superfluous since we include the arch.make file)
#
MPI_STRING=$(findstring -DMPI,$(EXTFLAGS))
ifneq (,$(MPI_STRING))
override FPPFLAGS += -DLAP -DSLAP -DCONV
else
override FPPFLAGS += -DLAP -DCONV
endif
#
#
ALLOBJ = MatrixSwitch.o

module: $(LIBNAME)
	cp $(LIBNAME) ../..
	@cp -p *.mod ../..

$(LIBNAME) : $(ALLOBJ)
	@echo "Incoming FPPFLAGS: " $(EXTFLAGS)
	@echo "Current FPPFLAGS: " $(FPPFLAGS)
	$(AR) $(ARFLAGS_EXTRA) cru $@ $(ALLOBJ)
	-$(RANLIB) $@

# Make sure that you include this in the "clean" rule
# in the top makefile. Otherwise these files will
# stay here and confuse the building logic, since the
# "copy" operation above is not a real target.
#
clean : 
	rm -f *.a *.o *.mod

