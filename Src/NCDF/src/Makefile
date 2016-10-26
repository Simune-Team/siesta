.SUFFIXES: .f90 .o .a

# The VPATH can be different, however the tree-structure
# has to be the same in this case.
VPATH?=$(shell pwd)

# Sadly the preprocessor for fortran is a bit of a pain.
# -C is needed to disregard C-comments (the operator // is a comment!?!)
# -nostdinc is needed to not include standard includes (which are C comments)
# -E only preprocess, do not try and compile
# -P do not add line markers
CPP = $(FC) -E -P -x c

ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)

# Setup the default things, in case they haven't been set.
AR         ?= ar
ARFLAGS    ?= cru
RANLIB     ?= ranlib


LIBFDICT_DEFAULT = ../fdict/libfdict.a
LIBFDICT ?= $(LIBFDICT_DEFAULT)

.PHONY: default
default: lib

# The different object files for this library.
OBJS = nf_ncdf.o

LIB = libncdf.a

# Ensures module are in include path
ifeq ($(LIBFDICT),$(LIBFDICT_DEFAULT))
$(OBJS): INC += -I$(VPATH)/../fdict
endif

.PHONY: lib
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

SED_DEL = 's/NEWLINE/\n/g;/^$$/d;/^\!.*&/d;\
s/[[:space:]]*\#\#[[:space:]]*\([^[:space:]]*\)/\1/g;\
s/[[:space:]]*\#\([^i][^[:space:]]*\)/"\1"/g;\
s:/ /://:g;\
s/"endif"/\n\#endif/g'
.PHONY: prep
prep:
	$(VPATH)/ncdf.sh
	$(CPP) $(CPPFLAGS) $(FPPFLAGS) -I. -I$(VPATH) $(VPATH)/nf_ncdf.F90 | sed -e $(SED_DEL) > nf_ncdf.f90 #2> /dev/null

.PHONY: clean
clean:
	-rm -f $(OBJS) $(LIB) *.o *.mod tmp.F90 nf_ncdf.f90
	-rm -f ncdf_interface.inc ncdf_funcs.inc

# Dependencies
nf_ncdf.f90: | prep
nf_ncdf.o: | prep
