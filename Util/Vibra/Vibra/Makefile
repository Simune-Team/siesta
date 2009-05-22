# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
# Makefile for Vibra package
#
# Pablo Ordejon, August 1998
# Modified by Alberto Garcia, February 2006
#
# Points to note, until we switch to a better building system:
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

.SUFFIXES:
.SUFFIXES: .f .o .a 
#
default: what copy_fdf fcbuild vibrator
#
# You should have compiled Siesta in the main Src directory (or
# in an architecture-specific object directory) before attempting 
# to compile the Vibra package.
#
OBJDIR=Src
include ../../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
# (Note: The following is needed only if INCFLAGS is not included
#  in the compilation rules in arch.make, but it is kept here in
#  case a crippled version of arch.make is used.)
#
# The fdf module file also needs to be present in the current directory
# Unfortunately the naming convention for the fdf module is compiler specific
# and so you might need to edit this. Should work for Intel and IBM compilers.
#
copy_fdf:
	cp -fp ../../../$(OBJDIR)/fdf.mod .
#
# Uncomment the following line for debugging support
#
FFLAGS=$(FFLAGS_DEBUG)
#
what:   
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_ARCH)"
	@echo
#
#
BUILD-SRCS= fcbuild.f chkdim.f recoor.f io.f
BUILD-OBJS=$(BUILD-SRCS:.f=.o)

VIBRA-SRCS= vibrator.f chkdim.f recoor.f hermdp.f paste.f klines.f outbands.f \
            reclat.f parse.f io.f
VIBRA-OBJS=$(VIBRA-SRCS:.f=.o)

#
# Routines using vibra.h
#
fcbuild.o vibrator.o:  vibra.h
#
FDF=../../../$(OBJDIR)/libfdf.a         
INCFLAGS:=$(INCFLAGS) -I../../../$(OBJDIR)           # To pick up fdf
#
fcbuild: $(BUILD-OBJS)
	$(FC) -o fcbuild \
	       $(LDFLAGS) $(BUILD-OBJS) $(FDF)

vibrator : $(VIBRA-OBJS) 
	$(FC) -o vibrator \
	       $(LDFLAGS) $(VIBRA-OBJS) $(FDF)
#
clean: 
	@echo "==> Cleaning object and executable files"
	rm -f fcbuild vibrator *.o *.mod



