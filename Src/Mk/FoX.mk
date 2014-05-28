# FoX helper snippets
#
ifeq ($(strip $(FOX_ROOT)),)
# Compile in source...
FOX_OBJ_ROOT=$(OBJDIR)/FoX

# FoX whatnot
#
# FoX build targets:
FoX_configured=$(FOX_OBJ_ROOT)/.config
FoX_built=$(FOX_OBJ_ROOT)/.FoX
#
# This is how we pick up modules and libraries for FoX:
FoX_FCFLAGS=`$(FOX_OBJ_ROOT)/FoX-config --fcflags`
FOX_LIBS=`$(FOX_OBJ_ROOT)/FoX-config --libs --wcml`
#
# And add them to global compiler flags:
INCFLAGS:=$(INCFLAGS) $(FoX_FCFLAGS)
#
# First, it needs to be configured. This may have been done
# by the main SIESTA configure, but in case not:
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
#
$(FoX_configured):
	(cd $(FOX_OBJ_ROOT); touch arch.make ; \
         CONFIGURE="$(VPATH)/FoX/configure"; \
         $$CONFIGURE VPATH="$(VPATH)/FoX" \
         FC="$(FC_SERIAL)" FCFLAGS="$(FFLAGS:$(IPO_FLAG)=)" \
         --enable-wcml $(DUMMY_FOX) || false )
#
# Note we have to make sure to use the same compiler as SIESTA,
# and pick up all the same FFLAGS, and also remember to maybe ask for a dummy library.
# Note also that the "false" clause will stop the 'make' process if the configuration fails.
#
# Then we want to make FoX itself. Like so:
$(FoX_built): $(FoX_configured)
	(cd $(FOX_OBJ_ROOT); $(MAKE) )

else

#
#  Pre-compiled library
#

FOX_CONFIG=$(FOX_ROOT)/bin/FoX-config
#
FOX_INCFLAGS=`$(FOX_CONFIG) --fcflags`
INCFLAGS:= $(INCFLAGS) $(FOX_INCFLAGS)
#
FOX_LIBS=`$(FOX_CONFIG) --libs`
#
endif