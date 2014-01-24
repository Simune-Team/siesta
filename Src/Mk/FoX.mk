# FoX helper snippets
#
FOX_CONFIG=$(FOX_ROOT)/bin/FoX-config
#
FOX_INCFLAGS=`$(FOX_CONFIG) --fcflags`
INCFLAGS:= $(INCFLAGS) $(FOX_INCFLAGS)
#
FOX_LIBS=`$(FOX_CONFIG) --libs`

