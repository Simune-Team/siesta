#
.SUFFIXES:
.SUFFIXES: .f .o .f90
#
# This can also be achieved with a "FORCE" empty target.      
# xml: FORCE
#	(cd ../../Src/xmlparser ; make)
# FORCE:
#
#
.PHONY: xml clean
#
default: xml2psf
#
include ../../Src/arch.make
#
FFLAGS=$(FFLAGS_DEBUG) -I../../Src/xmlparser
#
LIBS= ../../Src/xmlparser/libxmlparser.a
#
OBJS= sys.o io.o pseudopotential.o m_pseudo_types.o m_pseudo.o xml2psf_helper.o  xml2psf.o
#
xml2psf: xml $(OBJS)
	$(FC) -o xml2psf $(OBJS) $(LIBS)
#
xml:
	(cd ../../Src/xmlparser ; make)
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.mod xml2psf
#
# DO NOT DELETE THIS LINE - used by make depend
io.o: sys.o
m_pseudo.o: m_pseudo_types.o
xml2psf.o: m_pseudo.o m_pseudo_types.o pseudopotential.o
xml2psf_helper.o: m_pseudo_types.o pseudopotential.o
