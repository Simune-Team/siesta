#
#  Makefile for fdf example and fdf module.
#
.SUFFIXES: .f .F .o .a  .f90 .F90
include ../arch.make
#
module:  libwxml.a
OBJFILES= m_wxml_buffer.o m_wxml_dictionary.o m_wxml_elstack.o \
          m_wxml_text.o m_wxml_core.o flib_wxml.o
#
libwxml.a: $(OBJFILES)
	ar qv libwxml.a $(OBJFILES)
	-$(RANLIB) libwxml.a
	cp libwxml.a ..
	@cp -p *.*d ..
	@echo "wxml/work.pc" >> ../work.pcl
#
#         (Note trick for Intel Fortran Compiler)
#
clean:
	rm -f *.o *.*d *.a




