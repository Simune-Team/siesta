.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
module:  libwxml.a
include $(ARCH_MAKE)
#
#
OBJFILES= m_wxml_buffer.o m_wxml_dictionary.o m_wxml_elstack.o \
          m_wxml_text.o m_wxml_core.o flib_wxml.o
MODFILES=$(OBJFILES:.o=)

#
libwxml.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libwxml.a $(OBJFILES)
	-$(RANLIB) libwxml.a
	cp libwxml.a ..
	@cp -p *.*d ..
#
clean:
	rm -f *.o *.*d *.a
# DO NOT DELETE THIS LINE - used by make depend
flib_wxml.o: m_wxml_core.o m_wxml_text.o
m_wxml_core.o: m_wxml_buffer.o m_wxml_dictionary.o m_wxml_elstack.o
