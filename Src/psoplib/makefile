.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
#
default:  libpsop.a
module:  libpsop.a
#
# We allow for an external specification of the location 
# of the arch.make file. This must be an absolute path
#
#
PSOP_MK_DEFAULT=fortran.mk
PSOP_MK?=$(PSOP_MK_DEFAULT)
include $(PSOP_MK)
#
OBJFILES= m_psop.o m_kbgen.o m_localgen.o schrodinger.o arw.o \
          psop_params.o  flib_spline.o m_interpol.o
#
INCFLAGS=
#
libpsop.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libpsop.a $(OBJFILES) 
	-$(RANLIB) libpsop.a
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.*d *.a
#
# DO NOT DELETE THIS LINE - used by make depend
m_kbgen.o: psop_params.o schrodinger.o
m_localgen.o: flib_spline.o m_interpol.o
m_psop.o: m_kbgen.o m_localgen.o psop_params.o
schrodinger.o: psop_params.o
schrodinger_m.o: schrodinger.o
