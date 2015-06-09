.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: libxc-compat
#
default:  libncps.a
module:  libncps.a
#
# We allow for an external specification of the location 
# of the arch.make file. This must be an absolute path
#
#
ARCH_MAKE_DEFAULT=fortran.mk
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
VPATH=.
#
libxc-compat:
	@echo "Compiling the libxc-compat library..."
	(cd libxc-compat ; $(MAKE) \
            ARCH_MAKE=$(ARCH_MAKE)\
            VPATH=$(VPATH)/libxc-compat)
#
OBJFILES= m_ncps.o m_ncps_froyen_reader.o m_ncps_translators.o \
          m_ncps_froyen_ps_t.o m_ncps_reader.o	m_ncps_writers.o \
          flib_spline.o
#
INCFLAGS:= -I libxc-compat $(PSML_INCFLAGS)
#
libncps.a: libxc-compat $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libncps.a $(OBJFILES) libxc-compat/*.o
	-$(RANLIB) libncps.a
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.*d *.a
	(cd libxc-compat ; $(MAKE) \
                        "ARCH_MAKE=$(ARCH_MAKE)" clean)
#
# DO NOT DELETE THIS LINE - used by make depend
m_libxc_sxc_translation.o: m_libxc_list.o
m_ncps.o: m_ncps_froyen_ps_t.o m_ncps_reader.o m_ncps_translators.o
m_ncps.o: m_ncps_writers.o
m_ncps_froyen_reader.o: flib_spline.o m_ncps_froyen_ps_t.o
m_ncps_reader.o: m_ncps_froyen_ps_t.o m_ncps_froyen_reader.o
m_ncps_reader.o: m_ncps_translators.o
m_ncps_translators.o: m_ncps_froyen_ps_t.o
m_ncps_writers.o: m_ncps_froyen_ps_t.o
