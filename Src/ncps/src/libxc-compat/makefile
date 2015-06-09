.SUFFIXES: 
.SUFFIXES: .f .F .o .f90 .F90
#
ARCH_MAKE_DEFAULT=../fortran.mk
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
VPATH=.
#
OBJFILES= m_libxc_compat.o m_libxc_list.o m_libxc_sxc_translation.o \
          m_siestaxc_list.o
#
default: objs_and_modules
#
objs_and_modules: $(OBJFILES)
#
#------------------------------
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.mod
#
# DO NOT DELETE THIS LINE - used by make depend
m_libxc_compat.o: m_libxc_list.o m_libxc_sxc_translation.o
m_libxc_sxc_translation.o: m_libxc_list.o m_siestaxc_list.o
