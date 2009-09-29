#
#  Makefile for Sies2xcf set of tools
#
# --------------------------------- proper interface to Siesta building system
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
ARCH_MAKE=../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
LFLAGS=$(FFLAGS)
# --------------------------------- 
#
EXECS = eig2bxsf xv2xsf md2axsf rho2xsf vib2xsf
all: $(EXECS)
   eig_obj = eig2bxsf.o inver3.o opnout.o 
   xv_obj  = xv2xsf.o  opnout.o
   md_obj  = md2axsf.o makebox.o  fillbox.o inver3.o  hit.o \
             opnout.o  test_ani.o test_md.o wraxsf1.o wraxsf2.o
   rho_obj = rho2xsf.o read_xv.o  makebox.o fillbox.o inver3.o  hit.o \
             intpl04.o opnout.o 
   vib_obj = vib2xsf.o read_xv.o  makebox.o fillbox.o inver3.o  hit.o \
	     displa.o  opnout.o   read_ev.o itochar.o w_arrow.o w_movie.o

 eig2bxsf : $(eig_obj)
	$(FC) $(LFLAGS) -o $@ $(eig_obj)
  xv2xsf   : $(xv_obj)
	$(FC) $(LFLAGS) -o $@ $(xv_obj)
  md2axsf  : $(md_obj)
	$(FC) $(LFLAGS) -o $@ $(md_obj)
  rho2xsf  : $(rho_obj)
	$(FC) $(LFLAGS) -o $@ $(rho_obj)
  vib2xsf  : $(vib_obj)
	$(FC) $(LFLAGS) -o $@ $(vib_obj)

clean:
	@echo " Cleaning up"
	rm -f $(EXECS) *.o core

