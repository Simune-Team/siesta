#
# Makefile for tbtrans90 gengf90
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
#LINALG=linalg.a
#EXTRA=zgeev.o
#

default: what  tbtrans90 
#default: what  tbtrans90 gengf90 readTSHScdf
#
# System dependencies are dealt with by an include file. See Sys directory
# for examples and instructions. 
#

SIESTA_DIR=../../Src
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

SIESTA_SYS=${SIESTA_ARCH}

#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG) -g
#
what:
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo

#
# System dependent source:
#

SYSOBJ=$(SIESTA_DIR)/$(SYS).o
$(SYSOBJ): $(SIESTA_DIR)/$(SYS).f
	(cd $(SIESTA_DIR) ; $(FC) -c $(FFLAGS) $(SYS).f)

# Tbtrans:
MOD_TOBJS=  m_tbt_options.o m_tbt_gf.o m_tbt_kpts.o hamil2.o tsread2.o

SIESTA_OBJS= precision.o pxf.o  m_mpi_utils.o m_fdf_global.o paste.o \
             parallel.o sys.o io.o files.o reinit.o

TBTOBJS= $(SIESTA_OBJS) $(MOD_TOBJS) $(SYSOBJ) $(EXTRA) tbtrans.o transmission.o   \
         calc_green.o getsfe.o \
         mkRealContour.o  setupkham2.o TSiohs.o atompdos.o coop.o \
         csolve.o cdiag.o teigchana.o helpers.o OutputRegionData.o \
         Lowdin.o netcdfwrite.o \

#
# Interfaces to libraries
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
INCFLAGS=$(FDF_INCFLAGS)
#
MPI_MAKEFILE=$(VPATH)/MPI/Makefile
$(MPI_INTERFACE): 
	(cd MPI ; $(MAKE) -f $(MPI_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/MPI" \
                          "MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "FFLAGS=$(FFLAGS)" module_built)
#
LINALG_MAKEFILE=$(VPATH)/Libs/makefile
LINALG_INCFLAGS:=-I $(VPATH)/Libs $(INCFLAGS)
$(LINALG): 
	(cd Libs ; $(MAKE) -f $(LINALG_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/Libs" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(LINALG_INCFLAGS)" "FFLAGS=$(FFLAGS)" )
#
tbtrans90: $(TBTOBJS) $(FDF) $(SYSOBJ) $(NETCDF_INTERFACE) $(LINALG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o tbtrans90 \
              $(TBTOBJS) $(MPI_INTERFACE) $(NETCDF_INTERFACE)  $(FDF) $(LIBS) $(LINALG)
#	@(cd ../../bin ; rm -f tsScan ; ln -s tbtrans90 tsScan)


#gengf90: $(GENGFOBJS) $(FDF) $(SYSOBJ)
#	$(FC) -o gengf90 \
#	$(FFLAGS) $(LDFLAGS) $(GENGFOBJS)  $(FDF) $(LIBS)  


#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f *.o *.a *.mod
#	rm -f tbtrans90 gengf90 readTSHScdf	
	rm -f Sys/*.o
	rm -f fdf/*.o
	(cd fdf ; make clean)
	(cd MPI ; make clean)

readTSHScdf: readTSHScdf.F $(NETCDF_INTERFACE)
	$(FC) -o readTSHScdf readTSHScdf.F $(NETCDF_INTERFACE) $(NETCDF_LIBS)


#gengf90: $(GENGFOBJS) $(FDF) $(SYSOBJ)
#	$(FC) -o gengf90 \
#	$(FFLAGS) $(LDFLAGS) $(GENGFOBJS)  $(FDF) $(LIBS)  



#
# Dependencies
# 

$(TBTOBJS)      : $(FDF)
reinit.o        : $(MPI_INTERFACE)
sys.o           : $(MPI_INTERFACE) parallel.o
m_mpi_utils.o   : sys.o
m_tbt_gf.o      : m_tbt_kpts.o
m_tbt_kpts      : m_tbt_options.o $(MPI_INTERFACE)
TSprversion.o   : TSversion.h
setupkham2.o	: tsread2.o
tbtrans.o	: hamil2.o 
#gengf.o  setupcontour.o green4.o : fdf/fdfdefs.h
green4.o        : greensub.o
#gengf.o         : cont.o
m_tbt_gf.o       : tsread2.o
##atompdos.o      : fdf/fdfdefs.h
netcdfwrite.o   : $(NETCDF_INTERFACE)
OutputRegionData.o: Lowdin.F $(NETCDF_INTERFACE)
