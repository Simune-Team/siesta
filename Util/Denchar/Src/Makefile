#Makefile for denchar
#
.SUFFIXES: .f .F .o .a  .f90 .F90
MOD_DIR=Modules
#
#default: what denchar 
default: what denchar
#
include arch.make
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = atompla.o  \
        bessph.o chkdim.o colinear.o denchar.o dismin.o \
	dmna.o dot.o io.o iodm.o length.o matvect.o memory.o \
	neighb.o paste.o planed.o  prmem.o \
	radfft.o ranger.o readpla.o recipes.o \
        redata.o reinit.o rhoofr.o\
	version.o volcel.o wrout.o 
#
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" module)
#
# Routines using fdf calls.
#
atompla.o basis_io.o denchar.o iodm.o  rhoofr.o : $(FDF)
readpla.o reinit.o redata.o wrout.o : $(FDF)
#
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=precision.o parallel.o memoryinfo.o sys.o listsc.o  \
         atmparams.o atmfuncs.o ionew.o atm_types.o \
         radial.o spher_harm.o  basis_io.o basis_types.o \
         pseudopotential.o chemical.o xml.o
COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
$(COM_OBJS): precision.o parallel.o $(MPI_INTERFACE)
#
atmfuncs.o: atmfuncs.f atmparams.o atm_types.o radial.o spher_harm.o sys.o
atm_types.o: atm_types.f radial.o
basis_types.o: pseudopotential.o atmparams.o
basis_io.o: $(NETCDF_INTERFACE) basis_io.F atmfuncs.o ionew.o sys.o \
	atmparams.o atm_types.o radial.o bessph.o chemical.o xml.o \
	radfft.o recipes.o basis_types.o 
chkdim.o: sys.o
iodm.o : iodm.F 
ionew.o: ionew.F sys.o 
memory.o: memory.F memoryinfo.o
pseudopotential.o: sys.o
spher_harm.o: spher_harm.f sys.o
radial.o: xml.o
rhoofr.o denchar.o: rhoofr.f listsc.o
version.o: version.F 
atompla.o: sys.o
chemical.o: sys.o

# Interfaces to libraries
#
libmpi_f90.a: 
	@(cd MPI ; $(MAKE) )
libnetcdf_f90.a:
	@(cd NetCDF ; $(MAKE) )
netcdf.mod:
	@(cd NetCDF ; $(MAKE) )
#
# Libraries that need to be compiled
#
#
version:
	sed  "s/SIESTA_ARCH/${SIESTA_ARCH}/g" version.F > temp0.F
	sed  "s/FFLAGS/${FC} ${FFLAGS}/g" temp0.F > temp.F
	$(FC) -c $(FFLAGS) $(DEFS) temp.F
	@mv temp.o version.o
	@rm -f temp.F temp0.F 
	@echo

denchar: version $(MPI_INTERFACE) $(FDF) \
                $(COMP_LIBS) $(ALL_OBJS)
	$(FC) -o denchar \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(MPI_INTERFACE)\
               $(COMP_LIBS) $(LIBS) 
#

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f denchar *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)
	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi
	(cd NetCDF ; $(MAKE) clean )

