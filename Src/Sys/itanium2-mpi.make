SIESTA_ARCH=itanium2-mpi
#
#  itanium2 machine at Leioa with myrinet
#
FC=efc -Vaxlib
#
FC_ASIS=$(FC)
#
FFLAGS= -O2 -tpp2 -W0 
FFLAGS_DEBUG= -g -O0
LDFLAG=#-Vaxlib 
RANLIB=echo
LIBS=  
SYS=bsd
DEFS=-DMPI -DNODAT
MPI_INTERFACE=libmpi_f90.a
#
MPIROOT=/usr/local/mpich-1.2.5-10-IntelComp-7.1-GM-2.0
MPI_INCLUDE=$(MPIROOT)/include
#
#MPI_LIBS= $(MPIROOT)/lib/libmpich.a  $(MPIROOT)/lib/libpmpich.a  \
#         /usr/local/gm-2.0.6_Linux/binary/lib/libgm.a \
#          -L/opt/intel/mkl61/lib/64 -lpthread
#
MPI_LIBS=-L$(MPIROOT)/lib \
         -lmpichf90 -lpmpich -lmpich -lpmpich -lmpich \
         -L/usr/local/gm/binary/lib/ -L/usr/local/gm/lib/ \
         -lgm -lpthread -lPEPCF90
#
#LAPACK=-L/opt/intel/mkl61/lib/64 -llibmkl_lapack.a
#
LAPACK=-L/opt/intel/mkl61/lib/64 -lmkl_lapack64 -lmkl -lguide -lpthread
BLACS_LIBS=$(HOME)/lib/libblacs.a
SCALAPACK_LIBS=-L/home/wdpgaara/lib -lscalapack-arina.GM  \
                                    -lpblas-arina.GM \
                                    -ltools-arina.GM \
			            -lredist-arina.GM
#
LIBS=  $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK) \
       $(MPI_LIBS)  $(NETCDF_LIBS)
COMP_LIBS=
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#



