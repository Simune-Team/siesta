      program gen_basis

C Stand-alone program to generate the PAOs, KB projectors, Vlocal, etc
C for a set of atoms. 
C
C The input is exactly the same as that for SIESTA.
C This program generates .ion (optionally .ion.nc) files readable
C by SIESTA
C
C RESTRICTION: Just a single species allowed
C
      use fdf
      use ionew
      use precision
      use basis_types
      use basis_specs
      use basis_io
      use atmparams
      use old_atmfuncs
      use atom
C
C     It has to be MPI-aware since it calls routines down the line...
C
#ifdef MPI
      use mpi
#endif

      implicit none

      integer is, Node, Nodes

#ifdef MPI
      integer  MPIerror
#endif
      external transfer

c Reading input for the pseudopotentials and atomic orbitals 
C Initialise MPI and set processor number
#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node =  0
      Nodes = 1
#endif

      call io_setup     ! Initialize IOnode


      call fdf_init('stdin','gen-basis_out.fdf')

      call read_basis_specs
      call basis_specs_transfer

      nsmax = nsp                   ! For old_atmfuncs
      call allocate_old_arrays
      call clear_tables

      if (nsp .gt. 1) then
         write(6,'(//,a,//)') "JUST A SINGLE SPECIES ALLOWED..."
         STOP
      endif
      do is = 1,nsp
         call write_basis_specs(6,is)
         basp=>basis_parameters(is)
         call atom_main(iz(is), lmxkb(is), nkbl(0,is), 
     .           erefkb(1,0,is),lmxo(is), nzeta(0,1,is), rco(1,0,1,is), 
     .           lambda(1,0,1,is),
     .           atm_label(is), polorb(0,1,is), semic(is), nsemic(0,is),
     .           cnfigmx(0,is),charge(is), smass(is), basistype(is), is,
     $           rinn(0,1,is), vcte(0,1,is), basp)
      enddo 

      call prinput(nsp)

!        Create the new data structures

      call transfer

      call dump_basis_ascii
      call dump_basis_netcdf

      end program gen_basis
!
!     Use private copies of routines to avoid MPI issues
!
      subroutine die(str)
      character(len=*), intent(in), optional   :: str
      if (present(str)) write(6,'(a)') trim(str)
      write(6,'(a)') 'Stopping Program'
      stop
      end subroutine die

      subroutine bye(str)
      character(len=*), intent(in), optional   :: str
      if (present(str)) write(6,'(a)') trim(str)
      write(6,'(a)') 'Requested End of Run. Bye!!'
      stop
      end subroutine bye

      subroutine memory(Task, Type, NElements, CallingRoutine)
      integer  NElements
      character(len=1) Task, Type
      character(len=*) CallingRoutine
      end subroutine memory





