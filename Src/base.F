      program base

C Routine to initialize the Pseudopotentials and Atomic Orbitals.
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


      call fdf_init('stdin','base_out.fdf')

      call basis_read
      call basis_transfer
      call write_basis

      nsmax = nsp                   ! For old_atmfuncs
      call allocate_old_arrays
      call clear_tables

      do is = 1,nsp
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

      end program base
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





