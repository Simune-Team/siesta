      subroutine initatom

! Routine to initialize the Pseudopotentials and Atomic Orbitals.
! Substantially modified by Alberto Garcia (2000)
!
! The PAO and KB information can optionally be read from ASCII files
! (those produced by a standard run of Siesta or Base, but with the extension 
! renamed to '.input' instead of '.dump'), or from NetCDF files (if NetCDF
! is available). Note that there must be files for *all* species.
!
! This behavior is controlled by the FDF logicals 'user-basis' and
! 'user-basis-netcdf'.
!
! The old 'USER' basis type has been removed.
!
! This routine also outputs information about the basis specification
! determined by the routines in the 'basis_specs' modules. 
!


      use fdf
      use precision
      use basis_types
      use basis_specs
      use basis_io
      use old_atmfuncs
      use atom

      implicit none

C Internal variables ...................................................
      integer is

      logical user_basis, user_basis_netcdf
      
      external transfer

c Reading input for the pseudopotentials and atomic orbitals 

      write(6,'(/2a)') 
     .    'initatom: Reading input for the pseudopotentials ',
     .    'and atomic orbitals'

      user_basis = fdf_boolean('user-basis',.false.)
      user_basis_netcdf = fdf_boolean('user-basis-netcdf',.false.)

      if (user_basis_netcdf) then

         write(6,'(a)') 'Reading PAOs and KBs from NetCDF files...'
         call read_basis_netcdf

      else if (user_basis) then

         write(6,'(a)') 'Reading PAOs and KBs from ascii files...'
         call read_basis_ascii

      else
!
!     New routines in basis_specs and basis_types.
!
         call read_basis_specs
         call basis_specs_transfer

         nsmax = nsp             !! For old_atmfuncs
         call allocate_old_arrays
         call clear_tables

         do is = 1,nsp
            call write_basis_specs(6,is)
            basp=>basis_parameters(is)
            call atom_main( iz(is), lmxkb(is), nkbl(0,is), 
     .           erefkb(1,0,is),lmxo(is), nzeta(0,1,is), rco(1,0,1,is), 
     .           lambda(1,0,1,is),
     .           atm_label(is), polorb(0,1,is), semic(is), nsemic(0,is),
     .           cnfigmx(0,is),charge(is), smass(is), basistype(is), is,
     $           rinn(0,1,is), vcte(0,1,is),basp)
         enddo 

         call prinput(nsp)

!        Create the new data structures for atmfuncs.

         call transfer

      endif

      call dump_basis_ascii
      call dump_basis_netcdf
      call dump_basis_xml

      end subroutine initatom




