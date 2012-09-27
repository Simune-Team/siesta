subroutine ts_do_Green(tElec, HSFile, GFFile, GFTitle, &
     ElecValenceBandBot, ReUseGF, &
     nkpnt,kpoint,kweight, &
     NBufAt,NUsedAtoms,NA1,NA2, &
     ucell,xa,nua,NEn,contour,wGF,chem_shift,ZBulkDOS,nspin, &
     nua_GF,no_GF)

  use precision,  only : dp
  use parallel  , only : IONode
  use sys ,       only : die
#ifdef MPI
  use mpi_siesta, only : MPI_Comm_World
  use mpi_siesta, only : MPI_Bcast, MPI_Integer, MPI_Logical
#endif

  use m_ts_electrode, only : create_Green
  use m_ts_GF,        only : check_Green

  implicit none

! ***********************
! * INPUT variables     *
! ***********************
  character(len=1), intent(in) :: tElec   ! 'L' for Left electrode, 'R' for right
  character(len=*), intent(in) :: HSFile  ! The electrode TSHS file
  character(len=*), intent(in) :: GFFile  ! The electrode GF file to be saved to
  character(len=*), intent(in) :: GFTitle ! The title to be written in the GF file
  logical, intent(in)          :: ElecValenceBandBot ! Whether or not to calculate electrodes valence bandbottom
  logical, intent(in)          :: ReUseGF ! Should we re-use the GF files if they exists?    
  integer, intent(in)          :: nkpnt ! Number of k-points
  real(dp),dimension(3,nkpnt),intent(in) :: kpoint ! k-points
  real(dp),dimension(nkpnt),intent(in) :: kweight ! weights of kpoints
  integer, intent(in)            :: NBufAt,NA1,NA2 ! Buffer/Rep a1/Rep a2
  integer, intent(inout)         :: NUsedAtoms ! Needs update here
  integer, intent(in)            :: nua ! Full system count of atoms in unit cell
  real(dp), dimension(3,3)       :: ucell ! The unit cell of the CONTACT
  real(dp), intent(in)           :: xa(3,nua) ! Coordinates in the system for the TranSIESTA routine
  integer, intent(in)            :: nspin ! spin in system
  integer, intent(in)            :: NEn ! Number of energy points
  complex(dp), intent(in)        :: contour(NEn),wGF(NEn) !energy contours and weights for GF
  real(dp), intent(in)           :: chem_shift ! the Fermi-energy we REQUIRE the electrode

! ***********************
! * OUTPUT variables    *
! ***********************
  complex(dp), intent(out)       :: ZBulkDOS(NEn) ! DOS at energy points
  integer, intent(out)           :: nua_GF ! The number of atoms in the GF file
  integer, intent(out)           :: no_GF  ! The number of orbitals in the GF file

! ***********************
! * LOCAL variables     *
! ***********************
  integer :: uGF
  logical :: errorGF
#ifdef MPI
  integer :: MPIerror
#endif

#ifdef DEBUG
  call write_debug( 'PRE do_Green' )
#endif

! Create the GF file
  call create_Green(tElec,HSFile, GFFile, GFTitle, &
       ElecValenceBandBot, ReUseGF, &
       nkpnt,kpoint,kweight, &
       NBufAt,NUsedAtoms,NA1,NA2, &
       ucell,xa,nua,NEn,contour,wGF,chem_shift,ZBulkDOS,nspin)

!
!     Check that the Green's functions are correct!
!     This is needed as create_Green returns if user requests not to
!     overwrite an already existing file.
!     This check will read in the number of orbitals and atoms in the
!     electrode surface Green's function.
!

  errorGF = .false.
! Check the GF file
  if(IONode) then
     call io_assign(uGF)
     open(file=GFFile,unit=uGF,form='UNFORMATTED')

     call check_Green(uGF,chem_shift,ucell,nspin,nkpnt,kpoint, &
          kweight,NEn,contour,wGF,NA1,NA2,errorGF, &
          nua_GF,no_GF)

     call io_close(uGF)
  endif

#ifdef MPI
  call MPI_Bcast(nua_GF,1,MPI_integer,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(no_GF,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

!     Check the error in the GF file
#ifdef MPI
  call MPI_Bcast(errorGF,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
#endif
  if ( errorGF ) call die("Error in GFfile: "//trim(GFFile))

#ifdef DEBUG
  call write_debug( 'POS do_Green' )
#endif

end subroutine ts_do_Green
