module m_tbt_iotshs
!
! Routines that are used for Input  of files 
! This has been adapted to be used in TBtrans forms
!
!=============================================================================
! CONTAINS:
!          1) tbt_read_tshs

  implicit none

  public :: tbt_read_tshs

  private

contains

  subroutine tbt_read_tshs(HSfile, Gamma, no_s,no_u,nspin, &
       ucell, na_u, xa, lasto, &
       maxnh , numh , listhptr , listh , xij , indxuo, &
       H, S, Ef)

! *********************************************************************
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible
! Modified by M.Paulsson 2009 to:
! 1: To include information of which FC step for phonon calculations
! 2: To only save the overlap matrix if onlyS flag is set
!    (Used for e-ph coupling calculations)
! 3: File format changed to unify Copenhagen/Barcelona Transiesta vers.
! 4: Smaller files by writing arrays directly instead of element wise
! *************************** INPUT **********************************
! character(*) HSfile         : file name to read
! *************************** OUTPUT *********************************
! logical Gamma               : Where it a Gamma calculation
! integer no_s                : Number of basis orbitals per supercell
! integer no_u                : Number of basis orbitals per unit cell
! integer nspin               : Spin polarization (1 or 2)
! real*8  ucell(3,3)          : unit cell of TSHS file
! integer na_u                : Number of atoms per unit cell
! real*8  xa(3,na_u)          : Atomic coordinates in the unitcell
! integer lasto(0:na_u)       : The last orbital of the atom in the unit cell
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(no_u)          : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(no_u)      : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! integer indxuo(no_s)        : Index of orbitals in supercell
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  Ef                  : The Fermi level of the system

!
!  Modules
!
    use precision,     only : dp
    use fdf,           only : leqi
    use parallel,      only : IONode
    use sys,           only : die
    use m_tbt_kpoints, only : kscell, kdispl
    use m_ts_io,       only : ts_read_TSHS

#ifdef MPI
    use mpi_siesta
#endif

! **************************
! * INPUT variables        *
! **************************
    character(len=*),intent(in) :: HSfile    ! file name of TSHS

! **************************
! * OUTPUT variables       *
! **************************
    logical, intent(out)  :: Gamma           ! Gamma point calculation (i.e. one k-point)
    integer, intent(out)  :: no_s            ! # orbitals in supercell cell local to this processor
    integer, intent(out)  :: no_u            ! # orbitals in unit cell global
    integer, intent(out)  :: nspin           ! spins in the system
    real(dp), intent(out) :: ucell(3,3)      ! The unit cell that the transformation is happening in
    integer, intent(out)  :: na_u            ! The number of atoms in the unit cell
    real(dp), pointer, intent(out) :: xa(:,:)      ! The atomic coordinates in the unit cell
    integer, pointer, intent(out)  :: lasto(:)     ! The last orbital of each atom referenced in the unit cell
    integer, intent(out) :: maxnh           ! Maximum number of orbitals interacting
    integer, pointer, intent(out)  :: numh(:)      ! # of nonzero elements of each row of hamiltonian matrix
    integer, pointer, intent(out)  :: listhptr(:)  ! Pointer to each row (-1) of the hamiltonian matrix
    integer, pointer, intent(out)  :: listh(:)     ! Nonzero hamiltonian-matrix element column
!                                                    indexes for each matrix row
    real(dp), pointer, intent(out) :: xij(:,:)     ! Vectors between orbital centers (sparse)
!                                                    (not only gamma point)
    integer, pointer, intent(out)  :: indxuo(:)    ! Index of equivalent orbital in unit cell
!                                                    Unit cell orbitals must be the first in
!                                                    orbital lists, i.e. indxuo.le.no_l, with
!                                                    no_l the number of orbitals in unit cell
    real(dp), pointer, intent(out) :: H(:,:)        ! Hamiltonian in sparse form
    real(dp), pointer, intent(out) :: S(:)        ! Overlap in sparse form
    real(dp), intent(out) :: Ef              ! Fermi energy

! **************************
! * LOCAL variables        *
! **************************
    ! Loop counters
    integer :: i, j, fL

! Variables concerning the k-point sampling in the TranSIESTA run
    integer, dimension(3,3) :: kscell_file
    real(dp), dimension(3)  :: kdispl_file 
    integer, pointer :: iza(:)
    logical :: TSGamma, same, onlyS
    real(dp) :: Qtot, Temp
#ifdef MPI
    integer :: MPIerror
#endif

    fL = len_trim(HSfile)
    if ( leqi(HSfile(fL-4:fL),'.TSHS') ) then
       call ts_read_TSHS(HSfile,onlyS,Gamma,TSGamma, &
            ucell, na_u, no_u, no_u, no_s, maxnh, nspin, &
            kscell_file,kdispl_file, &
            xa, iza, lasto, &
            numh , listhptr , listh , xij , indxuo, &
            H, S, Ef, Qtot, Temp, i,j, &
            Bcast=.true.)
    else
       call die('Could not determine file format of: '//trim(HSfile))
    end if
    

! Check the k-point sampling
! In this case as we are in TBTrans we do not worry that they are not the same.
! In fact it can often be increased in accuracy while calculating the transmission.
! However, we write out that we have encountered a non-matching k-point sampling...
    if ( IONode ) then
       same = .true.
       do i = 1,3
          same = same .and. ( kdispl_file(i) == kdispl(i) )
       end do
       if ( .not. same ) then
          write(*,*) 'NOTICE: TSHS and TBTrans does not have same k-displacements.'
          write(*,'(a,2F8.4)') 'TSHS k displacements    :', (kdispl_file(j),j=1,2)
          write(*,'(a,2F8.4)') 'TBTrans k displacements :', (kdispl(j),j=1,2)
       end if
       same = .true.
       do j = 1,3
          do i = 1,3
             same = same .and. ( kscell_file(i,j) == kscell(i,j) )
          end do
       end do
       if ( .not. same ) then
          write(*,*) 'NOTICE: TSHS and TBTrans does not have same k-cell.'
          do i = 1,3
             write(*,'(a,/,3(2x,i3))') 'TSHS k cell    :', (kscell_file(i,j),j=1,3)
          end do
          do i = 1,3
             write(*,'(a,/,3(2x,i3))') 'TBTrans k cell :', (kscell(i,j),j=1,3)
          end do
       end if
    end if

  end subroutine tbt_read_TSHS

end module m_tbt_iotshs
