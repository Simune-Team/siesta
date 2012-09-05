module m_tbt_kpoints
!==============================================================================
! CONTAINS:
!          1) setup_ts_scf_kscell
!          2) setup_ts_kpoint_grid
!          3) ts_write_k_points
!          4) ts_iokp
!          5) ts_mkqgrid
!          6) ts_kpoint_convert

  USE precision, only : dp

  implicit none

  public
  save


  logical   :: siesta_Gamma
  logical   :: Gamma
  integer   :: nkpnt             ! Total number of k-points

  real(dp), pointer :: kweight(:) 
  real(dp), pointer :: kpoint(:,:)

  integer, dimension(3,3)  :: kscell = 0
  real(dp), dimension(3)   :: kdispl = 0.0_dp

  logical     :: spiral = .false.
  logical     :: firm_displ = .false.

  public :: setup_tbt_kscell, setup_tbt_kpoint_grid
  public :: tbt_write_k_points

contains


!-----------------------------------------------------------------------

  subroutine setup_tbt_kscell()

! ***************** INPUT **********************************************

!   The relevant fdf labels kgrid_Monkhorst_Pack and  tbt_kgrid_Monkhorst_Pack.
!   If both are present, tbt_kgrid_Monkhorst_Pack has priority. If none is
!   present, it will create the Gamma point.
!   Examples of fdf data specifications:
!     kgrid_cutoff  50. Bohr
!     %block kgrid_Monkhorst_Pack  # Defines kscell and kdispl
!     4  0  0   0.50               # (kscell(i,1),i=1,3), kdispl(1)
!     0  4  0   0.50               # (kscell(i,2),i=1,3), kdispl(2)
!     0  0  4   0.50               # (kscell(i,3),i=1,3), kdispl(3)
!     %endblock kgrid_Monkhorst_Pack
! **********************************************************************

!  Modules
    use fdf
    use sys, only : die

! Internal variables
    integer           i
    logical           mp_input

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline


    mp_input = fdf_block('TBT.kgrid_Monkhorst_Pack',bfdf)
    if ( .not. mp_input ) &
         mp_input = fdf_block('kgrid_Monkhorst_Pack',bfdf)
    if ( mp_input ) then
       do i = 1,3
          if (fdf_bline(bfdf,pline)) then
             kscell(1,i) = fdf_bintegers(pline,1)
             kscell(2,i) = fdf_bintegers(pline,2)
             kscell(3,i) = fdf_bintegers(pline,3)
             kdispl(i)   = fdf_breals(pline,1)
          else
             call die( 'setup_ts_scf_kscell: ERROR no data in' // &
                  bfdf%label // ' block' )
          endif
       enddo
       firm_displ = .true.
    else
       kscell(:,:) = 0
       kscell(1,1) = 1
       kscell(2,2) = 1
       kscell(3,3) = 1
       kdispl(:)   = 0.0_dp
    endif

    ! We need to do a "crude" initialization of the Gamma setting.
    ! This is because we can only create the k-grid after having
    ! recieved the unit-cell.
    siesta_Gamma = sum(kscell) == 3 .and. &
         sum(kdispl) == 0.0_dp


!     Modify the ts_kscell and ts_kdispl to obtain the 2D sampling
    kscell(1:3,3) = 0
    kscell(3,1:3) = 0
    kscell(3,3)   = 1
    kdispl(3)     = 0.0_dp

    ! We need to do a "crude" initialization of the Gamma setting.
    ! This is because we can only create the k-grid after having
    ! recieved the unit-cell.
    Gamma = sum(kscell) == 3 .and. &
         sum(kdispl) == 0.0_dp


  end subroutine setup_tbt_kscell

  subroutine setup_tbt_kpoint_grid( ucell )

! SIESTA Modules
    USE fdf, only       : fdf_defined
    USE m_find_kgrid, only : find_kgrid
    USE parallel, only  : IONode
    USE precision, only : dp       
#ifdef MPI
    USE mpi_siesta, only : MPI_Bcast, MPI_logical, MPI_Comm_World
#endif

! Local Variables
    real(dp) :: ucell(3,3), eff_kgrid_cutoff

#ifdef MPI
    integer :: MPIerror
#endif

    spiral = fdf_defined('SpinSpiral')

    call find_kgrid(ucell,kscell,kdispl,firm_displ, &
         (.not. spiral), &
         nkpnt,kpoint,kweight, eff_kgrid_cutoff)

    Gamma =  (nkpnt == 1 .and. &
         dot_product(kpoint(:,1),kpoint(:,1)) < 1.0e-20_dp)

    if (IONode) call tbt_write_k_points(ucell)

  end subroutine setup_tbt_kpoint_grid

  subroutine tbt_write_k_points(ucell)
    real(dp) :: ucell(3,3)
    real(dp) :: kpt(3)
    integer  :: ik, ix, i

    write(6,'(/a,i6)')  'TBtrans: k-grid: Number of Transport k-points =', nkpnt
    write(6,'(a)') 'TBtrans: k-grid: Supercell and displacements'
    write(6,'(a,3i4,3x,f8.3)') 'TBtrans: k-grid: ',        &
         (kscell(i,1),i=1,3), kdispl(1)
    write(6,'(a,3i4,3x,f8.3)') 'TBTrans: k-grid: ',        &
         (kscell(i,2),i=1,3), kdispl(2)
!      write(6,'(a,3i4,3x,f8.3)') 'TBTrans: k-grid: ',        &
!           (kscell(i,3),i=1,3), kdispl(3)
    write(6,*)


    write(6,'(/,a)') 'TBTrans: k-point coordinates (Bohr**-1) and weights:'
    do ik = 1 , nkpnt
       write(6,'(a,i4,3f12.6,3x,f12.6)') &
            'TBtrans: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik)
    end do

    write(6,'(/,a)') 'TBTrans: k-point coordinates (b**-1) and weights:'
    do ik = 1 , nkpnt
       call kpoint_convert(ucell,kpoint(:,ik),kpt,1)
       write(6,'(a,i4,3f12.6,3x,f12.6)') &
            'TBtrans: ', ik, (kpt(ix),ix=1,3), kweight(ik)
    end do

  end subroutine tbt_write_k_points

end module m_tbt_kpoints
