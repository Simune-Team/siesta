! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
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
             kdispl(i)   = mod(fdf_breals(pline,1), 1._dp)
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

    call tbt_find_kgrid(ucell,kscell,kdispl,firm_displ, &
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


  subroutine tbt_find_kgrid ( cell, kscell, displ, firm_displ, &
       nk, points, weight, eff_kgrid_cutoff )
! ***************** INPUT **********************************************
! real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
! integer kscell(3,3): Supercell reciprocal of k-grid unit cell
!                      scell(ix,i) = sum_j cell(ix,j)*kscell(j,i)
! ***************** INPUT/OUTPUT  **************************************
! real*8  displ(3)   : Grid origin in k-grid-vector coordinates:
!                      origin(ix) = sum_j gridk(ix,j)*displ(j)
! ***************** OUTPUT *********************************************
! integer nk           : Actual number of irreducible k-points
! real(dp) points(3:nk): Kpoints
! real(dp) weight(nk)  : weights
! real(dp) eff_kgrid_cutoff : actual equivalent kgrid cutoff 

!  Modules
    use precision,  only : dp
    use alloc,      only : re_alloc
    use units,      only : pi
    use m_minvec,   only : minvec
    use parallel,   only : Node
    
    implicit          none

! Passed variables
    integer, intent(in)      :: kscell(3,3)
    real(dp), intent(in)     :: cell(3,3)
    logical, intent(in)      :: firm_displ
    
    real(dp), intent(inout)  :: displ(3)
    
    integer, intent(out)     :: nk
    real(dp), intent(out)    :: eff_kgrid_cutoff
    real(dp), pointer        :: points(:,:)
    real(dp), pointer        :: weight(:)

    external                 ::  idiag, reclat

! Internal variables
    integer  :: i, ir, ix, j, i1, i2, i3
    integer  :: kdsc(3,3), maux(3,3,2), ml(3,3), mr(3,3)
    integer  :: ng(3), ni, nkr(3), nktot
    integer  :: proj(3,3), igmin(3), igmax(3)

    real(dp) :: d(3), dkg(3), dkx(3), dscell(3,3)
    real(dp) :: gridk(3,3), gscell(3,3)
    real(dp) :: scell(3,3), scmin(3,3), vmod, w1, wtot
    real(dp) :: ctransf(3,3)

    real(dp), parameter :: tiny = 1.d-12

    character(len=24) :: accum ! Used for Intel_12 bug

! Find total number of points (determinant of kscell)
    nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) + &
         kscell(2,1) * kscell(3,2) * kscell(1,3) + &
         kscell(3,1) * kscell(1,2) * kscell(2,3) - &
         kscell(1,1) * kscell(3,2) * kscell(2,3) - &
         kscell(2,1) * kscell(1,2) * kscell(3,3) - &
         kscell(3,1) * kscell(2,2) * kscell(1,3) )
    
! 
! Find k-grid supercell
!
    do i = 1,3
       do ix = 1,3
          scell(ix,i) = cell(ix,1) * kscell(1,i) + & 
               cell(ix,2) * kscell(2,i) + &
               cell(ix,3) * kscell(3,i)
       enddo
       vmod = sqrt( scell(1,i)**2 + scell(2,i)**2 + scell(3,i)**2 )
    enddo
    
! Find actual cutoff
    call minvec( scell, scmin, ctransf )
    eff_kgrid_cutoff = huge(1.0_dp)
    do i = 1,3
       vmod = sqrt( scmin(1,i)**2 + scmin(2,i)**2 + scmin(3,i)**2 )
       eff_kgrid_cutoff = min( eff_kgrid_cutoff, vmod/2.d0 )
    enddo
    
!     Equivalent supercell DA with the property that there exists
!     a primitive cell (pa') such that DA_i = N_i*pa'_i
!     (See Moreno and Soler)
!!    Direct route
!!    call DIGCEL( ucell, kSCELL, new_ucell, dscell, NSC, ISDIAG )

    call idiag( 3, kscell, kdsc, ml, mr, maux )
    proj(:,:) = 0  ! Possible sign changes
    do i = 1, 3
       proj(i,i) = 1
       if (kdsc(i,i) < 0) proj(i,i) = -1
    enddo
    kdsc = matmul(kdsc,proj)
    mr = matmul(mr,proj)
!
!     Set the displacements if not firm (i.e., specified by the
!     user). Even if firm, warn if a better choice is possible.
!
    do j = 1, 3
       if (mod(kdsc(j,j),2) .eq. 0) then
          if (firm_displ .and. displ(j) /= 0.5d0) then
             if (Node .eq. 0) &
                  write(6,"(a,i4,a,2f8.2)") &
                  "k-point displ. along", j, " input, could be: ", &
                  displ(j), 0.5d0
          else
             displ(j) = 0.5d0
          endif
       else
          if (firm_displ .and. displ(j) /= 0.0d0) then
             if (Node .eq. 0) &
                  write(6,"(a,i4,a,2f8.2)") &
                  "k-point displ. along", j, " input, could be: ", &
                  displ(j), 0.0d0
          else
             displ(j) = 0.0d0
          endif
       endif
    enddo
!
    dscell = matmul(scell,mr)
    
! Find k-grid unit vectors
    call reclat( dscell, gridk, 1 )
    
! Find grid origin in cartesian coordinates
    call reclat( scell, gscell, 1 )
    do ix = 1,3
       dkx(ix) = gscell(ix,1) * displ(1) + &
            gscell(ix,2) * displ(2) + &
            gscell(ix,3) * displ(3)
    enddo

! Find grid origin in gridk coordinates
    
    do i = 1,3
       dkg(i) = ( dkx(1) * dscell(1,i) + &
            dkx(2) * dscell(2,i) + &
            dkx(3) * dscell(3,i) ) / (2*pi)
    enddo

! Find total range of grid indexes
    do j = 1,3
       ng(j) = kdsc(j,j)
       igmin(j) = -( (ng(j)-1) / 2)
       igmax(j) = ng(j) / 2
    enddo
    
! Find number of points with time-reversal (inversion) symmetry,
! (if possible) after reflection on each alternative plane
    do j = 1,3
       ni = ng(j)
       if (abs(dkg(j)) .lt. tiny) then
          ni = ng(j)/2 + 1
       elseif (abs(dkg(j)-0.5d0) .lt. tiny) then
          ni = (ng(j)-1)/2 + 1
       endif
! To work around an Intel_12 compiler bug
       write(accum,"(3i8)") ni,nktot,kdsc(j,j)
       nkr(j) = ni * nktot / kdsc(j,j)
    enddo
    
! Select reflection plane
    ir = 3
    if (nkr(2) .lt. nkr(ir)) ir = 2
    if (nkr(1) .lt. nkr(ir)) ir = 1
    igmin(ir) = 0
    if (abs(dkg(ir)-0.5d0) .lt. tiny) &
         igmax(ir) = (ng(ir)-1)/2
    nk = nkr(ir)
    
    call re_alloc(points,1,3,1,nk,name='points', &
         routine='find_kgrid',copy=.false.)
    call re_alloc(weight,1,nk,name='weight', &
         routine='find_kgrid',copy=.false.)

! Find k points and weights
    w1 = 1.0d0 / nktot
    nk = 0
    do i3 = igmin(3),igmax(3)
       do i2 = igmin(2),igmax(2)
          do i1 = igmin(1),igmax(1)
             nk = nk + 1
             d(1) = i1 + dkg(1)
             d(2) = i2 + dkg(2)
             d(3) = i3 + dkg(3)
             if (d(1) .gt. 0.5d0*ng(1)+tiny) d(1) = d(1) - ng(1)
             if (d(2) .gt. 0.5d0*ng(2)+tiny) d(2) = d(2) - ng(2)
             if (d(3) .gt. 0.5d0*ng(3)+tiny) d(3) = d(3) - ng(3)
             do ix = 1,3
                points(ix,nk) = gridk(ix,1)*d(1) +  &
                     gridk(ix,2)*d(2) + &
                     gridk(ix,3)*d(3)
             enddo
             if ( abs(d(ir)) .lt. tiny .or. &
                  abs(d(ir)-0.5d0*ng(ir)) .lt. tiny) then
                weight(nk) = w1
             else
                weight(nk) = 2.0d0 * w1
             endif
          enddo
       enddo
    enddo

! Check that weight is normalised
    wtot = sum(weight(1:nk))
    if (abs(wtot-1.0d0) .gt. nk*tiny) then
       w1 = dble(nk)/wtot
       weight(1:nk) =  w1*weight(1:nk)
    endif
    
  end subroutine tbt_find_kgrid

end module m_tbt_kpoints
