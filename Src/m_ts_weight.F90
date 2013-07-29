!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! A module for doing different weighting schemes for the Transiesta voltage calculations
! *Notice* that for the Gamma-point they are all equivalent!

! In the following we describe the weights by:
!  w_L = the weight for the left correction term
!  w_R = the weight for the right correction term
!  w   = w_L / ( w_L + w_R )

! In order to accomodate several methods and precisions we introduce the following 3 methods
! 1)
!   do a single weight for the full correction term, i.e.:
!     w_L = (\sum_k \Delta_L(k))^2
!     w_R = (\sum_k \Delta_R(k))^2
!   this will assume that the k-point correction terms are correlated in some way

! 2)
!   do a single weight, assuming that each k-point correction term is uncorrelated
!     w_L = \sum_k \Delta_L^2(k)
!     w_R = \sum_k \Delta_R^2(k)
!   this will assume that the k-point correction terms are iid's (independent identically distributed)

! 3)
!   do a weight for each k, this is the same as 2), but we increase the k-point weight
!     w_L^k = \Delta_L^2(k)
!     w_R^k = \Delta_R^2(k)
!   this will actually save us more memory, as we don't need the full sum of the k-points.

module m_ts_weight

  use precision, only: dp

  implicit none
  
  private :: dp

  ! Method: 1)
  integer, parameter :: TS_W_CORRELATED = 1
  ! Method: 2)
  integer, parameter :: TS_W_UNCORRELATED = 2
  ! Method: 3)
  integer, parameter :: TS_W_K_UNCORRELATED = 3

  ! we default weight of uncorrelated as that should be 
  ! the most accurate
  integer, save :: TS_W_METHOD = TS_W_K_UNCORRELATED

contains


  subroutine weight_DM( &
       spDML   , spDMR,    &
       spDMneqL, spDMneqR, &
       spEDML  , spEDMR, nonEq_IsWeight)

#ifdef MPI
    use mpi_siesta
#endif
    use parallel,  only: IONode, Node, Nodes
    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData1D

    implicit none

! *********************
! * OUTPUT variables  *
! *********************
    ! Contour part of DM integration
    type(dSpData1D), intent(inout) :: spDML, spDMR
    ! Real-axis part of DM integration
    type(dSpData1D), intent(inout) :: spDMneqL, spDMneqR
    ! L-R estimates of EDM
    type(dSpData1D), intent(inout) :: spEDML, spEDMR
    ! Determine whether DMneqL already is the weight for the
    ! current weighting scheme
    logical, intent(in), optional :: nonEq_IsWeight

! *********************
! * LOCAL variables   *
! *********************
    real(dp) :: wL, wR
    
    ! arrays for looping in the sparsity pattern
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: DML(:), DMR(:)
    real(dp), pointer :: DMneqL(:), DMneqR(:)
    real(dp), pointer :: EDML(:), EDMR(:)
    integer,  pointer :: l_ncol(:)
    integer,  pointer :: l_ptr(:)
    integer,  pointer :: l_col(:)
    integer :: nr
    integer :: io, gio, jo, ind, j
    ! For error estimation
    integer  :: eM_i,eM_j
    real(dp) :: eM, DM, ee
    logical :: l_nonEq_IsWeight
#ifdef MPI
    integer :: MPIerror
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDM' )
#endif
    
    l_nonEq_IsWeight = .false.
    if ( present(nonEq_IsWeight) ) l_nonEq_IsWeight = nonEq_IsWeight

    ! TODO Enforce that sparsity is the same
    ! (however, we know that they are the same)
    sp => spar(spDML)
    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=nr)

    ! Obtain the values in the arrays...
    DML    => val(spDML)
    DMR    => val(spDMR)
    DMneqL => val(spDMneqL)
    DMneqR => val(spDMneqR)
    EDML   => val(spEDML)
    EDMR   => val(spEDMR)

    ! initialize the errors
    eM = 0._dp
    DM = 0._dp

    if ( l_nonEq_IsWeight ) then
    
       do io = 1 , nr
          ! We are in a buffer region...
          if ( l_ncol(io) == 0 ) cycle
          do j = 1 , l_ncol(io)

             ind = l_ptr(io) + j
             ! Retrieve the connecting orbital
             jo = l_col(ind)

             call get_weight(DMneqL(ind),DMneqR(ind),wL,wR)

             ! Do error estimation (capture before update)
             ee = DML(ind) - DMR(ind)

             ! Calculate density...
             DML(ind) = wL * DML(ind) + wR * DMR(ind)

             if ( select_error(ee,io,jo,eM ,eM_i ,eM_j ) ) then
                DM = DML(ind)
             end if

             ! EDML \equiv EDML + EDMneqR
             ! EDMR \equiv EDMR + EDMneqL
             EDML(ind) = wL * EDML(ind) + wR * EDMR(ind)

          end do
       end do

    else

       do io = 1 , nr
          ! We are in a buffer region...
          if ( l_ncol(io) == 0 ) cycle
          do j = 1 , l_ncol(io)

             ind = l_ptr(io) + j
             ! Retrieve the connecting orbital
             jo = l_col(ind)

             call get_weight(DMneqL(ind)**2,DMneqR(ind)**2,wL,wR)

             ! Do error estimation (capture before update)
             ee = DML(ind) + DMneqR(ind) - DMR(ind) - DMneqL(ind)

             DML(ind) = wL * (DML(ind) + DMneqR(ind)) &
                      + wR * (DMR(ind) + DMneqL(ind))

             if ( select_error(ee,io,jo,eM ,eM_i ,eM_j ) ) then
                DM = DML(ind)
             end if

             ! EDML \equiv EDML + EDMneqR
             ! EDMR \equiv EDMR + EDMneqL
             EDML(ind) = wL * EDML(ind) + wR * EDMR(ind)

          end do
       end do

    end if

#ifdef MPI
    ! remove pointer
    nullify(DMR,EDMR)
    allocate(DMR(4*Nodes),EDMR(4*Nodes))

    dit => dist(spDML)
    io = Node + 1
    ! Initialize
    DMR(:) = 0._dp
    DMR(io)         = eM
    DMR(io+Nodes)   = real(index_local_to_global(dit,eM_i,Node),dp)
    DMR(io+2*Nodes) = real(eM_j,dp)
    DMR(io+3*Nodes) = DM
    call MPI_Reduce(DMR(1),EDMR(1),4*Nodes, &
         MPI_Double_Precision, MPI_Sum, 0, MPI_Comm_World, MPIerror)
    if ( IONode ) then
       io = maxloc(EDMR(1:Nodes),1)
       eM   = EDMR(io)
       eM_i = nint(EDMR(io+Nodes))
       eM_j = nint(EDMR(io+2*Nodes))
       io = 3*Nodes + maxloc(EDMR(3*Nodes+1:4*Nodes),1)
       DM  = EDMR(io)
    endif
    deallocate(DMR,EDMR)
#endif

    call print_error_estimate(IONode,'ts: int. EE.:', &
         eM,eM_i,eM_j,DM)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS weightDM' )
#endif

  end subroutine weight_DM

  ! ************* Commonly used modules ****************

  ! Write out the error-estimate for the current iteration (or, k-point)
  subroutine print_error_estimate(IONode,a,eM,eM_i,eM_j,DM)
    logical, intent(in)  :: IONode
    character(len=*), intent(in) :: a
    real(dp), intent(in) :: eM, DM
    integer, intent(in)  :: eM_i,eM_j

    if ( IONode ) then
       write(*,'(a,'' DM_out('',i5,'','',i5,'')'',a,g10.5e1, &
            &'' , d_LR|('',i5,'','',i5,'')|'',a,g10.5e1)') &
            trim(a), &
            eM_i,eM_j,' = ',DM, &
            eM_i,eM_j,' = ',eM
    end if

  end subroutine print_error_estimate


  ! do simple weight calculation and return correct numbers
  pure subroutine get_weight(L,R,wL,wR)
    real(dp), intent(in) :: L,R
    real(dp), intent(out) :: wL,wR

    wL = L + R

    ! The weights
    if ( wL > 0._dp ) then
       wR = R / wL
       wL = L / wL
    else
       wL = 0.5_dp
       wR = 0.5_dp
    end if
  end subroutine get_weight

  ! do simple maximum choice
  function select_error(n_err,n_io,n_jo,err,io,jo) result(val)
    real(dp), intent(in) :: n_err
    integer, intent(in) :: n_io, n_jo
    real(dp), intent(inout) :: err
    integer, intent(inout) :: io, jo
    logical :: val
    val = abs(n_err) > abs(err)
    if ( val ) then
       err = n_err
       io  = n_io
       jo  = n_jo
    end if
  end function select_error

end module m_ts_weight
