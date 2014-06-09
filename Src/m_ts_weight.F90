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

  ! Generic methods for the k-weighting method
  ! Method: 1)
  integer, parameter :: TS_W_K_CORRELATED = 1
  ! Method: 2)
  integer, parameter :: TS_W_K_HALF_CORRELATED = 2
  ! Method: 3)
  integer, parameter :: TS_W_K_UNCORRELATED = 3

  ! General weighting method
  ! 1) orb-orb weighting
  integer, parameter :: TS_W_ORB_ORB = 1
  ! 2) atom-atom weighting using the trace of the nEq DM
  integer, parameter :: TS_W_TR_ATOM_ATOM = 2
  ! 3) atom-atom weighting using the sum of the nEq DM
  integer, parameter :: TS_W_SUM_ATOM_ATOM = 3
  ! 4) atom-orb weighting using the trace of the nEq DM
  integer, parameter :: TS_W_TR_ATOM_ORB = 4
  ! 5) atom-orb weighting using the sum of the nEq DM
  integer, parameter :: TS_W_SUM_ATOM_ORB = 5
  
  ! The general weighting can be:
  !   UNCORRELATED or CORRELATED
  ! The correlated is the default
  integer, parameter :: TS_W_CORRELATED = 100

  integer, save :: TS_W_METHOD = TS_W_ORB_ORB

  ! we default weight of uncorrelated as that should be 
  ! the most accurate
  integer, save :: TS_W_K_METHOD = TS_W_K_UNCORRELATED

contains


  subroutine weight_DM(N_Elec,Elecs,N_mu, na_u, lasto, &
       spDM, spDMneq, spEDM, &
       nonEq_IsWeight)
    
#ifdef MPI
    use mpi_siesta
#endif
    use parallel,  only: IONode, Node, Nodes
    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData2D

    use m_ts_electype
    use m_ts_contour_neq, only : N_nEq_ID, ID2mu
    use geom_helper, only : iaorb, ucorb

    implicit none

! *********************
! * OUTPUT variables  *
! *********************
    integer,            intent(in) :: N_Elec, N_mu
    type(Elec),         intent(in) :: Elecs(N_Elec)
    ! The last-orbital of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! Contour part of DM integration
    type(dSpData2D), intent(inout) :: spDM
    ! Real-axis part of DM integration
    type(dSpData2D), intent(inout) :: spDMneq
    ! Estimates of EDM
    type(dSpData2D), intent(inout) :: spEDM
    ! Determine whether DMneq already is the weight for the
    ! current weighting scheme
    logical, intent(in), optional :: nonEq_IsWeight

! *********************
! * LOCAL variables   *
! *********************
    real(dp) :: w(N_mu), neq(N_mu)
    
    ! arrays for looping in the sparsity pattern
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: DM(:,:), DMneq(:,:), EDM(:,:)
    integer,  pointer :: l_ncol(:)
    integer,  pointer :: l_ptr(:)
    integer,  pointer :: l_col(:)
    integer :: nr
    integer :: io, jo, ind, j
    integer :: mu_i, mu_j
    integer, allocatable :: ID_mu(:)
    ! For error estimation
    integer  :: eM_i,eM_j
    real(dp) :: eM, DMe, ee, ee_i, tmp
    logical :: l_nonEq_IsWeight, hasEDM, is_correlated, is_trace
    ! collecting the error contribution for each atom
    real(dp), allocatable :: atom_w(:,:), atom_neq(:,:)
    integer :: ng, lio, ia1, ia2, ia, TS_W
    
#ifdef MPI
    integer :: MPIerror
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDM' )
#endif
    
    l_nonEq_IsWeight = .false.
    if ( present(nonEq_IsWeight) ) l_nonEq_IsWeight = nonEq_IsWeight
    if ( l_nonEq_IsWeight .and. TS_W_METHOD /= TS_W_ORB_ORB ) then
       call die('Not implemented yet... Programming error')
    end if

    ! TODO Enforce that sparsity is the same
    ! (however, we know that they are the same)
    sp => spar(spDM)
    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=nr,nrows_g=ng)

    ! Obtain the values in the arrays...
    DM     => val(spDM)
    DMneq  => val(spDMneq)
    hasEDM = initialized(spEDM)
    if ( hasEDM ) EDM => val(spEDM)

    ! point to the orbital-distribution
    dit => dist(spDM)

    allocate(ID_mu(N_nEq_ID))
    do io = 1 , N_nEq_ID
       ID_mu(io) = ID2mu(io)
    end do

    ! initialize the errors
    eM  = 0._dp
    DMe = 0._dp

    if ( l_nonEq_IsWeight ) then
    
       do io = 1 , nr
          ! We are in a buffer region...
          if ( l_ncol(io) == 0 ) cycle
          do j = 1 , l_ncol(io)

             ind = l_ptr(io) + j
             ! Retrieve the connecting orbital
             jo = l_col(ind)

             call get_weight(N_Elec,N_mu,N_nEq_ID,ID_mu,DMneq(ind,:),w)

             ! Do error estimation (capture before update)
             ee = 0._dp
             do mu_i = 1 , N_mu - 1
                ee_i = DM(ind,mu_i)
                do mu_j = mu_i + 1 , N_mu
                   tmp = ee_i - DM(ind,mu_j)
                   if ( abs(tmp) > abs(ee) ) ee = tmp
                end do
             end do
             
             DM(ind,1)  = w(1) * ( DM(ind,1) )
             if ( hasEDM ) EDM(ind,1) = w(1) *  EDM(ind,1)
             do mu_i = 2 , N_mu
                DM(ind,1)  = DM(ind,1) + &
                     w(mu_i) * DM(ind,mu_i)
                if ( hasEDM ) EDM(ind,1) = EDM(ind,1) + &
                     w(mu_i) * EDM(ind,mu_i)
             end do

             if ( abs(ee) > abs(eM) ) then
                eM   = ee
                eM_i = io
                eM_j = jo
                DMe = DM(ind,1)
             end if

          end do
       end do

    else

       ! Is the data correlated
       is_correlated = TS_W_METHOD >= TS_W_CORRELATED
       TS_W = TS_W_METHOD
       if ( is_correlated ) TS_W = TS_W - TS_W_CORRELATED

       if ( TS_W /= TS_W_ORB_ORB ) then
          ! we are doing weighting per trace/sum of each atom

          ! this will not be that large an array... :)
          allocate(atom_neq(N_nEq_ID,na_u))
          atom_neq(:,:) = 0._dp

          ! determine whether it is trace or sum
          select case ( TS_W )
          case ( TS_W_SUM_ATOM_ATOM , TS_W_SUM_ATOM_ORB )
             is_trace = .false.
          case ( TS_W_TR_ATOM_ATOM , TS_W_TR_ATOM_ORB )
             is_trace = .true.
          case default
             call die('Error in weights: determine trace')
          end select

          do lio = 1 , nr

             ! We are in a buffer region...
             if ( l_ncol(lio) == 0 ) cycle
             io = index_local_to_global(dit,lio,Node)

             ia = iaorb(io,lasto) ! atom-index

             lio_connect: do j = 1 , l_ncol(lio)
                
                ind = l_ptr(lio) + j
                
                if ( .not. is_Trace ) then ! TS_W == TS_W_SUM_ATOM_?
                   
                   ia2 = iaorb(l_col(ind),lasto)
                   ! Only allow the same atom to contribute
                   if ( ia2 /= ia ) cycle lio_connect

                else ! TS_W == TS_W_TR_ATOM_?
                   
                   ! This is a SC sparsity pattern

                   ! Only allow the diagonal entry of
                   ! the density matrix
                   if ( l_col(ind) /= io ) cycle lio_connect

                end if

                if ( is_correlated ) then
                   atom_neq(:,ia) = atom_neq(:,ia) + DMneq(ind,:)
                else
                   atom_neq(:,ia) = atom_neq(:,ia) + DMneq(ind,:) ** 2
                end if

             end do lio_connect
          end do

#ifdef MPI
          allocate(atom_w(N_nEq_ID,na_u))
          ! We need to reduce the things
          call MPI_AllReduce(atom_neq(1,1),atom_w(1,1),N_nEq_ID*na_u, &
               MPI_Double_Precision, MPI_Sum, MPI_Comm_World, MPIerror)
          if ( is_correlated ) then
             atom_neq(:,:) = atom_w(:,:) ** 2
          else
             atom_neq(:,:) = atom_w(:,:)
          end if
          deallocate(atom_w)
#else
          if ( is_correlated ) then
             atom_neq(:,:) = atom_neq(:,:) ** 2
          end if
#endif

          ! in case of Bulk .and. DM_CrossTerms we
          ! can set the equivalent atom_w to 2 * maximum value
          ! This will force the nearest electrode to contribute the
          ! most. :)
          if ( N_Elec > 2 ) then
             call die('Check this code segment')
          end if
          tmp = maxval(atom_neq)
          allocate(atom_w(N_mu,na_u))
          do ia = 1 , na_u
             do io = 1 , N_Elec
                ! If we DO NOT use bulk electrodes we
                ! do have access to the diagonal correction
                ! contribution. Hence we only overwrite the electrode
                ! weight if Elec%Bulk
                if ( .not. Elecs(io)%Bulk ) cycle
                ! if we are not in the electrode we do not correct weight
                if ( .not. AtomInElec(Elecs(io),ia) ) cycle
                ! in case of sum with the off-diagonal terms we have to sum (otherwise it should be (:,ia) = tmp ; (mu%ID,ia) = 0._dp) 
                atom_neq(:,ia) = atom_neq(:,ia) + tmp
                atom_neq(Elecs(io)%mu%ID,ia) = atom_neq(Elecs(io)%mu%ID,ia) - tmp
             end do

             ! Calculate weights for this atom
             call get_weight(N_Elec,N_mu,N_nEq_ID,ID_mu, &
                  atom_neq(:,ia), atom_w(:,ia) )
             !if (node==0) &
             !     write(*,'(a,i2,2(tr1,g10.5))')'W: ', ia,atom_w(:,ia)

          end do

          ! clean up
          deallocate(atom_neq)

       end if

       do io = 1 , nr
          ! We are in a buffer region...
          if ( l_ncol(io) == 0 ) cycle

          ! Update the weight of the row-atom
          if ( TS_W /= TS_W_ORB_ORB ) then
             ! Calculate the weight of the atom corresponding to this
             ! orbital
             lio = index_local_to_global(dit,io,Node)
             ia1 = iaorb(lio,lasto)
          end if

          do j = 1 , l_ncol(io)

             ind = l_ptr(io) + j
             ! Retrieve the connecting orbital
             jo = l_col(ind)

             if ( TS_W == TS_W_ORB_ORB ) then

                ! Get the non-equilibrium contribution and the weight associated
                call get_neq_weight(N_Elec,N_mu,N_nEq_ID,ID_mu, &
                     DMneq(ind,:),neq,w)

             else ! we have weight per atom "somewhere"

                ! To compare the weights... For DEBUGging purposes...
                !call get_neq_weight(N_Elec,N_mu,N_nEq_ID,ID_mu, &
                !     DMneq(ind,:),neq,w)
                !write(*,'(a,i2,tr1,i2,4(tr1,g10.5))')'Wi: ', ia1,ia2,w
                
                ! Re-calculate the weight for special weighting...
                ia2 = iaorb(jo,lasto)

                ! [[ this is the geometric mean method...
                ! This should probably be used as the geometric mean
                ! retains the tendency of the data, however I am not 
                ! quite sure of its arguments...
                ! For test:
                ! A1:    [ 0.95  0.05 ]   # weight 1
                ! A2:    [ 0.5   0.5  ]   # weight 2
                ! GM-N:  [ 0.813 0.187]   # geometric mean of weights
                ! AM-N:  [ 0.725 0.275]   # arithmetic mean of weights

                w = sqrt(atom_w(:,ia1) * atom_w(:,ia2))
                   
                ! geometric mean does not retain normalization
                w = w / sum(w) 
                ! ]]

                ! [[ arithmetic mean
                ! the mean value between the atomic weights
                ! will be used as the actual weight
                ! w = (atom_w(:,ia1) + atom_w(:,ia2)) * 0.5_dp
                ! ]] 

                select case ( TS_W )                   
                case ( TS_W_TR_ATOM_ATOM , TS_W_SUM_ATOM_ATOM )

                   ! do nothing...

                case ( TS_W_TR_ATOM_ORB , TS_W_SUM_ATOM_ORB )
                   
                   ! this ensures that the atomic weights are
                   ! both taken into account, as well as the orb-orb
                   call get_weight(N_Elec,N_mu,N_nEq_ID,ID_mu, &
                        DMneq(ind,:) ** 2, neq(:) )

                   ! see above for arguments
                   w = sqrt(w(:) * neq(:))
                   w = w / sum(w) 

                case default
                   call die('Error in weights...')
                end select

                ! Get the non-equilibrium contribution
                call get_neq(N_mu,N_nEq_ID,ID_mu,DMneq(ind,:),neq)

                !write(*,'(a,i2,tr1,i2,4(tr1,g10.5))')'Wt: ', ia1,ia2,w,sum(w)
             end if
             
             ! Do error estimation (capture before update)
             ee = 0._dp
             do mu_i = 1 , N_mu - 1
                ee_i = DM(ind,mu_i) + neq(mu_i)
                do mu_j = mu_i + 1 , N_mu
                   tmp = ee_i - DM(ind,mu_j) - neq(mu_j)
                   if ( abs(tmp) > abs(ee) ) ee = tmp
                end do
             end do
             
             DM(ind,1) = w(1) * ( DM(ind,1) + neq(1) )
             if ( hasEDM ) EDM(ind,1) = w(1) * EDM(ind,1)
             do mu_i = 2 , N_mu
                DM(ind,1) = DM(ind,1) + &
                     w(mu_i) * ( DM(ind,mu_i) + neq(mu_i) )
                if ( hasEDM ) EDM(ind,1) = EDM(ind,1) + &
                     w(mu_i) * EDM(ind,mu_i)
             end do

             if ( abs(ee) > abs(eM) ) then
                eM   = ee
                eM_i = io
                eM_j = jo
                DMe = DM(ind,1)
             end if

          end do
       end do

       if ( TS_W /= TS_W_ORB_ORB ) then
          deallocate(atom_w)
       end if

    end if

    deallocate(ID_mu)

#ifdef MPI
    ! remove pointer
    nullify(DM,EDM)
    allocate(DM(Nodes,4),EDM(Nodes,4))

    io = Node + 1
    ! Initialize
    DM(:,:)  = 0._dp
    DM(io,1) = eM
    DM(io,2) = real(index_local_to_global(dit,eM_i,Node),dp)
    DM(io,3) = real(eM_j,dp)
    DM(io,4) = DMe
    call MPI_Reduce(DM(1,1),EDM(1,1),4*Nodes, &
         MPI_Double_Precision, MPI_Sum, 0, MPI_Comm_World, MPIerror)
    if ( IONode ) then
       io   = maxloc(abs(EDM(:,1)),1)
       eM   = EDM(io,1)
       eM_i = nint(EDM(io,2))
       eM_j = nint(EDM(io,3))
       DMe  = EDM(io,4)
    endif
    deallocate(DM,EDM)
#endif

    call print_error_estimate(IONode,'ts: int. EE.:', &
         eM,eM_i,eM_j,DMe)

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
       write(*,'(a,2(tr1,a,''('',i5,'','',i5,'')'',a,g10.5e1))') &
            trim(a), &
            'DM_out', eM_i,eM_j,' = ',DM, &
            ', d_ij', eM_i,eM_j,' = ',eM
    end if

  end subroutine print_error_estimate


  ! Calculate the weight for a given set of variances
  subroutine get_weight(N_El,N_mu,N_id,ID_mu,w_ID,w)
    use m_ts_contour_neq, only : ID2mu, ID2mult
    integer,  intent(in)  :: N_El, N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: w_ID(N_id)
    real(dp), intent(out) :: w(N_mu)
    integer :: mu_i, mu, ID
    real(dp) :: total, tmp

    !if ( any(w_ID < 0._dp) ) call die('get_weight: Error in code')

    ! TODO check that this is correct for several electrodes
    total  = 0._dp
    w(:)   = 0._dp
    do ID = 1 , N_id
       mu = ID_mu(ID)
       tmp = w_ID(ID) * ID2mult(ID)
       total = total + tmp
       do mu_i = 1 , N_mu
          if ( mu_i == mu ) cycle
          w(mu_i) = w(mu_i) + tmp
       end do
    end do

    total = total * ( N_El - 1 )
    
    if ( total > 0._dp ) then
       w(:) = w(:) / total
    else
       w = 1._dp / N_mu
    end if

  end subroutine get_weight
  
  ! Calculate both the non-equilibrium contribution
  ! and the weight associated with those points
  subroutine get_neq_weight(N_El,N_mu,N_id,ID_mu,neq_ID,neq,w)
    use m_ts_contour_neq, only : ID2mu, ID2mult
    integer,  intent(in)  :: N_El, N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: neq_ID(N_id)
    real(dp), intent(out) :: neq(N_mu)
    real(dp), intent(out) :: w(N_mu)
    integer :: mu_i, ID, mu
    real(dp) :: total, tmp, mult

    ! TODO check that this is correct for several electrodes
    total  = 0._dp
    neq(:) = 0._dp
    w(:)   = 0._dp
    do ID = 1 , N_id
       mu   = ID2mu(ID)
       mult = ID2mult(ID)
       neq(mu) = neq(mu) + neq_ID(ID) * mult
       tmp = neq_ID(ID) ** 2 * mult
       total = total + tmp
       do mu_i = 1 , N_mu
          if ( mu == mu_i ) cycle
          w(mu_i) = w(mu_i) + tmp
       end do
    end do

    total = total * ( N_El - 1 )

    if ( total > 0._dp ) then
       w(:) = w(:) / total
    else
       w = 1._dp / N_mu
    end if

  end subroutine get_neq_weight

  ! Calculate both the non-equilibrium contribution
  ! and the weight associated with those points
  subroutine get_neq(N_mu,N_id,ID_mu,neq_ID,neq)
    use m_ts_contour_neq, only : ID2mu, ID2mult
    integer,  intent(in)  :: N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: neq_ID(N_id)
    real(dp), intent(out) :: neq(N_mu)
    integer :: ID, mu

    ! TODO check that this is correct for several electrodes
    neq(:) = 0._dp
    do ID = 1 , N_id
       mu      = ID2mu(ID)
       neq(mu) = neq(mu) + neq_ID(ID) * ID2mult(ID)
    end do

  end subroutine get_neq

end module m_ts_weight
