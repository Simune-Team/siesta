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
!   do a weight for each k
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
  integer, parameter :: TS_W_K_UNCORRELATED = 2

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
  ! 6) Simple mean of contributions
  integer, parameter :: TS_W_MEAN = 6

  ! The general weighting can be:
  !   UNCORRELATED or CORRELATED
  ! The correlated is the default
  integer, parameter :: TS_W_CORRELATED = 100

  integer, save :: TS_W_METHOD = TS_W_ORB_ORB

  ! we default weight of uncorrelated as that should be 
  ! the most accurate
  integer, save :: TS_W_K_METHOD = TS_W_K_UNCORRELATED

contains

  subroutine read_ts_weight( )

    use fdf, only : fdf_get, leqi
    character(len=200) :: chars
    integer :: i

    ! Update the weight function
    chars = fdf_get('TS.Weight.k.Method','correlated')
    if ( leqi(chars,'correlated') ) then
       TS_W_K_METHOD = TS_W_K_CORRELATED
    else if ( leqi(chars,'uncorrelated') ) then
       TS_W_K_METHOD = TS_W_K_UNCORRELATED
    else
       call die('Could not determine flag TS.Weight.k.Method, &
            &please see manual.')
    end if

    ! The default weighting method is correlated if
    ! atom-atom is utilised
    TS_W_METHOD = TS_W_CORRELATED
    chars = fdf_get('TS.Weight.Method','orb-orb')
    ! first check whether we have correlated weighting
    i = index(chars,'+')
    if ( i > 0 ) then
       ! we do have something else
       if ( leqi(chars(1:i-1),'correlated') .or. &
            leqi(chars(1:i-1),'corr') ) then
          TS_W_METHOD = TS_W_CORRELATED
       else if ( leqi(chars(1:i-1),'uncorrelated') .or. &
            leqi(chars(1:i-1),'uncorr') ) then
          TS_W_METHOD = 0 ! non-correlated
       else
          call die('Unrecognized second option for TS.Weight.Method &
               &must be [[un]correlated+][orb-orb|tr-atom-atom|sum-atom-atom|mean]')
       end if
       chars = chars(i+1:)
    end if
    if ( leqi(chars,'orb-orb') ) then
       TS_W_METHOD = TS_W_ORB_ORB
       ! this does not make sense to make correlated, hence always assign
    else if ( leqi(chars,'tr-atom-atom') ) then
       TS_W_METHOD = TS_W_METHOD + TS_W_TR_ATOM_ATOM 
    else if ( leqi(chars,'tr-atom-orb') ) then
       TS_W_METHOD = TS_W_METHOD + TS_W_TR_ATOM_ORB
    else if ( leqi(chars,'sum-atom-atom') ) then
       TS_W_METHOD = TS_W_METHOD + TS_W_SUM_ATOM_ATOM 
    else if ( leqi(chars,'sum-atom-orb') ) then
       TS_W_METHOD = TS_W_METHOD + TS_W_SUM_ATOM_ORB
    else if ( leqi(chars,'mean') ) then
       TS_W_METHOD = TS_W_MEAN
    else
       call die('Unrecognized option for TS.Weight.Method &
            &must be [[un]correlated+|][orb-orb|tr-atom-[atom|orb]|sum-atom-[atom|orb]|mean]')
    end if

  end subroutine read_ts_weight

  subroutine weight_DM(N_Elec,Elecs,N_mu, na_u, lasto, &
       spDM, spDMneq, spEDM, n_s, sc_off)
    
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
    ! Number of supercells
    integer, intent(in) :: n_s
    ! the offsets
    real(dp), intent(in) :: sc_off(3,0:n_s-1)

! *********************
! * LOCAL variables   *
! *********************
    real(dp) :: w(N_mu), neq(N_mu)
    real(dp), parameter :: EPS = 1.e-4_dp
    
    ! arrays for looping in the sparsity pattern
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: DM(:,:), DMneq(:,:), EDM(:,:)
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: nr, n_nzs
    integer :: io, jo, ind, j, is
    integer :: mu_i, mu_j
    integer, allocatable :: ID_mu(:)
    ! For error estimation
    integer  :: eM_i, eM_j
    real(dp) :: eM, DMe, ee, tmp, m_err, e_f, ew
    logical :: hasEDM, is_correlated, is_trace
    ! collecting the error contribution for each atom
    real(dp), allocatable :: atom_w(:,:), atom_neq(:,:)
    integer :: ng, lio, ia1, ia2, ia, TS_W
    
#ifdef MPI
    integer :: MPIerror
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDM' )
#endif

    ! TODO Enforce that sparsity is the same
    ! (however, we know that they are the same)
    sp => spar(spDM)
    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=nr,nrows_g=ng,nnzs=n_nzs)

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
    ew  = 0._dp
    m_err = 0._dp

    ! Is the data correlated
    is_correlated = TS_W_METHOD >= TS_W_CORRELATED
    TS_W = TS_W_METHOD
    if ( is_correlated ) TS_W = TS_W - TS_W_CORRELATED

    if ( TS_W /= TS_W_ORB_ORB .and. TS_W /= TS_W_MEAN ) then
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

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ia,j,ind,ia2,is), &
!$OMP&reduction(+:atom_neq)
       do lio = 1 , nr

          ! We are in a buffer region...
          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)
          
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
                is = (l_col(ind)-1) / ng
                ! Check the unit-cell offset
                if ( sum(abs(sc_off(:,is))) > EPS ) cycle lio_connect
                
             end if

             if ( is_correlated ) then
                atom_neq(:,ia) = atom_neq(:,ia) + DMneq(ind,:)
             else
                atom_neq(:,ia) = atom_neq(:,ia) + DMneq(ind,:) ** 2
             end if
             
          end do lio_connect
          end if
       end do
!$OMP end parallel do
       
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

       ! in case of Bulk or DM_update /= update all we
       ! can set the equivalent atom_w to 2 * maximum value
       ! This will force the nearest electrode to contribute the
       ! most. :)
       tmp = maxval(atom_neq)
       allocate(atom_w(N_mu,na_u))
       l_atom: do ia = 1 , na_u
          do io = 1 , N_Elec
             ! If we DO NOT use bulk electrodes we
             ! do have access to the diagonal correction
             ! contribution. Hence we only overwrite the electrode
             ! weight if Elec%Bulk
             if ( (.not. Elecs(io)%Bulk) .and. Elecs(io)%DM_update /= 2 ) cycle
             ! if we are not in the electrode we do not correct weight
             if ( .not. AtomInElec(Elecs(io),ia) ) cycle
             ! in case of sum with the off-diagonal terms we have to sum (otherwise it should be (:,ia) = tmp ; (mu%ID,ia) = 0._dp) 
             atom_w(:,ia) = 0._dp
             atom_w(Elecs(io)%mu%ID,ia) = 1._dp

             !if (node==0) &
             !     write(*,'(a,i2,2(tr1,g10.5))')'W: ', ia,atom_w(:,ia)

             cycle l_atom
             
          end do

          ! Calculate weights for this atom
          call calc_weight(N_mu,N_nEq_ID,ID_mu, &
               atom_neq(:,ia), atom_w(:,ia) )
          !if (node==0) &
          !     write(*,'(a,i2,4(tr1,g10.5))')'W: ', ia,atom_w(:,ia),atom_neq(:,ia)

       end do l_atom

       ! clean up
       deallocate(atom_neq)
       
    end if

    if ( TS_W == TS_W_MEAN ) then
       ! The weight will always be divided
       w(:) = 1._dp / real(N_mu,dp)
    end if
    
! DM is accessed individually, so we will never have a data race
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ia1,ia2,j,ind,jo,neq), &
!$OMP&private(ee,e_f,mu_i,mu_j,tmp), &
!$OMP&firstprivate(w), reduction(+:m_err)
    do lio = 1 , nr
       ! We are in a buffer region...
       if ( l_ncol(lio) /= 0 ) then

       ! The global orbital
       io = index_local_to_global(dit,lio)
       
       ! Update the weight of the row-atom
       if ( TS_W /= TS_W_ORB_ORB ) then
          ! Calculate the weight of the atom corresponding to this
          ! orbital
          ia1 = iaorb(io,lasto)
       end if
       
       do j = 1 , l_ncol(lio)
          
          ind = l_ptr(lio) + j
          ! Retrieve the connecting orbital
          jo = l_col(ind)
          
          if ( TS_W == TS_W_ORB_ORB ) then

             ! Get the non-equilibrium contribution and the weight associated
             call calc_neq_weight(N_mu,N_nEq_ID,ID_mu, &
                  DMneq(ind,:),neq,w)

          else if ( TS_W == TS_W_MEAN ) then

             ! "w" already set
             ! Get the non-equilibrium contribution
             call calc_neq(N_mu,N_nEq_ID,ID_mu,DMneq(ind,:),neq)

          else ! we have weight per atom "somewhere"

             ! To compare the weights... For DEBUGging purposes...
             !call get_neq_weight(N_mu,N_nEq_ID,ID_mu, &
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
                call calc_weight(N_mu,N_nEq_ID,ID_mu, &
                     DMneq(ind,:) ** 2, neq(:) )
                
                ! see above for arguments
                w = sqrt(w(:) * neq(:))
                w = w / sum(w) 
                
             case default
                call die('Error in weights...')
             end select
             
             ! Get the non-equilibrium contribution
             call calc_neq(N_mu,N_nEq_ID,ID_mu,DMneq(ind,:),neq)
             
             !write(*,'(a,i2,tr1,i2,4(tr1,g10.5))')'Wt: ', ia1,ia2,w,sum(w)
          end if
             
#ifdef TRANSIESTA_WEIGHT_DEBUG
          if ( io == ucorb(jo,ng) .and. io == 28 ) then
             print '(2(a7,3(tr1,f10.5)))','Left',DM(ind,1),neq(1),w(1), &
                  'Right',DM(ind,2),neq(2),w(2)
          end if
#endif

          ! Calculate each contribution
          do mu_i = 1 , N_mu
             DM(ind,mu_i) = DM(ind,mu_i) + neq(mu_i)
          end do

          ! Do error estimation (capture before update)
          ee = 0._dp
          do mu_i = 1 , N_mu - 1
             do mu_j = mu_i + 1 , N_mu
                tmp = DM(ind,mu_i) - DM(ind,mu_j)
                ! Calculate sum of all errors
                m_err = m_err + tmp
                if ( abs(tmp) > abs(ee) ) ee = tmp
             end do
          end do
          
          ! Store for later estimation of the "final" error
          e_f = DM(ind,1)
          
          DM(ind,1) = w(1) * DM(ind,1)
          if ( hasEDM ) EDM(ind,1) = w(1) * EDM(ind,1)
          do mu_i = 2 , N_mu
             DM(ind,1) = DM(ind,1) + w(mu_i) * DM(ind,mu_i)
             if ( hasEDM ) &
                  EDM(ind,1) = EDM(ind,1) + w(mu_i) * EDM(ind,mu_i)
          end do

          ! Calculate error from estimated density
          e_f = e_f - DM(ind,1)
          do mu_i = 2 , N_mu
             tmp = DM(ind,mu_i) - DM(ind,1)
             if ( abs(tmp) > abs(e_f) ) e_f = tmp
          end do
          
          if ( abs(ee) > abs(eM) ) then
!$OMP critical
             eM   = ee
             eM_i = io
             eM_j = jo
             DMe = DM(ind,1)
             ew   = e_f
!$OMP end critical
          end if

       end do
       end if
    end do
!$OMP end parallel do

    ! Calculate mean of mean difference
    ! First calculate number of differences used
    io = 0
    do mu_i = 1, N_mu - 1
       do mu_j = mu_i + 1 , N_mu
          io = io + 1
       end do
    end do
    m_err = m_err / real(io,dp)
    
    if ( TS_W /= TS_W_ORB_ORB .and. TS_W /= TS_W_MEAN ) then
       deallocate(atom_w)
    end if

    deallocate(ID_mu)

#ifdef MPI
    if ( Nodes > 1 ) then
    ! remove pointer
    nullify(DM,EDM)
    allocate(DM(7,Nodes),EDM(7,Nodes))

    io = Node + 1
    ! Initialize
    DM(:,:)  = 0._dp
    DM(1,io) = eM
    DM(2,io) = real(eM_i,dp)
    DM(3,io) = real(eM_j,dp)
    DM(4,io) = DMe
    DM(5,io) = m_err
    DM(6,io) = n_nzs
    DM(7,io) = ew
    call MPI_Reduce(DM(1,1),EDM(1,1),7*Nodes, &
         MPI_Double_Precision, MPI_Sum, 0, MPI_Comm_World, MPIerror)
    if ( IONode ) then
       io   = maxloc(abs(EDM(1,:)),1)
       eM   = EDM(1,io)
       eM_i = nint(EDM(2,io))
       eM_j = nint(EDM(3,io))
       DMe  = EDM(4,io)
       ! Sum of errors
       m_err = sum( EDM(5,:) )
       n_nzs = nint( sum( EDM(6,:) ) )
       ew = EDM(7,io)
    end if
    deallocate(DM,EDM)
    end if
#endif

    ! Calculate mean
    m_err = m_err / real(n_nzs,dp)

    call print_error_estimate(IONode,'ts-EE:', &
         eM,ew,eM_i,eM_j,DMe,m_err)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS weightDM' )
#endif

  end subroutine weight_DM

  ! ************* Commonly used modules ****************

  ! Write out the error-estimate for the current iteration (or, k-point)
  subroutine print_error_estimate(IONode,a,eM,ew,eM_i,eM_j,DM,m_err)
    logical, intent(in)  :: IONode
    character(len=*), intent(in) :: a
    real(dp), intent(in) :: eM, ew, DM, m_err
    integer, intent(in)  :: eM_i,eM_j

    if ( IONode ) then
       write(*,'(a,tr1,a,i5,'','',i6,3(a,g10.5e1)&
            &,a,g11.5)') trim(a), &
            'ij(',eM_i,eM_j, '), DM = ',DM, &
            ', ew = ',ew, ', em = ',eM, &
            '. avg_m = ',m_err
    end if

  end subroutine print_error_estimate


  ! Calculate the theta values
  pure subroutine calc_theta(N_mu,N_id,ID_mu,w_ID,theta)
    integer,  intent(in)  :: N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: w_ID(N_id)
    real(dp), intent(out) :: theta(N_mu)
    integer :: ID

    ! TODO check that this is correct for several electrodes
    theta(:) = 0._dp
    do ID = 1 , N_id
       theta(ID_mu(ID)) = theta(ID_mu(ID)) + w_ID(ID)
    end do

  end subroutine calc_theta

  ! Calculate the theta values
  pure subroutine calc_weight(N_mu,N_id,ID_mu,w_ID,w)
    integer,  intent(in)  :: N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: w_ID(N_id)
    real(dp), intent(out) :: w(N_mu)
    real(dp) :: theta(N_mu), tmp
    integer :: i

    call calc_theta(N_mu,N_id,ID_mu,w_ID,theta)

    w(:) = product(theta)
    do i = 1 , N_mu
       if ( theta(i) > 0._dp ) then
          w(i) = w(i) / theta(i)
       else
          w(i) = 0._dp
       end if
    end do

    ! The denominator
    tmp = sum(w)
    if ( tmp == 0._dp ) then
       w = 1._dp / real(N_mu,dp)
    else
       w = w / tmp
    end if

  end subroutine calc_weight

  ! Calculate both the non-equilibrium contribution
  ! and the weight associated with those points
  pure subroutine calc_neq_weight(N_mu,N_id,ID_mu,neq_ID,neq,w)
    integer,  intent(in)  :: N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: neq_ID(N_id)
    real(dp), intent(out) :: neq(N_mu)
    real(dp), intent(out) :: w(N_mu)
    integer :: i
    real(dp) :: tmp

    ! TODO check that this is correct for several electrodes
    call calc_theta(N_mu,N_id,ID_mu,neq_ID**2,neq)
    w(:) = product(neq)
    do i = 1 , N_mu
       if ( neq(i) > 0._dp ) then
          w(i) = w(i) / neq(i)
       else
          w(i) = 0._dp
       end if
    end do

    ! The denominator
    tmp = sum(w)
    if ( tmp == 0._dp ) then
       w = 1._dp / real(N_mu,dp)
    else
       w = w / tmp
    end if

    neq(:) = 0._dp
    do i = 1 , N_id
       neq(ID_mu(i)) = neq(ID_mu(i)) + neq_ID(i)
    end do

  end subroutine calc_neq_weight

  ! Calculate both the non-equilibrium contribution
  ! and the weight associated with those points
  pure subroutine calc_neq(N_mu,N_id,ID_mu,neq_ID,neq)
    integer,  intent(in)  :: N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: neq_ID(N_id)
    real(dp), intent(out) :: neq(N_mu)
    integer :: ID

    ! TODO check that this is correct for several electrodes
    neq(:) = 0._dp
    do ID = 1 , N_id
       neq(ID_mu(ID)) = neq(ID_mu(ID)) + neq_ID(ID)
    end do

  end subroutine calc_neq

end module m_ts_weight
