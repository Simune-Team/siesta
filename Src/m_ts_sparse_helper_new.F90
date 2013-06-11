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

! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************


! In general many of the loops below can be followed by examining this loop:

!    ! This loop is across the local rows...
!    do lio = 1 , lnr
!
!       ! Quickly go past the empty regions... (we have nothing to update)
!       if ( l_ncol(lio) == 0 ) cycle
!
!       ! obtain the global index of the local orbital.
!       io = index_local_to_global(dit,lio,Node)
!
!       ! Quickly go past the empty regions... (we have nothing to update)
!       if ( up_ncol(io) == 0 ) cycle
!
!       ! Do a loop in the local sparsity pattern...
!       ! The local sparsity pattern is more "spread", hence
!       ! we do fewer operations by having this as an outer loop
!       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
!
!          jo = UCORB(l_col(lind),nr)
!
!          ! Now search the update region
!          ! This one must *per definition* have less elements.
!          ! Hence, we can exploit this, and find equivalent
!          ! super-cell orbitals.
!          rind = up_ptr(io)
!          ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
!          if ( ind <= rind ) cycle ! The element does not exist
!
!          ! Obtain the phase for the orbital ij
!          kx = k(1) * xij(1,pnt(lind)) + &
!               k(2) * xij(2,pnt(lind)) + &
!               k(3) * xij(3,pnt(lind))
!
!          ph = fact * cdexp(dcmplx(0._dp,-1._dp)*kx)
!
!          ! The integration is this:
!          ! \rho = e^{-i.k.R} [ \int (Gf^R-Gf^A)/2 dE + \int Gf^R\Gamma Gf^A dE ]
!
!          if ( non_Eq ) then
!
!             ! The integration is this:
!             ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
!             kx = aimag( ph*zDu(ind) )
!             dD(lind) = dD(lind) + kx
!             dE(lind) = dE(lind) + aimag( ph*zEu(ind) )
!
!          else
!
!             ! The fact that we have a SYMMETRIC
!             ! update region makes this *tricky* part easy...
!             rin  = up_ptr(jo)
!             rind = rin + SFIND(up_col(rin+1:rin+up_ncol(jo)),io)
!             ! We do a check, just to be sure...
!             if ( rind <= rin ) &
!                  call die('ERROR: Conjugated symmetrization point does not exist')
!
!             ! This integration is this:
!             ! \rho = e^{-i.k.R} \int (Gf^R-Gf^A)/2 dE
!             
!             dD(lind) = dD(lind) + aimag( ph*(zDu(ind) - conjg(zDu(rind))) )
!             dE(lind) = dE(lind) + aimag( ph*(zEu(ind) - conjg(zEu(rind))) )
!
!          end if
!             
!       end do
!    end do
    

module m_ts_sparse_helper_new

  use precision, only : dp

  implicit none
  
  private :: dp

contains

  ! ***
  ! The following scheme should be followed:
  !   add_*_DM routines are constructed to be able to handle different schemes
  !   They are called after each k-point and thus the arguments are the following:
  !     1. the local update sparsity pattern
  !     2. the global update sparsity pattern (suffix 'u')
  ! The k-point routine is constructed to handle three different methods of doing
  ! the weighting which is performed here.

  subroutine add_k_DM(dit,spDM,spEDM,spDMu,spEDMu,k,ipnt,n_nzs,xij,&
       non_Eq,spW)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_iSpData1D
    use class_dSpData1D
    use class_zSpData1D
    use intrinsic_missing, only : SFIND, UCORB => MODP
    use parallel, only : Node

! *********************
! * INPUT variables   *
! *********************
    type(OrbitalDistribution), intent(inout) :: dit
    ! The local integrated sparsity arrays
    type(dSpData1D), intent(inout) :: spDM, spEDM
    ! The current k-point global sparsity arrays
    type(zSpData1D), intent(inout) :: spDMu, spEDMu
    ! The k-point
    real(dp), intent(in) :: k(3)
    ! The pointer from xij -> spar(spDM).
    ! I.e. a pointer from the local update sparsity to the local sparsity
    ! (only needed to refrain from creating a duplicate xij array)
    type(iSpData1D), intent(inout) :: ipnt
    ! The orbitals distances (in the local SIESTA sparsity pattern)
    integer, intent(in) :: n_nzs
    real(dp), intent(in) :: xij(3,n_nzs)

    ! If the sparsity-weight is provided we will do this:
    ! DM = DM + DMu
    ! spW = DMu ** 2
    ! It MUST meen that we do TS_W_UNCORRELATED (see m_ts_weight)
    logical, intent(in), optional :: non_Eq
    type(dSpData1D), intent(inout), optional :: spW

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:)
    integer, pointer :: pnt(:)
    real(dp), pointer :: dD(:) , dE(:), dW(:)
    complex(dp), pointer :: zDu(:), zEu(:)
    integer :: lnr, lio, lind, io, ind, nr, jo
    integer :: rin, rind
    logical :: save_weight
    real(dp) :: kx
    complex(dp) :: ph

    l_s  => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD   => val(spDM)
    dE   => val(spEDM)

    up_s => spar(spDMu)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    zDu  => val(spDMu)
    zEu  => val(spEDMu)

    ! The pointer
    pnt  => val(ipnt)

    ! If the weight-array is clear, then save to that.
    save_weight = present(spW)
    if ( save_weight ) then
       ! We also need to capture the k-point weight
       dW => val(spW)
    end if
     
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    
    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected.')

    if ( save_weight .and. .not. non_Eq ) call die('Cannot save weight &
         &if not a non-equilibrium density')
    
    ! primary option (if we need to save the weight, then it must be 
    
    if ( save_weight ) then
       
       do lio = 1 , lnr

          if ( l_ncol(lio) == 0 ) cycle
          io = index_local_to_global(dit,lio,Node)
          if ( up_ncol(io) == 0 ) cycle

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)
             
             rind = up_ptr(io)
             ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
             if ( ind <= rind ) cycle ! The element does not exist

             kx = k(1) * xij(1,pnt(lind)) + &
                  k(2) * xij(2,pnt(lind)) + &
                  k(3) * xij(3,pnt(lind))

             ph = cdexp(dcmplx(0._dp,-1._dp)*kx)

             ! The integration is this:
             ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
             kx = aimag( ph*zDu(ind) )
             dD(lind) = dD(lind) + kx
             dE(lind) = dE(lind) + aimag( ph*zEu(ind) )

             ! Sum up the weight here
             dW(lind) = dW(lind) + kx ** 2

          end do
       end do

    else if ( non_Eq ) then

       do lio = 1 , lnr

          if ( l_ncol(lio) == 0 ) cycle
          io = index_local_to_global(dit,lio,Node)
          if ( up_ncol(io) == 0 ) cycle

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)
             
             rind = up_ptr(io)
             ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
             if ( ind <= rind ) cycle ! The element does not exist

             kx = k(1) * xij(1,pnt(lind)) + &
                  k(2) * xij(2,pnt(lind)) + &
                  k(3) * xij(3,pnt(lind))

             ph = cdexp(dcmplx(0._dp,-1._dp)*kx)

             ! The integration is this:
             ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
             dD(lind) = dD(lind) + aimag( ph*zDu(ind) )
             dE(lind) = dE(lind) + aimag( ph*zEu(ind) )

          end do
       end do

    else

       do lio = 1 , lnr

          if ( l_ncol(lio) == 0 ) cycle
          io = index_local_to_global(dit,lio,Node)
          if ( up_ncol(io) == 0 ) cycle

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)

             rind = up_ptr(io)
             ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
             if ( ind <= rind ) cycle ! The element does not exist

             kx = k(1) * xij(1,pnt(lind)) + &
                  k(2) * xij(2,pnt(lind)) + &
                  k(3) * xij(3,pnt(lind))

             ph = 0.5_dp * cdexp(dcmplx(0._dp,-1._dp)*kx)

             rin  = up_ptr(jo)
             rind = rin + SFIND(up_col(rin+1:rin+up_ncol(jo)),io)
             if ( rind <= rin ) &
                  call die('ERROR: Conjugated symmetrization point does not exist')

             ! This integration is this:
             ! \rho = e^{-i.k.R} \int (Gf^R-Gf^A)/2 dE
             dD(lind) = dD(lind) + aimag( ph*(zDu(ind) - conjg(zDu(rind))) )
             dE(lind) = dE(lind) + aimag( ph*(zEu(ind) - conjg(zEu(rind))) )

          end do
       end do

    end if

  end subroutine add_k_DM

  subroutine add_Gamma_DM(dit,spDM,spEDM,spDMu,spEDMu)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use m_ts_weight
    use intrinsic_missing, only : SFIND, UCORB => MODP
    use parallel, only : Node

! *********************
! * INPUT variables   *
! *********************
    type(OrbitalDistribution), intent(inout) :: dit
    ! The local integrated sparsity arrays
    type(dSpData1D), intent(inout) :: spDM, spEDM
    ! The current k-point global sparsity arrays
    type(dSpData1D), intent(inout) :: spDMu, spEDMu

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:)
    real(dp), pointer :: dD(:) , dE(:)
    real(dp), pointer :: dDu(:), dEu(:)
    integer :: lnr, lio, lind, io, ind, nr, jo
    integer :: rind

    l_s  => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD   => val(spDM)
    dE   => val(spEDM)

    up_s => spar(spDMu)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    dDu  => val(spDMu)
    dEu  => val(spEDMu)

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    
    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected.')

    do lio = 1 , lnr
       
       if ( l_ncol(lio) == 0 ) cycle
       io = index_local_to_global(dit,lio,Node)
       if ( up_ncol(io) == 0 ) cycle
       
       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
          jo = l_col(lind)
          
          rind = up_ptr(io)
          ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
          if ( ind <= rind ) cycle ! The element does not exist

          ! This integration is this:
          ! \rho = \int Re(Gf^R) dE
          dD(lind) = dDu(ind)
          dE(lind) = dEu(ind)

       end do
    end do

  end subroutine add_Gamma_DM

  subroutine update_DM(dit,sp,n_nzs,DM,EDM, spDM, spEDM, ipnt)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_iSpData1D
    use class_dSpData1D
    use parallel, only : Node

! *********************
! * INPUT variables   *
! *********************
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs), EDM(n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(dSpData1D), intent(inout) :: spDM, spEDM
    ! The pointer from xij -> spar(spDM).
    ! I.e. a pointer from the local update sparsity to the local sparsity
    ! (only needed to refrain from creating a duplicate xij array)
    type(iSpData1D), intent(inout), optional :: ipnt

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    type(OrbitalDistribution), pointer :: orb_dit
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
    integer, pointer :: pnt(:)
    real(dp), pointer :: dD(:), dE(:)
    integer :: lnr, uind, io, ind, nr, jo

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    dD => val(spDM)
    dE => val(spEDM)

    ! Check orbital distribution
    orb_dit => dist(spDM)
    if ( index_local_to_global(dit,1,Node) /= &
         index_local_to_global(orb_dit,1,Node) ) &
         call die('They are not the same sparsity pattern')
     
    if ( lnr /= nrows(s) ) call die('The sparsity format is not as &
         &expected.')

    if ( present(ipnt) ) then
       ! The pointer
       pnt  => val(ipnt)

       ! This loop is across the local rows...
       do io = 1 , lnr

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) == 0 ) cycle

          ! Do a loop in the local update sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)
             
             jo = lup_col(uind)
             ind = pnt(uind)

             DM(ind)  = DM(ind)  + dD(uind)
             EDM(ind) = EDM(ind) + dE(uind)
             
          end do
       end do

    else       

       ! This loop is across the local rows...
       do io = 1 , lnr

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) == 0 ) cycle

          ! Do a loop in the local update sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)

             jo = lup_col(uind)

             ! Now we loop across the local region
             ind = l_ptr(io)
             ind = l_ptr(io) + minloc(abs(l_col(ind+1:ind+l_ncol(io))-jo),1)
             if ( l_col(ind) /= jo ) then
                do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                   if ( l_col(ind) /= jo ) cycle
                   exit ! we have the correct ind-value
                end do
             end if

             ! We need to add in case of special weighting...
             DM(ind)  = DM(ind)  + dD(uind)
             EDM(ind) = EDM(ind) + dE(uind)

          end do
       end do
    end if

  end subroutine update_DM


#ifdef OLD

  subroutine update_DM_Bias(dit,sp,n_nzs,DM,EDM, spDM, spEDM)
      ! The DM and EDM equivalent matrices
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use parallel, only : Node
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs), EDM(n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(dSpData1D), intent(inout) :: spDM, spEDM

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    type(OrbitalDistribution), pointer :: orb_dit
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
    real(dp), pointer :: dD(:), dE(:)
    integer :: lnr, uind, io, ind, nr, jo

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    dD     => val(spDM)
    dE     => val(spEDM)

    ! Check orbital distribution
    orb_dit => dist(spDM)
    if ( index_local_to_global(dit,1,Node) /= &
         index_local_to_global(orb_dit,1,Node) ) &
         call die('They are not the same sparsity pattern')
     
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
     
    if ( lnr /= nrows(s) ) call die('The sparsity format is not as &
         &expected.')

    ! This loop is across the local rows...
    do io = 1 , lnr

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( lup_ncol(io) == 0 ) cycle

       ! Do a loop in the local update sparsity pattern...
       ! The local sparsity pattern is more "spread", hence
       ! we do fewer operations by having this as an outer loop
       do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)

          jo = lup_col(uind)

          ! Now we loop across the local region
          ind = l_ptr(io)
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
             if ( l_col(ind) /= jo ) cycle

             ! Probably we dont need to "add"
             ! We only have one k-point...
             DM(ind)  = DM(ind)  + dD(uind)
             EDM(ind) = EDM(ind) + dE(uind)
             exit
          end do

       end do
    end do
    
  end subroutine update_DM_Bias

  subroutine update_DM_new(dit,spDM,spEDM,spDMu,spEDMu)
      ! The DM and EDM equivalent matrices
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use intrinsic_missing, only : SFIND, UCORB => MODP
    use parallel, only : Node
    type(OrbitalDistribution), intent(inout) :: dit
    ! Updated sparsity arrays (they contain the full integration)
    type(dSpData1D), intent(inout) :: spDM, spEDM
    ! Updated sparsity arrays (they contain the current k-integration)
    type(dSpData1D), intent(inout) :: spDMu, spEDMu

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:)
    real(dp), pointer :: dD(:) , dE(:)
    real(dp), pointer :: dDu(:), dEu(:)
    integer :: lnr, lio, lind, io, ind, nr, ljo

    l_s => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD     => val(spDM)
    dE     => val(spEDM)

    up_s => spar(spDMu)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    dDu    => val(spDMu)
    dEu    => val(spEDMu)
     
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    
    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected.')
    
    ! This loop is across the local rows...
    do lio = 1 , lnr

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( l_ncol(lio) == 0 ) cycle

       ! obtain the global index of the local orbital.
       io = index_local_to_global(dit,lio,Node)

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( up_ncol(io) == 0 ) cycle

       ! Do a loop in the local sparsity pattern...
       ! The local sparsity pattern is more "spread", hence
       ! we do fewer operations by having this as an outer loop
       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ljo = UCORB(l_col(lind),nr)

          ! Now we loop across the update region
          ! This one must *per definition* have less elements.
          ! Hence, we can exploit this, and find equivalent
          ! super-cell orbitals.
          ! Ok, this is Gamma (but to be consistent)
          ind = up_ptr(io)
          ind = ind + SFIND(up_col(ind+1:ind+up_ncol(io)),ljo)
          if ( ind <= up_ptr(io) ) cycle
          call die('need to implement weight')

          ! Probably we dont need to "add"
          ! We only have one k-point...
          dD(lind) = dD(lind) + dDu(ind)
          dE(lind) = dE(lind) + dEu(ind)

       end do
    end do
    
  end subroutine update_DM_new


  subroutine update_zDM_new(dit,spDM,spEDM,spDMu,spEDMu,k,ipnt,n_nzs,xij, &
       spW)
      ! The DM and EDM equivalent matrices
    use class_OrbitalDistribution
    use class_Sparsity
    use class_iSpData1D
    use class_dSpData1D
    use class_zSpData1D
    use intrinsic_missing, only : SFIND, UCORB => MODP
    use parallel, only : Node
    type(OrbitalDistribution), intent(inout) :: dit
    ! Updated sparsity arrays (they contain the local full integration)
    type(dSpData1D), intent(inout) :: spDM, spEDM
    ! Updated sparsity arrays (they contain the current k-integration)
    type(zSpData1D), intent(inout) :: spDMu, spEDMu
    ! The k-point
    real(dp), intent(in) :: k(3)
    ! The pointer from xij -> spar(spDM).
    ! I.e. a pointer from the local update sparsity to the local sparsity
    type(iSpData1D), intent(inout) :: ipnt
    ! The orbitals distances
    integer, intent(in) :: n_nzs
    real(dp), intent(in) :: xij(3,n_nzs)
    ! Determine whether we have to also have the conjugated
    type(dSpData1D), intent(inout), optional :: spW

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:)
    integer, pointer :: pnt(:)
    real(dp), pointer :: dD(:) , dE(:), dW(:)
    complex(dp), pointer :: zDu(:), zEu(:)
    logical :: non_Eq
    integer :: lnr, lio, lind, io, ind, nr, jo
    integer :: rin, rind
    real(dp) :: fact, kx
    complex(dp) :: ph

    l_s => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD     => val(spDM)
    dE     => val(spEDM)

    up_s => spar(spDMu)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    zDu    => val(spDMu)
    zEu    => val(spEDMu)

    ! The pointer
    pnt    => val(ipnt)

    ! If the weight-array is clear, then save to that.
    non_Eq = present(spW)

    if ( non_Eq ) then
       fact = 1.0_dp
       ! We also need to capture the k-point weight
       dW => val(spW)
    else
       fact = 0.5_dp
    end if
     
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    
    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected.')
    
    ! This loop is across the local rows...
    do lio = 1 , lnr

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( l_ncol(lio) == 0 ) cycle

       ! obtain the global index of the local orbital.
       io = index_local_to_global(dit,lio,Node)

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( up_ncol(io) == 0 ) cycle

       ! Do a loop in the local sparsity pattern...
       ! The local sparsity pattern is more "spread", hence
       ! we do fewer operations by having this as an outer loop
       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          jo = UCORB(l_col(lind),nr)

          ! Now search the update region
          ! This one must *per definition* have less elements.
          ! Hence, we can exploit this, and find equivalent
          ! super-cell orbitals.
          rind = up_ptr(io)
          ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
          if ( ind <= rind ) cycle ! The element does not exist

          ! Obtain the phase for the orbital ij
          kx = k(1) * xij(1,pnt(lind)) + &
               k(2) * xij(2,pnt(lind)) + &
               k(3) * xij(3,pnt(lind))

          ph = fact * cdexp(dcmplx(0._dp,-1._dp)*kx)

          ! The integration is this:
          ! \rho = e^{-i.k.R} [ \int (Gf^R-Gf^A)/2 dE + \int Gf^R\Gamma Gf^A dE ]

          if ( non_Eq ) then

             ! The integration is this:
             ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
             kx = aimag( ph*zDu(ind) )
             dD(lind) = dD(lind) + kx
             dE(lind) = dE(lind) + aimag( ph*zEu(ind) )

             ! Sum up the weight here
             dW(lind) = dW(lind) + kx ** 2

          else

             ! The fact that we have a SYMMETRIC
             ! update region makes this *tricky* part easy...
             rin  = up_ptr(jo)
             rind = rin + SFIND(up_col(rin+1:rin+up_ncol(jo)),io)
             ! We do a check, just to be sure...
             if ( rind <= rin ) &
                  call die('ERROR: Conjugated symmetrization point does not exist')

             ! This integration is this:
             ! \rho = e^{-i.k.R} \int (Gf^R-Gf^A)/2 dE
             
             dD(lind) = dD(lind) + aimag( ph*(zDu(ind) - conjg(zDu(rind))) )
             dE(lind) = dE(lind) + aimag( ph*(zEu(ind) - conjg(zEu(rind))) )

          end if
             
       end do
    end do
    
  end subroutine update_zDM_new


! ##################################################################
! ##   Mixing the Density matrixes according to the smallest      ##
! ##    realspace integral                                        ##
! ##                                                              ##
! ##  Sparse version and corrected the error-estimation           ##
! ##                                                              ##
! ##################################################################
  subroutine weight_DM( &
       spDML, spDMR, spDMneqL, spDMneqR, &
       spEML, spEMR)
!  This routine find weight for the DM integral. On output
!  DML := w (DML+DMneqR) + (1-w) (DMR+DMneqL)
!  EML := w (EML+EMneqR) + (1-w) (EMR+EMneqL)
!     note that EML \equiv EML + EMneqR
!     note that EMR \equiv EMR + EMneqL
    use parallel,  only: IONode
    use class_Sparsity
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
    type(dSpData1D), intent(inout) :: spEML, spEMR

! *********************
! * LOCAL variables   *
! *********************
    real(dp) :: wL,wR,wSUM

    ! arrays for looping in the sparsity pattern
    type(Sparsity), pointer :: sp
    real(dp), pointer :: DML(:), DMR(:)
    real(dp), pointer :: DMneqL(:), DMneqR(:)
    real(dp), pointer :: EML(:), EMR(:)
    integer,  pointer :: l_ncol(:)
    integer,  pointer :: l_ptr(:)
    integer,  pointer :: l_col(:)
    integer :: nr
    integer :: io, jo, ind, j
    ! For error estimation
    integer  :: eM_i,eM_j,neM_i,neM_j
    real(dp) :: eM, neM, tmp

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDM' )
#endif

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
    EML    => val(spEML)
    EMR    => val(spEMR)

    ! initialize the errors
    eM  = -1._dp
    neM = -1._dp

    do io = 1 , nr
       ! We are in a buffer region...
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)

          ind = l_ptr(io) + j
          ! Retrieve the connecting orbital
          jo = l_col(ind)

          wL = DMneqL(ind)!*DMneqL(ind)
          wR = DMneqR(ind)!*DMneqR(ind)
          wSUM = wL + wR

          ! The weights
          if ( wSUM > 0._dp ) then
             wL = wL / wSUM
             wR = wR / wSUM
          else
             wL = 0.5_dp
             wR = 0.5_dp
          end if
          
          ! Do error estimation (capture before update)
          tmp = (DML(ind) - DMR(ind))**2
          !tmp = (DML(ind) + DMneqR(ind) - DMR(ind) - DMneqL(ind))**2

          DML(ind) = wL * DML(ind)  + wR * DMR(ind)
          !DML(ind) = wL * (DML(ind) + DMneqR(ind)) &
          !         + wR * (DMR(ind) + DMneqL(ind))
          ! EML \equiv EML + EMneqR
          ! EMR \equiv EMR + EMneqL
          EML(ind) = wL * EML(ind)  + wR * EMR(ind)

          ! this is absolute error
          if ( tmp > eM ) then
             eM   = tmp
             eM_i = io
             eM_j = jo
          end if
          ! this is normalized absolute error
          tmp = tmp * wL * wR
          if ( tmp > neM ) then
             neM   = tmp
             neM_i = io
             neM_j = jo
          end if

       end do
    end do

    call print_error_estimate(IONode,eM,eM_i,eM_j,neM,neM_i,neM_j)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS weightDM' )
#endif

  end subroutine weight_DM


  subroutine print_error_estimate(IONode,eM,eM_i,eM_j,neM,neM_i,neM_j)
    logical, intent(in) :: IONode
    real(dp), intent(in) :: eM, neM
    integer, intent(in) :: eM_i,eM_j, neM_i,neM_j

    if ( IONode ) then
       write(*,'(a,'' |('',i5,'','',i5,'')| = '',g10.5e1,&
            &'' , |('',i5,'','',i5,'')|~ = '',g10.5e1)') &
            'ts: int. EE.:',&
            eM_i,eM_j,sqrt(eM), &
            neM_i,neM_j,sqrt(neM)
    end if

  end subroutine print_error_estimate
#endif

end module m_ts_sparse_helper_new
