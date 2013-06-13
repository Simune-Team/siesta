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
!          ph = fact * cdexp(dcmplx(0._dp,-kx))
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
    

module m_ts_dm_update

  use precision, only : dp

  implicit none
  
  private :: dp

contains


  subroutine select_dE(cNEn,c, iPE, nspin, kw, Z, W, ZW)
    use precision, only: dp
    use units, only: Pi
    use parallel, only: Node, Nodes
    use m_ts_cctype
    integer, intent(in) :: cNEn
    type(ts_ccontour), intent(in) :: c(cNEn)
    integer, intent(in) :: iPE, nspin
    real(dp),intent(in) :: kw
    complex(dp), intent(out) :: Z, W, ZW

    integer :: iE

    ! obtain a valid energy point (truncate at NEn)
    iE = min(iPE,cNEn)

    ! save the current weight of the point
    ! This is where we include the factor-of-two for spin and
    ! and the (1/Pi) from DM = Im[G]/Pi
    ! Furthermore we include the weight of the k-point
    W = 1._dp/Pi*c(iE)%w * kw
    if ( nspin == 1 ) W = W * 2._dp

    ! save the contour energy point
    Z = c(iE)%c
    ! Save Z*W, used for E-arrays
    ZW = Z*W

  end subroutine select_dE


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

    if ( (.not. initialized(spDM)) .or. (.not. initialized(spDMu)) ) return

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

             ph = cdexp(dcmplx(0._dp,-kx))

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

             ph = cdexp(dcmplx(0._dp,-kx))

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

             ph = 0.5_dp * cdexp(dcmplx(0._dp,-kx))

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

    if ( (.not. initialized(spDM)) .or. (.not. initialized(spDMu)) ) return

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

  subroutine update_zDM(dit,sp,n_nzs,DM,EDM,xij,spDM, spEDM,k)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D

    use intrinsic_missing, only : SFIND, UCORB => MODP
    use parallel, only : Node
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs), EDM(n_nzs)
    ! The orbital distances
    real(dp), intent(in) :: xij(3,n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(zSpData1D), intent(inout) :: spDM, spEDM
    ! The k-point...
    real(dp), intent(in) :: k(3)

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
    complex(dp), pointer :: zD(:), zE(:)
    complex(dp) :: ph
    real(dp) :: kx
    integer :: lio, io, jo, ind, nr
    integer :: lnr, lind, rin, rind

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    zD     => val(spDM)
    zE     => val(spEDM)
     
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
     
    if ( nr /= nrows(s) ) call die('The sparsity format is not as &
         &expected.')
     
    ! This loop is across the local rows...
    do lio = 1 , lnr

       ! obtain the global index of the local orbital.
       io = index_local_to_global(dit,lio,Node)

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( lup_ncol(io) == 0 ) cycle

       ! Do a loop in the local sparsity pattern...
       ! The local sparsity pattern is more "spread", hence
       ! we do fewer operations by having this as an outer loop
       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          jo = UCORB(l_col(lind),nr)

          ! Now search the update region
          ! This one must *per definition* have less elements.
          ! Hence, we can exploit this, and find equivalent
          ! super-cell orbitals.
          rind = lup_ptr(io)
          ind = rind + SFIND(lup_col(rind+1:rind+lup_ncol(io)),jo)
          if ( ind <= rind ) cycle ! The element does not exist
          
          kx = k(1) * xij(1,lind) + &
               k(2) * xij(2,lind) + &
               k(3) * xij(3,lind)
          
          ! The fact that we have a SYMMETRIC
          ! update region makes this *tricky* part easy...
          rin  = lup_ptr(jo)
          ! TODO, this REQUIRES that lup_col(:) is sorted
          rind = rin + SFIND(lup_col(rin+1:rin+lup_ncol(jo)),io)
          ! We do a check, just to be sure...
          if ( rind <= rin ) &
               call die('ERROR: Conjugated symmetrization point does not exist')
          
          ! The integration is this:
          ! \rho = e^{-i.k.R} [ \int (Gf^R-Gf^A) dE + \int Gf^R\Gamma Gf^A dE ]
          ! NOTE that weightDMC removes the daggered Gf^R\Gamma Gf^A
          ph = 0.5_dp * cdexp(dcmplx(0._dp,-kx))
          
          DM(lind)  = DM(lind)  + aimag( ph*(zD(ind) - conjg(zD(rind))) )
          
          EDM(lind) = EDM(lind) + aimag( ph*(zE(ind) - conjg(zE(rind))) )

       end do
    end do

  end subroutine update_zDM


  subroutine init_DM(dit,sp,maxn,DM,EDM, up_sp)
          ! The DM and EDM equivalent matrices
    use class_OrbitalDistribution
    use class_Sparsity
    use intrinsic_missing, only : UCORB => MODP
    use parallel, only : Node
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: maxn
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(maxn), EDM(maxn)
    ! The updated sparsity arrays...
    type(Sparsity), intent(inout) :: up_sp

    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
    integer :: lnr, lio, lind, io, jo, ind, nr

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    call attach(up_sp, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
     
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
     
    if ( nr /= nrows(up_sp) ) call die('The sparsity format is not as &
         &expected.')
     
    ! This loop is across the local rows...
    do lio = 1 , lnr

       ! obtain the global index of the local orbital.
       io = index_local_to_global(dit,lio,Node)

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( lup_ncol(io) == 0 ) cycle

       ! Now we loop across the update region
       ! This one must *per definition* have less elements.
       ! Hence, we can exploit this, and find equivalent
       ! super-cell orbitals.
       ! Ok, this is Gamma (but to be consistent)
       do ind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)

          jo = lup_col(ind)

          ! Do a loop in the local sparsity pattern...
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             ! We know that the update region is in 
             ! UC-format. Hence we can compare directly, via
             ! the orbital index in the unit-cell.
             if ( UCORB(l_col(lind),nr) == jo ) then
                DM (lind) = 0._dp
                EDM(lind) = 0._dp
             end if
             
          end do
       end do
    end do
    
  end subroutine init_DM

end module m_ts_dm_update
