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

  ! ***
  ! The following scheme should be followed:
  !   add_*_DM routines are constructed to be able to handle different schemes
  !   They are called after each k-point and thus the arguments are the following:
  !     1. the local update sparsity pattern
  !     2. the global update sparsity pattern (suffix 'u')
  ! The k-point routine is constructed to handle three different methods of doing
  ! the weighting which is performed here.

  subroutine add_k_DM(spDM,spuDM,D_dim2, spEDM, spuEDM, E_dim2, &
       n_s,sc_off,k, non_Eq)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use class_zSpData2D
    use geom_helper, only : UCORB
    use intrinsic_missing, only : SFIND
    use parallel, only : Node

! *********************
! * INPUT variables   *
! *********************
    ! The local integrated sparsity arrays
    type(dSpData2D), intent(inout) :: spDM, spEDM
    ! The current k-point global sparsity arrays
    type(zSpData2D), intent(inout) :: spuDM, spuEDM
    ! current update region of last dimension
    integer, intent(in) :: D_dim2, E_dim2
    ! The k-point
    real(dp), intent(in) :: k(3)
    ! The supercell offsets
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    logical, intent(in), optional :: non_Eq

    ! Arrays needed for looping the sparsity
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:)
    real(dp), pointer :: dD(:,:) , dE(:,:)
    complex(dp), pointer :: zDu(:,:), zEu(:,:)
    integer :: lnr, lio, lind, io, ind, nr, jo
    integer :: rin, rind
    logical :: hasEDM
    complex(dp) :: ph

    if ( (.not. initialized(spDM)) .or. (.not. initialized(spuDM)) ) return

    hasEDM = initialized(spEDM)

    ! get the distribution
    dit => dist(spDM)

    l_s  => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD   => val(spDM)
    if ( hasEDM ) dE => val(spEDM)

    up_s => spar(spuDM)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    zDu  => val(spuDM)
    if ( hasEDM ) zEu => val(spuEDM)

    if ( size(zDu,2) < D_dim2 .or. size(dD,2) < D_dim2 ) then
       call die('add_k_DM: Error in code')
    end if

    if ( hasEDM ) then
       if ( size(zEu,2) < E_dim2 .or. size(dE,2) < E_dim2 ) then
          call die('add_k_DM: Error in code')
       end if
    end if

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    
    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected.')

    if ( non_Eq ) then

! No data race will occur, sparsity pattern only tranversed once
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lind,jo,rind,ind,ph)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio,Node)
          if ( up_ncol(io) /= 0 ) then

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)
             
             rind = up_ptr(io)
             ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
             if ( ind <= rind ) cycle ! The element does not exist
             
             jo = (l_col(lind)-1) / nr
             ph = cdexp(dcmplx(0._dp, - &
                  k(1) * sc_off(1,jo) - &
                  k(2) * sc_off(2,jo) - &
                  k(3) * sc_off(3,jo)))
             
             ! The integration is this:
             ! \rho = e^{-i.k.R} \int Gf^R\Gamma Gf^A dE
             dD(lind,1:D_dim2) = dD(lind,1:D_dim2) + aimag( ph*zDu(ind,1:D_dim2) )
             if ( hasEDM ) dE(lind,1:D_dim2) = dE(lind,1:D_dim2) + &
                  aimag( ph*zEu(ind,1:D_dim2) )

          end do

          end if
          end if

       end do
!$OMP end parallel do

    else

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lind,jo,rin,rind,ind,ph)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio,Node)
          if ( up_ncol(io) /= 0 ) then

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             jo = UCORB(l_col(lind),nr)

             rind = up_ptr(io)
             ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
             if ( ind <= rind ) cycle ! The element does not exist
             
             rin  = up_ptr(jo)
             rind = rin + SFIND(up_col(rin+1:rin+up_ncol(jo)),io)
             if ( rind <= rin ) &
                  call die('ERROR: Conjugated symmetrization point does not exist')

             jo = (l_col(lind)-1) / nr
             ph = cdexp(dcmplx(0._dp, - &
                  k(1) * sc_off(1,jo) - &
                  k(2) * sc_off(2,jo) - &
                  k(3) * sc_off(3,jo)))

             ! This integration is this:
             ! \rho = e^{-i.k.R} \int (Gf^R-Gf^A) dE
             dD(lind,1:D_dim2) = dD(lind,1:D_dim2) + &
                  aimag( ph*(zDu(ind,1:D_dim2) - conjg(zDu(rind,1:D_dim2))) )
             if ( hasEDM ) dE(lind,1:D_dim2) = dE(lind,1:D_dim2) + &
                  aimag( ph*(zEu(ind,1:D_dim2) - conjg(zEu(rind,1:D_dim2))) )

          end do

          end if
          end if

       end do
!$OMP end parallel do

    end if

  end subroutine add_k_DM

  subroutine add_Gamma_DM(spDM,spuDM,D_dim2,spEDM,spuEDM,E_dim2)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use geom_helper, only : UCORB
    use intrinsic_missing, only : SFIND
    use parallel, only : Node

! *********************
! * INPUT variables   *
! *********************
    ! The local integrated sparsity arrays
    type(dSpData2D), intent(inout) :: spDM
    ! The current Gamma-point global sparsity arrays
    type(dSpData2D), intent(inout) :: spuDM
    integer, intent(in) :: D_dim2
    ! The local integrated sparsity arrays
    type(dSpData2D), intent(inout) :: spEDM
    ! The current Gamma-point global sparsity arrays
    type(dSpData2D), intent(inout) :: spuEDM
    integer, intent(in) :: E_dim2

    ! Arrays needed for looping the sparsity
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: l_s, up_s
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: up_ncol(:), up_ptr(:), up_col(:)
    real(dp), pointer :: dD(:,:), dE(:,:)
    real(dp), pointer :: dDu(:,:), dEu(:,:)
    integer :: lnr, lio, lind, io, ind, nr, jo, rind
    logical :: hasEDM

    if ( (.not. initialized(spDM)) .or. (.not. initialized(spuDM)) ) return

    hasEDM = initialized(spEDM)

    ! get distribution
    dit  => dist(spDM)

    l_s  => spar(spDM)
    call attach(l_s ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    dD   => val(spDM)
    if ( hasEDM ) dE  => val(spEDM)

    up_s => spar(spuDM)
    call attach(up_s,n_col=up_ncol,list_ptr=up_ptr,list_col=up_col)
    dDu  => val(spuDM)
    if ( hasEDM ) dEu => val(spuEDM)

    if ( size(dDu,2) < D_dim2 .or. size(dD,2) < D_dim2 ) then
       call die('add_Gamma_DM: Error in code')
    end if

    if ( hasEDM ) then
       if ( size(dEu,2) < E_dim2 .or. size(dE,2) < E_dim2 ) then
          call die('add_Gamma_DM: Error in code')
       end if
    end if

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    
    if ( nr /= nrows(up_s) ) call die('The sparsity format is not as &
         &expected.')

    if ( hasEDM ) then

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lind,jo,rind,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio,Node)
          if ( up_ncol(io) /= 0 ) then

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             ! we might still have SIESTA-non-Gamma
             jo = ucorb(l_col(lind),nr)

             ! This sparsity pattern is in UC
             rind = up_ptr(io)
             ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
             if ( ind <= rind ) cycle ! The element does not exist

             dD(lind,1:D_dim2) = dD(lind,1:D_dim2) + dDu(ind,1:D_dim2)
             dE(lind,1:E_dim2) = dE(lind,1:E_dim2) + dEu(ind,1:E_dim2)
             
          end do

          end if
          end if

       end do
!$OMP end parallel do
    else
! No data race condition will ever be encountered
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lind,jo,rind,ind)
       do lio = 1 , lnr

          if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio,Node)
          if ( up_ncol(io) /= 0 ) then

          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             jo = ucorb(l_col(lind),nr)

             rind = up_ptr(io)
             ind = rind + SFIND(up_col(rind+1:rind+up_ncol(io)),jo)
             if ( ind <= rind ) cycle ! The element does not exist

             dD(lind,1:D_dim2) = dD(lind,1:D_dim2) + dDu(ind,1:D_dim2)
             
          end do

          end if
          end if

       end do
!$OMP end parallel do
    end if

  end subroutine add_Gamma_DM

  subroutine update_DM(dit,sp,n_nzs,DM, spDM, Ef, &
       EDM, spEDM, ipnt, UpSpGlobal)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use class_iSpData1D
    use parallel, only : Node
    use geom_helper, only : UCORB
    use intrinsic_missing, only : SFIND

! *********************
! * INPUT variables   *
! *********************
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(dSpData2D), intent(inout) :: spDM
    ! fermi-level, we shift the energy density matrix back
    real(dp), intent(in) :: Ef
    ! Sparse energy-DM-arrays (local)
    real(dp), intent(inout) :: EDM(n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(dSpData2D), intent(inout) :: spEDM
    ! I.e. a pointer from the local update sparsity to the local sparsity
    type(iSpData1D), intent(in), optional :: ipnt
    ! Whether the update sparsity pattern is a global update sparsity pattern
    logical, intent(in), optional :: UpSpGlobal

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
    integer, pointer :: pnt(:)
    real(dp), pointer :: dD(:,:), dE(:,:)
    integer :: lnr, nr, uind, lio, io, lind, ind, ljo, jo
    logical :: hasipnt, hasEDM, lUpSpGlobal

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s  => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    dD => val(spDM)

    hasEDM = initialized(spEDM) 

    if ( hasEDM ) then

       dE => val(spEDM)

       ! the actual size of the shift
       jo = nnzs(spDM)

       ! As we have shifted the fermi-level up to 0, we need to shift the
       ! energy-density matrix back
       call daxpy(jo,Ef,dD(1,1),1,dE(1,1),1)

    end if

    ! We have that the update sparsity pattern is in local
    ! form.
    ! this means that sp == s (besides the non-update objects)
    ! Hence we don't need to utilize index_local_to_global

    
    hasipnt = present(ipnt)
    if ( hasipnt ) hasipnt = initialized(ipnt)

    lUpSpGlobal = .false.
    if ( present(UpSpGlobal) ) lUpSpGlobal = UpSpGlobal
    
    if ( .not. lUpSpGlobal ) then

    if ( lnr /= nrows(s) ) &
         call die('The sparsity format is not as expected.')

    if ( hasipnt ) then

       ! The pointer
       pnt => val(ipnt)

       ! This loop is across the local rows...
! No data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(io,uind,ind)
       do io = 1 , lnr

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Do a loop in the local update sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)
                
             ind = pnt(uind)
             
             DM(ind) = DM(ind) + dD(uind,1)
             if ( hasEDM ) EDM(ind) = EDM(ind) + dE(uind,1)
             
          end do

          end if
       end do
!$OMP end parallel do
       
    else 

       ! This loop is across the local rows...
! no data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(io,uind,jo,ind)
       do io = 1 , lnr

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Do a loop in the local update sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do uind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)
                
             ! We are dealing with a non-UC sparsity pattern
             jo = lup_col(uind)

             ! Now we loop across the local region
             ind = l_ptr(io) + &
                  minloc(abs(l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io))-jo),1)
             if ( l_col(ind) /= jo ) then
                do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                   if ( l_col(ind) /= jo ) cycle
                   exit ! we have the correct ind-value
                end do
                if ( l_col(ind) /= jo ) cycle
             end if
             
             ! We need to add in case of special weighting...
             DM(ind) = DM(ind) + dD(uind,1)
             if ( hasEDM ) EDM(ind) = EDM(ind) + dE(uind,1)
             
          end do
          end if
       end do
!$OMP end parallel do
    end if

    else
       ! We have a global update sparsity pattern

              ! This is the global sparsity pattern
       ! i.e. we require to call index_local_to_global
       ! The global sparsity pattern is not in supercell format

       ! This loop is across the local rows...
! We will never have a data race here (it is on local sparsity pattern)
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lind,ljo,ind)
       do lio = 1 , lnr

          ! obtain the global index of the local orbital.
          io = index_local_to_global(dit,lio,Node)

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

          ! Do a loop in the local sparsity pattern...
          ! The local sparsity pattern is more "spread", hence
          ! we do fewer operations by having this as an outer loop
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             ! we need to compare with the global update sparsity
             ljo = UCORB(l_col(lind),nr)

             ! Now we loop across the update region
             ! This one must *per definition* have less elements.
             ! Hence, we can exploit this, and find equivalent
             ! super-cell orbitals.
             ! Ok, this is Gamma (but to be consistent)
             ind = lup_ptr(io)
             ind = ind + SFIND(lup_col(ind+1:ind+lup_ncol(io)),ljo)
             if ( ind <= lup_ptr(io) ) cycle
             
             ! We only have one k-point, yet in case of non-Gamma siesta
             DM(lind)  = DM(lind) + dD(ind,1)
             if ( hasEDM ) EDM(lind) = EDM(lind) + dE(ind,1)
             
          end do
          
          end if
       end do
!$OMP end parallel do

    end if

  end subroutine update_DM

  ! This routine will ONLY be called if .not. IsVolt,
  ! Hence we don't have any sparsity patterns with local sparsity patterns
  ! that is dealing with this routine (hence we do need the index_local_to_global)
  subroutine update_zDM(dit,sp,n_nzs,DM,spDM, Ef, &
       EDM,spEDM, k, n_s, sc_off)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData2D

    use geom_helper, only : UCORB
    use intrinsic_missing, only : SFIND
    use parallel, only : Node
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs), EDM(n_nzs)
    ! Updated sparsity arrays (they contain the current integration)
    type(zSpData2D), intent(inout) :: spDM, spEDM
    ! The fermi level
    real(dp), intent(in) :: Ef
    ! The k-point...
    real(dp), intent(in) :: k(3)
    ! The supercell offset
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
    complex(dp), pointer :: zD(:,:), zE(:,:)
    complex(dp) :: ph
    integer :: lio, io, jo, ind, nr
    integer :: lnr, lind, rin, rind
    logical :: hasEDM

    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    s => spar(spDM)
    call attach(s, n_col=lup_ncol,list_ptr=lup_ptr,list_col=lup_col)
    zD     => val(spDM)

    hasEDM = initialized(spEDM)

    if ( hasEDM ) then

       zE => val(spEDM)

       ! the actual size of the shift
       jo = nnzs(spDM)

       ! As we have shifted the fermi-level up to 0, we need to shift the
       ! energy-density matrix back
       ph = dcmplx(Ef,0._dp)
       call zaxpy(jo,ph,zD(1,1),1,zE(1,1),1)

    end if

    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern (but still the global sparsity pattern)

    if ( nr /= nrows(s) ) call die('The sparsity format is not as &
         &expected.')
     
    ! This loop is across the local rows...
! No data race will occur
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,lind,jo,rind,ind,rin,ph)
    do lio = 1 , lnr

       ! obtain the global index of the local orbital.
       io = index_local_to_global(dit,lio,Node)

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( lup_ncol(io) /= 0 ) then

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
          jo = (l_col(lind)-1) / nr
          ph = cdexp(dcmplx(0._dp, - &
               k(1) * sc_off(1,jo) - &
               k(2) * sc_off(2,jo) - &
               k(3) * sc_off(3,jo)))
      
          DM(lind) = DM(lind) + aimag( ph*(zD(ind,1) - conjg(zD(rind,1))) )
          if ( hasEDM ) &
               EDM(lind) = EDM(lind) + aimag( ph*(zE(ind,1) - conjg(zE(rind,1))) )

       end do

       end if
    end do
!$OMP end parallel do

  end subroutine update_zDM


  subroutine init_DM(dit,sp,n_nzs,DM,EDM, up_sp, Calc_Forces)
    ! The DM and EDM equivalent matrices
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB
    use parallel, only : Node
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! Size of the sparsity arrays
    integer, intent(in) :: n_nzs
    ! Sparse DM-arrays (local)
    real(dp), intent(inout) :: DM(n_nzs), EDM(n_nzs)
    ! The updated sparsity arrays...
    type(Sparsity), intent(inout) :: up_sp
    logical, intent(in) :: Calc_Forces

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

    if ( Calc_Forces ) then
     
       ! This loop is across the local rows...
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,jo,ind,lind)
       do lio = 1 , lnr

          ! obtain the global index of the local orbital.
          io = index_local_to_global(dit,lio,Node)

          ! Quickly go past the empty regions... (we have nothing to update)
          if ( lup_ncol(io) /= 0 ) then

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
          end if
       end do
!$OMP end parallel do

    else
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,jo,ind,lind)
       do lio = 1 , lnr
          io = index_local_to_global(dit,lio,Node)
          if ( lup_ncol(io) /= 0 ) then
          do ind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)
             jo = lup_col(ind)
             do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
                if ( UCORB(l_col(lind),nr) == jo ) &
                     DM(lind) = 0._dp
             end do
          end do
          end if
       end do
!$OMP end parallel do
    end if
    
  end subroutine init_DM

end module m_ts_dm_update
