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

module m_ts_sparse_helper

  use precision, only : dp

  implicit none
  
  private :: dp
  
contains

   ! Helper routine to create and distribute the sparse 
   ! k-point Hamiltonian.
  subroutine create_HS_kpt(dit,sp, &
       Ef, &
       no_BufL, no_BufR, &
       no_C_L, no_C_R, no_u, &
       maxnh, H, S, xij, SpArrH, SpArrS, k, &
       nwork, work)

    use parallel, only : Node
    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D

    use intrinsic_missing, only : SFIND
    use geom_helper, only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S and xij lives by
    type(Sparsity), intent(inout) :: sp
    ! Fermi-level
    real(dp), intent(in) :: Ef
    ! The number of orbitals we wish to cut-off at both ends
    integer, intent(in) :: no_BufL, no_BufR
    integer, intent(in) :: no_C_L, no_C_R, no_u
    ! The number of elements in the sparse arrays
    integer, intent(in) :: maxnh
    ! The hamiltonian and overlap sparse matrices 
    real(dp), intent(in) :: H(maxnh),S(maxnh)
    ! The orbital distance array
    real(dp), intent(in) :: xij(3,maxnh)
    ! The arrays we will save in...
    type(zSpData1D), intent(inout) :: SpArrH, SpArrS
    ! The k-point we will create
    real(dp), intent(in) :: k(3)
    ! we pass a work array
    integer, intent(in) :: nwork
    ! work-array
    complex(dp), intent(in out) :: work(nwork)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: k_ncol(:), k_ptr(:), k_col(:)
    complex(dp), pointer :: zH(:), zS(:)
    complex(dp) :: ph
    type(Sparsity), pointer :: sp_k
    integer :: no_l, lio, io, ind, jo, ind_k, kn
    integer :: no_max
     
    ! obtain the local number of rows and the global...
    no_l = nrows(sp)
    if ( no_u /= nrows_g(sp) ) then
       call die('Creating the k-&point matrix in &
            &transiesta went wrong. Please TODO...')
    end if

    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! obtain the full sparsity unit-cell
    sp_k   => spar    (SpArrH)
    call attach(sp_k, n_col=k_ncol,list_ptr=k_ptr,list_col=k_col)

    ! The boundary at the right buffer
    no_max = no_u - no_BufR
     
    call init_val(SpArrH)
    call init_val(SpArrS)
    ! obtain the value arrays...
    zH => val(SpArrH)
    zS => val(SpArrS)

    do lio = 1 , no_l

       ! obtain the global index of the orbital.
       io = index_local_to_global(dit,lio,Node)
       kn = k_ncol(io)
       ! if there is no contribution in this row
       if ( kn == 0 ) cycle

       ! The io orbitals are in the range [1;no_u_TS]
       ! This should be redundant as it is catched by kn==0
       if ( io <= no_BufL .or. no_max < io ) cycle

       ! Loop number of entries in the row... (index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)

          ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
          if ( jo <= no_BufL .or. no_max < jo ) cycle

          ! Do a check whether we have connections
          ! across the junction...
          ! This is the same as removing LEFT-RIGHT states..
          if ( io < no_C_L .and. no_C_R < jo ) cycle
          if ( jo < no_C_L .and. no_C_R < io ) cycle
           
          ! find the equivalent position in the sparsity pattern
          ! of the full unit cell
          ind_k = k_ptr(io)

          ! Notice that SFIND REQUIRES that the sparsity pattern
          ! is SORTED!
          ! Thus it will only work for UC sparsity patterns.
          ind_k = ind_k + SFIND(k_col(ind_k+1:ind_k+kn),jo)
! TODO, do we really need to check here?
!          if ( ind_k <= k_ptr(io) ) &
!               call die('Could not find k-point index')
          if ( ind_k <= k_ptr(io) ) cycle

          ph = cdexp(dcmplx(0._dp, &
               k(1) * xij(1,ind) + &
               k(2) * xij(2,ind) + &
               k(3) * xij(3,ind)) )
          
          zH(ind_k) = zH(ind_k) + H(ind) * ph
          
          zS(ind_k) = zS(ind_k) + S(ind) * ph

       end do

    end do
     
#ifdef MPI
    ! Note that zH => val(SpArrH)
    ! Note that zS => val(SpArrS)
    call AllReduce_zSpData1D(SpArrH,nwork,work)
    call AllReduce_zSpData1D(SpArrS,nwork,work)
#endif

    ! We symmetrize AND shift
    call symmetrize_HS_kpt(Ef,SpArrH,SpArrS)
     
    ! It could be argued that MPI reduction provides
    ! numeric fluctuations.
    ! However, the block-cyclic distribution ensures that
    ! there are no two elements accessed by two or more processors.
    ! This makes all non-local elements ZERO, and there should not
    ! be flucuations on adding ZEROS as they are *only* dependent
    ! on order of summation.

  end subroutine create_HS_kpt

  subroutine symmetrize_HS_kpt(Ef,SpArrH, SpArrS)
    use parallel, only : Node
    use class_Sparsity
    use class_zSpData1D
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: Ef
    ! The arrays we will save in... these are the entire TS-region sparsity
    type(zSpData1D), intent(inout) :: SpArrH, SpArrS

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    type(Sparsity), pointer :: s
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: zH(:), zS(:)
    integer :: nr, io, ind, jo, rin, rind

    s  => spar(SpArrH)
    call attach(s, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows_g=nr)
    zH => val(SpArrH)
    zS => val(SpArrS)

    ! This loop is across the local rows...
    do io = 1 , nr

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( l_ncol(io) == 0 ) cycle

       ! Now we loop across the update region
       ! This one must *per definition* have less elements.
       ! Hence, we can exploit this, and find equivalent
       ! super-cell orbitals.
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          jo = l_col(ind)

          ! As we symmetrize we do not need
          ! to cycle all points through two times...
          if ( jo < io ) cycle

          ! We will find the Hermitian part:
          ! The fact that we have a SYMMETRIC
          ! update region makes this *tricky* part easy...
          rin  = l_ptr(jo)
          ! TODO, this REQUIRES that l_col(:) is sorted
          rind = rin + SFIND(l_col(rin+1:rin+l_ncol(jo)),io)
          ! We do a check, just to be sure...
          if ( rind <= rin ) &
               call die('ERROR symmetrization orbital does not &
               &exist.')

          ! Symmetrize (notice that this is *transposed*)
          ! See prep_GF
          zS(rind)  = 0.5_dp * ( zS(ind) + dconjg(zS(rind)) )
          zH(rind)  = 0.5_dp * ( zH(ind) + dconjg(zH(rind)) ) &
               - Ef * zS(rind)
          
          zS(ind) = dconjg(zS(rind))
          zH(ind) = dconjg(zH(rind))

          if ( ind == rind ) then
             ! This is the diagonal matrix elements
             zS(ind) = dreal(zS(ind))
             zH(ind) = dreal(zH(ind))
          end if
                      
       end do
    end do

  end subroutine symmetrize_HS_kpt

  ! Helper routine to create and distribute the sparse 
  ! k-point Hamiltonian.
  subroutine create_HS_Gamma(dit,sp, &
       Ef, &
       no_BufL, no_BufR, &
       no_C_L, no_C_R, no_u, &
       maxnh, H, S, SpArrH, SpArrS, &
       nwork, work)

    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use parallel, only : Node
! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S and xij lives by
    type(Sparsity), intent(inout) :: sp
    ! Fermi-level
    real(dp), intent(in) :: Ef
    ! The number of orbitals we wish to cut-off at both ends
    integer, intent(in) :: no_BufL, no_BufR
    integer, intent(in) :: no_C_L, no_C_R, no_u
    ! The number of elements in the sparse arrays
    integer, intent(in) :: maxnh
    ! The hamiltonian and overlap sparse matrices 
    real(dp), intent(in) :: H(maxnh),S(maxnh)
    ! The arrays we will save in... these are the entire TS-region sparsity
    type(dSpData1D), intent(inout) :: SpArrH, SpArrS
    ! we pass a work array
    integer, intent(in) :: nwork
    ! work-array
    real(dp), intent(in out) :: work(nwork)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer  :: k_ncol(:), k_ptr(:), k_col(:)
    real(dp), pointer :: dH(:), dS(:)
    type(Sparsity), pointer :: sp_G
    integer :: no_l, lio, io, ind, jo, ind_k
    
    ! obtain the local number of rows and the global...
    no_l = nrows(sp)
    if ( no_u /= nrows_g(sp) ) then
       call die('Creating the k-&point matrix in &
            &transiesta went wrong. Please TODO...')
    end if

    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! obtain the full sparsity unit-cell
    sp_G   => spar(SpArrH)
    call attach(sp_G, n_col=k_ncol,list_ptr=k_ptr,list_col=k_col)
    
    ! initialize to 0
    call init_val(SpArrH)
    call init_val(SpArrS)
    ! obtain the value arrays...
    dH => val(SpArrH)
    dS => val(SpArrS)

    do lio = 1 , no_l

       ! obtain the global index of the orbital.
       io = index_local_to_global(dit,lio,Node)
       ! if there is no contribution in this row
       if ( k_ncol(io) == 0 ) cycle

       ! The io orbitals are in the range [1;no_u]
       if ( io <= no_BufL )       cycle
       if ( no_u - no_BufR < io ) cycle

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)

          ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
          if ( jo <= no_BufL .or. no_u - no_BufR < jo ) cycle

          ! Do a check whether we have connections
          ! across the junction...
          ! This is the same as removing LEFT-RIGHT states..
          if ( io < no_C_L .and. no_C_R < jo ) cycle
          if ( jo < no_C_L .and. no_C_R < io ) cycle
           
          ! find the equivalent position in the sparsity pattern
          ! of the full unit cell
          ind_k = k_ptr(io)
          ind_k = ind_k + SFIND(k_col(ind_k+1:ind_k+k_ncol(io)),jo)
! TODO, do we really need to check here
!          if ( ind_k <= k_ptr(io) ) &
!               call die('Could not find k-point index')
          if ( ind_k <= k_ptr(io) ) cycle
     
          ! Todo, as this is a Gamma-calculation
          ! we probably should NOT do 'dH = dH + H'
          ! rather 'dH = H'

          dH(ind_k) = dH(ind_k) + H(ind)
          dS(ind_k) = dS(ind_k) + S(ind)

       end do

    end do
     

#ifdef MPI
    ! Note that dH => val(SpArrH)
    ! Note that dS => val(SpArrS)
    call AllReduce_dSpData1D(SpArrH,nwork,work)
    call AllReduce_dSpData1D(SpArrS,nwork,work)
#endif

    ! We need to do symmetrization AFTER reduction as we need the full
    ! Hamiltonian before we can do anything
    call symmetrize_HS_Gamma(Ef,SpArrH,SpArrS)

    ! It could be argued that MPI reduction provides
    ! numeric fluctuations.
    ! However, the block-cyclic distribution ensures that
    ! there are no two elements accessed by two or more processors.
    ! This makes all non-local elements ZERO, and there should not
    ! be flucuations on adding ZEROS as they are *only* dependent
    ! on order of summation.

  end subroutine create_HS_Gamma

  subroutine symmetrize_HS_Gamma(Ef,SpArrH, SpArrS)
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    use class_Sparsity
    use class_dSpData1D
    use parallel, only : Node
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: Ef
    ! The arrays we will save in... these are the entire TS-region sparsity
    type(dSpData1D), intent(inout) :: SpArrH, SpArrS

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    type(Sparsity), pointer :: s
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: dH(:), dS(:)
    integer :: nr, io, ind, jo, rin, rind
    
    s    => spar(SpArrH)
    call attach(s, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows_g=nr)
    dH     => val(SpArrH)
    dS     => val(SpArrS)

    ! This loop is across the local rows...
    do io = 1 , nr

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( l_ncol(io) == 0 ) cycle

       ! Now we loop across the update region
       ! This one must *per definition* have less elements.
       ! Hence, we can exploit this, and find equivalent
       ! super-cell orbitals.
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          jo = l_col(ind)

          ! As we symmetrize we do not need
          ! to cycle all points through two times...
          if ( jo < io ) cycle

          ! We will find the Hermitian part:
          ! The fact that we have a SYMMETRIC
          ! update region makes this *tricky* part easy...
          rin  = l_ptr(jo)
          ! TODO, this REQUIRES that l_col(:) is sorted
          rind = rin + SFIND(l_col(rin+1:rin+l_ncol(jo)),io)
          ! We do a check, just to be sure...
          if ( rind <= rin ) &
               call die('ERROR symmetrization orbital does not &
               &exist.')

          ! Symmetrize
          dS(ind)  = 0.5_dp * ( dS(ind) + dS(rind) )
          dH(ind)  = 0.5_dp * ( dH(ind) + dH(rind) ) &
               - Ef * dS(ind)

          ! we have a real Matrix (so imaginary part is zero)
          dH(rind) = dH(ind)
          dS(rind) = dS(ind)
                      
       end do
    end do
    
  end subroutine symmetrize_HS_Gamma



! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************

#ifdef MPI
  subroutine AllReduce_zSpData1D(sp_arr,nwork,work)
    use mpi_siesta
    use class_zSpData1D
    type(zSpData1D), intent(inout) :: sp_arr
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    complex(dp), pointer :: arr(:)
    integer :: MPIerror, n_nzs
    n_nzs = nnzs(sp_arr)
    ! This should never happen, exactly due to the sparsity
    if ( n_nzs > nwork ) call die('Sparsity seems larger than &
         &work arrays, Transiesta????')
    arr => val(sp_arr)
    work(1:n_nzs) = arr(1:n_nzs)
    call MPI_AllReduce(work(1),arr(1),n_nzs, &
         MPI_Double_Complex, MPI_Sum, MPI_Comm_World, MPIerror)
  end subroutine AllReduce_zSpData1D

  subroutine AllReduce_dSpData1D(sp_arr,nwork,work)
    use mpi_siesta
    use class_dSpData1D
    type(dSpData1D), intent(inout) :: sp_arr
    integer, intent(in)     :: nwork
    real(dp), intent(inout) :: work(nwork)
    real(dp), pointer :: arr(:)
    integer :: MPIerror, n_nzs
    n_nzs = nnzs(sp_arr)
    ! This should never happen, exactly due to the sparsity
    if ( n_nzs > nwork ) call die('Sparsity seems larger than &
         &work arrays, Transiesta????')
    arr => val(sp_arr)
    work(1:n_nzs) = arr(1:n_nzs)
    call MPI_AllReduce(work(1),arr(1),n_nzs, &
         MPI_Double_Precision, MPI_Sum, MPI_Comm_World, MPIerror)
  end subroutine AllReduce_dSpData1D
#endif

  subroutine d_DM_EDM_Reduce_Shift(Ef,spDM, spEDM, nwork, work)
    use class_dSpData1D
    real(dp), intent(in) :: Ef
    type(dSpData1D), intent(inout) :: spDM, spEDM
    integer, intent(in) :: nwork
    real(dp), intent(inout) :: work(nwork)
    integer :: n
    real(dp), pointer :: DM(:), EDM(:)

#ifdef MPI
    call AllReduce_dSpData1D(spDM,nwork,work)
    call AllReduce_dSpData1D(spEDM,nwork,work)
#endif
    n   = nnzs(spDM)
    DM  => val(spDM)
    EDM => val(spEDM)
    call daxpy(n,Ef,DM,1,EDM,1)

  end subroutine d_DM_EDM_Reduce_Shift

  subroutine z_DM_EDM_Reduce_Shift(Ef,spDM, spEDM, nwork, work)
    use class_zSpData1D
    real(dp), intent(in) :: Ef
    type(zSpData1D), intent(inout) :: spDM, spEDM
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    integer :: n
    complex(dp) :: zEf
    complex(dp), pointer :: DM(:), EDM(:)

#ifdef MPI
    call AllReduce_zSpData1D(spDM,nwork,work)
    call AllReduce_zSpData1D(spEDM,nwork,work)
#endif
    zEf = cmplx(Ef,0._dp,dp)
    n   = nnzs(spDM)
    DM  => val(spDM)
    EDM => val(spEDM)
    call zaxpy(n,zEf,DM,1,EDM,1)

  end subroutine z_DM_EDM_Reduce_Shift
 
end module m_ts_sparse_helper
