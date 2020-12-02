!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************

! Creating a Hamiltonian/Overlap matrix at a specified k-point
! in a global UC sparsity pattern is enabled using these routines:

!  - create_HS
!    1. accepts a distributed matrix
!    2. requires the output matrix to be globalized in the sense
!       that the sparsity pattern is a UC sparsity pattern
!       NO dublicate entries. 
!       REQUIREMENT: each row MUST be sorted in column index
!  - create_U
!    1. accepts a globalized sparsity pattern of a matrix
!    2. Creates a matrix in Upper Triangular form
!       which directly can be inserted in LAPACK [dz]spgvd/[dz]spev
!       routines.
!  - create_Full
!    1. accepts a globalized sparsity pattern of a matrix
!    2. creates the full matrix (with all the zeroes) which can
!       be directly used for BLAS/LAPACK matrix operations
!    3. Symmetrization is not enforced by this routine

!  - AllReduce_SpData
!    1. Accepts a matrix which will be reduced in the
!       MPI_Comm_World communicator.
!    2. It takes a work-array as an additional argument
!       (however we should consider the MPI_INPLACE)

module m_ts_sparse_helper

  use precision, only : dp
  use m_ts_method, only : orb_type, TYP_BUFFER, TYP_DEVICE

  implicit none

  private

#ifdef MPI
  interface AllReduce_SpData
     module procedure AllReduce_dSpData1D
     module procedure AllReduce_zSpData1D
     module procedure AllReduce_dSpData2D
     module procedure AllReduce_zSpData2D
  end interface 
  public :: AllReduce_SpData
#endif

  interface create_HS
     module procedure create_HS_Gamma
     module procedure create_HS_kpt
  end interface 
  public :: create_HS

  interface create_U
     module procedure create_Gamma_U
     module procedure create_kpt_U
  end interface 
  public :: create_U

  interface create_Full
     module procedure create_Gamma_Full
     module procedure create_kpt_Full
  end interface 
  public :: create_Full

contains

   ! Helper routine to create and distribute the sparse 
   ! k-point Hamiltonian.
  subroutine create_HS_kpt(dit,sp, &
       Ef, &
       N_Elec, Elecs, no_u, n_s, &
       n_nzs, H, S, sc_off, SpArrH, SpArrS, k, &
       nwork, work)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D

    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    use m_ts_electype

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Fermi-level
    real(dp), intent(in) :: Ef
    ! The electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: no_u, n_s
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs
    ! The hamiltonian and overlap sparse matrices 
    real(dp), intent(in) :: H(n_nzs),S(n_nzs)
    ! The supercell offsets
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
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
    integer, pointer :: k_ncol(:), k_ptr(:), k_col(:), kp_col(:)
    complex(dp), pointer :: zH(:), zS(:)
    complex(dp) :: ph(0:n_s-1)
    type(Sparsity), pointer :: sp_k
    integer :: no_l, lio, io, io_T, ind, jo, jo_T
    integer :: ind_k, kp
    real(dp) :: E_Ef(0:N_Elec)
    logical :: Bulk(0:N_Elec)

    ! obtain the local number of rows and the global...
    no_l = nrows(sp)
    if ( no_u /= nrows_g(sp) ) then
       call die('Creating the k-&point matrix in &
            &transiesta went wrong. Please TODO...')
    end if

    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! obtain the full sparsity unit-cell
    sp_k => spar(SpArrH)
    call attach(sp_k, n_col=k_ncol,list_ptr=k_ptr,list_col=k_col)

    ! Pre-calculate phases
    do jo = 0 , n_s - 1
      ph(jo) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,jo)),dp))
    end do

    ! create the overlap electrode fermi-level
    ! Note that for bulk V_frac_CT will be set to 0.
    E_Ef(0) = Ef
    Bulk(0) = .false. ! value doesn't matter, this is to look it up
    do jo = 1, N_elec
      E_Ef(jo) = Ef - Elecs(jo)%V_frac_CT * Elecs(jo)%mu%mu
      Bulk(jo) = Elecs(jo)%bulk
    end do

    ! obtain the value arrays...
    zH => val(SpArrH)
    zS => val(SpArrS)

    zH(:) = cmplx(0._dp,0._dp,dp)
    zS(:) = cmplx(0._dp,0._dp,dp)

! No data race condition as each processor takes a separate row
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,io_T,kp,kp_col,ind,jo,jo_T,ind_k)
    do lio = 1 , no_l

       ! obtain the global index of the orbital.
      io = index_local_to_global(dit,lio)
      io_T = orb_type(io)
      ! if there is no contribution in this row
      if ( k_ncol(io) /= 0 ) then
        kp = k_ptr(io)
        kp_col => k_col(kp+1:kp+k_ncol(io))

#ifndef TS_NOCHECKS
       ! The io orbitals are in the range [1;no_u_TS]
       ! This should be redundant as it is catched by k_ncol(io)==0
       if ( io_T == TYP_BUFFER ) then
         call die('Error in code, please contact Nick Papior: nickpapior@gmail.com')
       end if
#endif

       ! Loop number of entries in the row... (index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

         ! as the local sparsity pattern is a super-cell pattern,
         ! we need to check the unit-cell orbital
         ! The unit-cell column index
         jo = UCORB(l_col(ind),no_u)
         jo_T = orb_type(jo)

         ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
         if ( jo_T == TYP_BUFFER ) cycle

         ! Do a check whether we have connections
         ! across the junction...
         ! This is the same as removing all electrode connections
         if ( io_T > 0 .and. jo_T > 0 .and. io_T /= jo_T ) cycle
         if ( io_T /= jo_T ) then
           ! we definitely have Elec -> device
           ! Choose the electrode fermi-level
           jo_T = max(io_T, jo_T)
         else if ( io_T == jo_T .and. Bulk(jo_T) ) then
           ! no need to shift since we have a bulk H/S
           jo_T = 0
         end if
           
         ! find the equivalent position in the sparsity pattern
         ! of the full unit cell
         ! Notice that SFIND REQUIRES that the sparsity pattern
         ! is SORTED!
         ! Thus it will only work for UC sparsity patterns.
         ind_k = kp + SFIND(kp_col,jo)
         if ( kp < ind_k ) then
           jo = (l_col(ind)-1) / no_u

           zH(ind_k) = zH(ind_k) + (H(ind) - E_Ef(jo_T)*S(ind)) * ph(jo)
           zS(ind_k) = zS(ind_k) + S(ind) * ph(jo)
         end if

       end do

      end if

    end do
!$OMP end parallel do
     
#ifdef MPI
    if ( dist_nodes(dit) > 1 ) then
       ! Note that zH => val(SpArrH)
       ! Note that zS => val(SpArrS)
       call AllReduce_SpData(SpArrH,nwork,work)
       call AllReduce_SpData(SpArrS,nwork,work)
    end if
#endif

    ! It could be argued that MPI reduction provides
    ! numeric fluctuations.
    ! However, the block-cyclic distribution ensures that
    ! there are no two elements accessed by two or more processors.
    ! This makes all non-local elements ZERO, and there should not
    ! be flucuations on adding ZEROS as they are *only* dependent
    ! on order of summation.

  end subroutine create_HS_kpt


  ! Helper routine to create and distribute the sparse 
  ! k-point Hamiltonian.
  subroutine create_HS_Gamma(dit,sp, &
       Ef, &
       N_Elec, Elecs, no_u, &
       n_nzs, H, S, SpArrH, SpArrS, &
       nwork, work)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D

    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    use m_ts_electype

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Fermi-level
    real(dp), intent(in) :: Ef
    ! The electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: no_u
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs
    ! The hamiltonian and overlap sparse matrices 
    real(dp), intent(in) :: H(n_nzs),S(n_nzs)
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
    integer, pointer  :: k_ncol(:), k_ptr(:), k_col(:), kp_col(:)
    real(dp), pointer :: dH(:), dS(:)
    type(Sparsity), pointer :: sp_G
    integer :: no_l, lio, io, io_T, ind, jo, jo_T
    integer :: ind_k, kp
    real(dp) :: E_Ef(0:N_Elec)
    logical :: Bulk(0:N_Elec)
    
    ! obtain the local number of rows and the global...
    no_l = nrows(sp)
    if ( no_u /= nrows_g(sp) ) then
       call die('Creating the k-&point matrix in &
            &transiesta went wrong. Please TODO...')
    end if

    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! obtain the full sparsity unit-cell
    sp_G => spar(SpArrH)
    call attach(sp_G, n_col=k_ncol,list_ptr=k_ptr,list_col=k_col)

    ! create the overlap electrode fermi-level
    ! Note that for bulk V_frac_CT will be set to 0.
    E_Ef(0) = Ef
    Bulk(0) = .true. ! necessary to force using the fermi-level
    do jo = 1, N_elec
      E_Ef(jo) = Ef - Elecs(jo)%V_frac_CT * Elecs(jo)%mu%mu
      Bulk(jo) = Elecs(jo)%bulk
    end do

    ! obtain the value arrays...
    dH => val(SpArrH)
    dS => val(SpArrS)

    dH(:) = 0._dp
    dS(:) = 0._dp

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,io_T,kp,kp_col,ind,jo,jo_T,ind_k)
    do lio = 1 , no_l

      ! obtain the global index of the orbital.
      io = index_local_to_global(dit,lio)
      io_T = orb_type(io)
      ! if there is no contribution in this row
      if ( k_ncol(io) /= 0 ) then
        kp = k_ptr(io)
        kp_col => k_col(kp+1:kp+k_ncol(io))

#ifndef TS_NOCHECKS
       ! The io orbitals are in the range [1;no_u]
       if ( io_T == TYP_BUFFER ) then
         call die('Error in code, please contact Nick Papior: nickpapior@gmail.com')
       end if
#endif

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

         ! as the local sparsity pattern is a super-cell pattern,
         ! we need to check the unit-cell orbital
         ! The unit-cell column index
         jo = UCORB(l_col(ind),no_u)
         jo_T = orb_type(jo)

         ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
         if ( jo_T == TYP_BUFFER ) cycle

         ! Do a check whether we have connections
         ! across the junction...
         ! This is the same as removing all electrode connections
         if ( io_T > 0 .and. jo_T > 0 .and. io_T /= jo_T ) cycle
         if ( io_T /= jo_T ) then
           ! we definitely have Elec -> device
           ! Choose the electrode fermi-level
           jo_T = max(io_T, jo_T)
         else if ( io_T == jo_T .and. Bulk(jo_T) ) then
           ! no need to shift since we have a bulk H/S
           jo_T = 0
         end if
         
         ! find the equivalent position in the sparsity pattern
         ! of the full unit cell
         ind_k = kp + SFIND(kp_col,jo)
         if ( kp < ind_k ) then
           dH(ind_k) = dH(ind_k) + H(ind) - E_Ef(jo_T) * S(ind)
           dS(ind_k) = dS(ind_k) + S(ind)
         end if
         
       end do

      end if

    end do
!$OMP end parallel do
     
#ifdef MPI
    ! Note that dH => val(SpArrH)
    ! Note that dS => val(SpArrS)
    call AllReduce_SpData(SpArrH,nwork,work)
    call AllReduce_SpData(SpArrS,nwork,work)
#endif

    ! It could be argued that MPI reduction provides
    ! numeric fluctuations.
    ! However, the block-cyclic distribution ensures that
    ! there are no two elements accessed by two or more processors.
    ! This makes all non-local elements ZERO, and there should not
    ! be flucuations on adding ZEROS as they are *only* dependent
    ! on order of summation.

  end subroutine create_HS_Gamma


  ! Helper routine to create and distribute an upper
  ! tri-angular matrix of the hamiltonian in a specified
  ! region.
  subroutine create_Gamma_U(dit,sp, &
       no, r, &
       n_nzs, A, A_UT)

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Number of orbitals that form this array
    integer, intent(in) :: no
    ! The region that we wish to create the UT matrix in
    type(tRgn), intent(in) :: r
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs
    real(dp), intent(in) :: A(n_nzs)
    ! The UT format matrix
    real(dp), intent(out) :: A_UT(no*(no+1)/2)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind, j, i, idx
    
    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    A_UT(:) = 0._dp

    ! Loop over region orbitals
!$OMP parallel do default(shared), private(i,io,lio,ind,j,idx)
    do i = 1 , no
       
       ! Global orbital
       io = r%r(i)
       ! obtain the global index of the orbital.
       lio = index_global_to_local(dit,io)
       if ( lio > 0 ) then
       
       ! if there is no contribution in this row
       if ( l_ncol(lio) > 0 ) then

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          idx = UCORB(l_col(ind),no_u)

          ! If the orbital is not in the region, we skip it
          j = rgn_pivot(r,idx)
          if ( j <= 0 ) cycle

          ! Calculate position
          ! For Gamma, we do not need the complex conjugate...
          if ( i > j ) then
             idx = j + (i -1)*i /2
!$OMP atomic
             A_UT(idx) = A_UT(idx) + 0.5_dp * A(ind)
          else if ( i < j ) then
             idx =  i + (j-1)*j/2
!$OMP atomic
             A_UT(idx) = A_UT(idx) + 0.5_dp * A(ind)
          else
             idx =  i + (j-1)*j/2
!$OMP atomic
             A_UT(idx) = A_UT(idx) +          A(ind)
          end if

       end do

       end if
       end if

    end do
!$OMP end parallel do
     
  end subroutine create_Gamma_U

  ! Helper routine to create and distribute an upper
  ! tri-angular matrix of the hamiltonian in a specified
  ! region.
  subroutine create_Gamma_Full(dit,sp, &
       no, r, &
       n_nzs, A, A_full)

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Number of orbitals that form this array
    integer, intent(in) :: no
    ! The region that we wish to create the UT matrix in
    type(tRgn), intent(in) :: r
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs
    real(dp), intent(in) :: A(n_nzs)
    ! The matrix
    real(dp), intent(out) :: A_full(no,no)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind, jo, i
    
    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    A_full(:,:) = 0._dp

    ! Loop over region orbitals
!$OMP parallel do default(shared), private(i,io,lio,ind,jo)
    do i = 1 , no
       
       ! Global orbital
       io = r%r(i)
       ! obtain the global index of the orbital.
       lio = index_global_to_local(dit,io)
       if ( lio > 0 ) then
       
       ! if there is no contribution in this row
       if ( l_ncol(lio) > 0 ) then

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)

          ! If the orbital is not in the region, we skip it
          jo = rgn_pivot(r,jo)
          if ( jo <= 0 ) cycle

          ! Calculate position
          ! For Gamma, we do not need the complex conjugate...
          A_full(jo,io) = A_full(jo,io) + A(ind)

       end do

       end if
       end if

    end do
!$OMP end parallel do
     
  end subroutine create_Gamma_Full

  ! Helper routine to create and distribute an upper
  ! tri-angular matrix of the hamiltonian in a specified
  ! region.
  subroutine create_kpt_U(dit,sp, &
       no, r, &
       n_nzs, n_s, A, sc_off, A_UT, k)

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Number of orbitals that form this array
    integer, intent(in) :: no
    ! The region that we wish to create the UT matrix in
    type(tRgn), intent(in) :: r
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs, n_s
    real(dp), intent(in) :: A(n_nzs), sc_off(3,0:n_s-1)
    ! The k-point we will create
    real(dp), intent(in) :: k(3)
    ! The UT format matrix
    complex(dp), intent(out) :: A_UT(no*(no+1)/2)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind, j, i, idx, is
    complex(dp) :: ph(0:n_s-1), p
    
    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
        nrows=no_l,nrows_g=no_u)

    ! Pre-calculate phases
    ! Note that we don't have a -, this is because we are populating correctly
    do is = 0 , n_s - 1
      ph(is) = exp(cmplx(0._dp, dot_product(k, sc_off(:,is)),dp)) * 0.5_dp
    end do

    A_UT(:) = cmplx(0._dp,0._dp,dp)

!$OMP parallel do default(shared), &
!$OMP&private(i,io,lio,ind,j,is,idx,p)
    do i = 1 , no
       
       ! Global orbital
       io = r%r(i)
       ! obtain the global index of the orbital.
       lio = index_global_to_local(dit,io)
       if ( lio > 0 ) then
       
       ! if there is no contribution in this row
       if ( l_ncol(lio) > 0 ) then

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          is = UCORB(l_col(ind),no_u)

          ! If the orbital is not in the region, we skip it
          j = rgn_pivot(r,is)
          if ( j <= 0 ) cycle

          is = (l_col(ind)-1)/no_u

          ! Calculate position
          ! Excerpt from dspgvd.f
          !          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
          !          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
          if ( i < j ) then ! U
            idx = i + (j-1)*j/2
            p = ph(is)
          else if ( j < i ) then ! L
            idx = j + (i-1)*i/2
            p = conjg(ph(is))
          else ! i == j same as U, but we don't double count these
            idx = i + (j-1)*j/2
            p = ph(is) * 2
          end if

!$OMP atomic
          A_UT(idx) = A_UT(idx) + p * A(ind)

       end do

       end if
       end if

    end do
!$OMP end parallel do
     
  end subroutine create_kpt_U

  ! Helper routine to create and distribute an upper
  ! tri-angular matrix of the hamiltonian in a specified
  ! region.
  subroutine create_kpt_full(dit,sp, &
       no, r, &
       n_nzs, n_s, A, sc_off, A_full, k)

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Number of orbitals that form this array
    integer, intent(in) :: no
    ! The region that we wish to create the UT matrix in
    type(tRgn), intent(in) :: r
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs, n_s
    real(dp), intent(in) :: A(n_nzs), sc_off(3,0:n_s-1)
    ! The k-point we will create
    real(dp), intent(in) :: k(3)
    ! The UT format matrix
    complex(dp), intent(out) :: A_full(no,no)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind, jo, i, is
    complex(dp) :: ph(0:n_s-1)
    
    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    ! Pre-calculate phases
    do is = 0 , n_s - 1
      ph(is) = exp(cmplx(0._dp, -dot_product(k, sc_off(:,is)), dp))
    end do

    A_full(:,:) = cmplx(0._dp,0._dp,dp)

    ! Loop over region orbitals
!$OMP parallel do default(shared), private(i,io,lio,ind,jo,is)
    do i = 1 , no
       
       ! Global orbital
       io = r%r(i)
       ! obtain the global index of the orbital.
       lio = index_global_to_local(dit,io)
       if ( lio > 0 ) then
       
       ! if there is no contribution in this row
       if ( l_ncol(lio) > 0 ) then

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)

          ! If the orbital is not in the region, we skip it
          jo = rgn_pivot(r,jo)
          if ( jo <= 0 ) cycle

          is = (l_col(ind)-1)/no_u

          A_full(jo,io) = A_full(jo,io) + ph(is) * A(ind)

       end do

       end if
       end if

    end do
!$OMP end parallel do
     
  end subroutine create_kpt_full


! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************

#ifdef MPI

  ! **** Double precision complex ****
  subroutine AllReduce_z1D(nnzs,arr,nwork,work)
    use mpi_siesta
    integer, intent(in) :: nnzs
    complex(dp), intent(inout) :: arr(nnzs)
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    integer :: MPIerror, i
    i = 0
    do while ( i + nwork <= nnzs )
      call zcopy(nwork,arr(i+1),1,work(1),1)
      call MPI_AllReduce(work(1),arr(i+1),nwork, &
          MPI_Double_Complex, MPI_Sum, MPI_Comm_World, MPIerror)
      i = i + nwork
    end do
    if ( i < nnzs ) then
      call zcopy(nnzs-i,arr(i+1),1,work(1),1)
      call MPI_AllReduce(work(1),arr(i+1),nnzs-i, &
          MPI_Double_Complex, MPI_Sum, MPI_Comm_World, MPIerror)
    end if
  end subroutine AllReduce_z1D

  subroutine AllReduce_zSpData1D(sp_arr,nwork,work)
    use mpi_siesta
    use class_zSpData1D
    type(zSpData1D), intent(inout) :: sp_arr
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    complex(dp), pointer :: arr(:)
    integer :: n_nzs
    n_nzs = nnzs(sp_arr)
    arr => val(sp_arr)
    call AllReduce_z1D(n_nzs,arr,nwork,work)
  end subroutine AllReduce_zSpData1D

  subroutine AllReduce_zSpData2D(sp_arr,nwork,work,dim2_count)
    use mpi_siesta
    use class_zSpData2D
    type(zSpData2D), intent(inout) :: sp_arr
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    integer, intent(in), optional :: dim2_count
    complex(dp), pointer :: arr(:,:)
    integer :: n_nzs, d2
    arr => val(sp_arr)
    d2 = size(arr,dim=2)
    if ( present(dim2_count) ) d2 = dim2_count
    n_nzs = nnzs(sp_arr) * d2
    call AllReduce_z1D(n_nzs,arr,nwork,work)
  end subroutine AllReduce_zSpData2D



  ! **** Double precision ****
  subroutine AllReduce_d1D(nnzs,arr,nwork,work)
    use mpi_siesta
    integer, intent(in) :: nnzs
    real(dp), intent(inout) :: arr(nnzs)
    integer, intent(in) :: nwork
    real(dp), intent(inout) :: work(nwork)
    integer :: MPIerror, i
    i = 0
    do while ( i + nwork <= nnzs )
      call dcopy(nwork,arr(i+1),1,work(1),1)
      call MPI_AllReduce(work(1),arr(i+1),nwork, &
          MPI_Double_Precision, MPI_Sum, MPI_Comm_World, MPIerror)
      i = i + nwork
    end do
    if ( i < nnzs ) then
      call dcopy(nnzs-i,arr(i+1),1,work(1),1)
      call MPI_AllReduce(work(1),arr(i+1),nnzs-i, &
          MPI_Double_Precision, MPI_Sum, MPI_Comm_World, MPIerror)
    end if
  end subroutine AllReduce_d1D
  
  subroutine AllReduce_dSpData1D(sp_arr,nwork,work)
    use mpi_siesta
    use class_dSpData1D
    type(dSpData1D), intent(inout) :: sp_arr
    integer, intent(in)     :: nwork
    real(dp), intent(inout) :: work(nwork)
    real(dp), pointer :: arr(:)
    integer :: n_nzs
    n_nzs = nnzs(sp_arr)
    arr => val(sp_arr)
    call AllReduce_d1D(n_nzs,arr,nwork,work)
  end subroutine AllReduce_dSpData1D

  subroutine AllReduce_dSpData2D(sp_arr,nwork,work,dim2_count)
    use mpi_siesta
    use class_dSpData2D
    type(dSpData2D), intent(inout) :: sp_arr
    integer, intent(in) :: nwork
    real(dp), intent(inout) :: work(nwork)
    integer, intent(in), optional :: dim2_count
    real(dp), pointer :: arr(:,:)
    integer :: n_nzs, d2
    arr => val(sp_arr)
    d2 = size(arr,dim=2)
    if ( present(dim2_count) ) d2 = dim2_count
    n_nzs = nnzs(sp_arr) * d2
    call AllReduce_d1D(n_nzs,arr,nwork,work)
  end subroutine AllReduce_dSpData2D

#endif

end module m_ts_sparse_helper

#ifdef MPI
subroutine my_full_G_reduce(sp_arr,nwork,work,dim2_count)
  use precision, only : dp
  use class_dSpData2D
  use m_ts_sparse_helper, only : AllReduce_SpData
  type(dSpData2D), intent(inout) :: sp_arr
  integer, intent(in)     :: nwork
  real(dp), intent(inout) :: work(nwork)
  integer, intent(in) :: dim2_count
  call AllReduce_SpData(sp_arr,nwork,work,dim2_count)
end subroutine my_full_G_reduce
#endif

