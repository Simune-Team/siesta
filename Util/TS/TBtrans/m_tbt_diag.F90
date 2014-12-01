
module m_tbt_diag

  use precision, only : dp
  use class_Sparsity
  use class_OrbitalDistribution
  use class_dSpData1D
  use class_dSpData2D

  use m_region

  implicit none

  private

  interface calc_sqrt_S
     module procedure calc_sqrt_S_Gamma
     module procedure calc_sqrt_S_kpt
  end interface calc_sqrt_S

  interface calc_Eig
     module procedure calc_Eig_Gamma
     module procedure calc_Eig_kpt
  end interface calc_Eig

  interface norm_Eigenstate
     module procedure norm_Eigenstate_Gamma
     module procedure norm_eigenstate_kpt
  end interface norm_Eigenstate

  public :: calc_sqrt_S
  public :: calc_Eig
  public :: norm_Eigenstate

contains

  subroutine calc_sqrt_S_Gamma(spS,orb,S_sq)
    use m_ts_sparse_helper, only : create_U
    type(dSpData1D), intent(inout) :: spS
    type(tRegion), intent(in) :: orb
    real(dp), intent(out) :: S_sq(orb%n,orb%n)

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: no, n_nzs, i, j, info
    real(dp), pointer :: S(:)
    real(dp), allocatable :: eig(:), v(:,:), S_UT(:), work(:)

    no = orb%n
    dit => dist(spS)
    sp => spar(spS)
    S => val(spS)
    n_nzs = size(S)

    allocate(eig(no))
    allocate(v(no,no))
    allocate(work(3*no)) ! real-work
    allocate(S_UT(no*no))

    call create_U(dit, sp, no, orb, n_nzs,S, S_UT)

    ! Diagonalize overlap matrix
    call dspev('V','U',no,S_UT,eig,v,no,work,info)
    if ( info /= 0 ) then
       write(*,*) 'INFO = ',info,' when diagonalizing Gamma overlap matrix'
    end if

    ! There are faster ways of doing this, but let's play safe for
    ! now. The overlap matrix is a positive definite matrix,
    ! hence all eigenvalues are positive.
    ! This "check" ensures that this is enforced.
    if ( any(eig < 0._dp) ) then
       write(*,'(3a,e12.5)')'tbtrans: Projection ',trim(orb%name), &
            ' is not completely positive definite, lowest eig of S: ', &
            minval(eig)
    end if
    where ( eig < 0.0_dp )
       eig = 0._dp
    elsewhere
       eig = dsqrt(eig)
    end where

    ! Calculate S^(1/2)
    ! Calculate v.sqrt(eig)
    j = 1
    do i = 1 , no
       S_UT(j:j+no-1) = v(:,i) * eig(i)
       j = j + no
    end do
    
    ! Calculate v.sqrt(eig).v^\dagger = S^(1/2)
    call dgemm('N','T',no,no,no,1._dp, &
         S_UT,no, v,no, &
         0._dp,S_sq,no)
    
    deallocate(eig,v,work,S_UT)

  end subroutine calc_sqrt_S_Gamma

  subroutine calc_Eig_Gamma(spH,spS,orb,eig,state)
    use m_ts_sparse_helper, only : create_U
    type(dSpData2D), intent(inout) :: spH
    type(dSpData1D), intent(inout) :: spS
    type(tRegion), intent(in) :: orb
    real(dp), intent(out) :: eig(orb%n), state(orb%n,orb%n)

    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: no, info
    integer :: i, j, n_nzs
    real(dp), pointer :: H(:,:), S(:)
    real(dp), allocatable :: H_UT(:), S_UT(:), work(:)
#ifdef DIVIDE_AND_CONQUER
    integer :: lwork, liwork
    integer, allocatable :: iwork(:)
#endif

    no = orb%n

    ! Re-create H_UT, S_UT in UT format
    allocate(H_UT(no*(no+1)/2),S_UT(no*(no+1)/2))

#ifdef DIVIDE_AND_CONQUER
    ! Do a work-size query
    call dspgvd(1,'V','U', no, H_UT, S_UT, eig, state, no, &
         S_UT(1), -1, liwork, -1, info)
    lwork = max(no,nint(S_UT(1)))
    allocate(work(lwork))
    allocate(iwork(liwork))
#endif

    dit => dist(spH)
    sp => spar(spH)
    S => val(spS)
    n_nzs = size(S)
    call create_U(dit, sp, no, orb, n_nzs, S, S_UT)
    H => val(spH)
    call create_U(dit, sp, no, orb, n_nzs, H(:,1), H_UT)

#ifdef DIVIDE_AND_CONQUER
    call dspgvd(1,'V','U', no, H_UT, S_UT, eig, state, no, &
         work, lwork, iwork, liwork, info)
#else
    allocate(work(3*no))
    call dspgv(1,'V','U', no, H_UT, S_UT, eig, state, no, &
         work, info)
#endif

    deallocate(H_UT,S_UT)
    if ( info /= 0 ) then
       write(*,'(a)')'Error in diagonalization of molecule, H,S'
#ifdef DIVIDE_AND_CONQUER
       write(*,'(a,i0)')'LAPACK (dspgvd) error message: ',info
#else
       write(*,'(a,i0)')'LAPACK (dspgv) error message: ',info
#endif
       call die('Error in Gamma diagonalization of molecule, H, S')
    end if
#ifdef DIVIDE_AND_CONQUER
    deallocate(iwork)
#endif

    ! Sort the eigen-values AND vectors (they most probably
    ! are sorted)
    do i = 1 , no - 1
       ! find minimun eigen-value in non-sorted region
       j = i + minloc(eig(i+1:no),1)
       ! check whether we should switch
       if ( eig(i) > eig(j) ) then
          ! switch wf_eigenvalue
          work(1)    = eig(j)
          eig(j)     = eig(i)
          eig(i)     = work(1)
          ! switch eigenvector
          work(1:no) = state(:,j)
          state(:,j) = state(:,i)
          state(:,i) = work(1:no)
       end if
    end do

    deallocate(work)

  end subroutine calc_Eig_Gamma

  subroutine norm_eigenstate_Gamma(no,state,S_sq)
    use intrinsic_missing, only : VNORM
    integer, intent(in) :: no
    real(dp), intent(inout) :: state(no,no)
    real(dp), intent(in) :: S_sq(no,no)
    real(dp) :: work(no)
    integer :: i
    
    do i = 1 , no

       ! Normalize eigenvectors and create orthogonal basis
       call dgemm('N','N',no,1,no,1._dp, &
            S_sq,no, state(1,i),no, &
            0._dp, work(1), no)
       state(:,i) = work(1:no) / VNORM(work(1:no))

    end do
    
  end subroutine norm_eigenstate_Gamma


  subroutine calc_sqrt_S_kpt(spS,nsc,sc_off,orb,S_sq,kpt)
    use m_ts_sparse_helper, only : create_U
    type(dSpData1D), intent(inout) :: spS
    integer, intent(in) :: nsc
    real(dp), intent(in) :: sc_off(3,nsc)
    type(tRegion), intent(in) :: orb
    complex(dp), intent(out) :: S_sq(orb%n,orb%n)
    real(dp), intent(in) :: kpt(3)

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: no, n_nzs, i, j, info
    real(dp), pointer :: S(:)
    real(dp), allocatable :: eig(:), rwork(:)
    complex(dp), allocatable :: v(:,:), S_UT(:), work(:)

    no = orb%n
    dit => dist(spS)
    sp => spar(spS)
    S => val(spS)
    n_nzs = size(S)

    allocate(eig(no))
    allocate(v(no,no))
    allocate(rwork(3*no)) ! real-work
    allocate(work(2*no))
    allocate(S_UT(no*no))
    
    call create_U(dit, sp, no, orb, n_nzs, nsc, S, sc_off, S_UT, kpt)

    ! Diagonalize overlap matrix
    call zhpev('V','U',no,S_UT,eig,v,no,work,rwork,info)
    if ( info /= 0 ) then
       write(*,*) 'INFO = ',info,' when diagonalizing kpt overlap matrix'
    end if

    ! There are faster ways of doing this, but let's play safe for
    ! now. The overlap matrix is a positive definite matrix,
    ! hence all eigenvalues are positive.
    ! This "check" ensures that this is enforced.
    if ( any(eig < 0._dp) ) then
       write(*,'(3a,e12.5)')'tbtrans: Projection ',trim(orb%name), &
            ' is not completely positive definite, lowest eig of S: ', &
            minval(eig)
    end if
    where ( eig < 0.0_dp )
       eig = 0._dp
    elsewhere
       eig = dsqrt(eig)
    end where

    ! Calculate S^(1/2)
    ! Calculate v.sqrt(eig)
    j = 1
    do i = 1 , no
       S_UT(j:j+no-1) = v(:,i) * eig(i)
       j = j + no
    end do
    
    ! Calculate v.sqrt(eig).v^\dagger = S^(1/2)
    call zgemm('N','C',no,no,no,dcmplx(1._dp,0._dp), &
         S_UT,no,v,no, &
         dcmplx(0._dp,0._dp),S_sq,no)

    deallocate(eig,v,work,rwork,S_UT)

  end subroutine calc_sqrt_S_kpt

  subroutine calc_Eig_kpt(spH,spS,nsc,sc_off,orb,eig,state,kpt)
    use m_ts_sparse_helper, only : create_U

    type(dSpData2D), intent(inout) :: spH
    type(dSpData1D), intent(inout) :: spS
    integer, intent(in) :: nsc
    real(dp), intent(in) :: sc_off(3,nsc)
    type(tRegion), intent(in) :: orb
    real(dp), intent(out) :: eig(orb%n)
    complex(dp), intent(out) :: state(orb%n,orb%n)
    real(dp), intent(in) :: kpt(3)

    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: no, info
    integer :: i, j, n_nzs
    real(dp), pointer :: H(:,:), S(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: H_UT(:), S_UT(:), work(:)
#ifdef DIVIDE_AND_CONQUER
    integer :: lwork, lrwork
    integer, allocatable :: iwork(:)
#endif

    no = orb%n

    ! Re-create H_UT, S_UT in UT format
    allocate(H_UT(no*(no+1)/2),S_UT(no*(no+1)/2))

#ifdef DIVIDE_AND_CONQUER
    ! Do a work-size query
    call zhpgvd(1,'V','U', no, H_UT, S_UT, eig, state, no, &
         S_UT(1), -1, eig, -1, liwork, -1, info)
    lwork = max(no,nint(real(S_UT(1),dp)))
    allocate(work(lwork))
    lrwork = nint(eig(1))
    allocate(rwork(lrwork))
    allocate(iwork(liwork))
#endif

    dit => dist(spH)
    sp => spar(spH)
    S => val(spS)
    n_nzs = size(S)
    call create_U(dit, sp, no, orb, n_nzs, nsc, S, sc_off,S_UT, kpt)
    H => val(spH)
    call create_U(dit, sp, no, orb, n_nzs, nsc, H(:,1), sc_off,H_UT, kpt)

#ifdef DIVIDE_AND_CONQUER
    call zhpgvd(1,'V','U', no, H_UT, S_UT, eig, state, no, &
         work, lwork, rwork, lrwork, iwork, liwork, info)
#else
    allocate(work(2*no),rwork(3*no))
    call zhpgv(1,'V','U', no, H_UT, S_UT, eig, state, no, &
         work, rwork, info)
#endif

    deallocate(H_UT,S_UT)
    if ( info /= 0 ) then
       write(*,'(a)')'Error in diagonalization of molecule, H,S'
#ifdef DIVIDE_AND_CONQUER
       write(*,'(a,i0)')'LAPACK (zhpgvd) error message: ',info
#else
       write(*,'(a,i0)')'LAPACK (zhpgv) error message: ',info
#endif
       call die('Error in k-point diagonalization of molecule, H, S')
    end if
#ifdef DIVIDE_AND_CONQUER
    deallocate(iwork)
#endif

    ! Sort the eigen-values AND vectors (they most probably
    ! are sorted)
    do i = 1 , no - 1
       ! find minimun eigen-value in non-sorted region
       j = i + minloc(eig(i+1:no),1)
       ! check whether we should switch
       if ( eig(i) > eig(j) ) then
          ! switch wf_eigenvalue
          rwork(1)   = eig(j)
          eig(j)     = eig(i)
          eig(i)     = rwork(1)
          ! switch eigenvector
          work(1:no) = state(:,j)
          state(:,j) = state(:,i)
          state(:,i) = work(1:no)
       end if
    end do
    
    deallocate(work,rwork)

  end subroutine calc_Eig_kpt

  subroutine norm_eigenstate_kpt(no,state,S_sq)
    use intrinsic_missing, only : VNORM
    integer, intent(in) :: no
    complex(dp), intent(inout) :: state(no,no)
    complex(dp), intent(in) :: S_sq(no,no)
    complex(dp) :: work(no)
    integer :: i
    
    do i = 1 , no

       ! Normalize eigenvectors and create orthogonal basis
       call zgemm('N','N',no,1,no,dcmplx(1._dp,0._dp), &
            S_sq,no, state(1,i),no, &
            dcmplx(0._dp,0._dp), work(1), no)
       state(:,i) = work(1:no) / VNORM(work(1:no))

    end do
    
  end subroutine norm_eigenstate_kpt
   
end module m_tbt_diag
