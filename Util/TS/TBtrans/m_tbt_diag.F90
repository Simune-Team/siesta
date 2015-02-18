
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

  public :: init_diag
  public :: calc_sqrt_S
  public :: calc_Eig
  public :: norm_Eigenstate

  logical, save :: use_DC = .false.

contains

  subroutine init_diag( )
    use fdf
    use parallel, only : Node

    ! Let the user decide whether to use divide and conquer
    ! or not...
    use_DC = fdf_get('TBT.DivideAndConquer',.false.)
    if ( Node == 0 ) then
       write(*,'(a,t53,''='',tr4,l1)')'tbt_options: Divide and conquer diagonalization',use_DC
    end if

  end subroutine init_diag

  subroutine calc_sqrt_S_Gamma(spS,orb,S_sq)
    use m_ts_sparse_helper, only : create_U
    type(dSpData1D), intent(inout) :: spS
    type(tRgn), intent(in) :: orb
    real(dp), intent(out) :: S_sq(orb%n,orb%n)

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: no, n_nzs, i, j, info
    real(dp), pointer :: S(:)
    real(dp), allocatable :: eig(:), v(:,:), S_UT(:), work(:)
    integer :: lwork, liwork
    integer, allocatable :: iwork(:)

    no  =  orb%n
    dit => dist(spS)
    sp  => spar(spS)
    S   => val( spS)
    n_nzs = size(S)

    allocate(eig(no))
    allocate(v(no,no))
    allocate(S_UT(no*no))
    lwork = 3*no
    if ( use_DC ) then
       ! work-size query
       call dspevd('V','U',no,S_UT,eig,v,no,eig,-1,liwork,-1,info)
       lwork = max(no,nint(eig(1)))
       allocate(iwork(liwork))
    end if
    allocate(work(lwork)) ! real-work

    call create_U(dit, sp, no, orb, n_nzs,S, S_UT)

    ! Diagonalize overlap matrix
    if ( use_DC ) then
       call dspevd('V','U',no,S_UT,eig,v,no,work,lwork,iwork,liwork,info)
       deallocate(iwork)
    else       
       call dspev('V','U',no,S_UT,eig,v,no,work,info)
    end if
    if ( info /= 0 ) then
       write(*,'(a)')'Error in diagonalization of molecule, S'
       if ( use_DC ) then
          write(*,'(a,i0)')'LAPACK (dspevd) error message: ',info
       else
          write(*,'(a,i0)')'LAPACK (dspev) error message: ',info
       end if
       call die('Error in Gamma diagonalization of molecule, S')
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
    type(tRgn), intent(in) :: orb
    real(dp), intent(out) :: eig(orb%n)
    real(dp), intent(out), target, optional :: state(orb%n,orb%n)

    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: no, info
    integer :: i, j, n_nzs
    real(dp), pointer :: H(:,:), S(:), lstate(:,:)
    real(dp), allocatable :: H_UT(:), S_UT(:), work(:)
    integer :: lwork, liwork
    integer, allocatable :: iwork(:)
    ! Character to denote whether only eigenvalue or +eigenstates
    character(len=1) :: NV

    no = orb%n
    
    if ( present(state) ) then
       NV = 'V'
       lstate => state(:,:)
       ! if ( size(lstate,1) /= orb%n )
       ! if ( size(lstate,2) /= orb%n )
    else
       NV = 'N'
       nullify(lstate) ; allocate(lstate(1,1))
    end if

    ! Re-create H_UT, S_UT in UT format
    allocate(H_UT(no*(no+1)/2),S_UT(no*(no+1)/2))

    lwork = 3 * no
    if ( use_DC ) then
       ! Do a work-size query
       call dspgvd(1,NV,'U', no, H_UT, S_UT, eig, lstate, no, &
            S_UT(1), -1, liwork, -1, info)
       lwork = max(no,nint(S_UT(1)))
       allocate(iwork(liwork))
    end if
    allocate(work(lwork))

    dit => dist(spH)
    sp  => spar(spH)
    S   => val( spS)
    n_nzs = size(S)
    call create_U(dit, sp, no, orb, n_nzs, S, S_UT)
    H   => val( spH)
    call create_U(dit, sp, no, orb, n_nzs, H(:,1), H_UT)

    if ( use_DC ) then
       call dspgvd(1,NV,'U', no, H_UT, S_UT, eig, lstate, no, &
            work, lwork, iwork, liwork, info)
       deallocate(iwork)
    else
       call dspgv(1,NV,'U', no, H_UT, S_UT, eig, lstate, no, &
            work, info)
    end if

    deallocate(H_UT,S_UT)

    if ( info /= 0 ) then
       write(*,'(a)')'Error in diagonalization of molecule, H, S'
       if ( use_DC ) then
          write(*,'(a,i0)')'LAPACK (dspgvd) error message: ',info
       else
          write(*,'(a,i0)')'LAPACK (dspgv) error message: ',info
       end if
       call die('Error in Gamma diagonalization of molecule, H, S')
    end if

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
          if ( NV == 'V' ) then
             ! switch eigenvector
             work(1:no) = lstate(:,j)
             lstate(:,j) = lstate(:,i)
             lstate(:,i) = work(1:no)
          end if
       end if
    end do

    deallocate(work)
    if ( NV == 'N' ) deallocate(lstate)

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
    type(tRgn), intent(in) :: orb
    complex(dp), intent(out) :: S_sq(orb%n,orb%n)
    real(dp), intent(in) :: kpt(3)

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: no, n_nzs, i, j, info
    real(dp), pointer :: S(:)
    real(dp), allocatable :: eig(:), rwork(:)
    complex(dp), allocatable :: v(:,:), S_UT(:), work(:)
    integer :: lrwork, lwork, liwork
    integer, allocatable :: iwork(:)

    no  = orb%n
    dit => dist(spS)
    sp  => spar(spS)
    S   => val( spS)
    n_nzs = size(S)

    allocate(eig(no))
    allocate(v(no,no))
    allocate(S_UT(no*no))
    lrwork = 3 * no
    lwork  = 2 * no
    if ( use_DC ) then
       ! work-size query
       call zhpevd('V','U',no,S_UT,eig,v,no,S_UT,-1, &
            eig,-1, liwork, -1, info)
       lwork  = nint(real(S_UT(1),dp))
       lrwork = nint(eig(1))
       allocate(iwork(liwork))
    end if
    allocate(rwork(lrwork)) ! real-work
    allocate(work(lwork))
    
    call create_U(dit, sp, no, orb, n_nzs, nsc, S, sc_off, S_UT, kpt)

    ! Diagonalize overlap matrix
    if ( use_DC ) then
       call zhpevd('V','U',no,S_UT,eig,v,no,work,lwork, &
            rwork, lrwork, iwork, liwork, info)
       deallocate(iwork)
    else
       call zhpev('V','U',no,S_UT,eig,v,no,work,rwork,info)
    end if
    if ( info /= 0 ) then
       write(*,'(a)')'Error in diagonalization of molecule, S'
       if ( use_DC ) then
          write(*,'(a,i0)')'LAPACK (zhpevd) error message: ',info
       else
          write(*,'(a,i0)')'LAPACK (zhpev) error message: ',info
       end if
       call die('Error in k-point diagonalization of molecule, S')
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
    type(tRgn), intent(in) :: orb
    real(dp), intent(out) :: eig(orb%n)
    complex(dp), intent(out), target, optional :: state(orb%n,orb%n)
    real(dp), intent(in) :: kpt(3)

    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: no, info
    integer :: i, j, n_nzs
    real(dp), pointer :: H(:,:), S(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: H_UT(:), S_UT(:), work(:)
    complex(dp), pointer :: lstate(:,:)
    integer :: lwork, lrwork, liwork
    integer, allocatable :: iwork(:)
    character(len=1) :: NV

    no = orb%n

    if ( present(state) ) then
       NV = 'V'
       lstate => state(:,:)
       ! if ( size(lstate,1) /= orb%n )
       ! if ( size(lstate,2) /= orb%n )
    else
       NV = 'N'
       nullify(lstate) ; allocate(lstate(1,1))
    end if

    ! Re-create H_UT, S_UT in UT format
    allocate(H_UT(no*(no+1)/2),S_UT(no*(no+1)/2))

    lwork  = 2 * no
    lrwork = 3 * no
    if ( use_DC ) then
       ! work-size query
       call zhpgvd(1,NV,'U', no, H_UT, S_UT, eig, lstate, no, &
            S_UT(1), -1, eig, -1, liwork, -1, info)
       lwork = max(no,nint(real(S_UT(1),dp)))
       lrwork = nint(eig(1))
       allocate(iwork(liwork))
    end if
    allocate(work(lwork))
    allocate(rwork(lrwork))

    dit => dist(spH)
    sp  => spar(spH)
    S   => val( spS)
    n_nzs = size(S)
    call create_U(dit, sp, no, orb, n_nzs, nsc, S, sc_off,S_UT, kpt)
    H   => val( spH)
    call create_U(dit, sp, no, orb, n_nzs, nsc, H(:,1), sc_off,H_UT, kpt)

    if ( use_DC ) then
       call zhpgvd(1,NV,'U', no, H_UT, S_UT, eig, lstate, no, &
            work, lwork, rwork, lrwork, iwork, liwork, info)
       deallocate(iwork)
    else
       call zhpgv(1,NV,'U', no, H_UT, S_UT, eig, lstate, no, &
         work, rwork, info)
    end if

    deallocate(H_UT,S_UT)

    if ( info /= 0 ) then
       write(*,'(a)')'Error in diagonalization of molecule, H,S'
       if ( use_DC ) then
          write(*,'(a,i0)')'LAPACK (zhpgvd) error message: ',info
       else
          write(*,'(a,i0)')'LAPACK (zhpgv) error message: ',info
       end if
       call die('Error in k-point diagonalization of molecule, H, S')
    end if

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
          if ( NV == 'V' ) then
             ! switch eigenvector
             work(1:no)  = lstate(:,j)
             lstate(:,j) = lstate(:,i)
             lstate(:,i) = work(1:no)
          end if
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
