! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_lowdin

  implicit none
  private

  public :: Lowdin_Trans

contains
! ##################################################################
! ##                Make Lowdin transformation                    ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################
  subroutine Lowdin_Trans(flag,n,mat,trans)
    use precision, only : dp

! ********************
! * INPUT variables  *
! ********************
    logical,intent(in) :: flag       ! .true. : T = S and on return T=S**(-1/2)
!                                    ! .false.: T = S**(-1/2) and is not changed
    integer,intent(in) :: n          ! Matrix size

! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(inout) :: mat(n,n)   ! on return: T.M.T where T=S**(-1/2)
    complex(dp), intent(inout) :: trans(n,n) ! or T=S or T=S**(-1/2)     

! ********************
! * LOCAL variables  *
! ********************
    complex(dp), dimension(n,n) :: X
    complex(dp) :: a,b
    integer :: i,j

    call timer('Lowdin',1)

    a = dcmplx(1.0_dp,0.0_dp)
    b = dcmplx(0.0_dp,0.0_dp)

    if (flag) then
       do j=1,n     
          do i=1,n
             X(i,j)=trans(i,j)
          end do
       end do
       call zinvsqrtM(n,X,trans)
    end if

    call zgemm('N','N',n,n,n,a,trans,n,mat,n,b,X,n)
    call zgemm('N','N',n,n,n,a,X,n,trans,n,b,mat,n)

    call timer('Lowdin',2)

  end subroutine Lowdin_Trans

!*******************************************************************       
! ##################################################################
! ##         Calculate square root of hermitian matrix            ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################
  subroutine zsqrtM(n,M,SQM)
    use precision, only : dp

! ********************
! * INPUT variables  *
! ********************
    integer,    intent(in)   :: n
    complex(dp), intent(in)  :: M(n,n)
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(out) :: SQM(n,n)
! ********************
! * LOCAL variables  *
! ********************
!     LAPACK diagononalization:
    complex(dp), dimension(:),   allocatable :: zwork, smat
    complex(dp), dimension(:,:), allocatable :: svect
    real(dp),    dimension(:),   allocatable :: rwork, seig
    integer     :: i,j,k, info

    allocate(zwork(2*n-1))
    allocate(smat(n*(n+1)/2))
    allocate(rwork(3*n-2))
    allocate(svect(n,n))
    allocate(seig(n))

! upper triangular form
    do j=1,n
       do i=1,n
          smat(i + ((j-1)*j/2)) = M(i,j)
       end do !i
    end do    !j

!     Diagonalize 
    call zhpev('V','U',n,smat,seig,svect,n,zwork,rwork,info)
    if (info.ne.0) then
       write(*,*) 'INFO = ',info,' when diagonalizing in sqrtM'
    end if
    deallocate(zwork)
!
!     There are faster ways of doing this, but let's play safe for
!     now.  Form the sqM matrix
!
    do i = 1,n
       if(seig(i) .lt. 0.0_dp) then
          seig(i) = 0.0_dp
       else
          seig(i) = dsqrt(seig(i))
       end if
    end do

    do i = 1,n
       do j = 1,i
          sqM(i,j) = dcmplx(0.0_dp,0.0_dp)
          do k = 1,n
             sqM(i,j) = sqM(i,j) + seig(k)*svect(i,k)*dconjg(svect(j,k))
          end do
          sqM(j,i) = dconjg(sqM(i,j))
       end do
    end do

    deallocate(rwork)
    deallocate(smat)
    deallocate(svect)
    deallocate(seig)

  end subroutine zsqrtM

! ##################################################################
! ##     Calculate inverse square root of hermitian matrix        ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################
  subroutine zinvsqrtM(n,M,SQM)
    use precision, only  : dp

! ********************
! * INPUT variables  *
! ********************
    integer,     intent(in)  :: n
    complex(dp), intent(in)  :: M(n,n)
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(out) :: SQM(n,n)

! ********************
! * LOCAL variables  *
! ********************
    real(dp) , parameter :: EPSnegl = 1e-10_dp

    complex(dp), dimension(:),   allocatable :: zwork, smat
    complex(dp), dimension(:,:), allocatable :: svect
    real(dp)   , dimension(:),   allocatable :: rwork, seig
    integer :: i,j,k, info

    allocate(zwork(2*n-1))
    allocate(smat(n*(n+1)/2))
    allocate(rwork(3*n-2))
    allocate(svect(n,n))
    allocate(seig(n))

! upper triangular form
    smat=dcmplx(0.0_dp,0.0_dp)
    do j=1,n
       do i=1,n
          smat(i + ((j-1)*j/2)) = M(i,j)
       end do
    end do

!     Diagonalize 
    call zhpev('V','U',n,smat,seig,svect,n,zwork,rwork,info)
    deallocate(zwork)

    if (info.ne.0) then
       write(*,*) 'INFO = ',info,' when diagonalizing in sqrtM'
    end if

!    if (PRINTALOT) then
!       write(jTS,'(a)') 'Lowdin orthogonalizaton: diagonalization'
!       write(jTS,'(a)') 'Eigenvalues of matrix :'
!       write(jTS,'(6f15.5)') (seig(i),i=1,n)
!       write(jTS,*)      
!    end if
!
!     There are faster ways of doing this, but let's play safe for
!     now.  Form the sqM matrix
!
    do i = 1,n
       if(seig(i) .lt. EPSnegl) then 
          seig(i)= 0.0_dp
       else 
          seig(i) = 1.0_dp/dsqrt(seig(i))
       end if
    end do

    do i = 1,n
       do j = 1,i
          SQM(i,j) = dcmplx(0.0_dp,0.0_dp)
          do k = 1,n
             SQM(i,j) = SQM(i,j) + seig(k)*svect(i,k)*dconjg(svect(j,k))
          end do
          SQM(j,i) = dconjg(SQM(i,j))
       end do
    end do

    deallocate(rwork)
    deallocate(smat)
    deallocate(svect)
    deallocate(seig)

  end subroutine zinvsqrtM

end module m_lowdin
