module m_svd
public :: solve_with_svd
CONTAINS
subroutine solve_with_svd(ain,c,rank_out)
  use precision, only: dp

  real(dp), intent(in) :: ain(:,:)
  real(dp), intent(out) :: c(:)
  integer, intent(out), optional :: rank_out

  real(dp), allocatable :: a(:,:), s(:), work(:), b(:)

  integer, parameter :: nb = 64
  integer :: n, lwork
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)

      REAL(DP) RCOND, RNORM
      INTEGER  ::    I, INFO, J, M, RANK, LDA

      REAL(DP) DNRM2
      EXTERNAL         DNRM2

      EXTERNAL         DGELSS

  n = size(ain,dim=1)
  lda = n
  m = n
  allocate(a(n,n))
  a = ain
  allocate(b(n),s(n))

  lwork = 3*n+nb*(n+m)
  allocate(work(lwork))

      b(:) = 0.0_dp
      b(n) = 1.0_dp 

!        Choose RCOND to reflect the relative accuracy of the input data

         RCOND = 1.0e-6_dp
         ! Singular values s_i < rcond*s_1 will be neglected for
         ! the estimation of the rank.

!        Solve the least squares problem min( norm2(b - Ax) ) for the x
!        of minimum norm.

         CALL DGELSS(M,N,1,A,LDA,B,M,S,RCOND,RANK,WORK,LWORK,INFO)
!
         IF (INFO.EQ.0) THEN

            c(1:n-1) = b(1:n-1)

            WRITE (NOUT,"(a,i3,e11.2)") &
              'Estimated rank of DIIS matrix, (tol)', RANK, RCOND
            WRITE (NOUT,"(a,7f11.4)") 'Singular values: ', &
                                      (S(I),I=1,N)
            if (present(rank_out)) then
               rank_out = rank
            endif

!           Compute and print estimate of the square root of the
!           residual sum of squares

            IF (RANK.EQ.N) THEN
               RNORM = DNRM2(M-N,B(N+1),1)
               !WRITE (NOUT,*)  'Square root of the residual sum of squares'
               !WRITE (NOUT,*) RNORM
            END IF
         ELSE
            WRITE (NOUT,*) 'The SVD algorithm failed to converge'
            print *, "info: ", info
            c(:) = 0.0_dp
            c(n-1) = 1.0_dp
            WRITE (NOUT,*) 'Re-using last item only...'
         END IF

         WRITE (NOUT,"(a,f12.8)") 'Sum of coeffs: ', sum(c(1:n-1))
         deallocate(a,s,work,b)

END subroutine solve_with_svd
end module m_svd
