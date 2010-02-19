      subroutine choleskyFact( n, lastn, AP, LL )
      use precision
      use sys
      implicit none
C     Input variables
      integer,   intent(in) :: n, lastn
      real(dp),  intent(in) :: AP((n*(n+1))/2)
C     Output variables
      real(dp), intent(out) :: LL((n*(n+1))/2)
C     Local variables
      integer               :: ii, jj, kk, ind1, ind2
      real(dp)              :: tt

      ind1 = (lastn*(lastn+1))/2
      do ii= lastn+1, n
        ind2 = 0
        do jj = 1, ii-1
          tt = AP(ind1+jj)
          do kk = 1, jj-1
            tt = tt - LL(ind2+kk)*LL(ind1+kk)
          enddo
          ind2 = ind2+jj
          LL(ind1+jj) = tt/LL(ind2)
        enddo

        tt = AP(ind1+ii)
        do jj= 1, ii-1
          tt = tt - LL(ind1+jj)*LL(ind1+jj)
        enddo
        if (tt.le.0.0)  call die( 'ERROR in choleskyFact' )
        ind1 = ind1+ii
        LL(ind1) = sqrt( tt )
      enddo
      end subroutine choleskyFact

      subroutine choleskySolve( n, LL, xx )
      use precision
      use sys
      implicit none
C     Input variables
      integer,   intent(in) :: n
      real(dp),  intent(in) :: LL((n*(n+1))/2)
C     Output variables
      real(dp), intent(out) :: xx(n)
C     Local variables
      integer               :: ii, jj, ind
      real(dp)              :: tt

!     xx = L^(-1)*xx
      ind = 0
      do ii= 1, n
        tt = xx(ii)
        do jj= 1, ii-1
          ind = ind + 1
          tt  = tt - xx(jj)*LL(ind)
        enddo
        ind = ind + 1
        xx(ii) = tt/LL(ind)
      enddo

!     xx = L^(-t)*xx
      ind = (n*(n+1))/2 + 1
      do ii= n, 1, -1
        ind = ind - 1
        xx(ii) = xx(ii)/LL(ind)
        do jj= ii-1, 1, -1
          ind = ind - 1
          xx(jj) = xx(jj) - xx(ii)*LL(ind)
        enddo
      enddo
      end subroutine choleskySolve
