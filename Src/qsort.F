      recursive subroutine myQsort( N, A, perm )
      implicit none
C     Input variables
      integer       :: N, A(N), perm(N)
C     Local variables
      integer       :: iq

!------------------------------------------------------------------------- BEGIN
      if (N.gt.1) then
        call Partition( N, A, perm, iq )
        call myQsort( iq-1, A, perm )
        call myQsort( N-iq+1, A(iq), perm(iq) )
      endif
!--------------------------------------------------------------------------- END
      end subroutine myQsort

      subroutine Partition( N, A, perm, marker )
C     Input variables
      integer       :: N, A(N),  perm(N), marker
C     Local variables
      integer       :: i, j, temp, x

!------------------------------------------------------------------------- BEGIN
      x = A(1)
      i = 0
      j = N + 1
      do
        j = j-1
        do
          if (A(j) >= x) exit
          j = j-1
        end do
        i = i+1
        do
          if (A(i) <= x) exit
          i = i+1
        end do
        if (i < j) then
C         exchange A(i) and A(j)
          temp    = A(i)
          A(i)    = A(j)
          A(j)    = temp
          temp    = perm(i)
          perm(i) = perm(j)
          perm(j) = temp
        elseif (i == j) then
          marker = i+1
          return
        else
          marker = i
          return
        endif
      end do
!--------------------------------------------------------------------------- END
      end subroutine Partition
