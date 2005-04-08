      SUBROUTINE CPUTIM (TIME)

      DOUBLE PRECISION TIME
      REAL TIMEreal
C
      call cpu_time(timereal)
      TIME = timereal
      END

      subroutine abort(str)
      use f90_unix, only: local_abort=>abort
      character(len=*), optional, intent(in) :: str

      call local_abort(str)
      end subroutine abort
