      SUBROUTINE CPUTIM (TIME)

      DOUBLE PRECISION TIME
      REAL TIMEreal
C
      call cpu_time(timereal)
      TIME = timereal
      END
