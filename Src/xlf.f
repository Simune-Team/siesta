      SUBROUTINE CPUTIM (TIME)

      DOUBLE PRECISION TIME
      REAL TIMEreal
C
      call cpu_time(timereal)
      TIME = timereal
      END
      
      subroutine flush(lun)
      integer, intent(in) :: lun
      external flush_
      call flush_(lun)
      end subroutine flush
