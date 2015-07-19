!--------------------------------------------------
! Stand-alone 'die' routine for use by libraries and
! low-level modules.

      subroutine die(str)

      character(len=*), intent(in)  :: str

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)

         call pxfflush(6)
         call pxfflush(0)
      call pxfabort()
      end subroutine die

      subroutine psml_die(str)
      character(len=*), intent(in)  :: str
      call die(str)
      end subroutine psml_die

      subroutine timer(str,iopt)
      character(len=*), intent(in)  :: str
      integer, intent(in)  :: iopt
      ! do nothing
      end subroutine timer

      
