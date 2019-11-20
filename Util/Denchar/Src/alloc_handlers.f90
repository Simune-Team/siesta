   ! Error message and integer code
   ! If 'code' is 0, this is the last call in a series
   ! (see below for usage)
   subroutine alloc_error_report(str,code)
     character(len=*), intent(in) :: str
     integer, intent(in)          :: code
   end subroutine alloc_error_report
   !
   ! Logger for memory events
   !
   subroutine alloc_memory_event(bytes,name)
     integer, intent(in)          :: bytes
     character(len=*), intent(in) :: name
   end subroutine alloc_memory_event
