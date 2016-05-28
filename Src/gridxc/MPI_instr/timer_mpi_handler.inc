! This file is to be included as the body
! of the 'timer_mpi' routine in all the
! instrumented modules that call it.
!
! Each library using the modules can
! have its own version of this little file
!
! This version is for gridxc
!
! Alberto Garcia, June 2015

  interface 
     subroutine gridxc_timer_mpi(prog,iopt)
       character(len=*), intent(in) :: prog
       integer, intent(in)          :: iopt
     end subroutine gridxc_timer_mpi
  end interface

  call gridxc_timer_mpi(prog,iopt)
