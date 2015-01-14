program test_psml

use m_psml, only: ps_t, ps_destroy, psml_reader   ! Clarify this
use m_psml, only: ps_NProjectors, ps_NPotentials

!use m_semicore_info, only: get_n_semicore_shells
use m_getopts

type(ps_t)   :: ps
integer      :: nsemic(0:3), n

      character(len=200) :: filename
      character(len=200) :: opt_arg, mflnm, ref_line
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels, iorb, ikb
!
!     Process options
!
      n_opts = 0
      do
         call getopts('h',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('h')
            call manual()
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: test_psml [ -h ] FILE"
             write(0,*) "Use -h option for manual"
             STOP
          end select
       enddo

       nargs = command_argument_count()
       nlabels = nargs - n_opts + 1
       if (nlabels /= 1)  then
          write(0,*) "Usage: test_psml [ -h ] FILE"
          write(0,*) "Use -h option for manual"
          STOP
       endif

       call get_command_argument(n_opts,value=filename,status=iostat)
       if (iostat /= 0) then
          STOP "Cannot get filename"
       endif

call ps_destroy(ps)
call psml_reader(filename,ps,debug=.true.)
print *, "Number of semilocal potentials: ", ps_NPotentials(ps)
print *, "Number of fully non-local potentials: ", ps_NProjectors(ps)

!call get_n_semicore_shells(ps,nsemic)

!print "(a,4i4)", "Number of semicore shells:", nsemic(0:3)

CONTAINS
subroutine manual()
end subroutine manual
end program test_psml

