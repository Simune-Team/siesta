program test_ncps

use m_ncps
use m_getopts

integer, parameter :: dp = selected_real_kind(10,100)

type(froyen_ps_t)   :: ps

      character(len=200) :: filename
      logical            :: debug, write_psf
      character(len=200) :: opt_arg, mflnm, ref_line
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels

      integer :: status
      
!
!     Process options
!
      n_opts = 0
      debug = .false.
      write_psf = .false.

      Rmax = 4.0_dp   ! default maximum radius for plots
      npts = 200      ! default number of points

      do
         call getopts('dR:n:f',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('d')
              debug = .true.
           case ('f')
              write_psf = .true.
           case ('R')
              read(opt_arg,*) Rmax
           case ('n')
              read(opt_arg,*) npts
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: test_ncps [ -d -R Rmax -n npts ] FILE"
             write(0,*) " -d for debug output"
             STOP
          end select
       enddo

       nargs = command_argument_count()
       nlabels = nargs - n_opts + 1
       if (nlabels /= 1)  then
          write(0,*) "Usage: test_ncps [ -d ] <PS file>"
          write(0,*) " -d for debug output"
          STOP
       endif

       call get_command_argument(n_opts,value=filename,status=iostat)
       if (iostat /= 0) then
          STOP "Cannot get filename"
       endif
!
       print "(a)", "Processing: " // trim(filename)
       call pseudo_read_from_file(filename,ps)
       if (write_psf) then
          print "(a)", "Writing psf version..."
          call pseudo_write_formatted(trim(filename)//".psf_out",ps,.true.)
       endif

end program test_ncps

