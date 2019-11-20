!> Extracts the semilocal potentials, the valence charge, and the core
!> charge (if present) from a PSML pseudopotential file, and outputs a
!> PSF pseudopotential file (in the Froyen classic format used in
!> Siesta).
!>
!> The program will output as much metadata as it is
!> possible to include in a PSF file.
!> To cover the variety of XC functionals allowed in PSML,
!> this program generates an extended psf file, with an
!> extra packed libxc code at the end of the first line.
!> This is backward compatible.

!> The output grid parameters can be selected in the command line. By
!> default, they will be those used by modern versions of ATOM.

program psml2psf

use m_psml
use m_getopts
use m_ncps

integer, parameter :: dp = selected_real_kind(10,100)

  logical  :: reparametrize
  real(dp) :: a
  real(dp) :: b
  real(dp) :: rmax, znuc
  integer  :: n_slpots
  integer, allocatable :: sl_idx(:)
  
        ! These are the "current" ATOM parameters
        ! There is another set turned on by the UCB_COMPAT flag in ATOM
        real(dp), parameter :: aa_def = 6.0_dp
        real(dp), parameter :: bb_def = 80.0_dp    ! UCB_COMPAT: 40.0
        real(dp), parameter :: rmax_def = 120.0_dp ! UCB_COMPAT: 80.0

character(len=200) :: filename, outfile

      logical            :: debug
      character(len=200) :: opt_arg
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels

  type(ps_t)           :: ps
  type(froyen_ps_t)    :: p

  a = 0.0_dp
  b = 0.0_dp
  Rmax = 0.0_dp
  
!     Process options
!
      n_opts = 0
      debug = .false.
      reparametrize = .false.
      outfile = "PSF_from_PSML"     ! default output name
      do
         call getopts('do:ra:b:R:',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('d')
              debug = .true.
           case ('r')
              reparametrize = .true.
           case ('R')
              read(opt_arg,*) Rmax
           case ('o')
              read(opt_arg,*) outfile
           case ('a')
              read(opt_arg,*) a
           case ('b')
              read(opt_arg,*) b
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: psml2psf [ -d -o Outfile -r -a a -b b -R Rmax ] FILE"
             write(0,*) " -d for debug output"
             write(0,*) " -r for reparametrization"
             write(0,*) " -R Rmax : maximum range of grid when reparametrizing"
             write(0,*) "         default: 120.0 au    "
             write(0,*) " -a a :  r(i) = b*(exp(ia)-1) "
             write(0,*) "         default: 1./80.      "
             write(0,*) " -b b :  r(i) = b*(exp(ia)-1) "
             write(0,*) "         default: exp(-6)/znuc"
             write(0,*) " -o Outfile  : output file name"
             STOP
          end select
       enddo

       nargs = command_argument_count()
      nlabels = nargs - n_opts + 1
       if (nlabels /= 1)  then
             write(0,*) "Usage: psml2psf [ -d -o Outfile -r -a a -b b -R Rmax ] FILE"
             write(0,*) " -d for debug output"
             write(0,*) " -r for reparametrization"
             write(0,*) " -R Rmax : maximum range of grid when reparametrizing"
             write(0,*) "         default: 120.0 au    "
             write(0,*) " -a a :  r(i) = b*(exp(ia)-1) "
             write(0,*) "         default: 1./80.      "
             write(0,*) " -b b :  r(i) = b*(exp(ia)-1) "
             write(0,*) "         default: exp(-6)/znuc au"
             write(0,*) " -o Outfile  : output file name"
          STOP
       endif

       call get_command_argument(n_opts,value=filename,status=iostat)
       if (iostat /= 0) then
          STOP "Cannot get filename"
       endif
!
  call ps_destroy(ps)
  call psml_reader(filename,ps)
  !
  ! Sanity check
  !
  call ps_Potential_Filter(ps,set=SET_ALL,number=n_slpots,indexes=sl_idx)
  if (n_slpots == 0) then
     write(0,*) "File " // trim(filename) // " does not contain semilocal potentials"
     !stop
  endif
  !
  if (reparametrize) then
     call ps_PseudoAtomSpec_Get(ps,atomic_number=znuc)
     ! Set defaults if not set by the user
     if (a == 0.0_dp) a = 1.0_dp / bb_def
     if (b == 0.0_dp) b = exp(-aa_def)/znuc
     if (Rmax == 0.0_dp) Rmax = rmax_def
  endif

  ! Uses the functionality of the ncps library
  call ncps_psml2froyen(ps,p,reparametrize,a,b,rmax)
  call pseudo_write_formatted(outfile,p)
  
end program psml2psf
