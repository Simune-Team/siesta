! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
program filter_test

use m_filter,  only: kcphi, filter
use m_getopts, only: getopts

implicit none

integer, parameter :: dp = selected_real_kind(14,100)
real(dp), parameter :: EkinTolDefault = 0.003_dp
real(dp), parameter :: Kcmax = 30.0_dp  ! 900 Ry


integer  :: norm_opt = 0
integer  :: l = -1
integer  :: nr, n_basis_functions, i
   
real(dp)  :: kmax, kmax_tol, etol, filterFactor = 1.0_dp
real(dp)  :: rin, fin

real(dp), allocatable, dimension(:)  :: r, f, fold, rf, rfold


character(len=200) :: opt_arg
character(len=10)  :: opt_name 
integer :: nargs, iostat, n_opts, nlabels

character(len=200) :: filein = "INP", fileout = "OUT"


etol = EkinTolDefault

!
!     Process options
!
      n_opts = 0
      do
         call getopts('l:t:f:o:n:k:h',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
         case ('l')  
            read(opt_arg,*) l
         case ('t')
            read(opt_arg,*) etol
         case ('k')
            read(opt_arg,*) filterFactor
         case ('n')  
            read(opt_arg,*) norm_opt
         case ('f')
            filein = opt_arg
         case ('o')
            fileout = opt_arg
         case ('h')
            call print_help()
         case ('?',':')
            write(0,*) "Invalid option: ", opt_arg(1:1)
            call print_help()
            STOP
         end select
      enddo

if ( l == -1) then
   write(0,*) " ** Need to provide l"
   call print_help()
   STOP
endif

write(0,*) "L, etol, norm_opt: ", l, etol, norm_opt
write(0,*) "filterFactor: ", filterFactor
write(0,*) "FILEIN, FILEOUT: ", trim(filein), "  ", trim(fileout)

open(unit=1,file=filein,form="formatted",status="old",action="read",position="rewind")
nr = 1
do
      read(1,fmt=*,iostat=iostat)  rin, fin
      if (iostat /= 0) then
         nr = nr -1
         exit
      endif
      nr = nr + 1
enddo


allocate (r(nr), f(nr), fold(nr), rf(nr), rfold(nr))

rewind(unit=1)
do i = 1, nr
      read(1,fmt=*,iostat=iostat)  r(i), fold(i)
      if (iostat /= 0) then
         STOP "bad read"
      endif
enddo
close(unit=1)

! Note that the storage convention in Siesta is f(r) / r**l
! It is "rf" what we have to work with, if we are using
! files taken from the basis output of Siesta (.ion files)
! as is the case in this test directory.

rfold(:) = fold(:) * r(:) ** l
rf(:) = rfold(:)

kmax_tol = kcPhi(l,nr,r,rf,etol)
kmax = kmax_tol/filterFactor
if (kmax_tol > Kcmax) then
   write(6,"(a,f10.2,a)") "kmax too big ==> cutoff: ", kmax**2, "Ry"
   kmax = Kcmax/filterFactor
endif
write(6,'(a,f8.2,a)')   "Filter: Cutoff (before filtering):",kmax**2,' Ry'
call filter(l,nr,r,rf,kmax,norm_opt,n_basis_functions)

kmax_tol = kcPhi(l,nr,r,rf,etol)
write(6,'(a,f8.2,a)') "Filter: Appropiate mesh cutoff is (at least):" , &
                      kmax_tol**2*filterFactor,' Ry'  ! Why the factor?

do i = 2, nr
   f(i) = rf(i) / r(i) ** l
enddo
f(1) = f(2)

open(unit=2,file=fileout,form="formatted",status="unknown",action="write",position="rewind")

write(2,"(a,i3,f10.6,i4)") "# L, etol, norm_opt: ", l, etol, norm_opt
write(2,"(a,f10.6)") "# filterFactor: ", filterFactor
write(2,"(a)") "#"
write(2,fmt="(a1,a15,4a16)")  "#",  "r", "fold", "f", "r**l * fold", "r**l * f"
do i = 1, nr
      write(2,fmt="(5f16.8)")  r(i), fold(i), f(i), rfold(i), rf(i)
enddo
close(unit=2)


deallocate (r, f, fold, rf, rfold)

CONTAINS

      subroutine print_help()
         print *, "Usage: test_filter [options]"
         print *, "Options:"
         print *, " "
         print *, " -l l        : Set angular momentum l  (Mandatory)"
         print *, " -t etol     : Set kinetic energy tolerance (default = ", EkinTolDefault, ")"
         print *, " -k filterF  : Filter factor (usually 1.0, 0.7 for orbitals)"
         print *, " -n norm_opt : Normalization option (0: nothing, 1: norm f, 2: norm f**2)"
         print *, " -f INFILE   : Read data from FILE (default INP)"
         print *, " -o OUTFILE  : Write data to OUTFILE (default OUT)"
         print *, " -h          : Print this help message"
         print *, " "
      end subroutine print_help

end program filter_test
