! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
module version_info

implicit none

integer, dimension(3), save  :: num_version = (/0,0,0/)
character(len=*), parameter :: version_str =  &
"TBTRANS_VERSION"
character(len=*), parameter :: siesta_arch= &
"SIESTA_ARCH"
character(len=*), parameter :: fflags= &
"FFLAGS"
character(len=*), parameter :: fppflags= &
"FPPFLAGS"

private
public :: num_version, version_str
public :: siesta_arch, fflags, fppflags

end module version_info
!================================================================

subroutine prversion

use version_info
implicit none

#ifdef TBT_PHONON
write(6,'(2a)') "PHtrans Version: ", trim(version_str)
#else
write(6,'(2a)') "TBtrans Version: ", trim(version_str)
#endif
write(6,'(2a)') 'Architecture  : ', trim(siesta_arch)
write(6,'(2a)') 'Compiler flags: ', trim(fflags)
write(6,'(2a)') 'PP flags      : ', trim(fppflags)

#ifdef MPI
write(6,'(a)') 'PARALLEL version'
#else
write(6,'(a)') 'SERIAL version'
#endif

!$OMP parallel
!$OMP master
!$    write(*,'(a)') 'THREADED version'
!$OMP end master
!$OMP end parallel

#ifdef USE_GEMM3M
write(6,'(a)') 'GEMM3M support'
#endif
#ifdef CDF
write(6,'(a)') 'NetCDF support'
#endif

end subroutine prversion
!----------------------------------------------------------

subroutine get_version(v)
  use version_info
  implicit none
  integer, intent(out)  :: v(3)
  v = num_version
end subroutine get_version

