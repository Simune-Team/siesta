! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine tt_eig(noD,tt,NEigCh,teig)
  use precision, only : dp
  use sys,       only : die

  implicit none

! ***********************
! * INPUT variables     *
! ***********************
  integer,     intent(in)    :: noD ! no. st. in contact region
  integer,     intent(in)    :: NEigch

! ***********************
! * OUTPUT variables    *
! ***********************
  complex(dp), intent(inout) :: tt(noD,noD) ! Eqv. to t^dagger.t
  real(dp),    intent(out)   :: TEig(NEigch)           ! channels

! ***********************
! * LOCAL variables     *
! ***********************
  integer                                   :: i
  integer                                   :: diagINFO     
  complex(dp), dimension (:),   allocatable :: worklap
  complex(dp), dimension (:,:), allocatable :: UR,UL
  complex(dp), dimension (:),   allocatable :: Eig
  real(dp),    dimension (:),   allocatable :: rwork
!=======================================================================
!     BEGIN *** Eigenchannel calculation:
! **  Diagonalizing t^dagger.t, eigenvectors in UR, eigenvalues in Eig
!=======================================================================

  call timer('tt_eig',1)

! We shift the eigenvalues to increase precision...
  do i = 1, noD
     tt(i,i) = tt(i,i) + dcmplx(1e-2_dp,1e-2_dp)
  end do

! define the arrays
  allocate(worklap(64*noD))
  call memory('A','Z',64*noD,'teigchan')
  allocate(rwork(4*noD))
  call memory('A','D',2*noD,'teigchan')
  allocate(UR(noD,noD))
  call memory('A','Z',noD*noD,'teigchan')
  allocate(UL(noD,noD))
  call memory('A','Z',noD*noD,'teigchan')
  allocate(Eig(noD))
  call memory('A','Z',noD,'teigchan')

! Calculate the eigenvalues
  call zgeev('N','N',noD,tt,noD,Eig, &
       UL,noD,UR,noD,worklap, &
       64*noD,rwork,diagINFO)

  if (diagINFO.lt.0) then
     diagINFO = -diagINFO
     write(*,*)'ERROR: DIAG. FAILED: element no. ',diagINFO, &
          ' had illegal value. '
     call die('teigchan: ERROR, DIAG. FAILED')
  end if
  if (diagINFO.gt.0) then
     write(*,*)'ERROR: DIAG. FAILED: only the ', diagINFO,'+1:', &
          noD,' elements converged.'
     call die('teigchan: ERROR, DIAG. FAILED')
  end if

  
! Reshift the eigenvalues...
  do i = 1 , noD
     Eig(i) = Eig(i) - dcmplx(1e-2_dp,1e-2_dp)
  end do

  if ( NEigch > noD ) then
     TEig = 0.0_dp
     call eigsort(eig,TEig,noD,noD)
  else if ( NEigch /= noD ) then
     call eigsort(eig,TEig,noD,NEigch)
  else
     ! Do not sort if NEigch == noD
     TEig(:) = REAL(eig(:),dp)
  end if

! make a memory free
  call memory('D','Z',size(worklap),'teigchan')
  deallocate(worklap)
  call memory('D','D',size(rwork),'teigchan')
  deallocate(rwork)
  call memory('D','Z',size(UR),'teigchan')
  deallocate(UR)
  call memory('D','Z',size(UL),'teigchan')
  deallocate(UL)
  call memory('D','Z',size(Eig),'teigchan')
  deallocate(Eig)

  call timer('tt_eig',2)

contains
  
  subroutine eigsort(eig, eigS, N, MaxN)
    use precision, only : dp

! ***********************
! * INPUT variables     *
! ***********************
    integer,    intent(in)  :: N, MaxN
    complex(dp),intent(inout)  :: eig(N)

! ***********************
! * OUTPUT variables    *
! ***********************
    real(dp),   intent(out) :: eigS(MaxN)

! ***********************
! * LOCAL variables     *
! ***********************
    complex(dp) :: ctmp
    integer     :: i,j

    do i=1,N-1
       do j=i+1,N
          if ((abs(eig(i))).lt.(abs(eig(j)))) then
             ctmp   = eig(j)
             eig(j) = eig(i)
             eig(i) = ctmp
          end if
       end do
    end do

    do i = 1 , MaxN
       eigS(i) = real(eig(i),dp)
    end do

  end subroutine eigsort

end subroutine tt_eig


