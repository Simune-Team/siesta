! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine ts_show_regions(ucell,na_u,xa,naBufL,naElecL,naElecR,naBufR,NA1L,NA2L,NA1R,NA2R)
  use precision, only : dp
  use units, only : Ang
  use parallel, only : IONode
  implicit none
! ********************
! * INPUT variables  *
! ********************
  real(dp), intent(in) :: ucell(3,3)
  integer, intent(in)  :: na_u
  real(dp), intent(in) :: xa(3,na_u)
  integer, intent(in)  :: naBufL, naElecL, naElecR, naBufR
  integer, intent(in)  :: NA1L, NA2L, NA1R, NA2R

! ********************
! * LOCAL variables  *
! ********************
  integer :: ia, i, curR, mid
  character(len=4) :: marker

  if ( .not. IONode ) return

  ! Initialize ia counter
  ia = 0

  write(*,'(/,a)') 'transiesta: Atomic coordinates and regions (Ang):'
  call out_REGION(ia,naBufL,'Left buffer','/')
  call out_REGION(ia,naElecL*NA1L*NA2L,'Left electrode','#')

  mid = (na_u - naBufL-naElecL*NA1L*NA2L-naElecR*NA1R*NA2R-naBufR+1) / 2
  do i = 1 , na_u - naBufL-naElecL*NA1L*NA2L-naElecR*NA1R*NA2R-naBufR
     ia = ia + 1
     if ( i == mid ) then
        write(*,'(tr1,3(tr2,f12.7),tr8,a)') &
             xa(:,ia)/Ang,'Device'
     else
        write(*,'(tr1,3(tr2,f12.7))') &
             xa(:,ia)/Ang
     end if
  end do

  call out_REGION(ia,naElecR*NA1R*NA2R,'Right electrode','#')
  call out_REGION(ia,naBufR,'Right buffer','/')

  write(*,*)

contains 

  subroutine out_REGION(ia,NA,name,marker)
    integer, intent(inout) :: ia
    integer, intent(in)    :: NA
    character(len=*), intent(in) :: name
    character(len=1), intent(in) :: marker
    integer :: i, mid
    
    if ( NA < 1 ) return
    mid = (NA+1) / 2
    write(*,'(a)') repeat(marker,46)
    do i = 1, NA
       ia = ia + 1
       if ( i == mid ) then
          write(*,'(a1,3(tr2,f12.7),tr2,a1,tr5,a)') &
               marker,xa(:,ia)/Ang,marker,trim(name)
       else
          write(*,'(a1,3(tr2,f12.7),tr2,a1)') &
               marker,xa(:,ia)/Ang,marker
       end if
    end do
    write(*,'(a)') repeat(marker,46)
  end subroutine out_REGION

end subroutine ts_show_regions


