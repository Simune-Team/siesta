subroutine ts_show_regions(ucell,na_u,xa,naBufL, &
     nElecs,Elecs,naBufR)
  use m_ts_electype
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
  integer, intent(in) :: nElecs
  type(Elec), intent(in) :: Elecs(nElecs)
  integer, intent(in)  :: naBufL, naBufR

! ********************
! * LOCAL variables  *
! ********************
  integer :: ia, i, mid, ia_mid
  logical :: printed_elec

  if ( .not. IONode ) return

  ! Initialize ia counter
  ia = 0

  ia_mid = (na_u - naBufL-sum(TotUsedAtoms(Elecs))-naBufR+1) / 2

  write(*,'(/,a)') 'transiesta: Atomic coordinates and regions (Ang):'
  call out_REGION(ia,naBufL,'Left buffer','/')
  do while ( ia < na_u - naBufR )
     printed_elec = .false.
     do i = 1 , nElecs
        if ( ia + 1 == Elecs(i)%idx_na ) then
           call out_REGION(ia,TotUsedAtoms(Elecs(i)), &
                trim(name(Elecs(i)))//' electrode','#')
           printed_elec = .true.
        end if
     end do
     if ( .not. printed_elec ) then
        ia = ia + 1
        ia_mid = ia_mid - 1
        if ( ia_mid == 0 ) then
           write(*,'(tr1,3(tr2,f12.7),tr8,a)') xa(:,ia)/Ang,'Device'
        else
           write(*,'(tr1,3(tr2,f12.7))') xa(:,ia)/Ang
        end if
     end if
  end do

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


