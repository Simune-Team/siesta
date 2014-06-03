subroutine ts_show_regions(ucell,na_u,xa,N_Elec,Elecs)
  use m_ts_electype
  use m_ts_method, only : na_Buf, TYP_BUFFER, atom_type
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
  integer, intent(in) :: N_Elec
  type(Elec), intent(in) :: Elecs(N_Elec)

! ********************
! * LOCAL variables  *
! ********************
  integer :: ia, i, ia_mid
  logical :: printed_elec

  if ( .not. IONode ) return

  ! mid-point of device
  ia_mid = (na_u - na_Buf-sum(TotUsedAtoms(Elecs))+1) / 2

  write(*,'(/,a)') 'transiesta: Atomic coordinates and regions (Ang):'
  ia = 0
  do while ( ia < na_u )
     if ( atom_type(ia+1) == TYP_BUFFER ) then
        i = ia + 1
        do while ( atom_type(i) == TYP_BUFFER )
           i = i + 1
        end do
        call out_REGION(ia,i-ia-1,'Buffer','/')
     end if
     printed_elec = .false.
     do i = 1 , N_Elec
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


