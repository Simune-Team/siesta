subroutine ts_show_regions(ucell,na_u,xa,N_Elec,Elecs)
  use precision, only : dp
  use units, only : Ang
  use parallel, only : IONode
  use fdf, only : fdf_get

  use m_ts_electype
  use m_ts_method, only : na_Buf, atom_type
  use m_ts_method, only : TYP_BUFFER, TYP_DEVICE
  use m_region
  
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

  if ( .not. IONode ) return

  if ( .not. fdf_get('TS.Atoms.Print',.false.) ) then
     write(*,'(/,a)') 'transiesta: Regions of atoms:'
     ia = 1
     i = 1
     do while ( ia < na_u )
        do while ( atom_type(ia) == atom_type(i) )
           ia = ia + 1
           if ( ia > na_u ) exit
        end do
        call print_rgn(i,ia-1)
        i = ia
     end do
     if ( ia == na_u ) call print_rgn(na_u,na_u)
     write(*,*) ! new-line
     return
  end if

  ! mid-point of device
  ia_mid = (na_u - na_Buf-sum(TotUsedAtoms(Elecs))+1) / 2

  write(*,'(/,a)') 'transiesta: Atomic coordinates and regions (Ang):'
  ia = 1
  do while ( ia <= na_u )
     select case ( atom_type(ia) )
     case ( TYP_BUFFER )
        i = ia
        do while ( i <= na_u )
           ! this assures that we do not have a memory leak
           if ( atom_type(i) /= TYP_BUFFER ) exit
           i = i + 1
        end do
        i = i - ia
        ! steps ia counter
        call out_REGION(ia,i,'Buffer','/')
     case ( TYP_DEVICE )
        ia_mid = ia_mid - 1
        if ( ia_mid == 0 ) then
           write(*,'(tr1,i5,tr2,3(tr2,f12.7),tr8,a)') ia,xa(:,ia)/Ang,'Device'
        else
           write(*,'(tr1,i5,tr2,3(tr2,f12.7))') ia,xa(:,ia)/Ang
        end if
        ia = ia + 1
     case default
        ! electrode position
        do i = 1 , N_Elec
           if ( ia == Elecs(i)%idx_a ) then
              call out_REGION(ia,TotUsedAtoms(Elecs(i)), &
                   trim(name(Elecs(i)))//' electrode','#')
              exit
           end if
        end do
     end select
  end do

  write(*,*)

contains 

  subroutine out_REGION(ia,NA,name,marker,first,last)
    integer, intent(inout) :: ia
    integer, intent(in)    :: NA
    character(len=*), intent(in) :: name
    character(len=1), intent(in) :: marker
    character(len=1), intent(in), optional :: first, last
    integer :: i, mid
    logical :: lfirst, llast
    lfirst = present(first)
    llast = present(last)

    if ( NA < 1 ) return
    mid = (NA+1) / 2
    write(*,'(tr7,a)') repeat(marker,46)
    do i = 1, NA
       if ( i == 1 .and. lfirst ) then
          write(*,'(tr1,i5,tr1,a1,3(tr2,f12.7),tr2,a1,tr1,a1)') &
               ia,marker,xa(:,ia)/Ang,marker,first
       else if ( i == mid ) then
          write(*,'(tr1,i5,tr1,a1,3(tr2,f12.7),tr2,a1,tr5,a)') &
               ia,marker,xa(:,ia)/Ang,marker,trim(name)
       else if ( i == NA .and. llast ) then
          write(*,'(tr1,i5,tr1,a1,3(tr2,f12.7),tr2,a1,tr1,a1)') &
               ia,marker,xa(:,ia)/Ang,marker,last
       else
          write(*,'(tr1,i5,tr1,a1,3(tr2,f12.7),tr2,a1)') &
               ia,marker,xa(:,ia)/Ang,marker
       end if
       ia = ia + 1
    end do
    write(*,'(tr7,a)') repeat(marker,46)

  end subroutine out_REGION

  subroutine print_rgn(ia1,ia2)
    integer, intent(in) :: ia1, ia2
    type(tRgn) :: rgn
    character(len=2) :: prefix_rgn

    call rgn_range(rgn,ia1,ia2)
    select case ( atom_type(ia1) )
    case ( TYP_BUFFER )
       rgn%name = 'Buffer'
       prefix_rgn = '//'
    case ( TYP_DEVICE ) 
       rgn%name = 'Device'
       prefix_rgn = '--'
    case default
       rgn%name = 'Elec.'//trim(Elecs(atom_type(ia1))%name)
       prefix_rgn = '##'
    end select
    
    call rgn_print(rgn,name=prefix_rgn,seq_max=12,indent = 3)

    call rgn_delete(rgn)
    
   end subroutine print_rgn

end subroutine ts_show_regions


