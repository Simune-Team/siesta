! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_tbt_out

  use parallel, only : Node, Nodes, IONode

  implicit none 

  ! Module for writing out data to files...
  ! It handles all cases of data, and makes the overview in TBTrans a lot easier.
  private
  save
  
  character(len=200), parameter :: &
       defEIGfmt = '(i5,2(tr1,e16.8))', &
       defTEIGfmt = '(f10.5,20000(tr1,e16.8))', &
       defDOSfmt = '(f10.5,tr1,e16.8)', &
       defTfmt = '(f10.5,3(tr1,e16.8))', &
       defCOOPfmt = '(2(tr1,i4),tr1,f10.6,4(tr1,e16.8))', &
       defCOOPLRfmt = '(i4,tr1,f10.6,4(tr1,e16.8))', &
       defAtomPDOSTotfmt = '(i4,tr1,f10.6,3(tr1,e16.8))', &
       defAtomPDOSOrbfmt = '(i4,tr1,f10.6,100(tr1,e16.8))' ! If orbitals on a single atom exceeds 100 EDIT HERE

  public :: create_file, out_NEWLINE
  public :: out_kpt_header
  public :: out_EIG, out_TEIG
  public :: out_DOS
  public :: out_Trans
  public :: out_COOP, out_COOPLR
  public :: out_REGION, out_DEVICE
  public :: out_AtomPDOS_Tot, out_AtomPDOS_Orb

contains

  subroutine create_file(slabel,basename,ispin,nspin,funit)
    character(len=*), intent(in) :: slabel, basename
    integer, intent(in) :: ispin, nspin
    integer, intent(out) :: funit
    ! Local file name
    character(len=200) :: fname
    character(len=200), external :: paste
    fname = paste(slabel,"."//trim(basename))
    if ( nspin > 1 ) then
       if ( ispin == 1 ) fname = paste(slabel,"_UP."//trim(basename))
       if ( ispin == 2 ) fname = paste(slabel,"_DN."//trim(basename))
    end if
    if ( IONode ) then
       call io_assign(funit)
       open(funit,file=fname, status='unknown', form='formatted')
    end if
  end subroutine create_file

  subroutine out_EIG(funit,N,eig,fmt)
    use precision, only : dp
    use units, only : eV
    integer, intent(in) :: funit
    integer, intent(in) :: N
    real(dp), intent(in) :: eig(N)
    character(len=*), intent(in), optional :: fmt
    integer :: i 
    do i = 1 , N
       if ( present(fmt) ) then
          if ( IONode ) write(funit,fmt)       i,eig(i)/eV
       else
          if ( IONode ) write(funit,defEIGfmt) i,eig(i)/eV
       end if
    end do
  end subroutine out_EIG

  subroutine out_TEIG(funit,E,N,Teig,fmt)
    use precision, only : dp
    use units, only : eV
    integer, intent(in) :: funit
    real(dp), intent(in) :: E
    integer, intent(in) :: N
    real(dp), intent(in) :: Teig(N)
    character(len=*), intent(in), optional :: fmt
    integer :: i 
    if ( N <= 0 ) return
    if ( present(fmt) ) then
       if ( IONode ) write(funit,fmt)        E/eV,(Teig(i),i=1,N)
    else
       if ( IONode ) write(funit,defTEIGfmt) E/eV,(Teig(i),i=1,N)
    end if
  end subroutine out_TEIG


  subroutine out_DOS(funit,N,E,DOS,fmt)
    use precision, only : dp
    use units, only : eV
    integer, intent(in) :: funit
    integer, intent(in) :: N
    complex(dp), intent(in) :: E(N), DOS(N)
    character(len=*), intent(in), optional :: fmt
    ! Local variables
    integer :: i 
    do i = 1 , N
       if ( present(fmt) ) then
          if ( IONode ) write(funit,fmt) real(E(i),dp)/eV,real(DOS(i),dp)
       else
          if ( IONode ) write(funit,defDOSfmt) real(E(i),dp)/eV,real(DOS(i),dp)
       end if
    end do
  end subroutine out_DOS


  subroutine out_kpt_header(funit,ikpt,kpt,wkpt,ucell,fmt)
    use precision, only : dp
    use units,     only : eV
    integer, intent(in) :: funit, ikpt
    real(dp), intent(in) :: kpt(3), wkpt, ucell(3,3)
    character(len=*), intent(in), optional :: fmt
    real(dp) :: ktmp(3)
    if ( IONode ) write(funit,'(a,i5)') '# k-point: ',ikpt
    call kpoint_convert(ucell,kpt,ktmp,1)
    if ( present(fmt) ) then
       if ( IONode ) write(funit,fmt) &
            '# k  = ',kpt ,'w= ',wkpt
       if ( IONode ) write(funit,fmt) &
            '# kb = ',ktmp,'w= ',wkpt
    else
       if ( IONode ) write(funit,'(a6,3(f10.6,'', ''),a,f10.6)') &
            '# k  = ',kpt ,'w= ',wkpt
       if ( IONode ) write(funit,'(a6,3(f10.6,'', ''),a,f10.6)') &
            '# kb = ',ktmp,'w= ',wkpt
    end if
  end subroutine out_kpt_header

  subroutine out_Trans(funit,E,T,DOS,PDOS,fmt)
    use precision, only : dp
    use units, only : eV
    integer, intent(in)     :: funit
    real(dp), intent(in) :: E, DOS, PDOS
    real(dp), intent(in)    :: T
    character(len=*), intent(in), optional :: fmt
    if ( present(fmt) ) then
       if ( IONode ) write(funit,fmt) &
            E/eV,T,DOS*eV,PDOS*eV
    else
       if ( IONode ) write(funit,defTfmt) &
            E/eV,T,DOS*eV,PDOS*eV
    end if
  end subroutine out_Trans


  subroutine out_COOP(funit,ia1,ia2,E,cT,cL,cR,fmt)
    use precision, only : dp
    use units, only : eV
    integer, intent(in)  :: funit
    integer, intent(in)  :: ia1,ia2
    real(dp), intent(in) :: E, cT, cL, cR
    character(len=*), intent(in), optional :: fmt
    if ( present(fmt) ) then
       if ( IONode ) write(funit,fmt) &
            ia1,ia2,E/eV,cT*eV,cL*eV,cR*eV
    else
       if ( IONode ) write(funit,defCOOPfmt) &
            ia1,ia2,E/eV,cT*eV,cL*eV,cR*eV
    end if
  end subroutine out_COOP

  subroutine out_COOPLR(funit,ia1,E,cLL,cL,cLR,fmt)
    use precision, only : dp
    use units, only : eV
    integer, intent(in)  :: funit
    integer, intent(in)  :: ia1
    real(dp), intent(in) :: E, cLL, cL, cLR
    character(len=*), intent(in), optional :: fmt
    if ( present(fmt) ) then
       if ( IONode ) write(funit,fmt) &
            ia1,E/eV,cLL*eV,cL*eV,cLR*eV
    else
       if ( IONode ) write(funit,defCOOPLRfmt) &
            ia1,E/eV,cLL*eV,cL*eV,cLR*eV
    end if
  end subroutine out_COOPLR

  subroutine out_NEWLINE(funit)
    integer, intent(in) :: funit
    if ( IONode ) write(funit,*) ""
  end subroutine out_NEWLINE

  subroutine out_REGION(funit,na_u,xa,ia1,ia2,name,marker)
    use precision, only : dp
    use units, only : Ang
    integer, intent(in)  :: funit, na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: ia1, ia2
    character(len=*), intent(in) :: name
    character(len=1), intent(in) :: marker
    integer :: i, mid
    if ( ia2 - ia1 == 0 ) return
    mid = (ia2 - ia1) / 2
    if ( IONode ) then
       write(funit,'(a)') repeat(marker,46)
       do i = ia1, ia2
          if ( i == ia1 + mid ) then
             write(funit,'(a1,3(tr2,f12.7),tr2,a1,tr5,a)') &
                  marker,xa(:,i)/Ang,marker,trim(name)
          else
             write(funit,'(a1,3(tr2,f12.7),tr2,a1)') &
                  marker,xa(:,i)/Ang,marker
          end if
       end do
       write(funit,'(a)') repeat(marker,46)
    end if
  end subroutine out_REGION

  subroutine out_DEVICE(funit,na_u,xa,ia1,ia2,isoa1,isoa2,marker)
    use precision, only : dp
    use units, only : Ang
    integer, intent(in)  :: funit, na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: ia1, ia2, isoa1, isoa2
    character(len=1), intent(in) :: marker
    integer :: i, mid
    mid = (ia2 - ia1) / 2
    if ( IONode ) then
       do i = ia1, ia2
          if ( i == ia1 + mid ) then
             if ( isoa1 <= i .and. i <= isoa2 ) then
                write(funit,'(tr1,3(tr2,f12.7),tr3,a1,tr4,a)') &
                     xa(:,i)/Ang,marker,'Device'
             else
                write(funit,'(tr1,3(tr2,f12.7),tr8,a)') &
                     xa(:,i)/Ang,'Device'
             end if
          else
             if ( isoa1 <= i .and. i <= isoa2 ) then
                write(funit,'(tr1,3(tr2,f12.7),tr3,a1)') &
                     xa(:,i)/Ang,marker
             else
                write(funit,'(tr1,3(tr2,f12.7))') &
                     xa(:,i)/Ang
             end if
          end if
       end do
    end if
  end subroutine out_DEVICE

  subroutine out_AtomPDOS_Tot(funit,ia,E,wE,Tot,Left,Right,fmt)
    use precision, only : dp
    use units, only : eV
    integer, intent(in)  :: funit
    integer, intent(in)  :: ia
    real(dp), intent(in) :: E, wE, Tot, Left, Right
    character(len=*), intent(in), optional :: fmt
    if ( wE == 0.0_dp ) return
    if ( present(fmt) ) then
       if ( IONode ) write(funit,fmt) &
            ia,E/eV,Tot*eV,Left*eV,Right*eV
    else
       if ( IONode ) write(funit,defAtomPDOSTOTfmt) &
            ia,E/eV,Tot*eV,Left*eV,Right*eV
    end if
  end subroutine out_AtomPDOS_Tot
  
  subroutine out_AtomPDOS_Orb(funit,ia,E,wE,Orb,no,fmt)
    use precision, only : dp
    use units, only : eV
    integer, intent(in)  :: funit
    integer, intent(in)  :: ia, no
    real(dp), intent(in) :: E, wE, Orb(no)
    character(len=*), intent(in), optional :: fmt
    integer :: i
    if ( wE == 0.0_dp ) return
    if ( present(fmt) ) then
       if ( IONode ) write(funit,fmt) &
            ia,E/eV,(Orb(i)*eV,i=1,no)
    else
       if ( IONode ) write(funit,defAtomPDOSOrbfmt) &
            ia,E/eV,(Orb(i)*eV,i=1,no)
    end if
  end subroutine out_AtomPDOS_Orb

end module m_tbt_out
