module m_tbt_out

  use parallel, only : Node, Nodes, IONode

  implicit none 

  ! Module for writing out data to files...
  ! It handles all cases of data, and makes the overview in TBTrans a lot easier.
  private
  save
  
  character(len=200) :: defEIGfmt = '(1i5,2e16.8)'
  character(len=200) :: defTEIGfmt = '(f10.5,20000e16.8)'
  character(len=200) :: defDOSfmt = '(3e16.8)'
  character(len=200) :: defTfmt = '(f10.5,3e16.8))'
  character(len=200) :: defCOOPfmt = '(2i4,5f10.6))'
  character(len=200) :: defCOOPLRfmt = '(1i4,5f10.6))'

  public :: create_file, out_NEWLINE
  public :: out_kpt_header
  public :: out_EIG, out_TEIG
  public :: out_DOS
  public :: out_Trans
  public :: out_COOP, out_COOPLR
  public :: out_REGION, out_DEVICE

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
       if ( ispin == 1 ) fname = paste(slabel,".UP."//trim(basename))
       if ( ispin == 2 ) fname = paste(slabel,".DN."//trim(basename))
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
    real(dp) :: totDOS
    integer :: i 
    totDOS = 0.0_dp
    do i = 1 , N
       totDOS = totDOS + real(DOS(i),dp)
       if ( present(fmt) ) then
          if ( IONode ) write(funit,fmt) real(E(i),dp)/eV,real(DOS(i),dp),totDOS
       else
          if ( IONode ) write(funit,defDOSfmt) real(E(i),dp)/eV,real(DOS(i),dp),totDOS
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
            real(E,dp)/eV,T,real(DOS,dp),real(PDOS,dp)
    else
       if ( IONode ) write(funit,defTfmt) &
            real(E,dp)/eV,T,real(DOS,dp),real(PDOS,dp)
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
            ia1,ia2,E/eV,cT,cL,cR
    else
       if ( IONode ) write(funit,defCOOPfmt) &
            ia1,ia2,E/eV,cT,cL,cR
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
            ia1,E/eV,cLL,cL,cLR
    else
       if ( IONode ) write(funit,defCOOPLRfmt) &
            ia1,E/eV,cLL,cL,cLR
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
                  marker,xa(:,i),marker,trim(name)
          else
             write(funit,'(a1,3(tr2,f12.7),tr2,a1)') &
                  marker,xa(:,i),marker
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
                write(funit,'(tr1,3(tr2,f12.7),tr3,a1,tr5,a)') &
                     xa(:,i),marker,'Device'
             else
                write(funit,'(tr1,3(tr2,f12.7),tr8,a)') &
                     xa(:,i),'Device'
             end if
          else
             if ( isoa1 <= i .and. i <= isoa2 ) then
                write(funit,'(tr1,3(tr2,f12.7),tr3,a1)') &
                     xa(:,i),marker
             else
                write(funit,'(tr1,3(tr2,f12.7))') &
                     xa(:,i)
             end if
          end if
       end do
    end if
  end subroutine out_DEVICE

end module m_tbt_out
