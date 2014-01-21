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
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com


module m_ts_chem_pot

  use precision, only : dp

  use m_ts_io_ctype, only : C_N_NAME_LEN
  implicit none

  real(dp), public, parameter :: mu_same = 1.e-8_dp
  integer,  public, parameter :: NAME_MU_LEN = 20

  type :: ts_mu
     ! name of the chemical potential
     character(len=NAME_MU_LEN) :: name = ' '
     ! ID
     integer  :: ID = 0
     ! the chemical potential
     real(dp) :: mu = 0._dp
     ! number of electrodes having this chemical potential
     integer  :: N_El = 0
     ! array of electrode indices (conforming with the Elecs-array)
     integer, pointer :: el(:) => null()
     ! We must have a container which determines the contour segments
     ! that gets attached to the chemical potential
     character(len=C_N_NAME_LEN), allocatable :: Eq_seg(:)
  end type ts_mu
  public :: ts_mu

  interface hasC
     module procedure hasCio
     module procedure hasCeq
  end interface hasC
  public :: hasC

  interface hasEl
     module procedure hasEl_i
  end interface hasEl
  public :: hasEl

  interface name
     module procedure name_
  end interface name
  public :: name

  public :: Eq_segs

  public :: chem_pot_add_Elec

  public :: fdf_nmu, fdf_mu
  private

contains

  function fdf_nmu(prefix,this_n) result(n)
    use fdf

    character(len=*), intent(in) :: prefix
    type(ts_mu), allocatable :: this_n(:)
    integer :: n

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    integer :: i
    
    logical :: found

    n = 0
    found = fdf_block(trim(prefix)//'.ChemPots',bfdf)
    if ( .not. found ) return

    ! first count the number of electrodes
    n = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       n = n + 1 
    end do

    allocate(this_n(n))

    ! rewind to read again
    call fdf_brewind(bfdf)

    n = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       n = n + 1 
       this_n(n)%Name = trim(fdf_bnames(pline,1))
       if ( n > 1 ) then
          ! Check that no name is the same
          do i = 1 , n - 1 
             if ( leqi(name(this_n(i)),name(this_n(n))) ) then
                call die('Chemical potential names must not be the same')
             end if
          end do
       end if
    end do

  end function fdf_nmu

  function fdf_mu(prefix,slabel,this) result(found)
    use fdf

    character(len=*), intent(in) :: prefix,slabel
    type(ts_mu), intent(inout) :: this
    logical :: found

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    logical :: info(2)

    character(len=200) :: ln

    found = fdf_block(trim(prefix)//'.ChemPot.'//trim(Name(this)),bfdf)
    if ( .not. found ) return

    info(:) = .false.

    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       
       ln = trim(fdf_bnames(pline,1))
       
       ! We select the input
       if ( leqi(ln,'chemical-shift') .or. &
            leqi(ln,'mu') ) then
          if ( fdf_bnvalues(pline) < 1 ) call die('Chemical-shift not supplied')
          if ( fdf_bnnames(pline) < 2 ) call die('Unit of chemical-shift not supplied')
          this%mu = fdf_bvalues(pline,1) * fdf_convfac(fdf_bnames(pline,2),'Ry')
          info(1) = .true.

       else if ( leqi(ln,'contour.eq') ) then
          ! we automatically make room for one pole contour
          call read_contour_names('Equilibrium',this%Eq_seg,fakes=1)
          info(2) = .true.

       else
          
          call die('Unrecognized option "'//trim(ln)//'" &
               &for chemical potential: '//trim(this%name))

       end if

    end do
    
    if ( .not. all(info) ) then
       write(*,*)'You need to supply at least:'
       write(*,*)' - chemical-shift'
       write(*,*)' - contour.eq'
       call die('You have not supplied all chemical potential information')
    end if

  contains
    
    subroutine read_contour_names(name,con,fakes)
      character(len=*), intent(in) :: name
      character(len=C_N_NAME_LEN), allocatable :: con(:)
      integer, intent(in), optional :: fakes
      integer :: i, empty

      character(len=200) :: ln

      if ( allocated(con) ) call die("Contour already found.")

      ! we need to read in the equilibrium contour
      ! skip to "begin"
      if ( .not. fdf_bline(bfdf,pline) ) &
           call die("Chemical potential block ended prematurely.")

      ! read in the begin ... end block
      ln = fdf_bnames(pline,1)
      if ( .not. leqi(ln,"begin") ) &
           call die(trim(name)//" contour errorneously formatted. &
           &First line *must* be begin!")

      ! Count lines
      i = 0
      empty = 0
      do 
         if ( .not. fdf_bline(bfdf,pline) ) &
              call die("Chemical potential block ended prematurely.")
         if ( fdf_bnnames(pline) < 1 ) then
            empty = empty + 1
         else
            ln = fdf_bnames(pline,1)
            if ( leqi(ln,'end') ) exit
            i = i + 1
         end if
      end do
      
      ! allocate names
      if ( present(fakes) ) then
         allocate(con(i+fakes))
         empty = empty - fakes
      else
         allocate(con(i))
      end if
      con = ' '
      do i = 0 , size(con) + empty
         if ( .not. fdf_bbackspace(bfdf) ) &
              call die("Backspacing too much...")
      end do
      i = 0
      do 
         if ( .not. fdf_bline(bfdf,pline) ) &
              call die("Chemical potential block ended prematurely.")
         if ( fdf_bnnames(pline) < 1 ) cycle
         ln = fdf_bnames(pline,1)
         if ( leqi(ln,'end') ) exit
         i = i + 1
         if ( len_trim(ln) > C_N_NAME_LEN ) then
            call die('Contour name: '//trim(ln)//' is too long, please use a &
                 shorter name.')
         end if
         con(i) = trim(ln)
      end do
      
    end subroutine read_contour_names
    
  end function fdf_mu

  elemental function Eq_segs(this) result(count)
    type(ts_mu), intent(in) :: this
    integer :: count
    count = size(this%Eq_seg)
  end function Eq_segs

  subroutine chem_pot_add_Elec(this,iEl) 
    type(ts_mu), intent(inout) :: this
    integer, intent(in) :: iEl
    integer, pointer :: tmp(:), tmp2(:)

    nullify(tmp)
    allocate(tmp(this%N_El+1))

    if ( this%N_El == 0 ) then
       this%N_El = 1
       this%el => tmp
       this%el(1) = iEl
    else if ( all(this%el /= iEl) ) then
       ! copy over
       tmp(1:this%N_El) = this%el(:)
       this%N_El = this%N_El + 1
       tmp(this%N_El) = iEl
       ! clean up
       tmp2 => this%el
       this%el => tmp
       deallocate(tmp2)
    end if
  end subroutine chem_pot_add_Elec

  elemental function hasEl_i(this,iEl) result(has)
    type(ts_mu), intent(in) :: this
    integer, intent(in) :: iEl
    logical :: has
    has = this%N_El > 0
    if ( has ) has = any(this%el == iEl)
  end function hasEl_i

  elemental function hasCio(this,c_io) result(has)
    use m_ts_io_ctype, only : ts_c_io
    type(ts_mu), intent(in) :: this
    type(ts_c_io), intent(in) :: c_io 
    logical :: has
    integer :: i
    do i = 1 , Eq_segs(this)
       has = this%Eq_seg(i) .eq. c_io%name
       if ( has ) return
    end do
  end function hasCio

  elemental function hasCeq(this,c) result(has)
    use m_ts_cctype, only : ts_cw
    type(ts_mu), intent(in) :: this
    type(ts_cw), intent(in) :: c 
    logical :: has
    integer :: i
    do i = 1 , Eq_segs(this)
       has = this%Eq_seg(i) .eq. c%c_io%name
       if ( has ) return
    end do
  end function hasCeq

  elemental function name_(this) result(name)
    type(ts_mu), intent(in) :: this
    character(len=NAME_MU_LEN) :: name
    name = this%name
  end function name_

end module m_ts_chem_pot
  
