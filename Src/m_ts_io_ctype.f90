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
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!

module m_ts_io_ctype
!
! Routines that are used to read in and print out the contour for integration of the GFs
! 
  use precision, only : dp

  implicit none

  integer, parameter :: C_N_NAME_LEN = 10

  integer, parameter :: c_N = 200
  type :: ts_c_io
     ! the name of the integration segment
     character(len=C_N) :: name = ' '
     ! The integral bounds
     real(dp) :: a, b
     ! (text-form)
     character(len=c_N) :: ca = ' ', cb = ' '
     ! The division of the grid
     real(dp) :: d
     ! (text-form)
     character(len=c_N) :: cd = ' '
     ! Number of points in the grid
     integer :: N
     ! (text-form)
     character(len=c_N) :: cN = ' '
     ! the integration method (g-legendre etc.)
     character(len=c_N) :: method = ' '
     ! the type of integration (eq|neq|t)
     character(len=4) :: type = ' '
     ! the integration part (circle|line|tail|pole)
     character(len=c_N) :: part = ' '
     ! the options attached to the integration method
     ! (this will take the form of a linked list)
     type(ts_c_opt_ll), pointer :: opt => null()
  end type ts_c_io

  type ts_c_opt_ll
     character(len=c_N) :: opt = ' '
     type(ts_c_opt_ll), pointer :: next => null()
  end type ts_c_opt_ll

  interface copy
     module procedure copy_
  end interface copy
  private :: copy_

contains

  subroutine copy_(from,to)
    type(ts_c_io), intent(in) :: from
    type(ts_c_io), intent(out) :: to
    to%name = from%name
    to%a = from%a
    to%b = from%b
    to%ca = from%ca
    to%cb = from%cb
    to%d = from%d
    to%cd = from%cd
    to%N = from%N
    to%cN = from%cN
    to%method = from%method
    to%type = from%type
    to%part = from%part
    to%opt => from%opt
  end subroutine copy_

  function fdf_nc_iotype(prefix,suffix) result(n)
    use fdf

    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: suffix
    integer :: n

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    integer :: i
    
    logical :: found

    n = 0
    found = fdf_block(trim(prefix)//'.Contours.'//trim(suffix),bfdf)
    if ( .not. found ) return

    ! first count the number of electrodes
    n = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       n = n + 1 
    end do

  end function fdf_nc_iotype

  function fdf_name_c_iotype(prefix,suffix,i) result(name)
    use fdf
    character(len=*), intent(in) :: prefix, suffix
    integer, intent(in) :: i
    character(len=C_N_NAME_LEN) :: name

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    integer :: n
    logical :: found

    name = ' '

    n = 0
    found = fdf_block(trim(prefix)//'.Contours.'//trim(suffix),bfdf)
    if ( .not. found ) return

    ! first count the number of electrodes
    n = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       n = n + 1 
       if ( n == i ) then
          name = fdf_bnames(pline,1)
          return
       end if
    end do

  end function fdf_name_c_iotype

  subroutine ts_read_contour_block(prefix,suffix,bName,c, kT, V) 

    use parallel, only : IONode
    use units, only : eV
    use fdf
    use parse, only : search_fun, characters

    character(len=*), intent(in) :: prefix, suffix
    character(len=C_N_NAME_LEN), intent(in) :: bName
    type(ts_c_io), intent(inout) :: c
    real(dp), intent(in) :: kT
    real(dp), intent(in), optional :: V

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()

    type(ts_c_opt_ll), pointer :: opt, opt_tmp

    character(len=200) :: g
    integer :: pNames, my_loop, cI, iV
    integer :: iS, iE
    logical :: read_method, read_points, read_int, error
    
    ! if the block does not exist, return
    if ( len_trim(suffix) ) then
       g = trim(prefix)//'.Contour.'//trim(suffix)
    else
       g = trim(prefix)//'.Contour'
    end if
    if ( .not. fdf_block(trim(g)//'.'//trim(bName),bfdf) ) return

    ! the contour has already been assigned values
    if ( len_trim(c%name) /= 0 ) return

    ! Read the information in block
    c%name = trim(bName)
    
    ! We must ensure the block be organized as this:
    !  part circle|line|tail
    !  from <a> to <b> 
    !    points <N> / separation <d>
    !      method <method>
    !       opt <opt-1>
    !       opt <opt-2>
    !       opt ...
    !       opt <opt-N>

    ! { "part <part>" circle|line|tail|pole
    if ( .not. move2names() ) then
       call die('Block: '//trim(bName)//'. &
            &Could not find part segment in contour.')
    end if

    ! get the 'part'
    iS = search_fun('part',pline)
    if ( iS < 0 ) iS = search_fun('p',pline)
    if ( iS < 0 ) then
       call die('Block: '//trim(bName)//' is not build correctly. &
            &part <part> line cannot find "part".')
    end if
    if ( fdf_bnnames(pline) < 2 ) then
       call die('Block: '//trim(bName)//' has not described the part properly. &
            &Must have part <part>')
    end if
    c%part = fdf_bnames(pline,2)
    if ( leqi(c%part,'circle') ) then
       c%part = 'circle' ! ensures capitalization!!!! DON'T EDIT
    else if ( leqi(c%part,'line') ) then
       c%part = 'line'
    else if ( leqi(c%part,'tail') ) then
       c%part = 'tail'
    else
       call die('Part of the contour could not be recognized as circle|line|tail')
    end if
   

    ! } "part

    ! { "from <a> to <b>"
    if ( .not. move2names() ) then
       call die('Block: '//trim(bName)//'. &
            &Could not find from <a> to <b> segment in contour')
    end if

    ! get the 'from'
    iS = search_fun('from',pline)
    if ( iS <= 0 ) then
       call die('Block: '//trim(bName)//' is not build correctly. &
            &from <a> to <b> line cannot find "from".')
    end if
    
    iE = search_fun('to',pline)
    if ( iE <= 0 ) then
       call die('Block: '//trim(bName)//' is not build correctly. &
            &from <a> to <b> line cannot find "to".')
    end if
    ! get the line as input
    call pline_E_parse(pline,iS,c%ca,val=c%a,V=V,kT=kT, &
         before=iE)
    if ( leqi(c%ca,'next') ) call die('Block: '//trim(bName)//' can not &
         &have a==next')

    ! get <b>
    iS = iE
    call pline_E_parse(pline,iS,c%cb,val=c%b,V=V,kT=kT)
    if ( leqi(c%cb,'previous') .or. & 
         leqi(c%cb,'prev') ) call die('Block: '//trim(bName)//' can not &
         &have b==previous')

    ! } "from"

    ! { "points <N>"
    if ( .not. move2names() ) then
       call die('Block: '//trim(bName)//'. &
            &Could not find points segment in contour')
    end if
    
    ! we now read the points or separation
    iS = search_fun('points',pline)
    if ( iS < 0 ) iS = search_fun('p',pline)
    iE = search_fun('delta',pline)
    if ( iE < 0 ) iE = search_fun('sep',pline)
    if ( iS < 0 .and. iE < 0 ) then
       call die('Block: '//trim(bName)//' is not build correctly. &
            &Could not decipher points/delta/separation')
    end if
    
    ! if we have points we simply read in the number
    if ( 0 <= iS ) then
       c%N  = fdf_bintegers(pline,1,after=iS) ! first integer
       c%cN = characters(pline,1,-1,after=iS)
    else ! we have a delta designation
       ! notice that we can actually use kT here
       call pline_E_parse(pline,iE,c%cd,val=c%d,kT=kT)
       if ( c%d <= 0._dp ) then
          call die('Block: '//trim(bName)//' is not valid. &
               &The dE designator is negative or zero.')
       end if
    end if

    ! } "points"

    ! { "method <method>"
    if ( .not. move2names() ) then
       call die('Block: '//trim(bName)//'. &
            &Could not find method <method> segment in contour')
    end if
    if ( fdf_bnnames(pline) < 2 ) then
       call die('Block: '//trim(bName)//' has not described the method properly. &
            &Must have method <method>')
    end if

    ! the method should be a one-name thing
    c%method = fdf_bnames(pline,2)

    ! } "method"

    ! { "opt <option>"
    do 
       ! if we don't find anything simply exit the optional reading
       if ( .not. move2names() ) exit
       
       ! we can allocate the opt
       nullify(opt)
       allocate(opt)
       
       if ( .not. associated(c%opt) ) then
          c%opt => opt
       else
          ! find the last element
          opt_tmp => c%opt
          do while ( associated(opt_tmp%next) ) 
             opt_tmp => opt_tmp%next
          end do
          opt_tmp%next => opt
       end if

       ! populate the option
       opt%opt = trim(characters(pline,1,-1,after=1))

    end do
    ! } "opt"

  contains 

    function move2names() result(found)
      logical :: found
      found = .false.
      do while ( fdf_bline(bfdf,pline) ) 
       
         ! if no names exist we can loop to the next line
         ! thus commenting out a contour works very good!
         if ( fdf_bnnames(pline) == 0 ) then
            cycle
         else
            found = .true.
            exit
         end if
      end do
      
    end function move2names

  end subroutine ts_read_contour_block

  function ts_c_bisphysical(pline,ind)
    use parse
    implicit none
!------------------------------------------------- Input Variables
    type(parsed_line), pointer        :: pline
    integer(ip), intent(in)           :: ind

!------------------------------------------------ Output Variables
    logical                           :: ts_c_bisphysical
    
    ts_c_bisphysical = match(pline,'vn',after=ind)

  end function ts_c_bisphysical
    
  function ts_c_bphysical(pline,ind,defunit)
    use fdf
    use parse, only : die
    
    implicit none
!------------------------------------------------- Input Variables
    type(parsed_line), pointer        :: pline
    integer,     intent(in)           :: ind
    character(len=*), intent(in)      :: defunit

!------------------------------------------------ Output Variables
    real(dp)                          :: ts_c_bphysical
    character(len=10)                 :: unitstr

    if (.not. ts_c_bisphysical(pline,ind) ) then
       call die('PARSE module: ts_c_bphysical', 'Not enough values and names in line',   &
            'm_ts_io_contour', -1)
    endif
    
    ! Label with value
    ts_c_bphysical = fdf_bvalues(pline, 1, ind)
    
    unitstr = fdf_bnames(pline, 1, ind+1)
    if (.not. leqi(unitstr, defunit)) &
         ts_c_bphysical = ts_c_bphysical * fdf_convfac(unitstr, defunit)
    
  end function ts_c_bphysical

  subroutine pline_E_parse(pline,after,c, &
       val,V,kT, &
       before)
    use fdf
    use parse, only : characters, ntokens
    use units, only : eV
    
    ! This routine will parse a pline from an index and collect the values from units
    ! into one value and return the value and the string which can build it again
    type(parsed_line), pointer :: pline
    integer, intent(in) :: after ! search physical values from after this index
    character(len=c_N), intent(out) :: c ! the string that can construct this physical value in the shortest sense
    ! value of the energy
    real(dp), intent(out), optional :: val
    real(dp), intent(in),  optional :: V, kT
    integer, intent(in), optional   :: before

    ! Local parameters
    real(dp) :: ret, tmp
    character(len=50) :: g, tmp_g
    logical :: add, get_val, has_V, has_kT
    integer :: i, j, stat, offset
    
    ! initialize the output string
    c = ' '
    ! whether we should collect the value
    get_val = present(val)
    if ( get_val ) then
       val = 0._dp  ! initialize
       has_V = present(V)
       has_kT = present(kT)
    end if

    i = after - 1
!    if ( present(before) ) then
!       write(*,*)'CHARS: ',trim(characters(pline,after+1,before-1))
!    else
!       write(*,*)'CHARS: ',trim(characters(pline,1,-1,after=after))
!    end if
    add = .true.

    do while ( ntokens(pline,after=i) > 1 ) 
       i = i + 1
       if ( present(before) ) then
          !print *,'Before: ',i,before
          if ( i >= before ) return
       end if
         
       ! first we need to figure out if we are dealing with 
       ! a single plus or negative sign...
       g = fdf_bnames(pline,1,after=i)
!       print *,'CHECK: ',trim(g),' ',trim(characters(pline,1,1,after=i))
       offset = 2
       if ( trim(characters(pline,1,1,after=i)) == trim(g) ) then
          offset = 1 ! for units the unit is at the index
       else if ( ntokens(pline,after=i) > 0 ) then
!          print *,'CHECK2: ',trim(g),' ',trim(characters(pline,1,1,after=i+1))
          if ( trim(characters(pline,1,1,after=i+1)) == trim(g) ) then
             offset = 0
          end if
       end if
!print *,'PE: ',trim(g),' offset:',offset
       if ( leqi(g,'+') ) then
          add = .true.
          c = trim(c)//' +'
          cycle
          
       else if ( leqi(g,'-') ) then
          add = .false.
          c = trim(c)//' -'
          cycle
          
       end if

       ! else we need to figure out if it is something which is not a specific unit...
         
       ! Parse compressed elements
       if ( leqi(g(1:1),'+') ) then
          add = .true.
          c = trim(c)//' +'
          g = g(2:)

       else if ( leqi(g(1:1),'-') ) then
          add = .false.
          c = trim(c)//' -'
          g = g(2:)

       end if

       ! we need to be sure to capture everything...
       if ( leqi(g,'prev') .or. &
            leqi(g,'previous') .or. &
            leqi(g,'next') ) then
          c = trim(g)
          return ! RETURN

       else if ( leqi(g,'inf') ) then
          ! note: we have already added sign
          c = trim(c)//' inf'

          if ( get_val ) then
             if ( add ) then
                val = huge(1._dp)
             else
                val = huge(-1._dp)
             end if
          end if

          ! you can't add/subtract anything meaningful to
          ! "inf" without getting inf again...
          return

       else if ( leqi(g(1:1),'v') ) then
          ! we still need to make sure that we can interpret the
          ! bias-fraction
          if ( leqi(g(2:2),'/') ) then
             ! read in the bias fraction
             read(g(3:),'(i9)',iostat=stat) j
             if ( stat /= 0 ) then
                call die('Fractional parameter chosen cannot be distinguished: '//trim(g))
             end if
             write(g,'(a,i0)') 'V/',j
          else
             j = 1
             g = 'V'
          end if

          ! get the V designation
          c = trim(c)//' '//trim(g)
          
          if ( get_val ) then
             if ( .not. has_V ) call die('You cannot request &
                  &V value in this segment')
             
             tmp = abs(V) / real(j,dp)
             !print*,'Found : '//trim(g)//' at value: ',tmp/eV

             if ( add ) then
                val = val + tmp
             else
                val = val - tmp
             end if

          end if

       else

          if ( offset > 1 ) then
             return ! we have too many things unrecognizable
          end if

          ! We are dealing with a proper unit...
          ! check that we have a physical quantity that can be described...
          if ( .not. ts_c_bisphysical(pline,i-offset) ) then
             !print *,'NOT PHY: ',characters(pline,1,2,after=i-offset)
             return ! if we can't find
          end if
          
          if ( leqi(g,'kt') .or. leqi(g,'kbt') ) then
             ! Read in the former value
             ! We should not remove pm as that will not be read in
             g = remove_pm(characters(pline,1,1,after=i-offset))
             c = trim(c)//' '//trim(g)//' kT'
             
             if ( get_val .and. .not. has_kT ) call die('You cannot request &
                  &kT value in this segment')
             ! get the value...
             tmp = fdf_bvalues(pline,1,after=i-offset) * kT

          else
             ! this will make the code break if the unit is wrongly assigned...
             tmp = ts_c_bphysical(pline,i-offset,'Ry')

             c = trim(c)//' '//trim(characters(pline,1,1,after=i-offset))//' '&
                  //trim(characters(pline,2,2,after=i-offset))

          end if

          if ( get_val ) then
             if ( add ) then
                val = val + tmp
             else
                if ( tmp > 0._dp ) then
                   val = val - tmp
                else
                   val = val + tmp
                end if
             end if
          end if


          ! skip both the value and unit...
          if ( offset == 0 ) then
             i = i + 1
          end if

       end if

       ! we must have the resetting here, otherwise single +- statements
       ! will not be caught
       add = .true.
    end do

  contains
    
    function remove_pm(s)
      character(len=*), intent(in) :: s
      character(len=len_trim(s)) :: remove_pm
      if ( leqi(s(1:1),'+') .or. &
           leqi(s(1:1),'-') ) then
         remove_pm = trim(s(2:))
      else
         remove_pm = trim(s)
      end if
    end function remove_pm

  end subroutine pline_E_parse

  subroutine ts_print_contour_block(prefix,c)

    use parallel, only : IONode

    character(len=*), intent(in) :: prefix
    type(ts_c_io), intent(in) :: c
    type(ts_c_opt_ll), pointer :: opt

    if ( .not. IONode ) return

    ! Start by writing out the block beginning
    write(*,'(a,a)') '%block ',trim(prefix)//trim(c%name)

    write(*,'(2a)') 'part ',trim(c%part)

    ! move to third column...
    write(*,'(t3)',advance='no')
      
    write(*,'(4a)') 'from ',trim(c%ca),' to ',trim(c%cb)

    if ( len_trim(c%cd) /= 0 ) then
       ! we have delta designation
       write(*,'(t7,a,tr1,a)') 'delta', trim(c%cd)

    else       
       ! Print the number of points...
       write(*,'(t7,a,tr1,i0)') 'points', c%N
      
    end if

    ! Print the method
    write(*,'(t9,a,tr1,a)') 'method', trim(c%method)

    opt => c%opt
    do while ( associated(opt) )
       write(*,'(t10,a,tr1,a)') 'opt', trim(opt%opt)
       opt => opt%next
    end do
    
    write(*,'(a,a)') '%endblock ',trim(prefix)//trim(c%name)
    
  end subroutine ts_print_contour_block
  
end module m_ts_io_ctype
