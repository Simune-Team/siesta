
module fdf_extra

  ! Module which extends the fdf routines
  ! by a little...
  use fdf
  use m_region
  
  implicit none
  
  private
  
  public :: fdf_bnext
  public :: fdf_brange
  
contains

  function fdf_bnext(bfdf,pline) result(has)
    type(block_fdf), intent(inout) :: bfdf
    type(parsed_line), pointer :: pline
    logical :: has
    do 
       has = fdf_bline(bfdf,pline)
       if ( .not. has ) return
       if ( fdf_bntokens(pline) > 0 ) return
    end do
  end function fdf_bnext

  subroutine fdf_brange(pline,r,low,high)
    type(parsed_line), pointer :: pline
    type(tRgn), intent(out) :: r
    ! Practically min/max for the region
    ! in question, however, it allows for 
    ! specifying negative numbers.
    integer, intent(in) :: low, high

    integer :: j, n, i, i1, i2, step
    integer, allocatable :: list(:)
    character(len=20) :: g

    ! We allocate a temporary list
    allocate(list(high-low+1))

    n = 0
    if ( fdf_bnnames(pline) == 1 ) then
       ! The input line is a simple list of integers
       do j = 1 , fdf_bnintegers(pline)
          i = fdf_bintegers(pline,j)
          n = n + 1 
          list(n) = correct(i,low,high)
       end do
    else if ( fdf_bnnames(pline) == 3 ) then
       g = fdf_bnames(pline,2)
       if ( .not. leqi(g,'from') ) then
          call die('Error in range block: &
               &from <int> to/plus/minus <int> is ill formatted')
       end if
       g = fdf_bnames(pline,3)
       if ( fdf_bnintegers(pline) < 2 ) then
          call die('Error in range block &
               &from <int> to/plus/minus <int> is ill formatted')
       end if
       i1 = fdf_bintegers(pline,1)
       ! The first index *must* be given 
       ! as an atom index
       i1 = correct(i1,low,high)
       i2 = fdf_bintegers(pline,2)
       if ( i1 < low ) then
          call die('Error in range block: &
               &from <int> is below lowest allowed value')
       end if
       ! Initialize step
       step = 1
       if ( fdf_bnintegers(pline) > 2 ) then
          step = abs(fdf_bnintegers(pline,3))
       end if
       if ( step == 0 ) call die('Stepping MUST be different from 0')
       if ( leqi(g,'to') ) then
          ! do nothing....
       else if ( leqi(g,'plus') ) then
          i2 = i1 + i2 - 1
       else if ( leqi(g,'minus') ) then
          step = - step
          i2 = i1 - i2 + 1
       else
          call die('Unrecognized designator of ending range, &
               [to, plus, minus] accepted.')
       end if
       i2 = correct(i2,low,high)
       !print *,i1,i2,step
       if ( high < i2 ) then
          call die('Error in range block: &
               &to <int> is above highest allowed value')
       end if
       if ( (i1 < i2 .and. step < 0) .or. &
            (i1 > i2 .and. step > 0) ) then
          call die('Block range is not consecutive')
       end if
       do i = i1 , i2 , step
          n = n + 1 
          list(n) = i
       end do
    else
       call die('Error in range block, input not recognized')
    end if

    do j = 1 , n
       if ( list(j) < low .or. high < list(j) ) then
          call die('Error in range block. Input is beyond range')
       end if
    end do

    call rgn_list(r,n,list)
    deallocate(list)
    
  contains

    function correct(in,low,high) result(out)
      integer, intent(in) :: in, low, high
      integer :: out
      if ( in < low ) then
         out = high + in + 1
      else if ( high < in ) then
         out = in - high + low - 1
      else
         out = in
      end if
    end function correct
          
  end subroutine fdf_brange

end module fdf_extra
