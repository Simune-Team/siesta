 logical function byteSizesOK()
  ! Users may call this function once to ensure that the kind parameters 
  !   the module defines are available with the current compiler. 
  ! We can't ensure that the two REAL kinds are actually four and 
  !   eight bytes long, but we can ensure that they are distinct. 
  ! Early Fortran 90 compilers would sometimes report incorrect results for 
  !   the bit_size intrinsic, but I haven't seen this in a long time. 

    ! Local variables
    integer (kind =  OneByteInt) :: One
    integer (kind =  TwoByteInt) :: Two
    integer (kind = FourByteInt) :: Four

    if (bit_size( One) == 8  .and. bit_size( Two) == 16 .and.  &
        bit_size(Four) == 32 .and.                             &
        FourByteReal > 0 .and. EightByteReal > 0 .and. &
        FourByteReal /= EightByteReal) then
      byteSizesOK = .true.
    else
      byteSizesOK = .false.
    end if
  end function byteSizesOK
