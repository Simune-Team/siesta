! Description:
!   Provide named kind parameters for use in declarations of real and integer 
!    variables with specific byte sizes (i.e. one, two, four, and eight byte
!    integers; four and eight byte reals). The parameters can then be used
!    in (KIND = XX) modifiers in declarations.
!   A single function (byteSizesOK()) is provided to ensure that the selected 
!    kind parameters are correct.
!  
! Input Parameters:
!   None.
!
! Output Parameters:
!   Public parameters, fixed at compile time:
!     OneByteInt, TwoByteInt, FourByteInt, EightByteInt
!                                     FourByteReal, EightByteRadl
!
! References and Credits:
!   Written by
!    Robert Pincus
!    Cooperative Institue for Meteorological Satellite Studies
!    University of Wisconsin - Madison
!    1225 W. Dayton St. 
!    Madison, Wisconsin 53706
!    Robert.Pincus@ssec.wisc.edu
!
! Design Notes:
!   Fortran 90 doesn't allow one to check the number of bytes in a real variable;
!     we check only that four byte and eight byte reals have different kind parameters. 
!

  integer, parameter, public ::   OneByteInt = selected_int_kind(2), &
                          TwoByteInt = selected_int_kind(4), &
                         FourByteInt = selected_int_kind(9), &
                        EightByteInt = selected_int_kind(18)

  integer, parameter, public ::                                          &
                        FourByteReal = selected_real_kind(P =  6, R =  37), &
                       EightByteReal = selected_real_kind(P = 15, R = 307)


