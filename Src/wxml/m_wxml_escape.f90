module m_wxml_escape

implicit none

private

public  :: escape_char_array
public  :: check_Name

character(len=*), parameter :: LETTER = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
character(len=*), parameter :: DIGIT  = '1234567890'
character(len=*), parameter :: ESCAPE = '&<'//'"'//"'"

CONTAINS

!This must accept input in the form of a character array,
!otherwise calculation of output length is not possible
!in a specification expression in f90.
function escape_char_array(s_a) result(s_out)
  character, dimension(:) :: s_a
  character(len=size(s_a) + 4*count(s_a=='&')  &
                          + 3*count(s_a=='<')  &
                          + 5*count(s_a=='"')  &
                          + 5*count(s_a=="'")) :: s_out

  integer :: i, i_out
  character(len=1) :: c
  
  i_out = 1
  do i = 1, size(s_a)
    c = s_a(i)
    select case (c)
      case ('&')
        s_out(i_out:i_out+4) = '&amp;'
        i_out=i_out+5
      case('<')
        s_out(i_out:i_out+3) = '&lt;'
        i_out=i_out+4
      case ('"')
        s_out(i_out:i_out+3) = '&apos;'
        i_out=i_out+6
      case ("'")
        s_out(i_out:i_out+3) = '&quot;'
        i_out=i_out+6
      case default
        s_out(i_out:i_out) = c
        i_out = i_out + 1
    end select
  enddo

end function escape_char_array


function check_Name(name) result(good)
 character(len=*), intent(in) :: name
 logical :: good
! Validates a string against the XML requirements for a NAME
! Is not fully compliant; ignores UTF issues.

 integer :: n, i

 good=.true.

 n = len(name)
 if (n == 0) then
   write(0,*) "tagName is an empty string."
   good=.false.
 elseif (good) then
   if (scan(name(1:1), LETTER//'_'//':') == 0) then
     write(0,*) "tagName must begin with a letter, underscore or colon: '", name(1:1), "'."
     good=.false.
   endif
 elseif (good) then
   do i = 1, n
     if (scan(name(i:i), LETTER//DIGIT//'.'//'-'//'_'//':') == 0) then
       write(0,*) "tagName contains a forbidden character: '", name(i:i), "'."
       good=.false.
       exit
     endif
   enddo
 elseif (good) then
   if (scan(name(1:1), 'Xx') == 1 .and. &
       scan(name(2:2), 'Mm') == 1 .and. &
       scan(name(3:3), 'Ll') == 1) then
     write(0,*) "tagName cannot start with the characters 'XML'."
     good=.false.
   endif
 endif
       
end function check_Name

end module m_wxml_escape
