module m_wxml_escape

use m_wxml_buffer

private

public  :: add_to_buffer_escaping_markup
public  :: check_Name

CONTAINS

subroutine add_to_buffer_escaping_markup(s,buf)
character(len=*), intent(in)      ::   s
type(buffer_t), intent(inout)     ::   buf

integer           :: len_s, i
character(len=1)  :: c

len_s = len(s)
i = 0
do 
 if (i==len_s) exit
 i = i + 1
 c = s(i:i)
 if (c == "<") then
    call add_to_buffer("&lt;",buf)
 else if (c == "&") then
    call add_to_buffer("&amp;",buf)
 else
    call add_to_buffer(c,buf)
 endif
enddo

end subroutine add_to_buffer_escaping_markup

subroutine check_Name(name)
 character(len=*), intent(in) :: name

 integer :: n

 n = len_trim(name)
 if (n == 0) then
    write(*,*) "tagName cannot be an empty string: '"
    stop
 endif

    if (scan(name(1:1), "1234567890") .ne. 0) then
       write(*,*) "tagName can not begin with a digit: '", name, "'"
       stop
    endif

    if (scan(name(1:1),       &
       "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_:") .eq. 0) &
               then          
       write(*,*)  &
       "tagName must begin with a letter, underscore or colon: '", &
               name, "'"
       stop
    endif

    if (scan(name, " ") .ne. 0) then
       write(*,*) "tagName can not contain whitespace: '", name, "'"
       stop
    end if

    if (name(1:1).eq."x" .or. name(1:1) .eq. "X") then
       if (n == 1) RETURN
       if (name(2:2).eq."m" .or. name(2:2) .eq. "M") then
          if (n == 2) RETURN
          if (name(3:3).eq."l" .or. name(3:3) .eq. "L") then
             write(*,*) "tagName can not start with the letters XML: '", &
                  name, "'"
             stop
          endif
       endif
    endif

end subroutine check_Name

end module m_wxml_escape
