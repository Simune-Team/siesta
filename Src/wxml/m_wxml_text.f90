module m_wxml_text

implicit none
!
integer, private, parameter ::  sp = selected_real_kind(6,30)
integer, private, parameter ::  dp = selected_real_kind(14,100)
!
private
public :: str

interface str
   module procedure  str_integer_fmt, str_integer, &
                     str_logical_fmt, str_logical, &
                     str_real_dp, str_real_sp
end interface

CONTAINS

      function str_integer_fmt(i,format) result(s)
      integer, intent(in)   :: i
      character(len=*), intent(in) :: format
      character(len=100)    :: s

      write(s,format) i
      s = adjustl(s)
      end function str_integer_fmt

      function str_integer(i) result(s)
        ! This will work correctly (return an appropriately-sized
        ! string) for integers i s.t. -99999999<=i<=999999999
        integer, intent(in) :: i
        character(len=int(merge(log10(real(max(abs(i),1)))+1, &
                                log10(real(max(abs(i),1)))+2, &
                          sign(i,1)>0))) :: s

        character(len=4) :: form
        
        write(form,'(a,i1,a)') '(i',len(s),')'
        write(s, form) i

      end function str_integer

      function str_logical_fmt(l,format) result(s)
      logical, intent(in)   :: l
      character(len=*), intent(in) :: format
      character(len=100)    :: s

      write(s,format) l
      s = adjustl(s)

      end function str_logical_fmt

      function str_logical(l) result(s)
        logical, intent(in)   :: l
        character(len=merge(4,5,l)) :: s
        
        if (l) then
          s='True'
        else
          s='False'
        endif
      end function str_logical

      function str_real_dp(x,format) result(s)
      real(kind=dp), intent(in)   :: x
      character(len=*), intent(in), optional  :: format
      character(len=100)    :: s

      if (present(format)) then
         write(s,format) x
      else
         write(s,"(g22.12)") x
      endif
      s = adjustl(s)
      end function str_real_dp

      function str_real_sp(x,format) result(s)
      real(kind=sp), intent(in)   :: x
      character(len=*), intent(in), optional  :: format
      character(len=100)    :: s

      if (present(format)) then
         write(s,format) x
      else
         write(s,"(g22.12)") x
      endif
      s = adjustl(s)
      end function str_real_sp


end module m_wxml_text
