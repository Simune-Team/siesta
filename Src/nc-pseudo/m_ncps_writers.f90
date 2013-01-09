module m_ncps_writers

  implicit none

  integer, parameter  :: dp = selected_real_kind(14)
      
  public :: pseudo_write_formatted

CONTAINS

!----
  subroutine pseudo_write_formatted(fname,p,print_gen_zval)
    use m_ncps_froyen_ps_t, only: froyen_ps_t

    character(len=*), intent(in)  :: fname
    type(froyen_ps_t), intent(in) :: p
    logical, intent(in), optional :: print_gen_zval

    integer io_ps, i, j

    call get_free_lun(io_ps)
    open(io_ps,file=fname,form='formatted',status='unknown', &
         action="write",position="rewind")

 8000   format(1x,i2)
 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010   format(1x,6a10,/,1x,a70)
 8015   format(1x,2i3,i5,4g20.12)
 8030   format(4(g20.12))
 8040   format(1x,a)

    write(io_ps,8005) p%name, p%icorr, p%irel, p%nicore
    write(io_ps,8010) (p%method(i),i=1,6), p%text
    if (present(print_gen_zval)) then
       if (print_gen_zval) then
          write(io_ps,8015) p%npotd, p%npotu, p%nr, &
               p%b, p%a, p%zval, p%gen_zval
       else
          write(io_ps,8015) p%npotd, p%npotu,    &
               p%nr, p%b, p%a, p%zval
       endif
    else
       write(io_ps,8015) p%npotd, p%npotu, p%nr, &
            p%b, p%a, p%zval
    endif

    write(io_ps,8040) "Radial grid follows"
    write(io_ps,8030) (p%r(j),j=2,p%nrval)

    do i=1,p%npotd
       write(io_ps,8040) "Down Pseudopotential follows (l on next line)"
       write(io_ps,8000) p%ldown(i)
       write(io_ps,8030) (force_underflow(p%vdown(i,j)), j=2,p%nrval)
    enddo

    do i=1,p%npotu
       write(io_ps,8040) "Up Pseudopotential follows (l on next line)"
       write(io_ps,8000) p%lup(i)
       write(io_ps,8030) (force_underflow(p%vup(i,j)), j=2,p%nrval)
    enddo

    write(io_ps,8040) "Core charge follows"
    write(io_ps,8030) (force_underflow(p%chcore(j)),j=2,p%nrval)
    write(io_ps,8040) "Valence charge follows"
    write(io_ps,8030) (force_underflow(p%chval(j)),j=2,p%nrval)

    close(io_ps)
  end subroutine pseudo_write_formatted
!--------
!
  function force_underflow(x) result(res)
    real(dp), intent(in) ::  x
    real(dp)             ::  res

!     Avoid very small numbers that might need a three-character
!     exponent field in formatted output
      
    if (abs(x) .lt. 1.0e-99_dp) then
       res = 0.0_dp
    else
       res = x
    endif

  end function force_underflow

  subroutine get_free_lun(lun)
    integer, intent(out) :: lun

    interface
       subroutine die(str)
         character(len=*), intent(in), optional :: str
       end subroutine die
    end interface

    logical :: used
    integer :: iostat

    do lun= 10,90
       inquire(unit=lun, opened=used, iostat=iostat)
       if (iostat .ne. 0) used = .true.
       if (.not. used) return  ! normal return with 'lun' value
    enddo
    call die("No luns available")

  end subroutine get_free_lun
      

end module m_ncps_writers



