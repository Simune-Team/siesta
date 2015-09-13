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

module m_ts_io_contour
!
! Routines that are used to read in and print out the contour for integration of the GFs
! 
  use precision, only : dp

  implicit none

  character(len=200), parameter :: OPT_N   = '(''ts: '',a)'
  character(len=200), parameter :: OPT_C   = '(''ts: '',a,t53,''=    '',a)'
  character(len=200), parameter :: OPT_CC  = '(''ts: '',a,t53,''=    '',a,tr2,a)'
  character(len=200), parameter :: OPT_F   = '(''ts: '',a,t53,''='',f10.4)'
  character(len=200), parameter :: OPT_INT = '(''ts: '',a,t53,''='',i5)'
  character(len=200), parameter :: OPT_F_U = '(''ts: '',a,t53,''='',f10.4,tr1,a)'
  character(len=200), parameter :: OPT_G_U = '(''ts: '',a,t53,''='',g11.4,tr1,a)'

contains

  subroutine write_e(str,val,unit)

    use units, only : eV, Kelvin

    character(len=*), intent(in) :: str
    real(dp), intent(in) :: val
    character(len=*), intent(in), optional :: unit

    ! This is our definition of infinity....
    if ( abs(val) > 10000._dp ) then
       if ( val < 0._dp ) then
          write(*,opt_c) trim(str),'-Infinity'
       else
          write(*,opt_c) trim(str),' Infinity'
       end if
    else if ( .not. present(unit) ) then
       write(*,opt_f_u) trim(str),val / eV,'eV'
    else
       if ( unit == 'eV' ) then
          write(*,opt_f_u) trim(str),val / eV,'eV'
       else if ( unit == 'Ry' ) then
          write(*,opt_f_u) trim(str),val,'Ry'
       else if ( unit == 'K' ) then
          write(*,opt_f_u) trim(str),val/Kelvin,'K'
       else
          call die('Programming error, unknown unit')
       end if
    end if

  end subroutine write_e

end module m_ts_io_contour
