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
module m_os
!
! Contains routines for changing file names, append file names and conforming
! with OS standards.
!
  implicit none

  private

  public :: file_ext
  public :: file_noext
  public :: basename
  public :: dirname
!  public :: cwd

contains

! Retrive the extension of the file (without the ".")
! This can only handle "." as a separator. 
! If "." is not found in the filename it will return and empty extension.
  pure function file_ext(filename)
    character(len=*), intent(in) :: filename
    character(len=512) :: file_ext
    integer :: i, fL
    fL = len(trim(filename))
    file_ext = ''
    do i = fL , 1 , -1
       if ( filename(i:i) == "/" .or. filename(i:i) == "\\" ) then !Handles dots in directory names
          return
       else if ( filename(i:i) == "." ) then
          file_ext = filename(i+1:fL)
          return
       end if
    end do
  end function file_ext

! Retrive the filename without the extension
! This can only handle "." as a separator. 
! If "." is not found in the filename it will return and empty extension.
! It relies on the function file_ext
  pure function file_noext(filename)
    character(len=*), intent(in) :: filename
    character(len=512) :: file_noext
    integer :: i, fL, fExtL
    fL = len(trim(filename))
    fExtL = len(trim(file_ext(filename)))
    if ( fExtL == 0 ) then
       file_noext = filename(1:fL)
    else
       file_noext = filename(1:fL-fExtL-1)
    end if
  end function file_noext

! Retrive the basename of the file, that is the name without the directory.
! This can only handle "/" or "\" as a directory separator.  (DO NOT NAME FILES with "/" or "\" in the name)
  pure function basename(filename)
    character(len=*), intent(in) :: filename
    character(len=512) :: basename
    integer :: i, fL
    fL = len(trim(filename))
    basename = filename(1:fL)
    do i = fL , 1 , -1
       if ( filename(i:i) == "/" .or. filename(i:i) == '\\' ) then
          basename = filename(i+1:fL)
          return
       end if
    end do
  end function basename

! Retrive the directoryname of the file, that is the directory without the name (with any potential directory divider).
! This can only handle "/" or "\" as a directory separator.  (DO NOT NAME FILES with "/" or "\" in the name)
  pure function dirname(filename)
    character(len=*), intent(in) :: filename
    character(len=512) :: dirname
    integer :: i, fL, bL
    fL = len(trim(filename))
    bL = len(trim(basename(filename)))
    if ( bL == fL ) then
       dirname = './'
    else
       dirname = filename(1:fL-bL)
    end if
  end function dirname

! Maybe this function needs a configure step to check for the compiler
! having the routine.
! This is ambiguous in that either "call getcwd" or "status = getcwd" 
! is both valid, however, using one of them disables the other...
!  function cwd()
!    character(len=512) :: cwd
!    integer :: status
!    call getcwd(cwd,status)
!    if ( status /= 0 ) then
!       cwd = './'
!    end if
!  end function cwd

end module m_os
