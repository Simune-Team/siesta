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
      module radial

      use precision
      use xml
      use flib_spline, only: spline, splint

      implicit none

      private

      public :: rad_alloc, rad_get, rad_setup_d2, rad_zero
      public :: radial_read_ascii, radial_dump_ascii
      public :: radial_dump_xml, reset_rad_func

      type, public :: rad_func
         integer          n
         double precision cutoff         
         double precision delta
         double precision, pointer :: f(:)=>null()   ! Actual data
         double precision, pointer :: d2(:)=>null()  ! 2nd derivative
      end type rad_func

      CONTAINS

      subroutine reset_rad_func( func )
      implicit none
      type(rad_func)   :: func

      func%n = 0
      nullify(func%f)
      nullify(func%d2)
      end subroutine reset_rad_func

      subroutine rad_alloc(func,n)
      use alloc, only: re_alloc
      implicit none
!
!     Sets the 'size' n of the arrays and allocates f and d2.
!
      type(rad_func), intent(inout)    :: func
      integer, intent(in)        :: n
      func%n = n
      nullify(func%f,func%d2)
      call re_alloc( func%f, 1, n, 'func%f', 'rad_alloc' )
      call re_alloc( func%d2, 1, n, 'func%d2', 'rad_alloc' )
!      allocate(func%f(n),func%d2(n))
      end subroutine rad_alloc

      subroutine rad_get(func,r,fr,dfdr)
      type(rad_func), intent(in) :: func
      real(dp), intent(in)         :: r
      real(dp), intent(out)        :: fr
      real(dp), intent(out)        :: dfdr

      if (func%n .eq. 0) then
          fr = 0._dp
          dfdr = 0._dp
       else
          call splint(func%delta,func%f,func%d2,func%n,r,fr,dfdr)
       endif
      
      end subroutine rad_get
!
!     Set up second derivative in a radial function
!
      subroutine rad_setup_d2(func,yp1,ypn)
      type(rad_func), intent(inout) :: func
      real(dp), intent(in)          :: yp1, ypn

      if (func%n .eq. 0) return
      call spline(func%delta,func%f,func%n,yp1,ypn,func%d2)
      
      end subroutine rad_setup_d2

      subroutine rad_zero(func)
      type(rad_func), intent(inout) :: func
      func%n      = 0
      end subroutine rad_zero
!
!     Do not use yet... interface in need of fuller specification
!
      function rad_rvals(func) result (r)
      use alloc, only: re_alloc
      implicit none
      real(dp), dimension(:), pointer :: r
      type(rad_func), intent(in) :: func

      integer i

      nullify(r)
      if (func%n .eq. 0) return
!      allocate(r(func%n))
      call re_alloc( r, 1, func%n, 'r', 'rad_rvals' )
      do i=1,func%n
        r(i) = func%delta *(i-1)
      enddo
      end function rad_rvals


      subroutine radial_read_ascii(op,lun,yp1,ypn)
      type(rad_func)    :: op 
      real(dp), intent(in)          :: yp1, ypn

      integer lun
      integer j, npts
      real(dp) dummy

      read(lun,*) npts, op%delta, op%cutoff
      call rad_alloc(op,npts)
      do j=1,npts
         read(lun,*) dummy, op%f(j)
      enddo
      call rad_setup_d2(op,yp1,ypn)
      end subroutine radial_read_ascii
!
!--------------------------------------------------------------------
      subroutine radial_dump_ascii(op,lun,header)
      type(rad_func)    :: op
      integer           :: lun
      logical, intent(in), optional :: header

      integer :: j
      logical :: print_header
!
!     The standard dump is to unit "lun"
!     and includes a header with npts, delta, and cutoff
!
      print_header = .true.
      if (present(header)) then
         print_header = header
      endif
!
      if (print_header) then
         write(lun,'(i4,2g22.12,a)') op%n,
     $        op%delta, op%cutoff, " # npts, delta, cutoff"
      endif
      do j=1,op%n
         write(lun,'(2g22.12)') (j-1)*op%delta, op%f(j)
      enddo

      end subroutine radial_dump_ascii
!--------------------------------------------------------------------
!
      subroutine radial_dump_xml(op,lun)

      type(rad_func)    :: op
      integer lun
      integer j

      write(lun,'(a)') '<radfunc>'
      call xml_dump_element(lun,'npts',str(op%n))
      call xml_dump_element(lun,'delta',str(op%delta))
      call xml_dump_element(lun,'cutoff',str(op%cutoff))
      write(lun,'(a)') '<data>'
      do j=1,op%n
         write(lun,'(2g22.12)') (j-1)*op%delta, op%f(j)
      enddo
      write(lun,'(a)') '</data>'
      write(lun,'(a)') '</radfunc>'
      end subroutine radial_dump_xml
!
      end module radial










