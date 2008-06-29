! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
program cdf_fft

!
! Computes the Fourier Transform of a grid function in netCDF format
!
!        The input file is called "GridFunc.nc"
!
!  The program gives as output the fourier components, scaled
!  by the cell volume in Bohr^3, and tagged by the corresponding
!  reciprocal lattice vectors.
!
!  If used on a density file, the first coefficient should give
!  the number of electrons in the unit cell.
!
!  Only the first spin component is processed in this prototype version.
!
!  Alberto Garcia, May 2008
!
use netcdf

implicit none

integer, parameter  :: dp = selected_real_kind(14,100)
integer, parameter  :: sp = selected_real_kind(6,20)

integer  :: ncid 
integer  :: xyz_id, step_id, abc_id, spin_id
integer  :: cell_id, n1_id, n2_id, n3_id, gridfunc_id

integer   ::    nspin  ! Number of spins 
integer   ::    n(3), nc(3), ig, jg, kg, i, j, k, ip
integer   ::    n1, n2, n3
integer   ::   ispin, iostat, ix, iy, iz, iret

real(dp)  ::    cell(3,3), cell_volume
real(sp), dimension(:), allocatable :: gridfunc
complex(sp), dimension(:), allocatable :: aa

external c3d_fft
!-----------------------------------------------------

call check( nf90_open('GridFunc.nc',NF90_NOWRITE,ncid))
       call check( nf90_inq_dimid(ncid,'spin',spin_id) )
       call check( nf90_inquire_dimension(ncid, dimid=spin_id, len=nspin) )

       call check( nf90_inq_dimid(ncid,'n1',n1_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n1_id, len=n(1)) )
       call check( nf90_inq_dimid(ncid,'n2',n2_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n2_id, len=n(2)) )
       call check( nf90_inq_dimid(ncid,'n3',n3_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n3_id, len=n(3)) )

       call check( nf90_inq_varid(ncid, "cell", cell_id) )
       call check( nf90_inq_varid(ncid, "gridfunc", gridfunc_id) )

       iret = nf90_get_var(ncid, cell_id, cell, start = (/1, 1 /), &
                        count = (/3, 3/) )
       call check(iret)
       cell_volume = volcel(cell)

   allocate(gridfunc(product(n(1:3))))
   allocate(aa(product(n(1:3))))
   ispin = 1
      iret = nf90_get_var(ncid, gridfunc_id, gridfunc, start = (/1, 1, 1, ispin /), &
           count = (/n(1), n(2), n(3), 1/) )
      call check(iret)

   n1 = n(1)
   n2 = n(2)
   n3 = n(3)

   nc(1:3) = (n(1:3) - 1) /2
   
   aa = cmplx(gridfunc)
   call c3d_fft(n1,n2,n3,aa,+1,n1,n2)
   aa = aa * cell_volume   ! To convert to number of electrons
   do ip = 1, size(aa)
      i = modulo(ip-1,n1) + 1
      j = modulo( (ip - i)/n1, n2) + 1
      k = (ip - i - (j-1)*n1)/(n1*n2)  + 1

      ig = i - n1 - 1
      if (ig < -nc(1)) ig = ig + n1

      jg = j - n2 - 1
      if (jg < -nc(2)) jg = jg + n2

      kg = k - n3 - 1
      if (kg < -nc(3)) kg = kg + n3

      print "(3i4,3x,2f10.6)", ig, jg, kg, aa(ip)
   enddo

        call check( nf90_close(ncid) )

CONTAINS

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) then
  print *, "netCDF error: " // NF90_STRERROR(code)
  STOP
endif
end subroutine check

FUNCTION VOLCEL( C )
real(dp) :: volcel

!  CALCULATES THE VOLUME OF THE UNIT CELL

      real(dp) ::  C(3,3)
      VOLCEL = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) + &
               ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) + &
               ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
      VOLCEL = ABS( VOLCEL )
END function volcel

end program cdf_fft