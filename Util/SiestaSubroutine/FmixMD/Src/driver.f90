program driver

! Driver to test fmixmd and fsiesta
! J.M.Soler. Nov.2003

  use sample_m

  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: maxa = 1000

  character(len=80)  :: fconv, ffast
  integer            :: i, ia, na
  real(dp)           :: cell(3,3), dtconv, dtfast, dtsamp, ma(maxa), &
                        t0, tf, va(3,maxa), xa(3,maxa)

! Read data
  read(5,*) ffast
  read(5,*) fconv
  read(5,*) t0, tf
  read(5,*) dtfast, dtconv, dtsamp
  do i = 1,3
    read(5,*) cell(:,i)
  end do
  read(5,*) na
  do ia = 1,na
    read(5,*) ma(ia), xa(:,ia), va(:,ia)
  end do

! Initialize sampling
  call sample('init')

! Perform mixed-force MD simulation
  call fmixmd( t0, tf, dtfast, dtconv, dtsamp, ffast, fconv, &
               na, ma(1:na), xa(:,1:na), va(:,1:na), cell )

! Print sampling results
  call sample('print')

end
