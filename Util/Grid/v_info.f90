program v_info
  use m_gridfunc

  implicit none

  type(gridfunc_t) :: gf
  character(len=256) :: fname

  integer, parameter :: dp = selected_real_kind(10,100)
  
  real(dp), parameter :: Ang    = 1.0_dp / 0.529177_dp
  real(dp) :: z_coord
  integer :: i

  real(grid_p), allocatable :: average(:,:)
  
  fname = "VH"
  call read_gridfunc(fname,gf)
  call get_planar_average(gf,3,average)
  if (monoclinic_z(gf%cell)) then
     print *, "#   z (Ang)     Value "
     do i = 1, size(average,dim=1)
        z_coord = gf%origin(3) + (i-1) * gf%cell(3,3) / gf%n(3)
        print *, z_coord/Ang, average(i,1)
     enddo
  else
     do i = 1, size(average,dim=1)
        print *, i, average(i,1)
     enddo
  endif
end program v_info
