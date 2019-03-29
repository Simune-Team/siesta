program v_info
  use m_gridfunc

  type(gridfunc_t) :: gf
  character(len=256) :: fname

  integer :: i
  
  real(grid_p), allocatable :: average(:,:)
  
  fname = "VH"
  call read_gridfunc(fname,gf)
  call get_planar_average(gf,3,average)
  do i = 1, size(average,dim=1)
     print *, i, average(i,1)
  enddo
end program v_info
