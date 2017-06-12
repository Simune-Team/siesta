subroutine  write_proj_psml(l,ikb,ekb,nrc,erefkb,dkbcos,nrval,rofi,proj)

  use local_xml, only: xf, rl, fval, nrl
  use xmlf90_wxml

  integer, parameter :: dp = selected_real_kind(10,100)
!
! Custom routine to process the information about each projector
! This version for psop -- psml generation
!
  integer, intent(in)  :: l
  integer, intent(in)  :: ikb
  integer, intent(in)  :: nrc
  real(dp), intent(in) :: ekb
  real(dp), intent(in) :: erefkb
  real(dp), intent(in) :: dkbcos
  integer, intent(in)  :: nrval
  real(dp), intent(in) :: rofi(:)
  real(dp), intent(in) :: proj(:)

  integer  :: ir
  real(dp) :: rc
  real(dp), allocatable :: f(:)
  character(len=1), dimension(0:4) ::  lsymb = (/'s','p','d','f','g'/)

  rc = rofi(nrc)

!  write(filename,"(a,i1,a,i1)")  "KBproj-rl+1.", l, ".", ikb
!  call file_out(nrc,rofi,proj,trim(filename)) 

  ! Restore factor of r**l (See in KBproj)
  ! for compliance with PSML standard (1.0 and above!)
  allocate(f(nrval))
  do ir = 1, nrval
     f(ir) = proj(ir) * rofi(ir)**l
  enddo

  call xml_NewElement(xf,"proj")
  call my_add_attribute(xf,"l",lsymb(l))
  call my_add_attribute(xf,"seq",str(ikb))
  call my_add_attribute(xf,"ekb",str(0.5_dp*ekb))
  call my_add_attribute(xf,"type","KB")
  call xml_NewElement(xf,"radfunc")
  call xml_NewElement(xf,"data")
  call xml_AddArray(xf, f(1:nrval))
  
  call xml_EndElement(xf,"data")
  call xml_EndElement(xf,"radfunc")
  call xml_EndElement(xf,"proj")

!!         Common block with the information about the  KB projectors
!          call comKB(is, a, b, rofi, proj,
!     .               l, ikb, rc(ikb,l), ekb(ikb,l), nrc)

  deallocate(f)

  CONTAINS

      subroutine my_add_attribute(xf,name,value)
      use xmlf90_wxml 
      type(xmlf_t), intent(inout)   :: xf
      character(len=*), intent(in)  :: name
      character(len=*), intent(in)  :: value

       call xml_AddAttribute(xf,name,trim(value))
      end subroutine my_add_attribute

end subroutine write_proj_psml
