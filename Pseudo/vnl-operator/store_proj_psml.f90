subroutine  store_proj_psml(l,ikb,ekb,nrc,erefkb,dkbcos,nrval,rofi,proj)

  use m_kb, only: nprojs, kbprojs

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

!  write(filename,"(a,i1,a,i1)")  "KBproj-rl+1.", l, ".", ikb
!  call file_out(nrc,rofi,proj,trim(filename)) 


  nprojs = nprojs + 1
  kbprojs(nprojs)%l = l
  kbprojs(nprojs)%nrval = nrval
  kbprojs(nprojs)%seq = ikb
  kbprojs(nprojs)%erefkb = erefkb
  kbprojs(nprojs)%ekb = ekb
  kbprojs(nprojs)%dkbcos = dkbcos
  kbprojs(nprojs)%l = l
  
  ! Restore factor of r**l (See in KBproj)
  ! for compliance with PSML standard (1.0 and above!)
  allocate(kbprojs(nprojs)%r(nrval))
  allocate(kbprojs(nprojs)%proj(nrval))
  do ir = 1, nrval
     kbprojs(nprojs)%r(ir) = rofi(ir)
     kbprojs(nprojs)%proj(ir) = proj(ir) * rofi(ir)**l
  enddo

  end subroutine store_proj_psml
