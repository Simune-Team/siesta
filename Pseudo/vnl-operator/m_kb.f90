module m_kb
  integer, parameter :: dp = selected_real_kind(10,100)
  
type, public :: kb_t
  integer  :: l
  integer  :: j
  integer  :: seq
  real(dp) :: ekb
  real(dp) :: erefkb
  real(dp) :: dkbcos
  integer  :: nrval
  real(dp), allocatable :: r(:)
  real(dp), allocatable  :: proj(:)
end type kb_t

! This is used as a counter...
integer   :: nprojs

! Array of structures to hold the KB information
type(kb_t), public, allocatable, target  :: kbprojs(:)

public :: init_kbproj_structs
public :: store_proj_psml

private

CONTAINS

  subroutine init_kbproj_structs(n_pjnl,nkbtot)
    integer, intent(in) :: n_pjnl
    integer, intent(in) :: nkbtot

    ! nprojs is used as a counter by store_proj_psml...

    allocate (kbprojs(n_pjnl))
  end subroutine init_kbproj_structs

  subroutine  store_proj_psml(l,lj_projs,jk,ikb,ekb,nrc,erefkb,dkbcos,nrval,rofi,proj)

  integer, parameter :: dp = selected_real_kind(10,100)
!
! Custom routine to process the information about each projector
! This version for psop -- psml generation
!
  integer, intent(in)  :: l
  logical, intent(in)  :: lj_projs
  integer, intent(in)  :: jk
  integer, intent(in)  :: ikb
  integer, intent(in)  :: nrc
  real(dp), intent(in) :: ekb
  real(dp), intent(in) :: erefkb
  real(dp), intent(in) :: dkbcos
  integer, intent(in)  :: nrval
  real(dp), intent(in) :: rofi(:)
  real(dp), intent(in) :: proj(:)

  integer  :: ir

  nprojs = nprojs + 1
  kbprojs(nprojs)%l = l
  if (lj_projs) then
     kbprojs(nprojs)%j = l + (2*jk-3)*0.5_dp ! l -/+ 1/2
     if (l == 0) kbprojs(nprojs)%j =  0.5_dp ! special case: only jk=1
  else
     kbprojs(nprojs)%j = 0
  endif     
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

end module m_kb
