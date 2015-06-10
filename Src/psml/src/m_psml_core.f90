!> Data structures to handle the PSML pseudopotential format.
!!
!> @author Alberto Garcia
!

module m_psml_core

use assoc_list, only: ps_annotation_t => assoc_list_t
use external_interfaces, only: die => psml_die

implicit none

private

character(len=3), parameter, public  :: PSML_TARGET_VERSION = "0.8"

!----------------------------------------------------------------
! Hardwired parameters (to be made dynamical in a later version)

! Possible maximum number of semilocal potentials:
! s, p, d, f, (g?)  -->  5(sr) + 4(so) or 5(up)+5(down)
integer, parameter, private    :: MAXN_POTS = 20

! Possible maximum number of projectors:
! It might depend on the type of 'so' conversion, and the
! multiplicity, and on whether [sr,so] and [lj] versions
! are stored...
integer, parameter, private    :: MAXN_PROJS = 96

! Maximum number of valence shells (including semicore):
integer, parameter, private    :: MAXN_SHELLS = 20

! Maximum number of pseudo-wavefunctions (as sl potentials):
integer, parameter, private    :: MAXN_WFNS = 10
!----------------------------------------------------------------

integer, parameter, private    :: dp = selected_real_kind(14)
!
!-----------------------------------------------------------

type, public :: provenance_t
        character(len=40)       :: creator = "-----"
        character(len=30)       :: date    = "-----"
end type provenance_t
!------
type, public :: header_t
        character(len=30)       :: atomic_label    !< generalized symbol
        real(kind=dp)           :: z  !< atomic number (might be non-integer)
        real(kind=dp)           :: zpseudo !< Z - ncore-electrons
        character(len=50)       :: flavor  !< pseudization method
        character(len=6)        :: relativity !< "no|scalar|dirac"
        logical                 :: polarized !< is spin_polarized?
        !
        character(len=4)        :: core_corrections !< are there NLCC's?
end type header_t
!------
type, public :: config_val_t
      integer                          :: nshells
      real(kind=dp)                    :: total_charge
      integer, dimension(MAXN_SHELLS)  :: n
      character(len=1), dimension(MAXN_SHELLS) :: l
      real(dp), dimension(MAXN_SHELLS) :: occ
      real(dp), dimension(MAXN_SHELLS) :: occ_up
      real(dp), dimension(MAXN_SHELLS) :: occ_down
end type config_val_t
!------
type, public :: xc_t
        integer                         :: n_functs_libxc = 0
        character(len=50), allocatable  :: libxc_name(:)
        integer, allocatable            :: libxc_id(:)
        real(dp), allocatable           :: libxc_weight(:)
        type(ps_annotation_t)           :: annotation
end type xc_t
!------
type, public :: grid_t
!
!     Note that the preferred option is to have explicit grid data.
!
      integer                        :: npts = 0
      real(dp), pointer              :: grid_data(:) => null()
      type(ps_annotation_t)          :: annotation 
end type grid_t      
!
type, public :: radfunc_t
      type(grid_t), pointer                   :: grid => null()
      real(kind=dp), dimension(:), pointer    :: data => null()
end type radfunc_t      
!
type, public :: semilocal_t
      integer                          :: npots = 0
      integer, dimension(MAXN_POTS)           :: n
      character(len=1), dimension(MAXN_POTS)  :: l
      real(dp), dimension(MAXN_POTS)          :: j
      integer, dimension(MAXN_POTS)           :: set
      character(len=40), dimension(MAXN_POTS) :: flavor
      real(dp), dimension(MAXN_POTS)          :: rc
      type(radfunc_t), dimension(MAXN_POTS)   :: V
end type semilocal_t

type, public :: psoperator_t
   !
   ! Optional provenance information
   !
      type(provenance_t)                        :: provenance
   !
   ! Optional private grid
   !
      type(grid_t), pointer                     :: grid => null()
   !
   ! Vlocal
   !
      type(radfunc_t)                          :: Vlocal
      character(len=40)                        :: vlocal_type
   ! Projectors
      integer                          :: nprojs = 0
      integer, dimension(MAXN_PROJS)           :: seq
      character(len=1), dimension(MAXN_PROJS)  :: l
      real(dp), dimension(MAXN_PROJS)          :: j
      integer, dimension(MAXN_PROJS)           :: set
      character(len=40), dimension(MAXN_PROJS) :: type
      real(dp), dimension(MAXN_PROJS)          :: ekb
      type(radfunc_t), dimension(MAXN_PROJS)   :: proj
end type psoperator_t

type, public :: pswfs_t
      integer                          :: npswfs = 0
      integer, dimension(MAXN_WFNS)           :: n
      character(len=1), dimension(MAXN_WFNS)  :: l
      real(dp), dimension(MAXN_WFNS)          :: j
      integer, dimension(MAXN_WFNS)           :: set
      type(radfunc_t), dimension(MAXN_WFNS)   :: Phi
end type pswfs_t

type, public :: valence_charge_t
      real(dp)        :: total_charge
      type(radfunc_t) :: rho_val
end type valence_charge_t

type, public :: core_charge_t
      integer         :: n_cont_derivs
      real(dp)        :: rcore
      type(radfunc_t) :: rho_core
end type core_charge_t


!> @brief Main derived type to hold the PSML information
!> @author Alberto Garcia
!> @todo It could be used also by psf and vps readers
type, public :: ps_t
      type(provenance_t)                 :: provenance
      type(header_t)                     :: header
      type(config_val_t)                 :: config_val
      type(xc_t)                         :: xc_info
      type(grid_t), pointer              :: global_grid => null()
      type(semilocal_t)                  :: semilocal
      type(psoperator_t)                 :: psoperator
      type(pswfs_t)                      :: pswfs
      !
      type(valence_charge_t)             :: valence_charge
      type(core_charge_t)                :: core_charge

   end type ps_t

   integer,  parameter, public   &
                          :: SET_SREL     =   1, &
                             SET_NONREL   =   2, &
                             SET_SO       =   4, &
                             SET_LJ       =   8, &
                             SET_UP       =  16, &
                             SET_DOWN     =  32, &
                             SET_SPINAVE  =  64, &
                             SET_SPINDIFF = 128, &         ! 2^7
                             SET_USER1    = 256, &         ! 2^8
                             SET_USER2    = 512            ! 2^9

   integer, parameter, public    :: SET_ALL =  2**10 -1

 public  :: ps_destroy
 public  :: str_of_set

 public  :: setcode_of_string    ! utility function, not for client normal use

 CONTAINS
!> @brief Cleans the ps object
!> @author Alberto Garcia
!> @date March-July 2014
subroutine ps_destroy(ps)
type(ps_t), intent(inout)     :: ps

integer :: i

call destroy_xc(ps%xc_info)

!
! Note that freshly declared objects must have
! npots = 0 and npswfs = 0 !
!
do i = 1, ps%semilocal%npots
   call destroy_radfunc(ps%semilocal%V(i))
enddo
!
do i = 1, ps%psoperator%nprojs
   call destroy_radfunc(ps%psoperator%proj(i))
enddo

call destroy_radfunc(ps%psoperator%vlocal)
if (associated(ps%psoperator%grid)) then
   call destroy_grid(ps%psoperator%grid)
endif
!
do i = 1, ps%pswfs%npswfs
   call destroy_radfunc(ps%pswfs%Phi(i))
enddo
!
call destroy_radfunc(ps%valence_charge%rho_val)
call destroy_radfunc(ps%core_charge%rho_core)
!
if (associated(ps%global_grid)) then
   call destroy_grid(ps%global_grid)
endif

end subroutine ps_destroy

subroutine destroy_radfunc(rp)
type(radfunc_t) :: rp

if (associated(rp%grid)) then
   call destroy_grid(rp%grid)
endif
if (associated(rp%data)) then
   deallocate(rp%data)
   rp%data => null()
endif
end subroutine destroy_radfunc

subroutine destroy_grid(gp)
type(grid_t), pointer :: gp

if (associated(gp%grid_data)) then
   deallocate(gp%grid_data)
   gp%grid_data => null()
endif
gp => null()
end subroutine destroy_grid

subroutine destroy_xc(xp)
type(xc_t), intent(inout) :: xp

if (allocated(xp%libxc_name)) deallocate(xp%libxc_name)
if (allocated(xp%libxc_id)) deallocate(xp%libxc_id)
if (allocated(xp%libxc_weight)) deallocate(xp%libxc_weight)

end subroutine destroy_xc

function setcode_of_string(str) result(code)
       character(len=*), intent(in) :: str
       integer                      :: code

       select case (trim(str))
       case ("non_relativistic")
          code = SET_NONREL
       case ("scalar_relativistic")
          code = SET_SREL
       case ("spin_orbit")
          code = SET_SO
       case ("lj")
          code = SET_LJ
       case ("spin_up")
          code = SET_UP
       case ("spin_down")
          code = SET_DOWN
       case ("spin_average")
          code = SET_SPINAVE
       case ("spin_difference")
          code = SET_SPINDIFF
       case ("user_extension1")
          code = SET_USER1
       case ("user_extension2")
          code = SET_USER2
       case ("all","any")
          code = SET_ALL
       case default
          call die("Wrong set string: "//trim(str))
       end select

end function setcode_of_string

function str_of_set(code) result(str)
       integer, intent(in)          :: code
       character(len=20)            :: str

       character(len=100) :: msg

       select case (code)
       case (SET_NONREL)
          str ="non_relativistic"
       case (SET_SREL)
          str ="scalar_relativistic"
       case (SET_SO)
          str ="spin_orbit"
       case (SET_LJ)
          str ="lj"
       case (SET_UP)
          str ="spin_up"
       case (SET_DOWN)
          str ="spin_down"
       case (SET_SPINAVE)
          str ="spin_average"
       case (SET_SPINDIFF)
          str ="spin_difference"
       case (SET_USER1)
          str ="user_extension1"
       case (SET_USER2)
          str ="user_extension2"
       case (SET_ALL)
          str ="all"
       case default
          write(msg,"(a,i4)") "Wrong set code: ", code
          call die(msg)
       end select

end function str_of_set

end module m_psml_core






