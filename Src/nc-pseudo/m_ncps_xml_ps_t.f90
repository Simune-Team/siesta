module m_ncps_xml_ps_t
!
! Data structures to match the XML pseudopotential format
!
integer, parameter, private    :: MAXN_POTS = 8
integer, parameter, private    :: dp = selected_real_kind(14)
!
public  :: dump_pseudo
!
!-----------------------------------------------------------
type, public :: grid_t
!
!     It should be possible to represent both log and linear
!     grids with a few parameters here.
!
      character(len=20)              :: type
      real(kind=dp)                  :: scale
      real(kind=dp)                  :: step 
      integer                        :: npts 
end type grid_t      
!
type, public :: radfunc_t
      type(grid_t)                            :: grid
      real(kind=dp), dimension(:), pointer    :: data => null()
end type radfunc_t      
      
type, public :: vps_t
      character(len=1)               :: l
      integer                        :: n
      integer                        :: spin
      real(kind=dp)                  :: occupation
      real(kind=dp)                  :: cutoff
      type(radfunc_t)                :: V
end type vps_t

type, public :: pswf_t
      character(len=1)               :: l
      integer                        :: n
      integer                        :: spin
      type(radfunc_t)                :: V
end type pswf_t

type, public :: header_t
        character(len=2)        :: symbol
        real(kind=dp)           :: zval
        real(kind=dp)           :: gen_zval  ! Valence charge at generation
        character(len=10)       :: creator
        character(len=10)       :: date
        character(len=40)       :: flavor
        logical                 :: relativistic
        logical                 :: polarized
        character(len=30)       :: xcfunctionaltype
        character(len=30)       :: xcfunctionalparametrization
        character(len=4)        :: core_corrections
end type header_t

type, public :: xml_ps_t
      type(header_t)                     :: header
      integer                            :: npots 
      integer                            :: npswfs
      integer                            :: npots_down
      integer                            :: npots_up 
      type(vps_t), dimension(MAXN_POTS)  :: pot
      type(pswf_t), dimension(MAXN_POTS) :: pswf
      type(radfunc_t)                    :: core_charge
      type(radfunc_t)                    :: valence_charge
   end type xml_ps_t


CONTAINS !===============================================

subroutine dump_pseudo(pseudo,lun)
  type(xml_ps_t), intent(in), target   :: pseudo
  integer, intent(in)                  :: lun

integer  :: i
type(vps_t), pointer :: pp
type(pswf_t), pointer :: pw
type(radfunc_t), pointer :: rp

write(lun,*) "---PSEUDO data:"

do i = 1, pseudo%npots
      pp =>  pseudo%pot(i)
      rp =>  pseudo%pot(i)%V
      write(lun,*) "VPS ", i, " angular momentum: ", pp%l
      write(lun,*) "                 n: ", pp%n
      write(lun,*) "                 occupation: ", pp%occupation
      write(lun,*) "                 cutoff: ", pp%cutoff
      write(lun,*) "                 spin: ", pp%spin
      write(lun,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
enddo
do i = 1, pseudo%npswfs
      pw =>  pseudo%pswf(i)
      rp =>  pseudo%pswf(i)%V
      write(lun,*) "PSWF ", i, " angular momentum: ", pw%l
      write(lun,*) "                 n: ", pw%n
      write(lun,*) "                 spin: ", pw%spin
      write(lun,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
enddo
rp => pseudo%valence_charge
write(lun,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
rp => pseudo%core_charge
write(lun,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)

end subroutine dump_pseudo

end module m_ncps_xml_ps_t




















