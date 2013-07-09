!
! See module TrialOrbitalClass after this
!
module swan_com

use precision,         only: dp
implicit none

interface debug
  module procedure debug_integer_array,debug_complex_array,          &
   debug_real_array,debug_string,debug_integer,debug_complex,        &
   debug_real, debug_boolean_array,debug_boolean,debug_real_2,       &
   debug_complex_2
end interface debug

interface str
  module procedure str_char_integer,str_char_dp,str_char_dp3,        &
    str_char_boolean,str_char_complex
end interface str

contains
!
! Auxilliary communication subroutines--------------------------------
!
subroutine hail()
  print *,"swan: Hi! There is swan here."
end subroutine

subroutine echo(what)
  character(len=*)      :: what
  write(6,*) "swan: ",trim(what)
end subroutine

subroutine debug_boolean_array(msg,what)
  logical, intent(in),dimension(:)   :: what
  character(len=*),intent(in)        :: msg
  call echo("DeBuG:"// msg)
  print *,what
end subroutine debug_boolean_array

subroutine debug_boolean(msg,what)
  logical, intent(in)                :: what
  character(len=*),intent(in)        :: msg
  call echo(str("DeBuG:"// msg,what))
end subroutine debug_boolean

subroutine debug_integer_array(msg,what)
  character(len=*),intent(in)        :: msg
  integer,intent(in),dimension(:)    :: what
  call echo("DeBuG:"// msg)
  write(*,*) what
end subroutine debug_integer_array

subroutine debug_complex_2(msg,what)
  character(len=*),intent(in)        :: msg
  complex(dp),intent(in),dimension(:,:):: what
  call echo("DeBuG:"// msg)
  write(*,*) what
end subroutine debug_complex_2

subroutine debug_complex_array(msg,what)
  character(len=*),intent(in)        :: msg
  complex(dp),intent(in),dimension(:):: what
  call echo("DeBuG:"// msg)
  write(*,*) what
end subroutine debug_complex_array

subroutine debug_real_array(msg,what)
  character(len=*),intent(in)        :: msg
  real(dp),intent(in),dimension(:)   :: what
  call echo("DeBuG:"// msg)
  write(*,*) what
end subroutine debug_real_array

subroutine debug_integer(msg,what)
  character(len=*),intent(in)        :: msg
  integer,intent(in)                 :: what
  call echo("DeBuG:"// str(msg,what))
end subroutine debug_integer

subroutine debug_complex(msg,what)
  character(len=*),intent(in)        :: msg
  complex(dp),intent(in)             :: what
  call echo("DeBuG:"// str(msg,what))
end subroutine debug_complex

subroutine debug_real(msg,what)
  character(len=*),intent(in)        :: msg
  real(dp),intent(in)                :: what
  call echo("DeBuG:"// str(msg,what))
end subroutine debug_real

subroutine debug_string(what)
  character(len=*)      :: what
  call echo("DeBuG: " // what)
end subroutine debug_string

subroutine debug_real_2(msg,what)
  character(len=*),intent(in)      :: msg
  real(dp),intent(in)              :: what(:,:)
  call echo("DeBuG:"// msg)
  write(*,*) what
end subroutine

!
! Data concatenation: string + something -> string
!
function str_char_boolean(msg,what)
  character(len=*),intent(in):: msg
  logical,intent(in)         :: what
  character(len=len(msg)+8)  :: str_char_boolean
  write(str_char_boolean,fmt="(a,' ',l)") msg,what
end function str_char_boolean

function str_char_integer(msg,what)
  character(len=*)           :: msg
  character(len=len(msg)+6)  :: str_char_integer
  integer                    :: what
  write(str_char_integer,fmt="(a,' ',i5)") msg,what
end function str_char_integer

function str_char_complex(msg,what)
  character(len=*)           :: msg
  character(len=len(msg)+6)  :: str_char_complex
  complex(dp)                :: what
  write(str_char_complex,fmt="(a,' ',2g15.7)") msg,what
end function str_char_complex

function str_char_dp(msg,what)
  character(len=*)           :: msg
  character(len=len(msg)+16) :: str_char_dp
  real(dp)                   :: what
  write(str_char_dp,fmt="(a,' ',g15.7)") msg,what
end function str_char_dp

function str_char_dp3(msg,what)
  character(len=*)           :: msg
  character(len=len(msg)+16*3) :: str_char_dp3
  real(dp)                   :: what(:)
  write(str_char_dp3,fmt="(a,' ',3g15.7)") msg,what
end function str_char_dp3

end module swan_com
!
! WARNING: NEVER USE DIVISION WITH _DP
! ie 3_dp/2_dp = 1_dp
!    3d0/2d0   = 1.5d0
!
module TrialOrbitalClass
use swan_com
use sys,             only:die
use precision,       only:dp
use units,           only:pi,Ang

implicit none

type TrialOrbital
  real(dp),dimension(3) :: centre ! Orbital centre; Bohr
  real(dp),dimension(3) :: zaxis  ! Angular momentum z-axis,unit.vec.
  real(dp),dimension(3) :: xaxis  ! Angular momentum x-axis
  real(dp),dimension(3) :: yaxis  ! Angular momentum y-axis
  real(dp)              :: zovera ! z/a, diffusivity, spread: Bohr^-1
  integer               :: r      ! Radial quantum number
  integer               :: l      ! Angular momentum
  integer               :: mr     ! z-projection quantum number
  real(dp)              :: rCut   ! Siesta's cut-off radius: Bohr
  integer               :: lMax   ! Maximum total angular momentum
end type

! Cut-off radii in units of \alpha^-1
!real(dp),parameter,dimension(3),private ::                           &
!                           cutoffs = (/5.50_dp,15.86_dp,30.65_dp/)
!real(dp),parameter,private :: T = 0.001_dp
real(dp),parameter,dimension(3),private ::                           &
                           cutoffs = (/6.934_dp,18.87_dp,35.44_dp/)
! Squared norm tolerance governing the cut-off radii
real(dp),parameter,private :: T = 0.0001_dp

!! A static instance of the TrialOrbital class
!type(TrialOrbital),private :: static

contains

!
! Vyvojova poznamka: preco nezamenit jednotky uz v readNnkp()?
!
function construct(  &
!
! This function should be the only way of creating the TrialOrbital
! instance
  latVec,  &  ! Lattice vectors in Bohrs
  centre,  &  ! Orbital centre as given in .nnkp, ie ScaledByLatVec
  zaxis,   &  ! Angular momentum z-axis
  xaxis,   &  ! Angular momentum x-axis
  zovera,  &  ! z/a, diffusivity, spread; Units: Ang^-1
  r,       &  ! Radial quantum number
  l,       &  ! Angular momentum
  mr )     &  ! z-projection quantum number
 result(orb)
!
! Arguments:
!
  real(dp),dimension(3,3),intent(in) :: latVec
  real(dp),dimension(3),intent(in) :: centre,zaxis,xaxis
  real(dp),intent(in)              :: zovera
  integer,intent(in)               :: r,l,mr
! iterators
  integer                          :: i,j
! result
  type(TrialOrbital)               :: orb

! Now convert "centre" from ScaledByLatticeVectors to Bohrs
  do i=1,3
    orb%centre(i) = 0_dp
    do j=1,3
      orb%centre(i) = centre(j)*latVec(j,i) + orb%centre(i)
    enddo
  enddo
  orb%zaxis  = zaxis
  orb%xaxis  = xaxis
  orb%zovera = zovera*Ang ! To Bohr^-1
  orb%r      = r
  orb%l      = l
  orb%mr     = mr

! Compute the y-axis
  orb%yaxis = vectorProduct(zaxis,xaxis)
  orb%yaxis = orb%yaxis/                                 &
    sqrt(DOT_PRODUCT(orb%yaxis,orb%yaxis))
  select case(l)
    case(0:3)
      orb%lMax = l
    case(-3:-1)
      orb%lMax = 1 ! spN hybrids
    case(-5:-4)
      orb%lMax = 2 ! spdN hybrids
    case default
      call die(str("Invalid l in TrialOrbital.construct() :",l))
  end select
! Further checks
  if (mr.lt.1) call die(str("Invalid mr in TrialOrbital.construct() :",mr))
  if (r.gt.3.or.r.lt.1) then
    call die(str("Invalid r in TrialOrbital.construct() :",r))
  elseif (l.ge.0.and.mr.gt.2*l+1) then
    call die(str("Invalid mr in TrialOrbital.construct() :",mr))
  elseif (l.lt.0.and.mr.gt.1-l) then
    call die(str("Invalid mr in TrialOrbital.construct() :",mr))
  endif
! Cut-off initialization
  call debug("cutoffs(r)",cutoffs(r))
  orb%rCut = cutoffs(r)/orb%zovera ! in Bohr
  call debug("orb%rCut",orb%rCut)

  contains
  function vectorProduct(a,b)
    real(dp),intent(in),dimension(3) :: a,b
    real(dp),dimension(3)            :: vectorProduct
    vectorProduct(1) = a(2)*b(3)-a(3)*b(2)
    vectorProduct(2) = a(3)*b(1)-a(1)*b(3)
    vectorProduct(3) = a(1)*b(2)-a(2)*b(1)
  end function vectorProduct
end function construct

!
!<-----------------------WAVE FUNCTIONS----------------------------->
!
real(dp) function getTrialWaveFunction(orbital,atPoint)
!
! Trial orbital as defined in Wannier90.
! Yields the function value at a given point relative to its centre.
! It contains the library of Wannier90 trial orbitals in the form
! of statement functions.
!
! Argument: coordinates of the point in Bohrs
! Output unit: Bohr^{-3/2}
!
  real(dp),intent(in),dimension(3):: atPoint
  type(TrialOrbital),intent(in)   :: orbital

  real(dp),dimension(3)           :: arg
! Angular-dependent factor of the wave-function
  real(dp)                        :: angular 
!
! Statement functions: declaration
!
  real(dp)                        :: z,y,x,sphere,s,px,py,pz,rr
  real(dp)                        :: dz2,dxz,dyz,dxy,dx2y2
  real(dp)    :: fz3,fxz2,fyz2,fzx2y2,fxyz,fxx23y2,fy3x2y2
  real(dp)                        :: sp_1,sp_2,sp2_1,sp2_2,sp2_3
  real(dp)                        :: sp3_1,sp3_2,sp3_3,sp3_4
  real(dp)                        :: sp3d_1,sp3d_2,sp3d_3,sp3d_4,sp3d_5
  real(dp)    :: sp3d2_1,sp3d2_2,sp3d2_3,sp3d2_4,sp3d2_5,sp3d2_6
  integer                         :: rank
  real(dp)                        :: R1,R2,R3

!
! Inverse square roots
!
  real(dp),parameter              :: rs2 = 1_dp/sqrt(2d0)
  real(dp),parameter              :: rs3 = 1_dp/sqrt(3d0)
  real(dp),parameter              :: rs6 = 1_dp/sqrt(6d0)
  real(dp),parameter              :: rs12 = 1_dp/sqrt(12d0)
!
! Constants depending on l and power of z
!
  real(dp),parameter              :: l0norm = 1_dp/sqrt(4_dp*pi)
  real(dp),parameter              :: l1norm = sqrt(3_dp*1d0)*l0norm
  real(dp),parameter              :: l2z2norm = sqrt((1d0*5_dp/16_dp)/pi)
  real(dp),parameter              :: l2z1norm = sqrt((1d0*15_dp/4_dp)/pi)
  real(dp),parameter              :: l2z0norm = l2z1norm*0.5_dp
  real(dp),parameter              :: l3z3norm = sqrt(7_dp/pi)/4_dp
  real(dp),parameter              :: l3z2norm = sqrt(1d0*21_dp/2_dp/pi)/4_dp
  real(dp),parameter              :: l3z1norm = sqrt(105_dp/pi)/4_dp
  real(dp),parameter              :: l3z0norm = sqrt(1d0*35_dp/2_dp/pi)/4_dp
!
! Real Spherical Harmonics
!
! To get a dimensionless spherical harmonics
sphere(x,y,z,rank) = sqrt(x**2+y**2+z**2)**rank
          s(x,y,z) = l0norm
         px(x,y,z) = l1norm*x/sphere(x,y,z,1)
         py(x,y,z) = l1norm*y/sphere(x,y,z,1)
         pz(x,y,z) = l1norm*z/sphere(x,y,z,1)

        dz2(x,y,z) = l2z2norm*(2_dp*z**2-x**2-y**2)/sphere(x,y,z,2)
        dxz(x,y,z) = l2z1norm*z*x/sphere(x,y,z,2)
        dyz(x,y,z) = l2z1norm*z*y/sphere(x,y,z,2)
      dx2y2(x,y,z) = l2z0norm*(x**2-y**2)/sphere(x,y,z,2)
        dxy(x,y,z) = l2z0norm*x*y*2_dp/sphere(x,y,z,2)

        fz3(x,y,z) = l3z3norm*(2_dp*z**2-3_dp*x**2-3_dp*y**2)*z/sphere(x,y,z,3)
       fxz2(x,y,z) = l3z2norm*(4_dp*z**2-x**2-y**2)*x/sphere(x,y,z,3)
       fyz2(x,y,z) = l3z2norm*(4_dp*z**2-x**2-y**2)*y/sphere(x,y,z,3)
     fzx2y2(x,y,z) = l3z1norm*z*(x**2-y**2)/sphere(x,y,z,3)
       fxyz(x,y,z) = l3z1norm*z*x*y*2_dp/sphere(x,y,z,3)
    fxx23y2(x,y,z) = l3z0norm*(x**2-3_dp*y**2)*x/sphere(x,y,z,3)
    fy3x2y2(x,y,z) = l3z0norm*(3_dp*x**2-y**2)*y/sphere(x,y,z,3)
!
! Hybrids
!
       sp_1(x,y,z) = rs2*s(x,y,z)+rs2*px(x,y,z)
       sp_2(x,y,z) = rs2*s(x,y,z)-rs2*px(x,y,z)

      sp2_1(x,y,z) = rs3*s(x,y,z)-rs6*px(x,y,z)+rs2*py(x,y,z)
      sp2_2(x,y,z) = rs3*s(x,y,z)-rs6*px(x,y,z)-rs2*py(x,y,z)
      sp2_3(x,y,z) = rs3*s(x,y,z)+rs6*px(x,y,z)*2_dp
      
      sp3_1(x,y,z) = 0.5_dp*(s(x,y,z)+px(x,y,z)+py(x,y,z)+pz(x,y,z))
      sp3_2(x,y,z) = 0.5_dp*(s(x,y,z)+px(x,y,z)-py(x,y,z)-pz(x,y,z))
      sp3_3(x,y,z) = 0.5_dp*(s(x,y,z)-px(x,y,z)+py(x,y,z)-pz(x,y,z))
      sp3_4(x,y,z) = 0.5_dp*(s(x,y,z)-px(x,y,z)-py(x,y,z)+pz(x,y,z))

     sp3d_1(x,y,z) = rs3*s(x,y,z)-rs6*px(x,y,z)+rs2*py(x,y,z)
     sp3d_2(x,y,z) = rs3*s(x,y,z)-rs6*px(x,y,z)-rs2*py(x,y,z)
     sp3d_3(x,y,z) = rs3*s(x,y,z)+rs6*px(x,y,z)*2_dp
     sp3d_4(x,y,z) = rs2*pz(x,y,z) + rs2*dz2(x,y,z)
     sp3d_5(x,y,z) = -rs2*pz(x,y,z) + rs2*dz2(x,y,z)

    sp3d2_1(x,y,z) = rs6*s(x,y,z)-rs2*px(x,y,z)-rs12*dz2(x,y,z)+dx2y2(x,y,z)/2d0
    sp3d2_2(x,y,z) = rs6*s(x,y,z)+rs2*px(x,y,z)-rs12*dz2(x,y,z)+dx2y2(x,y,z)/2d0
    sp3d2_3(x,y,z) = rs6*s(x,y,z)-rs2*px(x,y,z)-rs12*dz2(x,y,z)-dx2y2(x,y,z)/2d0
    sp3d2_4(x,y,z) = rs6*s(x,y,z)+rs2*px(x,y,z)-rs12*dz2(x,y,z)-dx2y2(x,y,z)/2d0
    sp3d2_5(x,y,z) = rs6*s(x,y,z)-rs2*pz(x,y,z)+rs3*dz2(x,y,z)
    sp3d2_6(x,y,z) = rs6*s(x,y,z)+rs2*pz(x,y,z)+rs3*dz2(x,y,z)
! Radial part
!
  R1(rr) = 2_dp*orbital%zovera**(3d0/2d0)*                          &
            exp(-orbital%zovera*rr)
  R2(rr) = 0.5_dp/sqrt(2d0)*orbital%zovera**(3d0/2d0)*              &
            (2_dp-orbital%zovera*rr)*exp(-orbital%zovera*rr/2d0)
  R3(rr) = sqrt(4d0/27d0)*orbital%zovera**(3d0/2d0)*                &
            (1_dp-2_dp*orbital%zovera*rr/3d0 +                          &
          2_dp*orbital%zovera**2*rr**2/27d0)*exp(-orbital%zovera*rr/3d0)

!
! Executables                 !----------->
!
 !arg = atPoint - orbital%centre Recall that this is done in phiatm()

  arg = atPoint
  x = DOT_PRODUCT(orbital%xaxis,arg)
  y = DOT_PRODUCT(orbital%yaxis,arg)
  z = DOT_PRODUCT(orbital%zaxis,arg)
  rr = sphere(x,y,z,1)
!
! If out of the rCut sphere then vanish
!
  if (rr.gt.orbital%rCut) then
    getTrialWaveFunction = 0_dp
    return
  endif
!
! Decipher arguments: radial part
!
  select case(orbital%r)
    case(1)
      getTrialWaveFunction = R1(rr)
    case(2)
      getTrialWaveFunction = R2(rr)
    case(3)
      getTrialWaveFunction = R3(rr)
  end select
! Renormalize the wave function, since we cut it at R_c
  getTrialWaveFunction = getTrialWaveFunction/sqrt(1_dp-T)
! Decipher arguments: angular one
!
  if (rr.eq.0_dp) then
! Special treatment of the origin
    select case(orbital%l)
      case (0)
        angular = l0norm
      case (2)
        if(orbital%mr.eq.1) then
! l=2,mr=1 doesn't vanish!
!          angular = -l2z2norm  ... this is the limit in the z=0 plane
! but other limiting value is +2*l2z2norm along the z axis
! so there is a cusp and we take it's average: zero
        else
          angular = 0_dp
        endif
! Hybrids are combinations of dz2 and s limits
      case (-1)
        angular = l0norm/sqrt(2d0)
      case (-2)
        angular = l0norm/sqrt(3d0)
      case (-3)
        angular = l0norm/2_dp
      case (-4)
        if (orbital%mr.lt.4) then
          angular = l0norm/sqrt(3d0)
        else
          !angular = -l2z2norm/sqrt(2d0)
        endif
      case(-5)
        if (orbital%mr.lt.5) then
          angular = l0norm/sqrt(6d0)!+l2z2norm/sqrt(12d0)
        else
          angular = l0norm/sqrt(6d0)!-l2z2norm/sqrt(3d0)
        endif
      case default
        angular = 0_dp
    end select
  else 
!
! rr.not equal.0
!
  select case(orbital%l)
    case(0)
        angular = s(x,y,z)
    case(1)
      select case(orbital%mr)
        case(1)
          angular = pz(x,y,z)
        case(2)
          angular = px(x,y,z)
        case(3)
          angular = py(x,y,z)
      end select
    case(2)
      select case(orbital%mr)
        case(1)
          angular = dz2(x,y,z)
        case(2)
          angular = dxz(x,y,z)
        case(3)
          angular = dyz(x,y,z)
        case(4)
          angular = dx2y2(x,y,z)
        case(5)
          angular = dxy(x,y,z)
      end select
    case(3)
      select case(orbital%mr)
        case(1)
          angular = fz3(x,y,z)
        case(2)
          angular = fxz2(x,y,z)
        case(3)
          angular = fyz2(x,y,z)
        case(4)
          angular = fzx2y2(x,y,z)
        case(5)
          angular = fxyz(x,y,z)
        case(6)
          angular = fxx23y2(x,y,z)
        case(7)
          angular = fy3x2y2(x,y,z)
      end select
!
! Hybrids
!
    case(-1)
      select case(orbital%mr)
        case(1)
          angular = sp_1(x,y,z)
        case(2)
          angular = sp_2(x,y,z)
      end select
    case(-2)
      select case(orbital%mr)
        case(1)
          angular = sp2_1(x,y,z)
        case(2)
          angular = sp2_2(x,y,z)
        case(3)
          angular = sp2_3(x,y,z)
      end select
    case(-3)
      select case(orbital%mr)
        case(1)
          angular = sp3_1(x,y,z)
        case(2)
          angular = sp3_2(x,y,z)
        case(3)
          angular = sp3_3(x,y,z)
        case(4)
          angular = sp3_4(x,y,z)
      end select
    case(-4)
      select case(orbital%mr)
        case(1)
          angular = sp3d_1(x,y,z)
        case(2)
          angular = sp3d_2(x,y,z)
        case(3)
          angular = sp3d_3(x,y,z)
        case(4)
          angular = sp3d_4(x,y,z)
        case(5)
          angular = sp3d_5(x,y,z)
      end select
    case(-5)
      select case(orbital%mr)
        case(1)
          angular = sp3d2_1(x,y,z)
        case(2)
          angular = sp3d2_2(x,y,z)
        case(3)
          angular = sp3d2_3(x,y,z)
        case(4)
          angular = sp3d2_4(x,y,z)
        case(5)
          angular = sp3d2_5(x,y,z)
        case(6)
          angular = sp3d2_6(x,y,z)
      end select
  end select
  endif
  !print *,x,getTrialWaveFunction,orbital%zovera**(3_dp/2_dp)
  !print *,3_dp/2_dp,1.5_dp,3_dp/2d0,3d0/2_dp
  getTrialWaveFunction = getTrialWaveFunction*angular
end function getTrialWaveFunction

real(dp) function getTrialRCut(ofWhat)
  type(TrialOrbital),intent(in)    :: ofWhat
  getTrialRCut = ofWhat%rCut
end function

integer function getTrialLMax(ofWhat)
  type(TrialOrbital),intent(in)    :: ofWhat
  getTrialLMax = ofWhat%lMax
end function

function getTrialCentre(ofWhat)
  real(dp),dimension(3)            :: getTrialCentre
  type(TrialOrbital),intent(in)    :: ofWhat
  getTrialCentre = ofWhat%centre
end function getTrialCentre


subroutine print_trialOrb(what)
  type(TrialOrbital),intent(in)    :: what
  write(*,fmt='(a,3f8.3,a)') "centre =",what%centre,"Bohr"
  write(*,fmt='(a,3f8.3)') "zaxis  =",what%zaxis
  write(*,fmt='(a,3f8.3)') "xaxis  =",what%xaxis
  write(*,fmt='(a,3f8.3)') "yaxis  =",what%yaxis
  write(*,fmt='(a,1f8.3,a)') "zovera =",what%zovera,"/Bohr"
  write(*,*) "r      =",what%r
  write(*,*) "mr     =",what%mr
  write(*,*) "l      =",what%l
  write(*,'(a,1f8.3,a)') "rCut   =",what%rCut,"Bohr"
  write(*,*) "lMax   =",what%lMax
end subroutine print_trialOrb

end module TrialOrbitalClass

!/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|
