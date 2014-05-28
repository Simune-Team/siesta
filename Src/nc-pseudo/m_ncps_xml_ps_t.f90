module m_ncps_xml_ps_t
!
! Data structures and functions to handle the XML pseudopotential format.
!
! The PSML library will eventually have three sections:
!
! A. Definition of data structures.
! B. Accessors for pseudopotential information and evaluation
! C. Parsing helpers
!
! In order for this module will hide completely the implementation of
! the parsing of the XML file and the later querying of the data
! structures, the psxml_t type should be made private. But then all the
! accessors and the parsing helpers should reside in this module too, making
! it quite big.
!
! Currently the parsing helpers are in 'm_ncps_parsing_helpers'
! This scattering will be solved with the final packaging of the library.
!
implicit none

integer, parameter, private    :: MAXN_POTS = 8
integer, parameter, private    :: dp = selected_real_kind(14)
!
!-----------------------------------------------------------
type, public :: grid_t
!
!     It should be possible to represent both log and linear
!     grids with a few parameters here.
!     Note that the preferred option is to have explicit grid data.
!
      character(len=20)              :: type
      real(kind=dp)                  :: scale
      real(kind=dp)                  :: step 
      integer                        :: npts = 0
      real(dp), pointer              :: grid_data(:) => null()
end type grid_t      
!
type, public :: radfunc_t
      type(grid_t), pointer                   :: grid => null()
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
        real(kind=dp)           :: z         ! atomic number (might be non-integer)
        real(kind=dp)           :: zval
        real(kind=dp)           :: gen_zval  ! Valence charge at generation
        character(len=10)       :: creator
        character(len=10)       :: date
        character(len=40)       :: flavor
        logical                 :: relativistic
        logical                 :: polarized
        character(len=80)       :: xc_libxc_string ! LibXC string
        integer                 :: xc_libxc_code   ! LibXC numerical code
        ! Generator's own terminology for XC 
        character(len=30)       :: xcfunctionaltype
        character(len=30)       :: xcfunctionalparametrization
        !
        character(len=4)        :: core_corrections
        logical                 :: rV
end type header_t

type, public :: xml_ps_t
      type(header_t)                     :: header
      integer                            :: npots 
      integer                            :: npswfs
      integer                            :: npots_down
      integer                            :: npots_up 
      type(grid_t), pointer              :: global_grid => null()
      type(vps_t), dimension(MAXN_POTS)  :: pot
      type(pswf_t), dimension(MAXN_POTS) :: pswf
      type(radfunc_t)                    :: core_charge
      type(radfunc_t)                    :: valence_charge
   end type xml_ps_t

!
public  :: xml_ps_destroy
public  :: dump_pseudo           ! For initial debugging
!

! Accessor list
public :: psxmlEvaluateValenceCharge
public :: psxmlEvaluateCoreCharge
public :: psxmlAtomicSymbol
public :: psxmlAtomicNumber
public :: psxmlCreator
public :: psxmlDate
public :: psxmlPseudoFlavor
public :: psxmlPseudoZval
public :: psxmlGenerationZval
public :: psxmlXCLibXCCode
public :: psxmlXCLibXCString
public :: psxmlXCFunctional
public :: psxmlXCFunctionalType
public :: psxmlIsRelativistic
public :: psxmlIsSpinPolarized
public :: psxmlHasCoreCorrections
public :: psxmlHasGlobalLogGrid
public :: psxmlGridNpoints
public :: psxmlLogGridStep
public :: psxmlLogGridScale
public :: psxmlGridRmax
public :: psxmlPotentialsUp
public :: psxmlPotentialsDown
public :: psxmlPotAngMomentum
public :: psxmlOccupation
public :: psxmlGenerationCutoff
public :: psxmlPrincipalN
public :: psxmlEvaluatePotential
public :: psxmlEvaluatePseudoOrbital

CONTAINS !===============================================

subroutine xml_ps_destroy(p)
type(xml_ps_t), pointer :: p

integer :: i

if (.not. associated(p)) RETURN

if (associated(p%global_grid)) then
   call destroy_grid(p%global_grid)
endif
do i = 1, p%npots
   call destroy_radfunc(p%pot(i)%V)
enddo
do i = 1, p%npswfs
   call destroy_radfunc(p%pswf(i)%V)
enddo
call destroy_radfunc(p%valence_charge)
call destroy_radfunc(p%core_charge)

deallocate(p)
p => null()

end subroutine xml_ps_destroy

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

!========================================================

function psxmlEvaluateValenceCharge(p,r,debug) result(val)
type(xml_ps_t), intent(in) :: p
real(dp), intent(in)       :: r
logical, intent(in), optional :: debug
real(dp)                   :: val

val = eval_radfunc(p%valence_charge,r,debug)
end function psxmlEvaluateValenceCharge

function psxmlEvaluateCoreCharge(p,r,debug) result(val)
type(xml_ps_t), intent(in) :: p
real(dp), intent(in)       :: r
logical, intent(in), optional :: debug
real(dp)                   :: val

val = eval_radfunc(p%core_charge,r,debug)
end function psxmlEvaluateCoreCharge
!
function psxmlAtomicSymbol(psxml) result(name)
type(xml_ps_t), intent(in) :: psxml
character(len=2) :: name
name = psxml%header%symbol
end function psxmlAtomicSymbol
!
function psxmlAtomicNumber(psxml) result(z)
type(xml_ps_t), intent(in) :: psxml
real(dp) :: z
 z = psxml%header%z
end function psxmlAtomicNumber
!
function psxmlCreator(psxml) result(name)
type(xml_ps_t), intent(in) :: psxml
character(len=len_trim(psxml%header%creator)) :: name
name = trim(psxml%header%creator)
end function psxmlCreator
!
function psxmlDate(psxml) result(str)
type(xml_ps_t), intent(in) :: psxml
character(len=len_trim(psxml%header%date)) :: str
str = trim(psxml%header%date)
end function psxmlDate
!
function psxmlPseudoFlavor(psxml) result(str)
type(xml_ps_t), intent(in) :: psxml
character(len=len_trim(psxml%header%flavor)) :: str
str = trim(psxml%header%flavor)
end function psxmlPseudoFlavor
!
function psxmlPseudoZval(psxml) result(zval)
type(xml_ps_t), intent(in) :: psxml
real(dp)                   :: zval
zval = psxml%header%zval
end function psxmlPseudoZval
!
function psxmlGenerationZval(psxml) result(zval)
type(xml_ps_t), intent(in) :: psxml
real(dp)                   :: zval
zval = psxml%header%gen_zval
end function psxmlGenerationZval
!
function psxmlXCLibXCString(psxml) result(xc_string)
type(xml_ps_t), intent(in) :: psxml
character(len=len_trim(psxml%header%xc_libxc_string)) :: xc_string
xc_string = trim(psxml%header%xc_libxc_string)
end function psxmlXCLibXCString
!
function psxmlXCLibXCCode(psxml) result(xc_code)
type(xml_ps_t), intent(in) :: psxml
integer  :: xc_code
xc_code = psxml%header%xc_libxc_code
end function psxmlXCLibXCCode
!
function psxmlXCFunctional(psxml) result(xc_string)
type(xml_ps_t), intent(in) :: psxml
character(len=len_trim(psxml%header%xcfunctionalparametrization)) :: xc_string
xc_string = trim(psxml%header%xcfunctionalparametrization)
end function psxmlXCFunctional
!
function psxmlXCFunctionalType(psxml) result(xc_type)
type(xml_ps_t), intent(in) :: psxml
character(len=len_trim(psxml%header%xcfunctionaltype)) :: xc_type
xc_type = trim(psxml%header%xcfunctionaltype)
end function psxmlXCFunctionalType
!
function psxmlIsRelativistic(psxml) result(rel)
type(xml_ps_t), intent(in) :: psxml
logical                    :: rel
rel = psxml%header%relativistic
end function psxmlIsRelativistic
!
function psxmlIsSpinPolarized(psxml) result(pol)
type(xml_ps_t), intent(in) :: psxml
logical                    :: pol
pol = psxml%header%polarized
end function psxmlIsSpinPolarized
!
function psxmlHasCoreCorrections(psxml) result(cc)
type(xml_ps_t), intent(in) :: psxml
logical                    :: cc
cc = (psxml%header%core_corrections == "yes")
end function psxmlHasCoreCorrections
!
function psxmlHasGlobalLogGrid(psxml) result(log_grid)
type(xml_ps_t), intent(in) :: psxml
logical                    :: log_grid
! to be fixed: use a global grid...
if (.not.(associated(psxml%global_grid))) then
   log_grid = .false.
else
   log_grid = (psxml%global_grid%type == "log")
endif
end function psxmlHasGlobalLogGrid

function psxmlGridNpoints(psxml) result(npts)
type(xml_ps_t), intent(in) :: psxml
integer                    :: npts

if (.not.(associated(psxml%global_grid))) then
   call die("Do not have global grid to get npts...")
endif
npts = psxml%global_grid%npts
end function psxmlGridNpoints
!
function psxmlLogGridStep(psxml) result(step)
type(xml_ps_t), intent(in) :: psxml
real(dp)                   :: step

if (.not. psxmlHasGlobalLogGrid(psxml)) then
   call die("Do not have global log grid to get step...")
endif
step = psxml%global_grid%step
end function psxmlLogGridStep
!
function psxmlLogGridScale(psxml) result(scale)
type(xml_ps_t), intent(in) :: psxml
real(dp)                   :: scale

if (.not. psxmlHasGlobalLogGrid(psxml)) then
   call die("Do not have global log grid to get scale...")
endif
scale = psxml%global_grid%scale
end function psxmlLogGridScale
!
function psxmlGridRmax(psxml) result(rmax)
type(xml_ps_t), intent(in) :: psxml
real(dp)                   :: rmax

integer :: npts

! We should make sure that there is a global grid...
! Perhaps this should be the maximum rmax among all 
! possible grids...

npts = psxml%global_grid%npts
rmax = psxml%global_grid%grid_data(npts)
end function psxmlGridRmax
!
function psxmlPotentialsUp(psxml) result(n)
type(xml_ps_t), intent(in) :: psxml
integer                    :: n
n = psxml%npots_up
end function psxmlPotentialsUp
!
function psxmlPotentialsDown(psxml) result(n)
type(xml_ps_t), intent(in) :: psxml
integer                    :: n
n = psxml%npots_down
end function psxmlPotentialsDown
!
function psxmlPotAngMomentum(psxml,ud,i) result(l)
type(xml_ps_t), intent(in) :: psxml
character, intent(in)      :: ud
integer,   intent(in)      :: i
integer                    :: l

character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)
integer :: ndown
character(len=1) :: str

ndown = psxmlPotentialsDown(psxml)

select case(ud)
   case ( "u", "U")
      if (i > psxmlPotentialsUp(psxml)) then
         call die("attempt to get l from non-existing Up potential")
      endif
      str = psxml%pot(ndown+i)%l
   case ( "d", "D")
      if (i > ndown) then
         call die("attempt to get l from non-existing Down potential")
      endif
      str = psxml%pot(i)%l
end select
!
do l = 0,4
   if (str == sym(l)) RETURN
enddo
call die("Wrong l symbol in potential")

end function psxmlPotAngMomentum
!
function psxmlOccupation(psxml,ud,i) result(zo)
type(xml_ps_t), intent(in) :: psxml
character, intent(in)      :: ud
integer,   intent(in)      :: i
real(dp)                   :: zo

integer :: ndown

ndown = psxmlPotentialsDown(psxml)

select case(ud)
   case ( "u", "U")
      if (i > psxmlPotentialsUp(psxml)) then
         call die("attempt to get occupation from non-existing Up potential")
      endif
      zo = psxml%pot(ndown+i)%occupation
   case ( "d", "D")
      if (i > ndown) then
         call die("attempt to get l from non-existing Down potential")
      endif
      zo = psxml%pot(i)%occupation
end select
end function psxmlOccupation
!
function psxmlGenerationCutoff(psxml,ud,i) result(rc)
type(xml_ps_t), intent(in) :: psxml
character, intent(in)      :: ud
integer,   intent(in)      :: i
real(dp)                   :: rc

integer :: ndown

ndown = psxmlPotentialsDown(psxml)

select case(ud)
   case ( "u", "U")
      if (i > psxmlPotentialsUp(psxml)) then
         call die("attempt to get cutoff from non-existing Up potential")
      endif
      rc = psxml%pot(ndown+i)%cutoff
   case ( "d", "D")
      if (i > ndown) then
         call die("attempt to get cutoff from non-existing Down potential")
      endif
      rc = psxml%pot(i)%cutoff
end select
end function psxmlGenerationCutoff
!
function psxmlPrincipalN(psxml,ud,i) result(n)
type(xml_ps_t), intent(in) :: psxml
character, intent(in)      :: ud
integer,   intent(in)      :: i
integer                    :: n

integer :: ndown

ndown = psxmlPotentialsDown(psxml)

select case(ud)
   case ( "u", "U")
      if (i > psxmlPotentialsUp(psxml)) then
         call die("attempt to get n from non-existing Up potential")
      endif
      n = psxml%pot(ndown+i)%n
   case ( "d", "D")
      if (i > ndown) then
         call die("attempt to get n from non-existing Down potential")
      endif
      n = psxml%pot(i)%n
end select
end function psxmlPrincipalN
!
function psxmlEvaluatePotential(psxml,ud,i,r,debug) result(val)
type(xml_ps_t), intent(in) :: psxml
character, intent(in)      :: ud
integer,   intent(in)      :: i
real(dp),  intent(in)      :: r
logical, intent(in), optional :: debug
real(dp)                   :: val


integer :: ndown
real(dp), parameter :: tiny = 1.0e-8_dp
real(dp) :: reff
logical :: hasRfactor

ndown = psxml%npots_down
hasRfactor = psxml%header%rV

reff = r
if (r == 0.0_dp) then
   reff = tiny
endif

select case(ud)
   case ( "u", "U")
      if (i > psxmlPotentialsUp(psxml)) then
         call die("attempt to evaluate non-existing Up potential")
      endif
      val = eval_radfunc(psxml%pot(ndown+i)%V,reff)
   case ( "d", "D")
      if (i > ndown) then
         call die("attempt to get n from non-existing Down potential")
      endif
      val = eval_radfunc(psxml%pot(i)%V,reff)
end select

if (hasRfactor)   val = val / reff

end function psxmlEvaluatePotential
!
function psxmlEvaluatePseudoOrbital(psxml,ud,i,r,debug) result(val)
type(xml_ps_t), intent(in) :: psxml
character, intent(in)      :: ud
integer,   intent(in)      :: i
real(dp),  intent(in)      :: r
logical, intent(in), optional :: debug
real(dp)                   :: val


integer :: ndown

ndown = psxml%npots_down

select case(ud)
   case ( "u", "U")
      if (i > psxmlPotentialsUp(psxml)) then
         call die("attempt to evaluate non-existing Up pswf")
      endif
      val = eval_radfunc(psxml%pswf(ndown+i)%V,r)
   case ( "d", "D")
      if (i > ndown) then
         call die("attempt to get n from non-existing Down pswf")
      endif
      val = eval_radfunc(psxml%pswf(i)%V,r)
end select

end function psxmlEvaluatePseudoOrbital

!====================================================
function eval_radfunc(f,r,debug) result(val)
type(radfunc_t), intent(in) :: f
real(dp), intent(in)      :: r
real(dp)                  :: val
logical, intent(in), optional :: debug

logical :: remove_r
real(dp), pointer :: x(:) => null(), y(:) => null()

x => f%grid%grid_data(:)
y => f%data(:)

call interpolate(x,y,r,val,debug)

end function eval_radfunc

subroutine interpolate(x,y,r,val,debug)
real(dp), intent(in) :: x(:), y(:)
real(dp), intent(in) :: r
real(dp), intent(out):: val
logical, intent(in), optional :: debug

integer, save :: i0
integer :: npts, nmin, nmax, nn
integer, parameter :: npoint = 2  ! basis for interpolation order
real(dp)  :: dy

npts = size(x)
if (size(y) /= npts) call die("x and y not conformable in interpolate")

call hunt(x,npts,r,i0)
nmin=max(1,i0-npoint)
nmax=min(npts,i0+npoint)
nn=nmax-nmin+1
call polint(x(nmin:),y(nmin:),nn,r,val,dy)
if (present(debug)) then
   if (debug) then
      print "(a,2g20.10,i4,2g20.10)", &
          "r ,r-x(i0), i0, val, d: ", r, r-x(i0), i0, val, (val-y(i0))
   endif
endif
end subroutine interpolate

      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL(dp), intent(in) ::  x,xx(n)

      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
    END SUBROUTINE hunt
!  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.



subroutine dump_pseudo(pseudo,lun)
  type(xml_ps_t), intent(in), target   :: pseudo
  integer, intent(in)                  :: lun

integer  :: i
type(vps_t), pointer :: pp
type(pswf_t), pointer :: pw
type(radfunc_t), pointer :: rp

real(dp), parameter :: rsmall = 1.e-3_dp

write(lun,*) "---PSEUDO data:"

      if (associated(pseudo%global_grid)) then
         write(lun,*) "global grid data: ",  &
           pseudo%global_grid%npts
      endif
         
do i = 1, pseudo%npots
      pp =>  pseudo%pot(i)
      rp =>  pseudo%pot(i)%V
      write(lun,*) "VPS ", i, " angular momentum: ", pp%l
      write(lun,*) "                 n: ", pp%n
      write(lun,*) "                 occupation: ", pp%occupation
      write(lun,*) "                 cutoff: ", pp%cutoff
      write(lun,*) "                 spin: ", pp%spin
      write(lun,*) "grid data: ", rp%grid%npts, rp%data(1)
      write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)/rsmall
enddo
do i = 1, pseudo%npswfs
      pw =>  pseudo%pswf(i)
      rp =>  pseudo%pswf(i)%V
      write(lun,*) "PSWF ", i, " angular momentum: ", pw%l
      write(lun,*) "                 n: ", pw%n
      write(lun,*) "                 spin: ", pw%spin
      write(lun,*) "grid data: ", rp%grid%npts, rp%data(1)
      write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)
enddo


rp => pseudo%valence_charge
write(lun,*) "grid data: ", rp%grid%npts, rp%data(1)
write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)

if (associated(pseudo%core_charge%data)) then
   write(lun,*) "grid data: ", rp%grid%npts, rp%data(1)
   write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)
endif

write(lun,*) "---------------------- Done with Dump_pseudo"

end subroutine dump_pseudo

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY) 
!*****************************************************************
! Polinomic interpolation. Modified and adapted to double 
! precision from same routine of Numerical Recipes.
! D. Sanchez-Portal, Oct. 1996
!*****************************************************************
! Input:
!   real*8  XA(N) : x values of the function y(x) to interpolate
!   real*8  YA(N) : y values of the function y(x) to interpolate
!   integer N     : Number of data points
!   real*8  X     : x value at which the interpolation is desired
! Output:
!   real*8  Y     : interpolated value of y(x) at X
!   real*8  DY    : accuracy estimate
!*****************************************************************

      IMPLICIT NONE
      INTEGER  :: N
      REAL(dp) :: XA(N),YA(N), X, Y, DY

      INTEGER  :: I, M, NS
      REAL(dp) :: C(N), D(N), DEN, DIF, DIFT, HO, HP, W
      REAL(dp), PARAMETER :: ZERO=0.D0

      NS=1
      DIF=ABS(X-XA(1))
      DO I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
      END DO ! I
      Y=YA(NS)
      NS=NS-1
      DO M=1,N-1
        DO I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF (DEN.EQ.ZERO) call die('polint: ERROR. Two XAs are equal')
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
        END DO ! I
        IF (2*NS.LT.N-M) THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
      END DO ! M

      END SUBROUTINE POLINT

!
      FUNCTION at_number(SYMBOL) result(z)

! Given the atomic symbol, it returns the atomic number
! Based on code by J. Soler

      character(len=2), intent(in)    :: SYMBOL  ! Atomic symbol
      integer                         :: Z       ! Atomic number

      integer, parameter  :: NZ=103
      character(len=2), parameter :: NAME(NZ) =  &
               (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
                 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
                 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
                 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
                 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
                 'Md','No','Lr'/)

     do z = 1, NZ
        if (SYMBOL == NAME(Z)) then
           RETURN
        endif
     enddo
     call die("Cannot find atomic number for " // symbol)
        
   end FUNCTION at_number

end module m_ncps_xml_ps_t




















