!> @brief Data structures and functions to handle
!> the PSML pseudopotential format.
!> @author Alberto Garcia
!
!> @detail The PSML library will eventually have three sections:
!>
!> A. Definition of data structures.
!> B. Accessors for pseudopotential information and evaluation
!> C. Parsing helpers

module m_ncps_xml_ps_t
! 
! In order for this module will hide completely the implementation of
! the parsing of the XML file and the later querying of the data
! structures, the ps_t type should be made private. But then all the
! accessors and the parsing helpers should reside in this module too, making
! it quite big.
!
! Currently the parsing helpers are in 'm_ncps_parsing_helpers'
! This scattering will be solved with the final packaging of the library.
!
implicit none

integer, parameter, private    :: MAXN_POTS = 8
integer, parameter, private    :: MAXN_SHELLS = 8
integer, parameter, private    :: MAXN_WFNS = 8
integer, parameter, private    :: dp = selected_real_kind(14)
!
!-----------------------------------------------------------
type, public :: provenance_t
        character(len=40)       :: creator
        character(len=30)       :: date
end type provenance_t
!------
type, public :: header_t
        character(len=30)       :: atomic_label    ! generalized symbol
        real(kind=dp)           :: z  ! atomic number (might be non-integer)
        real(kind=dp)           :: zpseudo ! Z - ncore-electrons
        character(len=50)       :: flavor
        logical                 :: relativistic
        logical                 :: polarized
        character(len=50)       :: xc_libxc_exchange ! LibXC string
        character(len=50)       :: xc_libxc_correlation ! LibXC string
        ! Generator's own terminology for XC 
        ! Could be something like "type:LDA//authors:PZ"
        character(len=80)       :: xc_functional
        !
        character(len=4)        :: core_corrections
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
type, public :: grid_t
!
!     Note that the preferred option is to have explicit grid data.
!
      integer                        :: npts = 0
      real(dp), pointer              :: grid_data(:) => null()
      character(len=80)              :: annotation
end type grid_t      
!
type, public :: radfunc_t
      type(grid_t), pointer                   :: grid => null()
      real(kind=dp), dimension(:), pointer    :: data => null()
end type radfunc_t      
!
type, public :: semilocal_t
      integer                          :: npots
      integer                          :: npots_major
      integer                          :: npots_minor
      integer, dimension(MAXN_POTS)           :: n
      character(len=1), dimension(MAXN_POTS)  :: l
      character(len=5), dimension(MAXN_POTS)  :: set
      character(len=40), dimension(MAXN_POTS) :: flavor
      real(dp), dimension(MAXN_POTS)          :: rc
      type(radfunc_t), dimension(MAXN_POTS)   :: V
      !
      ! indexes
      integer, dimension(MAXN_POTS)           :: major
      integer, dimension(MAXN_POTS)           :: minor
end type semilocal_t

type, public :: pswfs_t
      character(len=40)                :: format
      integer                          :: npswfs
      integer                          :: npswfs_major
      integer                          :: npswfs_minor
      integer, dimension(MAXN_WFNS)           :: n
      character(len=1), dimension(MAXN_WFNS)  :: l
      character(len=5), dimension(MAXN_WFNS)  :: set
      type(radfunc_t), dimension(MAXN_WFNS)   :: Phi
      !
      ! indexes
      integer, dimension(MAXN_WFNS)           :: major
      integer, dimension(MAXN_WFNS)           :: minor
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
      type(grid_t), pointer              :: global_grid => null()
      type(semilocal_t)                  :: semilocal
      type(pswfs_t)                      :: pswfs
      !
      type(valence_charge_t)             :: valence_charge
      type(core_charge_t)                :: core_charge

   end type ps_t

!
public  :: ps_destroy
!public  :: dump_pseudo           ! For initial debugging
!

! Accessor list
public :: ps_AtomicSymbol
public :: ps_AtomicLabel
public :: ps_AtomicNumber
public :: ps_ZPseudo
public :: ps_GenerationZval

public :: ps_Creator
public :: ps_Date

public :: ps_PseudoFlavor
public :: ps_XCLibXCExchange
public :: ps_XCLibXCCorrelation
public :: ps_XCFunctional

public :: ps_IsRelativistic
public :: ps_IsSpinPolarized
public :: ps_HasCoreCorrections

public :: ps_NValenceShells
public :: ps_ValenceShellL
public :: ps_ValenceShellN
public :: ps_ValenceShellOccupation

public :: ps_NPotentials
public :: ps_PotentialL
public :: ps_PotentialN
public :: ps_PotentialRc

public :: ps_NPseudoWfs
public :: ps_PseudoWfL
public :: ps_PseudoWfN

public :: ps_EvaluatePotential
public :: ps_EvaluatePseudoWf
public :: ps_EvaluateValenceCharge
public :: ps_EvaluateCoreCharge

CONTAINS !===============================================
!> @brief Cleans and deallocates the ps object
!> @author Alberto Garcia
!> @date March-July 2014
!> @detail ps is a pointer, and not a "value". This is confusing
subroutine ps_destroy(ps)
type(ps_t), pointer :: ps

integer :: i

if (.not. associated(ps)) RETURN

if (associated(ps%global_grid)) then
   call destroy_grid(ps%global_grid)
endif
do i = 1, ps%semilocal%npots
   call destroy_radfunc(ps%semilocal%V(i))
enddo
do i = 1, ps%pswfs%npswfs
   call destroy_radfunc(ps%pswfs%Phi(i))
enddo
call destroy_radfunc(ps%valence_charge%rho_val)
call destroy_radfunc(ps%core_charge%rho_core)

deallocate(ps)
ps => null()

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

!========================================================

!> @brief Returns the number of non-empty valence shells
!> in the ps generation configuration
!> @author Alberto Garcia
!> @date 2014
!> @param ps is a handle to the psml information
function ps_NValenceShells(ps) result(nshells)
type(ps_t), intent(in) :: ps
integer                :: nshells

nshells = ps%config_val%nshells

end function ps_NValenceShells

!> @brief Returns the angular momentum of the i'th valence shell
!> in the ps generation configuration
!> @author Alberto Garcia
!> @date 2014
!> @param ps is a handle to the psml information
!> @param i is the index of the shell
!> @note i should be within range
function ps_ValenceShellL(ps,i) result(l)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
integer                :: l

character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)
character(len=1) :: str

if (i > ps%config_val%nshells) then
   call die("Attempt to get l for non-existing shell")
endif
str = ps%config_val%l(i)
!
do l = 0,4
   if (str == sym(l)) RETURN
enddo
call die("Wrong l symbol in valence shell")

end function ps_ValenceShellL

!> @brief Returns the principal quantum number of the i'th valence shell
!> in the ps generation configuration
!> @author Alberto Garcia
!> @date 2014
!> @param ps is a handle to the psml information
!> @param i is the index of the shell
!> @note i should be within range
function ps_ValenceShellN(ps,i) result(n)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
integer                :: n

if (i > ps%config_val%nshells) then
   call die("Attempt to get n for non-existing shell")
endif
n = ps%config_val%n(i)

end function ps_ValenceShellN

!> @brief Returns the occupation of the i'th valence shell
!> in the ps generation configuration
!> @author Alberto Garcia
!> @date 2014
!> @param ps is a handle to the psml information
!> @param i is the index of the shell
!> @param channel is an optional parameter for spin-polarized 
!> calculations ("u" or "d"). 
!> @note i should be within range
!> @note If "channel" is present, the occupation returned
!> corresponds to the given channel
function ps_ValenceShellOccupation(ps,i,channel) result(occ)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
character(len=1), intent(in), optional :: channel
real(dp)                :: occ

if (i > ps%config_val%nshells) then
   call die("Attempt to get occupation for non-existing shell")
endif
if (present(channel)) then
   if (ps_IsSpinPolarized(ps)) then
      if (channel == "u") then
         occ = ps%config_val%occ_up(i)
      else if (channel == "d") then
         occ = ps%config_val%occ_down(i)
      else
         call die("Wrong channel in ValShellOccupation")
      endif
   else
      call die("Cannot speficy channel in ValShellOccupation")
   endif
else
   occ = ps%config_val%occ(i)
endif

end function ps_ValenceShellOccupation
!
!-------------------------------------------------------
function ps_EvaluateValenceCharge(ps,r,debug) result(val)
type(ps_t), intent(in) :: ps
real(dp), intent(in)       :: r
logical, intent(in), optional :: debug
real(dp)                   :: val

if (r > max_range(ps%valence_charge%rho_val)) then
   val = 0.0_dp
else
   val = eval_radfunc(ps%valence_charge%rho_val,r,debug)
endif
end function ps_EvaluateValenceCharge

function ps_EvaluateCoreCharge(ps,r,debug) result(val)
type(ps_t), intent(in) :: ps
real(dp), intent(in)       :: r
logical, intent(in), optional :: debug
real(dp)                   :: val

if (r > max_range(ps%core_charge%rho_core)) then
   val = 0.0_dp
else
   val = eval_radfunc(ps%core_charge%rho_core,r,debug)
endif
end function ps_EvaluateCoreCharge
!
function ps_AtomicSymbol(ps) result(name)
type(ps_t), intent(in) :: ps
character(len=2) :: name
name = ps%header%atomic_label(1:2)
end function ps_AtomicSymbol
!
function ps_AtomicLabel(ps) result(name)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%header%atomic_label)) :: name
name = trim(ps%header%atomic_label)
end function ps_AtomicLabel
!
function ps_AtomicNumber(ps) result(z)
type(ps_t), intent(in) :: ps
real(dp) :: z
 z = ps%header%z
end function ps_AtomicNumber
!
function ps_Creator(ps) result(name)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%provenance%creator)) :: name
name = trim(ps%provenance%creator)
end function ps_Creator
!
function ps_Date(ps) result(str)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%provenance%date)) :: str
str = trim(ps%provenance%date)
end function ps_Date
!
!**AG** Generalize to accept an index argument for potential
function ps_PseudoFlavor(ps) result(str)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%header%flavor)) :: str
str = trim(ps%header%flavor)
end function ps_PseudoFlavor
!
function ps_ZPseudo(ps) result(zpseudo)
type(ps_t), intent(in) :: ps
real(dp)                   :: zpseudo
zpseudo = ps%header%zpseudo
end function ps_ZPseudo
!
function ps_GenerationZval(ps) result(zval)
type(ps_t), intent(in) :: ps
real(dp)                   :: zval
zval = ps%config_val%total_charge
end function ps_GenerationZval
!
function ps_XCLibXCExchange(ps) result(xc_string)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%header%xc_libxc_exchange)) :: xc_string
xc_string = trim(ps%header%xc_libxc_exchange)
end function ps_XCLibXCExchange
!
function ps_XCLibXCCorrelation(ps) result(xc_string)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%header%xc_libxc_correlation)) :: xc_string
xc_string = trim(ps%header%xc_libxc_correlation)
end function ps_XCLibXCCorrelation
!
function ps_XCFunctional(ps) result(xc_string)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%header%xc_functional)) :: xc_string
xc_string = trim(ps%header%xc_functional)
end function ps_XCFunctional
!
function ps_IsRelativistic(ps) result(rel)
type(ps_t), intent(in) :: ps
logical                    :: rel
rel = ps%header%relativistic
end function ps_IsRelativistic
!
function ps_IsSpinPolarized(ps) result(pol)
type(ps_t), intent(in) :: ps
logical                    :: pol
pol = ps%header%polarized
end function ps_IsSpinPolarized
!
function ps_HasCoreCorrections(ps) result(cc)
type(ps_t), intent(in) :: ps
logical                    :: cc
cc = (ps%header%core_corrections == "yes")
end function ps_HasCoreCorrections
!
function ps_NPotentials(ps,set) result(n)
type(ps_t), intent(in) :: ps
character(len=5), intent(in), optional :: set
integer                    :: n

logical major

major = .true.

if (present(set)) then
   if (set == "minor") then
      major = .false.
   endif
endif

if (major) then
   n = ps%semilocal%npots_major
else
   n = ps%semilocal%npots_minor
endif
   
end function ps_NPotentials
!
function ps_PotentialL(ps,i,set) result(l)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
character(len=5), intent(in), optional :: set
integer                    :: l

character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)
character(len=1) :: str

integer :: idx

idx = ps_GetPotentialIndex(ps,i,set)
str = ps%semilocal%l(idx)
!
do l = 0,4
   if (str == sym(l)) RETURN
enddo
call die("Wrong l symbol in potential")

end function ps_PotentialL
!
function ps_PotentialRc(ps,i,set) result(rc)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
character(len=5), intent(in), optional :: set
real(dp)                   :: rc

integer :: idx

idx = ps_GetPotentialIndex(ps,i,set)
rc = ps%semilocal%rc(idx)

end function ps_PotentialRc
!
function ps_PotentialN(ps,i,set) result(n)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
character(len=5), intent(in), optional :: set
integer                    :: n

integer :: idx

idx = ps_GetPotentialIndex(ps,i,set)
n = ps%semilocal%n(idx)

end function ps_PotentialN
!
function ps_EvaluatePotential(ps,i,r,set,debug) result(val)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
real(dp),  intent(in)      :: r
character(len=5), intent(in), optional :: set
logical, intent(in), optional :: debug
real(dp)                   :: val

integer :: idx

idx = ps_GetPotentialIndex(ps,i,set)
if (r> max_range(ps%semilocal%V(idx))) then
   val = - ps_ZPseudo(ps)/r
else
   val = eval_radfunc(ps%semilocal%V(idx),r, debug)
endif

end function ps_EvaluatePotential

function ps_GetPotentialIndex(ps,i,set) result(idx)
type(ps_t), intent(in)                 :: ps
integer,   intent(in)                  :: i
character(len=5), intent(in), optional :: set
integer                                :: idx

logical major

if (i > ps_NPotentials(ps,set)) then
     call die("attempt to get index for non-existing potential")
endif

major = .true.
if (present(set)) then
   if (set == "minor") then
      major = .false.
   endif
endif

if (major) then
   idx = ps%semilocal%major(i)
else
   idx = ps%semilocal%minor(i)
endif
end function ps_GetPotentialIndex

function ps_NPseudoWfs(ps,set) result(n)
type(ps_t), intent(in) :: ps
character(len=5), intent(in), optional :: set
integer                    :: n

logical major

major = .true.

if (present(set)) then
   if (set == "minor") then
      major = .false.
   endif
endif

if (major) then
   n = ps%pswfs%npswfs_major
else
   n = ps%pswfs%npswfs_minor
endif
   
end function ps_NPseudoWfs

function ps_GetPswfIndex(ps,i,set) result(idx)
type(ps_t), intent(in)                 :: ps
integer,   intent(in)                  :: i
character(len=5), intent(in), optional :: set
integer                                :: idx

logical major

if (i > ps_NPseudoWfs(ps,set)) then
     call die("attempt to get index for non-existing pswf")
endif

major = .true.
if (present(set)) then
   if (set == "minor") then
      major = .false.
   endif
endif

if (major) then
   idx = ps%pswfs%major(i)
else
   idx = ps%pswfs%minor(i)
endif
end function ps_GetPswfIndex
!
function ps_PseudoWfL(ps,i,set) result(l)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
character(len=5), intent(in), optional :: set
integer                    :: l

character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)
character(len=1) :: str

integer :: idx

idx = ps_GetPswfIndex(ps,i,set)
str = ps%pswfs%l(idx)
!
do l = 0,4
   if (str == sym(l)) RETURN
enddo
call die("Wrong l symbol in pswf")

end function ps_PseudoWfL
!
function ps_PseudoWfN(ps,i,set) result(n)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
character(len=5), intent(in), optional :: set
integer                    :: n

integer :: idx

idx = ps_GetPswfIndex(ps,i,set)
n = ps%pswfs%n(idx)

end function ps_PseudoWfN
!
function ps_EvaluatePseudoWf(ps,i,r,set,debug) result(val)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
real(dp),  intent(in)      :: r
character(len=5), intent(in), optional :: set
logical, intent(in), optional :: debug
real(dp)                   :: val

integer :: idx

idx = ps_GetPswfIndex(ps,i,set)
if (r> max_range(ps%pswfs%Phi(idx))) then
   val = 0.0_dp
else
   val = eval_radfunc(ps%pswfs%Phi(idx),r,debug)
endif

end function ps_EvaluatePseudoWf

!====================================================
!> @brief Maximum radius in a radfunc's grid
function max_range(f) result(range)
type(radfunc_t), intent(in) :: f
real(dp)                  :: range

integer :: npts

npts = f%grid%npts
range = f%grid%grid_data(npts)
end function max_range
!----------
function eval_radfunc(f,r,debug) result(val)
type(radfunc_t), intent(in) :: f
real(dp), intent(in)      :: r
real(dp)                  :: val
logical, intent(in), optional :: debug

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



!!$subroutine dump_pseudo(pseudo,lun)
!!$  type(ps_t), intent(in), target   :: pseudo
!!$  integer, intent(in)                  :: lun
!!$
!!$integer  :: i
!!$type(vps_t), pointer :: pp
!!$type(pswf_t), pointer :: pw
!!$type(radfunc_t), pointer :: rp
!!$
!!$real(dp), parameter :: rsmall = 1.e-3_dp
!!$
!!$write(lun,*) "---PSEUDO data:"
!!$
!!$      if (associated(pseudo%global_grid)) then
!!$         write(lun,*) "global grid data: ",  &
!!$           pseudo%global_grid%npts
!!$      endif
!!$         
!!$do i = 1, pseudo%semilocal%npots
!!$      pp =>  pseudo%semilocal
!!$      rp =>  pseudo%semilocal%V(i)
!!$      write(lun,*) "VPS ", i, " angular momentum: ", pp%l
!!$      write(lun,*) "                 n: ", pp%n
!!$      write(lun,*) "                 occupation: ", pp%occupation
!!$      write(lun,*) "                 rc: ", pp%rc
!!$      write(lun,*) "                 spin: ", pp%spin
!!$      write(lun,*) "grid data: ", rp%grid%npts, rp%data(1)
!!$      write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)/rsmall
!!$enddo
!!$do i = 1, pseudo%npswfs
!!$      pw =>  pseudo%pswf(i)
!!$      rp =>  pseudo%pswf(i)%V
!!$      write(lun,*) "PSWF ", i, " angular momentum: ", pw%l
!!$      write(lun,*) "                 n: ", pw%n
!!$      write(lun,*) "                 spin: ", pw%spin
!!$      write(lun,*) "grid data: ", rp%grid%npts, rp%data(1)
!!$      write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)
!!$enddo
!!$
!!$
!!$rp => pseudo%valence_charge
!!$write(lun,*) "grid data: ", rp%grid%npts, rp%data(1)
!!$write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)
!!$
!!$if (associated(pseudo%core_charge%data)) then
!!$   write(lun,*) "grid data: ", rp%grid%npts, rp%data(1)
!!$   write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)
!!$endif
!!$
!!$write(lun,*) "---------------------- Done with Dump_pseudo"
!!$
!!$end subroutine dump_pseudo

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




















