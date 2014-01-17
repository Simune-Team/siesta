module m_ncps_xml_ps_t
!
! Data structures to match the XML pseudopotential format
!
integer, parameter, private    :: MAXN_POTS = 8
integer, parameter, private    :: dp = selected_real_kind(14)
!
public  :: dump_pseudo, eval_radfunc
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
      type(grid_t), pointer              :: global_grid => null()
      type(vps_t), dimension(MAXN_POTS)  :: pot
      type(pswf_t), dimension(MAXN_POTS) :: pswf
      type(radfunc_t)                    :: core_charge
      type(radfunc_t)                    :: valence_charge
   end type xml_ps_t


CONTAINS !===============================================

function eval_radfunc(f,r,remove_rfactor) result(val)
type(radfunc_t), intent(in) :: f
real(dp), intent(in)      :: r
real(dp)                  :: val
logical, intent(in), optional :: remove_rfactor

logical :: remove_r
real(dp), pointer :: x(:) => null(), y(:) => null()

x => f%grid%grid_data(:)
y => f%data(:)

remove_r = .false.
if (present(remove_rfactor)) then
   remove_r = remove_rfactor
endif

call interpolate(x,y,r,val)

end function eval_radfunc

subroutine interpolate(x,y,r,val)
real(dp), intent(in) :: x(:), y(:)
real(dp), intent(in) :: r
real(dp), intent(out):: val

integer, save :: i0
integer :: npts, nmin, nmax, nn
integer, parameter :: npoint = 2  ! interpolation order
real(dp)  :: dy

npts = size(x)
if (size(y) /= npts) call die("x and y not conformable in interpolate")

call hunt(x,npts,r,i0)
nmin=max(1,i0-npoint)
nmax=min(npts,i0+npoint)
nn=nmax-nmin+1
call polint(x(nmin:),y(nmin:),nn,r,val,dy)

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
           pseudo%global_grid%npts, pseudo%global_grid%scale
      endif
         
do i = 1, pseudo%npots
      pp =>  pseudo%pot(i)
      rp =>  pseudo%pot(i)%V
      write(lun,*) "VPS ", i, " angular momentum: ", pp%l
      write(lun,*) "                 n: ", pp%n
      write(lun,*) "                 occupation: ", pp%occupation
      write(lun,*) "                 cutoff: ", pp%cutoff
      write(lun,*) "                 spin: ", pp%spin
      write(lun,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
      write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)/rsmall
enddo
do i = 1, pseudo%npswfs
      pw =>  pseudo%pswf(i)
      rp =>  pseudo%pswf(i)%V
      write(lun,*) "PSWF ", i, " angular momentum: ", pw%l
      write(lun,*) "                 n: ", pw%n
      write(lun,*) "                 spin: ", pw%spin
      write(lun,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
      write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)
enddo
rp => pseudo%valence_charge
write(lun,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)
rp => pseudo%core_charge
write(lun,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
write(lun,*) "value at r=0: ", eval_radfunc(rp,rsmall)

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

end module m_ncps_xml_ps_t




















