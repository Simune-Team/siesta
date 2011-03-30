module m_radfunc_registry
!
! This module provides a general interface to the
! whole set of radial functions needed in Siesta
! 
! The only restriction is that the representation
! of the functions conforms to the "radfunc" convention.
!
! Sketch of usage (see examples in the code):
!
! Each type of function is "registered", and a global
! index is obtained. The user is responsible for the correct
! bookeeping of the indexes (see routine register_rfs for an
! example appropriate for the orbitals, KB projectors, and Vna
! in Siesta, and routine "overlap").
!
! Once registered, a function can be evaluated, its cutoff and 
! angular momentum obtained, etc.  So far only the functions
! needed by matel are implemented: rcut, lofio, phiatm (renamed
! as "evaluate").
!
! The registry contains some meta-data for each function (its
! kind as a string, and the original indexes handled by the
! user at registration time).
!
  use precision, only: dp
  use radial, only: rad_func, rad_get
  use spher_harm, only:  rlylm
  use sys, only: die

  implicit none

  private

  type id_t
     character(len=10) :: func_type = "null"
     integer           :: n_indexes = 0
     integer           :: indexes(4) = (/-1,-1,-1,-1/)
  end type id_t

!
! Each function is represented by the radfunc record (for its data)
! and a meta-data section, which includes l, m, and id.
!
  type ext_radfunc_t
     type(rad_func), pointer :: func => null()
     integer        :: l
     integer        :: m
     type(id_t)     :: id
  end type ext_radfunc_t

!
! There should be a mechanism to make sure of the size
!
  integer, parameter :: nmax_funcs = 500
  type(ext_radfunc_t), dimension(nmax_funcs)  :: rf_pool

  integer, private    :: nfuncs = 0

  integer, parameter               :: max_l = 5
  integer, parameter               :: max_ilm = (max_l+1)*(max_l+1)

!  public :: get_number_of_funcs
  public :: register_in_rf_pool
  public :: evaluate, rcut, lofio
  public :: evaluate_x, evaluate_y, evaluate_z

  CONTAINS

    function get_number_of_funcs() result(n)
      integer :: n

      n = nfuncs
    end function get_number_of_funcs
      
    function valid(gindex) result (ok)
      integer, intent(in)    :: gindex
      logical                :: ok

      ok = (gindex >= 0 .AND. gindex <= nfuncs)
    end function valid

!
!   This is the main entry to the registry
!
    subroutine register_in_rf_pool(func,l,m,func_type,indexes,gindex)
     type(rad_func), pointer :: func             ! function data
     integer, intent(in)     :: l, m           
     character(len=*), intent(in) :: func_type   ! mnemonic kind
     integer, intent(in)          :: indexes(:)  ! legacy indexes
     integer, intent(out)    :: gindex           ! global index

     integer :: n_indexes 
     logical, external :: io_node

     n_indexes = size(indexes)

     nfuncs = nfuncs + 1
     if (nfuncs > nmax_funcs) call die("Overflow in registry")
     gindex = nfuncs
     
     rf_pool(gindex)%func => func
     rf_pool(gindex)%l = l
     rf_pool(gindex)%m = m

     rf_pool(gindex)%id%func_type = func_type
     rf_pool(gindex)%id%n_indexes = n_indexes
     rf_pool(gindex)%id%indexes(1:n_indexes) = indexes(:)

     if (io_node()) then
        print *, "Registered: ", trim(func_type), &
                              indexes(1:n_indexes), func%cutoff
     endif

   end subroutine register_in_rf_pool

    function is_vna(gindex) 
      integer, intent(in)    :: gindex
      logical                :: is_vna

      is_vna = (rf_pool(gindex)%id%func_type(1:3) == "vna")

    end function is_vna

   function rcut(gindex) result(cutoff)
     integer, intent(in)    :: gindex
     real(dp)               :: cutoff

     if (valid(gindex)) then
        cutoff = rf_pool(gindex)%func%cutoff
     else
        call die("Invalid gindex")
     endif
   end function rcut

   function lofio(gindex) result(l)
     integer, intent(in)    :: gindex
     integer                :: l

     if (valid(gindex)) then
        l = rf_pool(gindex)%l
     else
        call die("Invalid gindex")
     endif
   end function lofio

   subroutine evaluate(gindex,r,f,grad)
     integer, intent(in)    :: gindex
     real(dp), intent(in)   :: r(3)
     real(dp), intent(out)  :: f
     real(dp), intent(out)  :: grad(3)

     type(rad_func), pointer :: func
     real(dp) rmod, phir, dphidr
     real(dp) rly(max_ilm), grly(3,max_ilm)
     integer i, l, m, ilm

     real(dp), parameter     :: tiny = 1.e-20_dp

     if (valid(gindex)) then
        
        func => rf_pool(gindex)%func
        l = rf_pool(gindex)%l
        m = rf_pool(gindex)%m

        ! the addition of "tiny" avoids division by 0
        !
        rmod = sqrt(sum(r*r)) + tiny
        
        if(rmod > func%cutoff) then

           f = 0.0_dp
           grad(1:3) = 0.0_dp

        else

           call rad_get(func,rmod,phir,dphidr)

           if (is_vna(gindex)) then
              f=phir
              grad(1:3)=dphidr*r(1:3)/rmod
           else
              ilm = l*l + l + m + 1
              call rlylm( l, r, rly, grly )
              f = phir * rly(ilm)
              do i = 1,3
                 grad(i)=dphidr*rly(ilm)*r(i)/rmod+phir*grly(i,ilm)
              enddo
           endif
        endif
     else
        call die("Invalid gindex")
     endif
   end subroutine evaluate

      SUBROUTINE evaluate_x(gindex,r,xphi,grxphi)
!      Calculates x*phiatm and its gradient

      integer, intent(in)   :: gindex
      real(dp), intent(in)  :: r(3)
      real(dp), intent(out) :: xphi, grxphi(3)

      real(dp) phi, grphi(3), x

      call evaluate(gindex,r,phi,grphi)
      x = r(1)
      xphi = x * phi
      grxphi(1) = x * grphi(1) + phi
      grxphi(2) = x * grphi(2)
      grxphi(3) = x * grphi(3)
      END SUBROUTINE evaluate_x

      SUBROUTINE evaluate_y(gindex,r,yphi,gryphi)
!     Calculates y*phiatm and its gradient

      integer, intent(in)   :: gindex
      real(dp), intent(in)  :: r(3)
      real(dp), intent(out) :: yphi, gryphi(3)

      real(dp) phi, grphi(3), y

      call evaluate(gindex,r,phi,grphi)
      y = r(2)
      yphi = y * phi
      gryphi(1) = y * grphi(1)
      gryphi(2) = y * grphi(2) + phi
      gryphi(3) = y * grphi(3)
      END SUBROUTINE evaluate_y

      SUBROUTINE evaluate_z(gindex,r,zphi,grzphi)
!     Calculates z*phiatm and its gradient

      integer, intent(in)   :: gindex
      real(dp), intent(in)  :: r(3)
      real(dp), intent(out) :: zphi, grzphi(3)

      real(dp) phi, grphi(3), z
      call evaluate(gindex,r,phi,grphi)
      z = r(3)
      zphi = z * phi
      grzphi(1) = z * grphi(1)
      grzphi(2) = z * grphi(2)
      grzphi(3) = z * grphi(3) + phi
      END SUBROUTINE evaluate_z

 end module m_radfunc_registry
   
