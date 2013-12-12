!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!

module m_ts_contour
!
! Routines that are used to setup the contour for integration of the GFs
! 
! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_cctype
  use m_ts_io_contour
  use precision, only : dp

  implicit none

  public :: read_contour_options
  public :: print_contour_options
  public :: print_contour_block
  public :: io_contour
  public :: sort_contour

  integer, public :: PNEn, NEn, NEn_eq,N_transport
  type(ts_ccontour), public, pointer :: contour(:)!, contour_Eq(:), contour_nEq(:)

  public :: contour_Eq, contour_EqL, contour_EqR, contour_nEq, contour_Transport
  private

contains

  function contour_Eq() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour_EqL()
  end function contour_Eq

  function contour_EqL() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(1:NEn_eq)
  end function contour_EqL

  function contour_EqR() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(NEn_eq+1:NEn_eq*2)
  end function contour_EqR

  function contour_nEq() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(NEn_eq*2+1:NEn-N_transport)
  end function contour_nEq

  function contour_Transport() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(NEn-N_transport+1:NEn)
  end function contour_Transport

  subroutine io_contour(IsVolt, Elecs,slabel,suffix)
    use m_ts_electype
    use m_ts_contour_eq, only : io_contour_eq
    use m_ts_contour_neq, only : io_contour_neq
    logical, intent(in) :: IsVolt
    type(Elec), intent(in) :: Elecs(:)
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

    call io_contour_eq(Elecs,slabel,suffix)

    if ( IsVolt ) then
       call io_contour_neq(slabel,suffix)
    end if

  end subroutine io_contour

  subroutine read_contour_options(Elecs, kT, IsVolt, Volt)

    use m_ts_electype
    use m_ts_contour_eq, only : read_contour_eq_options
    use m_ts_contour_neq, only : read_contour_neq_options

    type(Elec), intent(inout) :: Elecs(:)
    logical, intent(in) :: IsVolt
    real(dp), intent(in) :: kT, Volt

    call read_contour_eq_options(Elecs,kT,Volt)

    if ( IsVolt ) then
       call read_contour_neq_options(Elecs,kT,Volt)
    end if

  end subroutine read_contour_options

  subroutine print_contour_options(prefix,IsVolt)

    use m_ts_contour_eq, only : print_contour_eq_options
    use m_ts_contour_neq, only : print_contour_neq_options

    character(len=*), intent(in) :: prefix
    logical, intent(in) :: IsVolt

    call print_contour_eq_options(prefix)

    if ( IsVolt ) then
       call print_contour_neq_options(prefix)
    end if

  end subroutine print_contour_options

  subroutine print_contour_block(prefix, IsVolt)

    use parallel, only : IONode
    use m_ts_contour_eq, only : print_contour_eq_block
    use m_ts_contour_neq, only : print_contour_neq_block

    character(len=*), intent(in) :: prefix
    logical, intent(in) :: IsVolt

    ! Print out the contour blocks
    if ( IONode ) then
       write(*,'(/,a)') 'transiesta: contour input as perceived:'
       write(*,'(/,a)') 'transiesta: Equilibrium contours'
    end if
    
    call print_contour_eq_block(prefix)

    if ( IsVolt ) then
       if ( IONode ) then
          write(*,'(/,a)') 'transiesta: non-Equilibrium contours'
       end if
       call print_contour_neq_block(prefix)
    end if
    if ( IONode ) write(*,'(/)')

  end subroutine print_contour_block

  ! In order to retain the best numerical accuracy we introduce
  ! an ordering of the energy contour by weights.
  ! It is the only reasonable thing to access as the functional is
  ! non-deterministic.
  ! An example of the importance of this sorting:
  ! Take +200 voltage energy points. 
  ! The weights will be something like this:
  ! 1e-5,...,1e-4,...,1e-3,...,1e-4,...,1e-5
  ! Which means that the summation up till the half works great.
  ! But when we reach the last weight we could be in the situation where:
  !  1._dp + 1e-12_dp which will limit the accuracy obtained for that energy point.

  ! In order to circumvent this we simply sort by weights.
  subroutine sort_contour(NC,ce,cw)
    integer, intent(in) :: NC
    complex(dp), intent(inout) :: ce(NC), cw(NC)
    
    ! Local variables
    real(dp) :: tmp
    complex(dp) :: ctmp
    integer :: i,j
    ! As we will only do this once we dont need a fancy sorting
    ! algorithm...
    do i = 1 , NC - 1
       tmp = real(cw(i)*conjg(cw(i)),dp)
       do j = i+1, NC
          if ( tmp > real(cw(j)*conjg(cw(j))) ) then
             ctmp  = ce(i)
             ce(i) = ce(j)
             ce(j) = ctmp
             ctmp  = cw(i)
             cw(i) = cw(j)
             cw(j) = ctmp
          end if
       end do
    end do

  end subroutine sort_contour

end module m_ts_contour
