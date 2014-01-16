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
! Nick Papior Andersen, 2013, nickpapior@gmail.com
!

module m_ts_contour_eq

  use precision, only : dp
!
! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_electype

  use m_ts_cctype
  use m_ts_chem_pot
  use m_ts_io_contour
  use m_ts_io_ctype
  use m_ts_aux

  implicit none

  ! equilibrium contour IO-segments
  integer, save, public :: N_Eq
  type(ts_c_io), pointer, save, public :: Eq_io(:) => null()
  type(ts_cw)  , pointer, save, public :: Eq_c(:) => null()

  ! The contour specific variables
  real(dp), save, public :: Eq_Eta

  ! We need to retain the information about the contour here.
  ! It provides an easier overview as there are quite a few constants governing the
  ! methods.
  integer, save, public :: N_poles

  ! The contours for the equilibrium density are attributed a fruitful discussion with
  ! Hans Skriver. Previously the routine names reflected his contribution.
  ! However, for programming clarity we have employed a different naming scheme.
  ! In order to retain the contributions it is encouraged to keep this sentiment for his
  ! contributions.

  ! For further attributions see the original paper by Brandbyge, et. al, 2002: DOI: 10.1103/PhysRevB.65.165401

  ! this is heavily linked with CONTOUR_NEQ from m_ts_contour_neq
  integer, parameter, public :: CONTOUR_EQ = 1

  public :: read_contour_eq_options
  public :: print_contour_eq_options
  public :: print_contour_eq_block
  public :: io_contour_eq
  public :: N_Eq_E
  public :: Eq_E
  public :: c2weight_eq
  public :: ID2idx

  interface ID2idx
     module procedure ID2idx_cw
     module procedure ID2idx_cidx
  end interface

  private

contains

  subroutine read_contour_eq_options(N_mu, mus, kT, Volt)

    use units, only : eV
    use parallel, only : IONode, Nodes, operator(.parcount.)
    use fdf

    integer, intent(in)        :: N_mu
    type(ts_mu), intent(inout) :: mus(N_mu)
    real(dp), intent(in)       :: kT, Volt

    integer :: i,j,k
    integer :: different_poles, N
    character(len=C_N_NAME_LEN), allocatable :: tmp(:), nContours(:)
    integer :: cur, next, prev
    integer, allocatable :: ID(:)
    logical :: isCircle
    logical :: isTail

    call fdf_obsolete('TS.ComplexContour.NPoles')

    ! broadening
    Eq_Eta = fdf_get('TS.Contours.Eq.Eta',0._dp,'Ry')
    if ( Eq_Eta < 0._dp) call die('ERROR: Eq_Eta < 0, we do not allow &
         &for using the advanced Greens function, please correct.')
    
    ! We only allow the user to either use the old input format, or the new
    ! per-electrode input

    ! Read in the generic things about the contours...
    N_poles = fdf_get('TS.Contours.Eq.Pole.N',6)

    do i = 1 , N_mu
       write(mus(i)%Eq_seg(Eq_segs(mus(i))),'(a,i0)') 'pole',i
    end do

    ! Count the number of contour segments
    N = minval(Eq_segs(mus))
    if ( N == 0 ) then
       call die('You have not assigned any contours for the chemical &
            &potentials.')
    end if
    N = sum(Eq_segs(mus))

    ! collect all equilibrium names
    allocate(tmp(N))
    tmp(:) = ' '
    N = 0
    do i = 1 , N_mu
       do j = 1 , Eq_segs(mus(i)) - 1
          N = N + 1
          tmp(N) = mus(i)%Eq_seg(j)
       end do
    end do
    ! we need to have the pole contours as the last ones
    ! to be able to easily differ them from the read in values
    do i = 1 , N_mu
       j = Eq_segs(mus(i))
       N = N + 1
       tmp(N) = mus(i)%Eq_seg(j)
    end do

    ! find all unique equilibrium names
    N = 1
    uniq_names: do i = 2 , size(tmp)
       j = 0
       do 
          j = j + 1
          if ( i <= j ) exit
          if ( tmp(i) .eq. tmp(j) ) cycle uniq_names
       end do
       N = N + 1
    end do uniq_names

    allocate(nContours(N)) ! unique names
    nContours = ' '

    ! populate unique equilibrium names
    N = 0
    uniq_names_pop: do i = 1 , size(tmp)
       do j = 1 , N
          if ( tmp(i) .eq. nContours(j) ) cycle uniq_names_pop
       end do
       N = N + 1
       if ( N > size(nContours) ) call die('Unique contour setup went wrong, contact devs')
       nContours(N) = tmp(i)
    end do uniq_names_pop
    if ( N /= size(nContours) ) then
       call die('ERROR: We have not populated enough contours')
    end if
    deallocate(tmp)

    ! Allocate space for the equilibrium contours
    nullify(Eq_io,Eq_c)
    N_Eq = N
    allocate(Eq_io(N_Eq)) ! the equilibrium contour io container
    allocate(Eq_c(N_Eq))  ! the equilibrium contour container
    Eq_io(:)%type = 'eq'  ! set the equilibrium tag

    ! Attach all the contours
    do i = 1 , N
       Eq_c(i)%c_io => Eq_io(i)
    end do

    ! read in the equilibrium contours (the last will be the poles)
    do i = 1 , N - N_mu

       print *,'reading ',nContours(i) ! TODO DELETE
       ! read in the contour
       call ts_read_contour_block('TS','',nContours(i),Eq_io(i), kT, Volt)

    end do

    ! We here create the "fake" pole contours
    j = 1
    do i = N - N_mu + 1, N
       ! assign name to the Eq_io
       Eq_io(i)%name   = nContours(i)
       Eq_io(i)%N      = N_poles
       Eq_io(i)%part   = 'pole'
       Eq_io(i)%method = 'residual'
       Eq_io(i)%a = mus(j)%mu
       Eq_io(i)%b = mus(j)%mu
       Eq_io(i)%d = mus(j)%mu
       j = j + 1
           
    end do

    ! fix the read-in constants
    do i = 1 , N_mu

       j = 1
       cur = get_c_io_index(mus(i)%Eq_seg(j))
       isCircle = Eq_io(cur)%type == 'circle'
       isTail   = Eq_io(cur)%type == 'tail'
       if ( Eq_segs(mus(i)) > 2 ) then
          next = get_c_io_index(mus(i)%Eq_seg(j+1))
          call ts_fix_contour( Eq_io(cur), next=Eq_io(next) )
       end if
          
       ! we should not check the pole (hence minus 2)
       do j = 2 , Eq_segs(mus(i)) - 2
          prev = get_c_io_index(mus(i)%Eq_seg(j-1))
          cur  = get_c_io_index(mus(i)%Eq_seg(j  ))
          next = get_c_io_index(mus(i)%Eq_seg(j+1))
          call ts_fix_contour( Eq_io(cur), &
               prev=Eq_io(prev), next=Eq_io(next) )
          call consecutive_types(Eq_io(cur),isCircle,isTail)
       end do

       ! don't check the pole
       j = Eq_segs(mus(i)) - 1
       cur = get_c_io_index(mus(i)%Eq_seg(j))
       if ( Eq_segs(mus(i)) > 2 ) then
          prev = get_c_io_index(mus(i)%Eq_seg(j-1))
          call ts_fix_contour( Eq_io(cur), prev=Eq_io(prev) )
       end if
       call consecutive_types(Eq_io(cur),isCircle,isTail)

    end do

    ! Allocate sizes of the contours
    do i = 1 , N_Eq
       ! number of different weights for each energy-point
       k = count(hasC(mus,Eq_c(i)))
       if ( k == 0 ) then
          call die('No chemical potentials has been assigned this contour: '// &
               trim(Eq_c(i)%c_io%name))
       end if

       ! the id's referring to the chemical potential in the
       ! mus array.
       allocate(Eq_c(i)%ID(k))
       k = 0
       do j = 1 , N_mu
          if ( hasC(mus(j),Eq_c(i)) ) then
             k = k + 1
             Eq_c(i)%ID(k) = j
          end if
       end do

       ! Allocate the different weights and initialize
       allocate(Eq_c(i)%c(Eq_c(i)%c_io%N))
       Eq_c(i)%c = 0._dp
       allocate(Eq_c(i)%w(k,Eq_c(i)%c_io%N))
       Eq_c(i)%w = 0._dp

    end do
    
    do i = 1 , N_mu

       call setup_Eq_contour(mus(i),N_poles,kT,Eq_Eta)
       
    end do

    write(*,*) 'TODO correct empty cycles'

  contains
    
    subroutine consecutive_types(c_io,isCircle,isTail)
      type(ts_c_io), intent(in) :: c_io
      logical, intent(inout) :: isCircle, isTail

      if ( c_io%type == 'circle' .and. .not. isCircle ) then
         call die('The circle contour must be the first contour type &
              &as well as connected to other circle contours.')
      end if
      isCircle = c_io%type == 'circle'

      if ( c_io%type /= 'tail' .and. isTail ) then
         call die('The tail contour must be the last contour &
              &as well as connected to other tail contours.')
      end if
      isTail   = c_io%type == 'tail'

    end subroutine consecutive_types

  end subroutine read_contour_eq_options


  ! This routine assures that we have setup all the 
  ! equilibrium contours for the passed electrode
  subroutine setup_Eq_contour(mu,N_poles,kT,Eta)
    type(ts_mu), intent(in) :: mu
    integer, intent(in) :: N_poles
    real(dp), intent(in) :: kT, Eta

    ! Local variables
    integer :: i, idx
    real(dp) :: a, b, R, cR, lift

    if ( Eq_segs(mu) == 0 ) &
         call die('Chemical potential: '//trim(Name(mu))//' has &
         &no equilibrium contours.')

    ! retrieve circle bounds for the electrode
    call mu_circle_bounds(mu,a,b)

    ! Calculate the circle entries
    call calc_Circle(a,b,N_poles,kT,Eta,R,cR,lift)

    do i = 1 , Eq_segs(mu)
       
       idx = get_c_io_index(mu%Eq_seg(i))

       if ( Eq_c(idx)%c_io%part == 'circle' ) then

          call contour_Circle(Eq_c(idx),mu,kT,R,cR,Eta)

       else if ( Eq_c(idx)%c_io%part == 'line' ) then

          call contour_line(Eq_c(idx),mu,kT,lift)

       else if ( Eq_c(idx)%c_io%part == 'tail' ) then

          call contour_tail(Eq_c(idx),mu,kT,lift)

       else if ( Eq_c(idx)%c_io%part == 'pole' ) then

          ! the poles all have the same weight (Pi*kT*2)
          call contour_poles(Eq_c(idx),Eq_c(idx)%c_io%d,kT,Eta)

       else
          
          call die('Unrecognized contour type for the &
               &equilibrium part.')

       end if

    end do
    
! right-hand side of the center
    write(*,*)'TODO insert checks for contour'
    
  contains

    subroutine calc_Circle(a,b,N_poles,kT,Eta,R,cR,lift)
      use units, only : Pi
      real(dp), intent(in)  :: a,b ! the circle bounds
      integer, intent(in) :: N_poles ! number of poles the 'b' point is lifted
      real(dp), intent(in) :: kT, Eta ! the imaginary part lifted
      real(dp), intent(out) :: R, cR, lift ! the radius, center(real)
      
      ! local variables
      real(dp) :: alpha

      lift = Pi * kT * (2._dp*(N_poles-1)+1._dp)
      do while ( lift < Eta )
         lift = lift + 2._dp * Pi * kT
      end do
      lift = lift + Pi * kT
      ! this means that we place the line contour right in the middle of two poles!

      cR = b - a
      ! the angle between the axis and the line from the start
      ! of the circle to the end of the circle contour
      alpha = datan( (lift - Eta) / cR )
    
      ! the radius can be calculated using two triangles in the circle
      ! there is no need to use the cosine-relations
      R = .5_dp * cR / cos(alpha) ** 2

      ! the real-axis center
      cR = a + R

    end subroutine calc_Circle

    ! retrieve the bounds of the circle contour
    ! this allows to split the circle integral into as many parts as needed
    subroutine mu_circle_bounds(mu,a,b)
      type(ts_mu), intent(in) :: mu
      real(dp), intent(out) :: a,b
      if ( Eq_segs(mu) < 1 ) call die('Error Eq_seg CB')
      idx = get_c_io_index(mu%Eq_seg(1))
      a = Eq_c(idx)%c_io%a
      b = Eq_c(idx)%c_io%b
      do i = 1 , Eq_segs(mu)
         idx = get_c_io_index(mu%Eq_seg(i))
         if ( Eq_c(idx)%c_io%type == 'eq' .and. Eq_c(idx)%c_io%part == 'circle' ) then
            b = Eq_c(idx)%c_io%b
         end if
      end do
    end subroutine Mu_circle_bounds

  end subroutine setup_Eq_contour

  subroutine ID2idx_cw(c,ID,idx)
    type(ts_cw), intent(in) :: c
    integer,     intent(in) :: ID
    integer,    intent(out) :: idx
    do idx = 1 , size(c%ID)
       if ( c%ID(idx) == ID ) return
    end do
    idx = -1
  end subroutine ID2idx_cw

  subroutine ID2idx_cidx(c,ID,idx)
    type(ts_c_idx), intent(in) :: c
    integer, intent(in) :: ID
    integer, intent(out) :: idx
    if ( c%idx(1) /= CONTOUR_EQ ) call die('Could not locate ID')
    call ID2idx_cw(Eq_c(c%idx(2)),ID,idx)
  end subroutine ID2idx_cidx

  pure subroutine c2weight_eq(c,idx,k,W,ZW)
    type(ts_c_idx), intent(in) :: c
    integer, intent(in) :: idx
    real(dp), intent(in)  :: k
    complex(dp), intent(out) :: W,ZW
    
    if ( c%idx(1) /= CONTOUR_EQ ) then
       W = 0._dp
       ZW = 0._dp
       return
    end if

    W  = k * Eq_c(c%idx(2))%w(idx,c%idx(3))
    ! TODO full write out of the energy-density matrix calculation (suspect that 
    ! the fault lies in this line)!!
    ZW = W * Eq_c(c%idx(2))%c(c%idx(3))
    
  end subroutine c2weight_eq

  subroutine contour_Circle(c,mu,kT,R,cR,Eta)
    use m_integrate
    use m_gauss_quad
    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    real(dp), intent(in) :: kT, R, cR, Eta

    ! local variables
    character(len=c_N) :: tmpC
    logical :: set_c
    complex(dp) :: ztmp
    integer :: i, idx
    real(dp) :: a,b, tmp
    real(dp), allocatable :: ce(:), cw(:)

    if ( c%c_io%part /= 'circle' ) &
         call die('Contour is not a line')
   
    ! notice the switch (the circle has increasing contour
    ! from right to left)
    call calc_angle(c%c_io%a,R,cR,a)
    call calc_angle(c%c_io%b,R,cR,b)
    if ( b > a ) then
       call die('Contour equilibrium circle went wrong')
    end if

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

    select case ( method(c%c_io) ) 
    case ( CC_G_LEGENDRE ) 

       if ( c_io_has_opt(c%c_io,'right') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call Gauss_Legendre_Rec(2*c%c_io%N,0,2._dp*a-b,b,ce,cw)

          do i = 1 ,c%c_io%N
             ce(i) = ce(i+c%c_io%N)
             cw(i) = 2._dp * cw(i+c%c_io%N)
          end do

       else if ( c_io_has_opt(c%c_io,'left') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call Gauss_Legendre_Rec(2*c%c_io%N,0,a,2._dp*b-a,ce,cw)

          cw = 2._dp * cw

       else

          call Gauss_Legendre_Rec(c%c_io%N,0,a,b,ce,cw)

       end if

    case ( CC_TANH_SINH )

       ! we should also gain an option for this
       if ( c_io_has_opt(c%c_io,'precision') ) then
          tmpC = c_io_get_opt(c%c_io,'precision')
          read(tmpC,'(g20.10)') tmp
       else
          tmp = 2.e-2_dp * abs(b-a) / real(c%c_io%N,dp)
          write(tmpC,'(g20.10)') tmp
          call c_io_add_opt(c%c_io,'precision',tmpC)
       end if

       if ( c_io_has_opt(c%c_io,'right') ) then
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call TanhSinh_Exact(c%c_io%N*2,ce,cw,2._dp*a-b,b, p=tmp)

          do i = 1 ,c%c_io%N
             ce(i) = ce(i+c%c_io%N)
             cw(i) = 2._dp * cw(i+c%c_io%N)
          end do

       else if ( c_io_has_opt(c%c_io,'left') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call TanhSinh_Exact(c%c_io%N*2,ce,cw,a,2._dp*b-a, p=tmp)

          cw = 2._dp * cw

       else

          call TanhSinh_Exact(c%c_io%N,ce,cw,a,b, p=tmp)

       end if

    case ( CC_BOOLE_MIX )

       call Booles_Simpson_38_3_rule(c%c_io%N, ce, cw, a, b)

    case default
       call die('Unknown method for the circle integral, please correct')
    end select

    ! I know this is "bad" practice, however, zero is a well-defined quantity in FPU
    set_c = sum(abs(c%c(:))) == 0._dp

    ! get the index in the ID array (same index in w-array)
    call ID2idx(c,mu%ID,idx)

    do i = 1 , c%c_io%N
       ztmp = R * cdexp(dcmplx(0._dp,ce(i)))

       if ( set_c ) then
          c%c(i) = dcmplx(cR,Eta) + ztmp
       else
          if ( abs(c%c(i) - (dcmplx(cR,Eta) + ztmp) ) > 1.e-10_dp ) then
             print *,c%c(i),dcmplx(cR,Eta) + ztmp
             call die('contours does not match')
          end if
       end if

       ! Factor i, comes from Ed\theta=dE=iR e^{i\theta}
       c%w(idx,i) = cw(i) * nf((dcmplx(cR,0._dp)+ztmp-mu%mu)/kT) * dcmplx(0._dp,1._dp) * ztmp

    end do

    deallocate(ce,cw)
    
  contains

    subroutine calc_angle(a,R,cR,angle)
      use units, only : Pi
      real(dp), intent(in) :: a,R,cR
      real(dp), intent(out) :: angle
      real(dp) :: E

      E = abs(cR - a)

      ! the integration angles
      if ( E <  1.e-8_dp ) then
         ! we are at the apex of the circle
         angle = Pi * .5_dp
      else if ( abs(E - R) <  1.e-8_dp ) then
         ! we are at the radius of the circle
         angle = Pi
      else
         angle = asin(E/R)
         ! correct for left/right side of radius
         if ( a < cR ) then
            angle = angle + Pi * .5_dp
         else if ( a > cR ) then
            angle = .5_dp * Pi - angle
         end if
      end if
    end subroutine calc_angle

  end subroutine contour_Circle

  subroutine contour_line(c,mu,kT,Eta)
    use m_integrate
    use m_gauss_quad

    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    real(dp), intent(in) :: kT, Eta

    ! local variables
    character(len=c_N) :: tmpC
    logical :: set_c
    integer :: i, idx
    real(dp) :: a,b, tmp
    real(dp), allocatable :: ce(:), cw(:)

    if ( c%c_io%part /= 'line' ) &
         call die('Contour is not a line')

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b
    
    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

    select case ( method(c%c_io) )
    case ( CC_MID )
       
       call Mid_Rule(c%c_io%N,ce,cw,a,b)
       
    case ( CC_SIMP_MIX )
       
       call Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)
       
    case ( CC_BOOLE_MIX )
       
       call Booles_Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)       
       
    case ( CC_G_LEGENDRE ) 

       if ( c_io_has_opt(c%c_io,'right') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call Gauss_Legendre_Rec(c%c_io%N*2,0,2._dp*a-b,b,ce,cw)

          do i = 1 ,c%c_io%N
             ce(i) = ce(i+c%c_io%N)
             cw(i) = 2._dp * cw(i+c%c_io%N)
          end do

       else if ( c_io_has_opt(c%c_io,'left') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call Gauss_Legendre_Rec(c%c_io%N*2,0,a,2._dp*b-a,ce,cw)

          cw = 2._dp * cw

       else

          call Gauss_Legendre_Rec(c%c_io%N,0,a,b,ce,cw)

       end if
       
    case ( CC_TANH_SINH ) 

       ! we should also gain an option for this
       if ( c_io_has_opt(c%c_io,'precision') ) then
          tmpC = c_io_get_opt(c%c_io,'precision')
          read(tmpC,'(g20.10)') tmp
       else
          tmp = 2.e-2_dp * abs(b-a) / real(c%c_io%N,dp)
          write(tmpC,'(g20.10)') tmp
          call c_io_add_opt(c%c_io,'precision',tmpC)
       end if

       if ( c_io_has_opt(c%c_io,'right') ) then
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call TanhSinh_Exact(c%c_io%N*2,ce,cw,2._dp*a-b,b, p=tmp)

          do i = 1 ,c%c_io%N
             ce(i) = ce(i+c%c_io%N)
             cw(i) = 2._dp * cw(i+c%c_io%N)
          end do

       else if ( c_io_has_opt(c%c_io,'left') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call TanhSinh_Exact(c%c_io%N*2,ce,cw,a,2._dp*b-a, p=tmp)

          cw = 2._dp * cw

       else

          call TanhSinh_Exact(c%c_io%N,ce,cw,a,b, p=tmp)

       end if
       
       
    case default
       !write(*,*) 'Method for contour ',trim(c%c_io%name),' could not be deciphered: ', &
       !     c%c_io%method
       call die('Could not determine the line-integral')
    end select

    ! I know this is "bad" practice, however, zero is a well-defined quantity.
    set_c = sum(abs(c%c(:))) == 0._dp

    ! get the index in the ID array (same index in w-array)
    call ID2idx(c,mu%ID,idx)

    do i = 1 , c%c_io%N
       if ( set_c ) then
          c%c(i) = dcmplx(ce(i),Eta)
       else
          if ( abs(c%c(i) - dcmplx(ce(i),Eta)) > 1.e-10_dp ) then
             call die('contour_line: Error on contour match')
          end if
       end if

       !TODO check that this is correct, we do it with respect to the chemical potential
       c%w(idx,i) = cw(i) * nf((ce(i)-mu%mu)/kT)

    end do

    deallocate(ce,cw)
    
  end subroutine contour_line

  subroutine contour_tail(c,mu,kT,Eta)
    use m_gauss_fermi_inf
    use m_gauss_fermi_30
    use m_gauss_fermi_28
    use m_gauss_fermi_26
    use m_gauss_fermi_24
    use m_gauss_fermi_22
    use m_gauss_fermi_20
    use m_gauss_fermi_19
    use m_gauss_fermi_18
    use m_gauss_fermi_17
    use units, only : eV

    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    real(dp), intent(in) :: kT, Eta

    ! local variables
    integer :: idx, offset, infinity
    real(dp) :: a,b
    real(dp), allocatable :: ce(:), cw(:)

    if ( c%c_io%part /= 'tail' ) &
         call die('Contour is not a tail contour')

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b

    write(*,*) 'TODO check the contours for the gaussian quadrature'

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))
    
    select case ( method(c%c_io) )
    case ( CC_G_NF_MIN:CC_G_NF_MAX )

       offset = nint((c%c_io%a-mu%mu)/kT)

       if ( abs(offset * kT - (c%c_io%a-mu%mu)) > 1.e-7_dp ) then
          call die('The integer value of the kT offset for the &
               &Gauss-Fermi integral is not valid, please check input')
       end if
       if ( c%c_io%b < 0._dp ) then
          if ( c%c_io%b < -30.5_dp * kT ) then
             infinity = huge(1)
          else
             infinity = nint(abs(c%c_io%b-mu%mu)/kT)
          end if
       else
          if ( c%c_io%b > 30.5_dp * kT ) then
             infinity = huge(1)
          else
             infinity = nint(abs(c%c_io%b-mu%mu)/kT)
          end if
       end if
       if ( infinity > 30 ) infinity = huge(1)
       
       ! calculate the offset from the energy chemical potential tail
       select case ( infinity )
       case ( huge(1) )
          call GaussFermi_inf(offset,c%c_io%N,ce,cw)
       case ( 30 )
          call GaussFermi_30(offset,c%c_io%N,ce,cw)
       case ( 28 )
          call GaussFermi_28(offset,c%c_io%N,ce,cw)
       case ( 26 ) 
          call GaussFermi_26(offset,c%c_io%N,ce,cw)
       case ( 24 ) 
          call GaussFermi_24(offset,c%c_io%N,ce,cw)
       case ( 22 )
          call GaussFermi_22(offset,c%c_io%N,ce,cw)
       case ( 20 )
          call GaussFermi_20(offset,c%c_io%N,ce,cw)
       case ( 19 )
          call GaussFermi_19(offset,c%c_io%N,ce,cw)
       case ( 18 )
          call GaussFermi_18(offset,c%c_io%N,ce,cw)
       case ( 17 )
          call GaussFermi_17(offset,c%c_io%N,ce,cw)
       case default
          call die('Unknown tail integral ending')
       end select

       ! we might as well correct the method

       ! the Gauss-Fermi quadrature is wrt. E'->E/kT
       ce = ce * kT + mu%mu
       cw = cw * kT

       call ID2idx(c,mu%ID,idx)

       ! move over the weights and the contour values
       c%c        = dcmplx(ce,Eta)
       c%w(idx,:) = cw

    case default

       ! we revert so that we can actually use the line-integral
       c%c_io%part = 'line'

       call contour_line(c,mu,kT,Eta)

    end select

    deallocate(ce,cw)

  end subroutine contour_tail


  ! The residuals of the fermi-function at a real-energy
  subroutine contour_poles(c, E, kT, Eta)
    
    use precision, only : dp
    use units, only : Pi

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_cw), intent(inout) :: c

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E    ! at energy
    real(dp), intent(in) :: kT   ! temperature in Ry
    real(dp), intent(in) :: Eta  ! The lift of the equilibrium contour

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: i

    ! all pole-weights have the same weight (negative due to contour method)
    c%w(:,:) =  - dcmplx(0._dp, Pi * kT * 2._dp)
    ! Residuals
    do i = 1 , c%c_io%N
       c%c(i) = dcmplx(E , Pi * kT * (2._dp*(i-1)+1._dp))
    end do

    ! correct the imaginary 
    do while ( aimag(c%c(1)) < Eta ) 
       c%c(:) = c%c(:) + dcmplx(0._dp,2._dp*Pi*kT)
    end do

  end subroutine contour_poles

  function Eq_E(id,step) result(c)
    integer, intent(in) :: id
    integer, intent(in), optional :: step
    type(ts_c_idx) :: c ! the configuration of the energy-segment
    integer :: lstep, i, PN
    lstep = 1
    if ( present(step) ) lstep = step
    PN = N_Eq_E()
    if ( id <= PN ) then
       c = get_c(id)
       return
    end if
    c = get_c(-1)
    i = MOD(PN,lstep)
    if ( i /= 0 ) PN = PN + lstep - i
    if ( id <= PN ) then
       c%exist = .true.
       c%fake  = .true.
    end if
  end function Eq_E

  function get_c(id) result(c)
    integer, intent(in) :: id
    type(ts_c_idx) :: c
    integer :: i,j,iE
    c%exist = .false.
    c%fake  = .false.
    c%e     = dcmplx(0._dp,0._dp)
    c%idx   = 0
    if ( id < 1 ) return

    iE = 0
    do j = 1 , N_Eq ! number of contours
       if ( iE + Eq_c(j)%c_io%N < id ) then
          iE = iE + Eq_c(j)%c_io%N
          cycle
       end if
       i = id - iE
       if ( i <= Eq_c(j)%c_io%N ) then
          c%exist = .true.
          c%e     = Eq_c(j)%c(i)
          c%idx(1) = 1 ! designates the equilibrium contours
          c%idx(2) = j ! designates the index of the equilibrium contour
          c%idx(3) = i ! is the index of the equilibrium contour
          return
       end if
    end do

  end function get_c

  function N_Eq_E() result(N)
    integer :: N, i
    N = 0
    do i = 1 , N_Eq
       N = N + size(Eq_c(i)%c)
    end do
  end function N_Eq_E

  subroutine print_contour_eq_block(prefix)
    character(len=*), intent(in) :: prefix
    integer :: i

    do i = 1 , N_Eq
       if ( Eq_io(i)%part /= 'pole' ) then
          call ts_print_contour_block(trim(prefix)//'.Contour.',Eq_io(i))
       end if
    end do
  end subroutine print_contour_eq_block

  subroutine print_contour_eq_options(prefix)

    use parallel, only : IONode
    use units, only : eV

    use m_ts_io_contour

    character(len=*), intent(in) :: prefix
    character(len=200) :: chars
    integer :: i
    type(ts_c_opt_ll), pointer :: opt

    if ( .not. IONode ) return
    
    write(*,opt_n) ' ----------------- Contour ----------------- '

    write(*,opt_n) '           >> Residual contour << '
    do i = 1 , N_eq
       if ( eq_io(i)%part == 'pole' ) then
          chars = trim(eq_io(i)%part)
          call write_e('Pole chemical potential',eq_io(i)%d)
          write(*,opt_int) '  Pole points',eq_io(i)%N
       end if
    end do

    write(*,opt_n) '          >> Equilibrium contour << '
    write(*,opt_g_u) 'Equilibrium Greens function Eta',Eq_Eta/eV,'eV'
    do i = 1 , N_eq
       if ( eq_io(i)%part /= 'pole' ) then
          chars = '  '//trim(eq_io(i)%part)
          write(*,opt_c) 'Contour name',trim(prefix)//'.Contour.'//trim(eq_io(i)%name)
          call write_e(trim(chars)//' contour E_min',eq_io(i)%a)
          call write_e(trim(chars)//' contour E_max',eq_io(i)%b)
          write(*,opt_int) trim(chars)//' contour points',eq_io(i)%N
          write(*,opt_c) trim(chars)//' contour method', &
               trim(longmethod2str(eq_io(i)))
          opt => eq_io(i)%opt
          do while ( associated(opt) )
             if ( len_trim(opt%val) > 0 ) then
                write(*,opt_cc) '   Option for contour method',trim(opt%opt),trim(opt%val)
             else
                write(*,opt_c)  '   Option for contour method',trim(opt%opt)
             end if
             opt => opt%next
          end do
       end if
    end do
    
  end subroutine print_contour_eq_options

  function get_c_io_index(Name) result(idx)
    character(len=*), intent(in) :: Name
    integer :: idx
    do idx = 1 , N_Eq
       if ( trim(name) == trim(Eq_io(idx)%name) ) return
    end do
    idx = 0
  end function get_c_io_index

  subroutine io_contour_eq(mus,slabel,suffix)
    type(ts_mu), intent(in) :: mus(:)
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

    integer :: i

    do i = 1 , size(mus)
       call io_contour_eq_mu(mus(i),slabel,suffix)
    end do

  end subroutine io_contour_eq

  subroutine io_contour_eq_mu(mu,slabel,suffix)
    use parallel, only : IONode
    use fdf, only : leqi
    use units, only : eV
    type(ts_mu), intent(in) :: mu
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix
    
! *********************
! * LOCAL variables   *
! *********************
    character(len=200) :: fname
    integer :: i, j, unit, idx
    type(ts_c_idx) :: cidx
    
    if ( .not. IONode ) return
    
    if ( present(suffix) ) then
       fname = trim(slabel)//trim(suffix)
    else
       fname = trim(slabel)//'.TSCCEQ-'//trim(Name(mu))
    end if

    call io_assign( unit )
    open( unit, file=fname, status='unknown' )
    write(unit,'(a)') '# Contour path for the equilibrium contour segment.'
    write(unit,'(a)') '# This segment belongs to the chemical potential: '//trim(Name(mu))
    write(unit,'(a)') '# It has the chemical potential:'
    write(unit,'(a,tr1,f10.5,tr1,a)') '#',mu%mu/eV,'eV'
    write(unit,'(a,a12,3(tr1,a13))') '#','Re(c) [eV]','Im(c) [eV]','Re(w)','Im(w)'

    cidx%idx(1) = CONTOUR_EQ
    do i = 1 , Eq_segs(mu)

       cidx%idx(2) = c_ioname2idx(mu%Eq_seg(i))
       if ( cidx%idx(2) < 1 ) call die('io_contour_eq_mu: Error in code, C-ID')
       call ID2idx(Eq_c(cidx%idx(2)),mu%ID,idx)
       if ( idx < 1 ) call die('io_contour_eq_mu: Error in code')

       !print *,trim(Name(mu)),mu%Eq_seg(i)

       ! we have now retrieved the electrode index in the contour part
       ! write it out
       call io_contour_c(unit,cidx,idx)

    end do

    call io_close( unit )

  contains
    
    function c_ioname2idx(cName) result(idx)
      character(len=*), intent(in) :: cName
      integer :: idx
      do idx = 1 , N_Eq
         if ( Eq_c(idx)%c_io%Name == cName ) return
      end do
      idx = -1
    end function c_ioname2idx

  end subroutine io_contour_eq_mu

! Write out the contour to a contour file
  subroutine io_contour_c(unit,cidx,idx)
    use parallel, only : IONode
    use units, only : eV
    integer, intent(in) :: unit
    type(ts_c_idx), intent(inout) :: cidx
    integer, intent(in) :: idx

! *********************
! * LOCAL variables   *
! *********************
    integer :: i
    type(ts_cw), pointer :: c
    complex(dp) :: W, ZW

    if ( .not. IONode ) return
    c => Eq_c(cidx%idx(2))

    do i = 1 , size(c%c)
       cidx%e      = c%c(i)
       cidx%idx(3) = i
       call c2weight_eq(cidx,idx,1._dp,W,ZW)
       write(unit,'(4(e13.6,tr1))') c%c(i)/eV, W
    end do
    
  end subroutine io_contour_c

end module m_ts_contour_eq
