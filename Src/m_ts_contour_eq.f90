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

module m_ts_contour_eq

  use precision, only : dp
!
! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_electype

  use m_ts_cctype
  use m_ts_io_contour
  use m_ts_io_ctype
  use m_ts_aux

  implicit none

  ! equilibrium contour IO-segments
  integer, save, public :: N_eq
  type(ts_c_io), pointer, save, public :: Eq_io(:) => null()
  type(ts_eq_c), pointer, save, public :: Eq_c(:) => null()

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

  public :: read_contour_eq_options
  public :: print_contour_eq_options
  public :: print_contour_eq_block
  public :: io_contour_eq

  private

contains

  subroutine read_contour_eq_options(Elecs, kT, Volt)

    use units, only : eV
    use parallel, only : IONode, Nodes, operator(.parcount.)
    use fdf

    type(Elec), intent(inout) :: Elecs(:)
    real(dp), intent(in) :: kT, Volt

    integer :: i,j,k
    integer :: different_poles, N
    character(len=C_N_NAME_LEN), allocatable :: tmp(:), nContours(:)
    integer :: cur, next, prev

    real(dp), parameter :: mu_same = 1.e-10_dp

    call fdf_obsolete('TS.ComplexContour.NPoles')

    ! broadening
    Eq_Eta = fdf_get('TS.Contours.Eq.Eta',0._dp,'Ry')
    if ( Eq_Eta < 0._dp) call die('ERROR: Eq_Eta < 0, we do not allow &
         &for using the advanced Greens function, please correct.')
    
    ! We only allow the user to either use the old input format, or the new
    ! per-electrode input

    ! Read in the generic things about the contours...
    N_poles = fdf_get('TS.Contours.Eq.Pole.N',6)

    ! figure out how many different poles we need for the "fake" c_io types which contains the poles
    different_poles = 1 ! this is the last electrodes mu
    do i = 1 , size(Elecs) - 1
       ! we only count the different mu-levels for the last electrode
       if ( count(abs(Elecs(i+1:)%mu - Elecs(i)%mu) < mu_same) == 0 ) then
          different_poles = different_poles + 1
       end if
       ! assign names of the pole contour to the electrodes
       k = 1
       if ( i == 1 ) then
          write(Elecs(i)%Eq_seg(Eq_segs(Elecs(i))),'(a,i0)') 'pole',k
       else
          do j = i + 1 , size(Elecs)
             if ( abs(Elecs(j)%mu - Elecs(i)%mu) > mu_same ) then
                k = k + 1
             else
                write(Elecs(i)%Eq_seg(Eq_segs(Elecs(i))),'(a,i0)') 'pole',k
                exit
             end if
          end do
       end if
          
    end do
    ! Create the last pole
    i = size(Elecs)
    k = k + 1
    write(Elecs(i)%Eq_seg(Eq_segs(Elecs(i))),'(a,i0)') 'pole',k

    ! Count the number of contour segments
    N = 0
    do i = 1 , size(Elecs)
       N = N + Eq_segs(Elecs(i))
    end do
    if ( N == 0 ) then
       call die('You have not assigned any contours for the electrodes.')
    end if

    ! collect all equilibrium names
    allocate(tmp(N))
    tmp(:) = ' '
    N = 0
    do i = 1 , size(Elecs)
       do j = 1 , Eq_segs(Elecs(i)) - 1
          N = N + 1
          tmp(N) = Elecs(i)%Eq_seg(j)
       end do
    end do
    do i = 1 , size(Elecs)
       j = Eq_segs(Elecs(i))
       N = N + 1
       tmp(N) = Elecs(i)%Eq_seg(j)
    end do


    ! find all unique equilibrium names
    N = 1
    uniq_names: do i = 2 , size(tmp)
       j = 0
       do 
          j = j + 1
          if ( i <= j ) exit
          if ( tmp(i) == tmp(j) ) cycle uniq_names
       end do
       N = N + 1
    end do uniq_names

    allocate(nContours(N)) ! unique names
    nContours = ' '

    ! populate unique equilibrium names
    N = 0
    uniq_names_pop: do i = 1 , size(tmp)
       do j = 1 , N
          if ( tmp(i) == nContours(j) ) cycle uniq_names_pop
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
    allocate(Eq_io(N_Eq))     ! the equilibrium contour io container
    allocate(Eq_c(N_Eq))      ! the equilibrium contour container
    Eq_io(:)%type = 'eq' ! set the equilibrium tag

    ! Attach all the contours
    do i = 1 , N
       Eq_c(i)%c_io => Eq_io(i)
    end do

    ! read in the equilibrium contours
    do i = 1 , N - different_poles

       print *,'reading ',nContours(i)
       ! read in the contour
       call ts_read_contour_block('TS','',nContours(i),Eq_io(i), kT, Volt)

    end do

    ! We here create the "fake" pole contours
    j = 0
    k = 0
    do i = N - different_poles + 1, N
       ! assign name to the Eq_io
       Eq_io(i)%name = nContours(i)
       Eq_io(i)%N = N_poles
       Eq_io(i)%part = 'pole'
       Eq_io(i)%method = 'residual'
       ! find the next mu level
       do 
          j = j + 1
          if ( j > different_poles ) call die('Pole setup creation went wrong, contact devs')
          ! we only count the different mu-levels for the last electrode
          if ( j == different_poles ) then
             k = k + 1
             Eq_io(i)%a = Elecs(j)%mu
             Eq_io(i)%b = Elecs(j)%mu
             Eq_io(i)%d = Elecs(j)%mu
             exit
          else if ( count(abs(Elecs(j+1:)%mu - Elecs(j)%mu) < mu_same) == 0 ) then
             k = k + 1
             Eq_io(i)%a = Elecs(j)%mu
             Eq_io(i)%b = Elecs(j)%mu
             Eq_io(i)%d = Elecs(j)%mu
             exit
          end if
       end do
           
    end do

    ! check that we have assigned all poles
    if ( k < different_poles ) then
       call die('Pole setup went wrong, please contact the developers')
    end if

    ! fix the read-in constants
    do i = 1 , size(Elecs)

       j = 1
       cur = get_c_io_index(Elecs(i)%Eq_seg(j))
       if ( Eq_segs(Elecs(i)) > 2 ) then
          next = get_c_io_index(Elecs(i)%Eq_seg(j+1))
          call ts_fix_contour( Eq_io(cur), next=Eq_io(next) )
       end if
          
       ! we should not check the pole (hence minus 2)
       do j = 2 , Eq_segs(Elecs(i)) - 2
          prev = get_c_io_index(Elecs(i)%Eq_seg(j-1))
          cur  = get_c_io_index(Elecs(i)%Eq_seg(j))
          next = get_c_io_index(Elecs(i)%Eq_seg(j+1))
          call ts_fix_contour( Eq_io(cur), &
               prev=Eq_io(prev), next=Eq_io(next) )
       end do

       ! don't check the pole
       j = Eq_segs(Elecs(i)) - 1
       cur = get_c_io_index(Elecs(i)%Eq_seg(j))
       if ( Eq_segs(Elecs(i)) > 2 ) then
          prev = get_c_io_index(Elecs(i)%Eq_seg(j-1))
          call ts_fix_contour( Eq_io(cur), prev=Eq_io(prev) )
       end if

       write(*,*)'TODO need to check the circle-line-tail sequence'

    end do

    ! the best way to setup the contours is to
    ! set them up electrode wise...
    ! This is required if the circle is deformed
    j = 0

    ! Assign names to the contours
    do i = 1 , N_Eq

       call assign_Electrode(Elecs,Eq_c(i))

       if ( Eq_c(i)%Ne == 0 ) then
          call die('No electrodes has been assigned this contour: '// &
               trim(Eq_c(i)%c_io%name))
       end if

    end do

    do i = 1 , size(Elecs)

       write(*,*)'TODO check that the energy relative to the contours are mu'
       call setup_Eq_contour(Elecs(i),N_poles,kT,Eq_Eta)

    end do

    write(*,*) 'TODO correct empty cycles'

  end subroutine read_contour_eq_options


  ! Routine for creating the contour
  subroutine assign_Electrode(Elecs,c)
    
    use precision, only : dp
    use fdf, only : leqi

! **********************
! * INPUT variables    *
! **********************
    type(Elec), intent(in) :: Elecs(:)
    type(ts_eq_c), intent(inout) :: c ! The contour

    integer :: i,j, N
    integer :: NE
    
    NE = size(Elecs)

    ! first we find all the electrodes that needs to be attached
    if ( leqi(c%c_io%part,'pole') ) then ! this must be an equilibrium contour
       ! If it is a pole we need to compare against the chemical potential
       N = count(abs(Elecs(:)%mu - c%c_io%d) < 1.e-10_dp )
       c%Ne = N
       allocate(c%elec(c%Ne))
       j = 0
       do i = 1 , NE
          if ( abs(Elecs(i)%mu - c%c_io%d) < 1.e-10_dp ) then
             j = j + 1
             c%elec(j) = trim(Name(Elecs(i)))
          end if
       end do
       
    else
       ! we search in the electrodes equilibrium contour segments
       N = 0
       do i = 1 , NE
          N = N + count(Elecs(i)%Eq_seg(:) == c%c_io%name)
       end do
       c%Ne = N
       allocate(c%elec(c%Ne))
       j = 0
       do i = 1 , NE
          if ( count(Elecs(i)%Eq_seg(:) == c%c_io%name) > 0 ) then
             j = j + 1
             c%elec(j) = trim(Name(Elecs(i)))
          end if
       end do
    end if

    ! Allocate the electrode weights and contours
    allocate(c%c(c%c_io%N),c%w(c%Ne,c%c_io%N))

    c%c = 0._dp
    c%w = 0._dp

  end subroutine assign_Electrode

  ! This routine assures that we have setup all the 
  ! equilibrium contours for the passed electrode
  subroutine setup_Eq_contour(El,N_poles,kT,Eta)
    type(Elec), intent(in) :: El
    integer, intent(in) :: N_poles
    real(dp), intent(in) :: kT,Eta

    ! Local variables
    integer :: i, idx
    real(dp) :: a,b,R,cR,cI,lift

    if ( Eq_segs(El) == 0 ) &
         call die('Electrode: '//trim(Name(El))//' has no equilibrium contours.')

    ! retrieve circle bounds for the electrode
    call Elec_circle_bounds(El,a,b)

    ! Calculate the circle entries
    call calc_Circle(a,b,N_poles,kT,Eta,R,cR,lift)

    do i = 1 , Eq_segs(El)
       
       idx = get_c_io_index(El%Eq_seg(i))

       call ts_print_contour_block('Contour.',Eq_c(idx)%c_io)
       
       if ( Eq_c(idx)%c_io%part == 'circle' ) then

          call contour_Circle(Eq_c(idx),El,kT,R,cR,Eta)

       else if ( Eq_c(idx)%c_io%part == 'line' ) then

          call contour_line(Eq_c(idx),El,kT,lift)

       else if ( Eq_c(idx)%c_io%part == 'tail' ) then

          call contour_tail(Eq_c(idx),El,kT,lift)

       else if ( Eq_c(idx)%c_io%part == 'pole' ) then

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

      lift = Pi*kT*(2._dp*(N_poles-1)+1._dp)
      do while ( lift < Eta )
         lift = lift + 2._dp * Pi * kT
      end do
      lift = lift + Pi * kT
      ! this means that we place the line contour right in the middle of two poles!

      ! the angle between point 'a' and 'gamma'
      alpha = datan(lift - Eta) / (b-a)
    
      ! the radius can be calculated using two triangles in the circle
      ! there is no need to use the cosine-relations
      R = .5_dp * (b-a)/cos(alpha)**2
    
      ! the real-axis center
      cR = a + R

    end subroutine calc_Circle

    ! retrieve the bounds of the circle contour
    ! this allows to split the circle integral into as many parts as needed
    subroutine Elec_circle_bounds(El,a,b)
      type(Elec), intent(in) :: El
      real(dp), intent(out) :: a,b
      if ( Eq_segs(El) < 1 ) call die('Error Eq_seg CB')
      idx = get_c_io_index(El%Eq_seg(1))
      a = Eq_c(idx)%c_io%a
      b = Eq_c(idx)%c_io%b
      do i = 1 , Eq_segs(El)
         idx = get_c_io_index(El%Eq_seg(i))
         if ( Eq_c(idx)%c_io%type == 'eq' .and. Eq_c(idx)%c_io%part == 'circle' ) then
            b = Eq_c(idx)%c_io%b
         end if
      end do
    end subroutine Elec_circle_bounds

  end subroutine setup_Eq_contour

  subroutine contour_Circle(c,El,kT,R,cR,Eta)
    use m_integrate
    use m_gauss_quad
    type(ts_eq_c), intent(inout) :: c
    type(Elec), intent(in) :: El
    real(dp), intent(in) :: kT, R, cR, Eta

    ! local variables
    logical :: set_c
    complex(dp) :: ztmp
    integer :: i, j, iE
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

    ! we know that the contour has already been allocated
    iE = 0
    do 
       iE = iE + 1
       if ( iE > c%Ne ) call die('Error on finding the electrode')
       if ( c%elec(iE) == Name(El) ) exit
    end do

    write(*,*) 'TODO correct for right/left quadrature contours'
    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

    select case ( method(c%c_io) ) 
    case ( CC_G_LEGENDRE ) 

       call Gauss_Legendre_Rec(c%c_io%N,0,a,b,ce,cw)

    case ( CC_TANH_SINH )

       tmp = 1.8e-2_dp * abs(b-a) / real(c%c_io%N,dp)
       call TanhSinh_Exact(c%c_io%N, ce, cw, a, b, p=tmp)

    case ( CC_BOOLE_MIX )

       call Booles_Simpson_38_3_rule(c%c_io%N, ce, cw, a, b)

    case default
       call die('Unknown method for the circle integral, please correct')
    end select

    ! I know this is "bad" practice, however, zero is a well-defined quantity in FPU
    set_c = sum(abs(c%c(:))) == 0._dp

    do i = 1 , c%c_io%N
       ! we need to reverse the order due to the contours creation
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
       c%w(iE,i) = cw(i) * nf((dcmplx(cR,0._dp)+ztmp-El%mu)/kT) * dcmplx(0._dp,1._dp) * ztmp

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

  subroutine contour_line(c,El,kT,Eta)
    use m_integrate
    use m_gauss_quad

    type(ts_eq_c), intent(inout) :: c
    type(Elec), intent(in) :: El
    real(dp), intent(in) :: kT, Eta

    ! local variables
    logical :: set_c
    integer :: i, iE
    real(dp) :: a,b, tmp, loffset
    real(dp), allocatable :: ce(:), cw(:)

    if ( c%c_io%part /= 'line' ) &
         call die('Contour is not a line')

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b
    
    ! the offset in energy... Could be the chemical potential
    iE = 0
    do
       iE = iE + 1
       if ( iE > c%Ne ) call die('Error on finding the electrode')
       if ( c%elec(iE) .eq. El ) exit
    end do

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

    write(*,*) 'TODO correct TANH-sinh option gathering',trim(c%c_io%name),c%c_io%N

    select case ( method(c%c_io) )
    case ( CC_MID )
       
       call Mid_Rule(c%c_io%N,ce,cw,a,b)
       
    case ( CC_SIMP_MIX )
       
       call Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)
       
    case ( CC_BOOLE_MIX )
       
       call Booles_Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)       
       
    case ( CC_G_LEGENDRE ) 
       
       call Gauss_Legendre_Rec(c%c_io%N,0,a,b,ce,cw)
       
    case ( CC_TANH_SINH ) 
       
       tmp = 2.e-2_dp * abs(b-a) / real(c%c_io%N,dp)
       call TanhSinh_Exact(c%c_io%N,ce,cw,a,b, p=tmp)
       
    case default
       !write(*,*) 'Method for contour ',trim(c%c_io%name),' could not be deciphered: ', &
       !     c%c_io%method
       call die('Could not determine the line-integral')
    end select

    ! I know this is "bad" practice, however, zero is a well-defined quantity.
    set_c = sum(abs(c%c(:))) == 0._dp

    do i = 1 , c%c_io%N
       if ( set_c ) then
          c%c(i) = dcmplx(ce(i),Eta)
       else
          if ( abs(c%c(i) - dcmplx(ce(i),Eta)) > 1.e-10_dp ) then
             call die('contours does not match')
          end if
       end if

       !TODO check that this is correct, we do it with respect to the chemical potential
       c%w(iE,i) = cw(i) * nf((ce(i)-El%mu)/kT)

    end do

    deallocate(ce,cw)
    
  end subroutine contour_line

  subroutine contour_tail(c,El,kT,Eta)
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

    type(ts_eq_c), intent(inout) :: c
    type(Elec), intent(in) :: El
    real(dp), intent(in) :: kT, Eta

    ! local variables
    integer :: i, iE, offset, infinity
    real(dp) :: a,b
    real(dp), allocatable :: ce(:), cw(:)

    if ( c%c_io%part /= 'tail' ) &
         call die('Contour is not a tail contour')

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b

    ! Find the electrode position for the equilibrium contour
    iE = 0
    do
       iE = iE + 1
       if ( iE > c%Ne ) call die('Error on finding the electrode')
       if ( c%elec(iE) .eq. El ) exit
    end do

    write(*,*) 'TODO check the contours for the gaussian quadrature'

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))
    
    select case ( method(c%c_io) )
    case ( CC_G_NF_MIN:CC_G_NF_MAX )

       offset = nint((c%c_io%a-El%mu)/kT)

       if ( abs(offset * kT - (c%c_io%a-El%mu)) > 1.e-7_dp ) then
          call die('The integer value of the kT offset for the &
               &Gauss-Fermi integral is not valid, please check input')
       end if
       if ( c%c_io%b < 0._dp ) then
          if ( c%c_io%b < -30.5_dp * kT ) then
             infinity = huge(1)
          else
             infinity = nint(abs(c%c_io%b-El%mu)/kT)
          end if
       else
          if ( c%c_io%b > 30.5_dp * kT ) then
             infinity = huge(1)
          else
             infinity = nint(abs(c%c_io%b-El%mu)/kT)
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
       ce = ce * kT + El%mu
       cw = cw * kT

       ! move over the weights and the contour values
       c%c       = dcmplx(ce,Eta)
       c%w(iE,:) = cw

    case default

       ! we revert so that we can actually use the line-integral
       c%c_io%part = 'line'

       call contour_line(c,El,kT,Eta)

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
    type(ts_eq_c), intent(inout) :: c

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
             write(*,opt_c) '   Option for contour method', trim(opt%opt)
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

  subroutine io_contour_eq(Elecs,slabel,suffix)
    type(Elec), intent(in) :: Elecs(:)
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

    integer :: i

    do i = 1 , size(Elecs)
       call io_contour_eq_el(Elecs(i),slabel,suffix)
    end do

  end subroutine io_contour_eq

  subroutine io_contour_eq_el(El,slabel,suffix)
    use parallel, only : IONode
    use fdf, only : leqi
    type(Elec), intent(in) :: El
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix
    
! *********************
! * LOCAL variables   *
! *********************
    character(len=200) :: fname
    integer :: i, j, idx, unit
    
    if ( .not. IONode ) return
    
    if ( present(suffix) ) then
       fname = trim(slabel)//trim(suffix)
    else
       fname = trim(slabel)//'.TSCCEQ-'//trim(Name(El))
    end if

    call io_assign( unit )
    open( unit, file=fname, status='unknown' )
    write(unit,'(a)') '# Contour path for the equilibrium &
         &contour segment of electrode: '//trim(Name(El))
    write(unit,'(a,a12,3(tr1,a13))') '#','Re(c)[eV]','Im(c)[eV]','Re(w)','Im(w)'

    do i = 1 , Eq_segs(El)
       ! We put them out in contour sequence
       idx = 0
       do j = 1 , N_Eq
          if ( leqi(El%Eq_seg(i),Eq_c(j)%c_io%name) ) then
             idx = elec_idx(Eq_c(j),El)
             exit
          end if
       end do
       if ( idx == 0 ) &
            call die('Could not find the index of the electrode')
       
       ! we have now retrieved the electrode index in the contour part
       ! write it out
       call io_contour_c(unit,Eq_c(j),idx)

    end do

    call io_close( unit )

  end subroutine io_contour_eq_el
       
! Write out the contour to a contour file
  subroutine io_contour_c(unit,c,iE)
    use parallel, only : IONode
    use units, only : eV
    integer, intent(in) :: unit
    type(ts_eq_c), intent(in) :: c
    integer, intent(in) :: iE

! *********************
! * LOCAL variables   *
! *********************
    integer :: i

    if ( .not. IONode ) return

    do i = 1 , size(c%c)
       write(unit,'(4(e13.6,tr1))') c%c(i)/eV,c%w(iE,i)
    end do
    
  end subroutine io_contour_c

end module m_ts_contour_eq
