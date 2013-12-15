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

module m_ts_contour_neq

  use precision, only : dp

! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_electype

  use m_ts_cctype
  use m_ts_io_contour
  use m_ts_io_ctype
  use m_ts_aux

  implicit none

  ! The non-equilibrium density integration are attributed discussions with
  ! Antti-Pekka Jauho. 

  ! For further attributions see the original paper by Brandbyge, et. al, 2002: DOI: 10.1103/PhysRevB.65.165401

  ! Contour path
  integer, save, public :: N_nEq, N_nEq_tail
  type(ts_c_io),  pointer, save, public :: nEq_io(:) => null(), nEq_tail_io(:) => null()
  type(ts_nEq_c), pointer, save, public :: nEq_c(:) => null(), nEq_tail_c(:) => null()
  ! this is the actual tail integral read in from the options.
  ! it is merely a placeholder before we revert to the nEq_tail array
  integer, save, private :: N_tail
  type(ts_c_io),  pointer, save, private :: tail_io(:) => null()
  type(ts_nEq_c), pointer, save, private :: tail_c(:) => null()

  ! type to contain the information about each contour element.
  type :: ts_nEq_seg
     type(Elec_p), allocatable :: Els1(:), Els2(:)
     ! the corresponding chemical potentials of the electrodes
     real(dp) :: mu1, mu2
     ! the indices for the nEq_io contours so that we don't look them up each time
     integer, allocatable :: io(:)
     ! the indices for the tail_io contours so that we don't look them up each time
     ! we always have two tails in any one segment
     integer :: tail_io(2) = 0
  end type ts_nEq_seg
  type(ts_nEq_seg), pointer :: nEq_segs(:) => null()

  ! The contour specific variables
  real(dp), save, public :: nEq_Eta


  public :: read_contour_neq_options
  public :: print_contour_neq_options
  public :: print_contour_neq_block
  public :: io_contour_neq
  public :: N_nEq_E, N_nEq_window_E, N_nEq_tail_E
  public :: nEq_E

  private

contains

  subroutine read_contour_neq_options(Elecs, kT, Volt)

    use units, only : eV
    use parallel, only : IONode, Nodes, operator(.parcount.)
    use fdf
    use m_ts_electype

    ! only to gain access to the chemical shifts
    type(Elec), intent(in), target :: Elecs(:)
    real(dp), intent(in) :: kT, Volt
    
    integer :: i, j, k, N, diff_mu, left, right
    real(dp), parameter :: mu_same = 1.e-10_dp
    real(dp), allocatable :: mus(:)
    real(dp) :: tmp
    integer, allocatable :: mus_tail(:)
    logical, allocatable :: mus_used(:)

    call fdf_obsolete('TS.biasContour.Eta')

    ! broadening
    nEq_Eta = fdf_get('TS.Contours.nEq.Eta',0.000001_dp*eV,'Ry')
    if ( nEq_Eta <= 0._dp ) call die('ERROR: nEq_Eta <= 0, we do not allow &
         &for using the advanced Greens function, please correct.')

    ! We only allow the user to either use the old input format, or the new
    ! per-electrode input


    ! Bias-window setup
    call my_setup('Bias.Window',N_nEq,nEq_c,nEq_io)

    ! Here we setup the tail integral
    ! TODO consider only doing the tails at V/2 and -V/2
    ! TODO if the dE is small enough then we might not need the quadrature for the real-axis
    call my_setup('Bias.Tail',N_tail,tail_c,tail_io)
    if ( N_tail > 1 ) &
         call die('You can only use one tail integral')

    ! We need to create a contour on the real axis
    ! which is split in each of the electrodes fermi levels
    ! First we find the number of different fermi levels
    diff_mu = 1
    do i = 1 , size(Elecs) - 1
       ! we only count the different mu-levels for the last electrode
       if ( count(abs(Elecs(i+1:)%mu - Elecs(i)%mu) < mu_same) == 0 ) then
          diff_mu = diff_mu + 1
       end if
    end do
    allocate(mus(diff_mu))
    diff_mu = 1
    do i = 1 , size(Elecs) - 1
       ! we only count the different mu-levels for the last electrode
       if ( count(abs(Elecs(i+1:)%mu - Elecs(i)%mu) < mu_same) == 0 ) then
          mus(diff_mu) = Elecs(i)%mu
          diff_mu = diff_mu + 1
       end if
    end do
    mus(diff_mu) = Elecs(size(Elecs))%mu
    ! sort the potentials so that we have the chemical potentials
    ! following the non-equilibrium real-axis contour (i.e. in ascending order)
    do i = 1 , diff_mu - 1
       do j = i + 1 , diff_mu
          if ( mus(j) < mus(i) ) then
             tmp = mus(j)
             mus(j) = mus(i)
             mus(i) = tmp
          end if
       end do
    end do

    ! Create all the different tail segments
    ! We have one tail in both ends and two tails for each middle segment
    ! We also check whether the tail integral fits in any of the contours
    N_nEq_tail = 2 + (diff_mu - 2) * 2
    allocate(nEq_tail_io(N_nEq_tail))
    allocate(nEq_tail_c(N_nEq_tail))
    allocate(mus_tail(diff_mu))
    mus_tail = 0
    j = 1
    do i = 1 , diff_mu
       ! this should check that we are at an edge chemical potential
       if ( count(mus(:) - mus(i) < mu_same) == 1 .and. &
            count(mus(i) - mus(:) < mu_same) == diff_mu ) then
          mus_tail(i) = fits_left(mus(i),tail_io)

          call copy(tail_io(1),nEq_tail_io(j)) ! TODO N_tail > 1
          call assign_set_E(nEq_tail_c(j),nEq_tail_io(j), &
               tail_c(1),tail_io(:),mus(i)) ! TODO N_tail > 1
          j = j + 1
       else if ( count(mus(i) - mus(:) < mu_same) == 1 .and. &
            count(mus(:) - mus(i) < mu_same) == diff_mu ) then
          mus_tail(i) = fits_right(mus(i),tail_io)

          call copy(tail_io(1),nEq_tail_io(j)) ! TODO N_tail > 1
          call assign_set_E(nEq_tail_c(j),nEq_tail_io(j), & 
               tail_c(1), tail_io(:), mus(i))
          j = j + 1
       else
          ! we are in a middle segment
          mus_tail(i) = min(fits_left(mus(i),tail_io),fits_right(mus(i),tail_io))

          ! we need both the left and right tail
          call copy(tail_io(1),nEq_tail_io(j)) ! TODO N_tail > 1
          call assign_set_E(nEq_tail_c(j),nEq_tail_io(j), & 
               tail_c(1), tail_io(:), mus(i),reverse=.true.)
          j = j + 1

          call copy(tail_io(1),nEq_tail_io(j)) ! TODO N_tail > 1
          call assign_set_E(nEq_tail_c(j),nEq_tail_io(j), & 
               tail_c(1), tail_io(:), mus(i), reverse=.false.)
          j = j + 1
       end if
    end do
    ! check that a tail can be placed at all segments
    if ( any(mus_tail == 0) ) then
       call die('No real axis contour tail fits with your chemical potentials')
    end if

    do j = 1 , N_nEq_tail
       call io_contour_c(6,nEq_tail_c(j))
    end do

    ! Allocate all the segments, this comes from simple permutation rules
    ! 2 chemical potentials => 1 segment
    ! 3 chemical potentials => 3 segments
    ! 4 chemical potentials => 6 segments, etc.
    i = 0 ! count
    do j = diff_mu - 1 , 1 , -1
       i = i + j
    end do
    allocate(nEq_segs(i))

    ! populate the segments by their chemical potentials
    allocate(mus_used(size(Elecs)))
    diff_mu = 1
    do i = 1 , size(mus) - 1
       mus_used = abs(Elecs(:)%mu - mus(i)) < mu_same
       do j = 0 , size(mus) - i - 1
          nEq_segs(diff_mu+j)%mu1 = mus(i)
          call allocate_and_assign_electrodes(nEq_segs(diff_mu+j), &
               mus_used,.true.)
       end do
       do j = i + 1 , size(mus)
          mus_used = abs(Elecs(:)%mu - mus(j)) < mu_same
          nEq_segs(diff_mu)%mu2 = mus(j)
          call allocate_and_assign_electrodes(nEq_segs(diff_mu), &
               mus_used,.false.)
          ! we need to find the tail contours that fit
          left  = fits_left(nEq_segs(diff_mu)%mu1,tail_io)
          right = fits_right(nEq_segs(diff_mu)%mu2,tail_io)
          !print *,nEq_segs(diff_mu)%mu1,left,right,nEq_segs(diff_mu)%mu2
          if ( left == 0 .or. right == 0 ) &
               call die('Something went wrong with the segment tails')
          if ( left > right ) &
               call die('The contours have not been sorted properly, please &
               &contact the developers')
          ! Allocate the pointers to the non-equilibrium contours
          N = abs(left-right) + 1
          allocate(nEq_segs(diff_mu)%io(N))
          do k = left , right
             nEq_segs(diff_mu)%io(k-left+1) = k
          end do
          nEq_segs(diff_mu)%tail_io(:) = 0
          do k = 1 , N_nEq_tail
             if ( abs(nEq_io(nEq_segs(diff_mu)%io(1))%a - nEq_tail_io(k)%b) < mu_same ) then
                nEq_segs(diff_mu)%tail_io(1) = k
             end if
             if ( abs(nEq_io(nEq_segs(diff_mu)%io(N))%b - nEq_tail_io(k)%a) < mu_same ) then
                nEq_segs(diff_mu)%tail_io(2) = k
             end if
          end do
          if ( any(nEq_segs(diff_mu)%tail_io == 0) ) then
             call die('Could not find all tails')
          end if

          diff_mu = diff_mu + 1
       end do
    end do

    write(*,*) 'Mu1',nEq_segs(:)%mu1
    write(*,*) 'Mu2',nEq_segs(:)%mu2

    deallocate(mus_used,mus,mus_tail)
    write(*,*) 'TODO check that the bias window stops at every \mu and that &
         &a equivalent electrode has that \mu'
    write(*,*) 'TODO correct empty cycles, i.e. if two line contours are neighbours &
         &then we have overlying energy points...'

  contains 

    subroutine assign_set_E(c,c_io,c_from,c_io_from,mu,reverse)
      type(ts_nEq_c), intent(inout) :: c
      type(ts_c_io), target :: c_io
      type(ts_c_io) :: c_io_from(:)
      type(ts_nEq_c), intent(in) :: c_from
      real(dp), intent(in) :: mu
      logical, intent(in), optional :: reverse
      logical :: lreverse
      integer :: i, j
      ! assign
      c%c_io => c_io

      lreverse = mu < 0._dp
      if ( present(reverse) ) lreverse = reverse

      ! correct the end points
      if ( lreverse ) then
         ! we are the lower tail (so reverse it)
         c_io%a = -c_io_from(size(c_io_from))%b + mu
         c_io%b = -c_io_from(1)%a               + mu
      else
         c_io%a =  c_io_from(1)%a               + mu
         c_io%b =  c_io_from(size(c_io_from))%b + mu
      end if

      ! create the contours
      allocate(c%c(c_io%N),c%w(c_io%N))
      if ( lreverse ) then
         do j = c_io%N , 1 , -1
            i = c_io%N - j + 1
            c%c(i) = dcmplx(-dreal(c_from%c(j)),dimag(c_from%c(j))) + mu
            c%w(i) = c_from%w(j)
         end do
      else
         do i = 1 , c_io%N
            c%c(i) = c_from%c(i) + mu
            c%w(i) = c_from%w(i)
         end do
      end if
    end subroutine assign_set_E

    function fits_left(mu,tail_io) result(i)
      real(dp), intent(in) :: mu
      type(ts_c_io), intent(in) :: tail_io(:)
      integer :: i
      do i = 1 , N_nEq
         if ( abs(mu - nEq_io(i)%a - tail_io(1)%a) < mu_same ) then
            return
         end if
      end do
      i = 0
    end function fits_left

    function fits_right(mu,tail_io) result(i)
      real(dp), intent(in) :: mu
      type(ts_c_io), intent(in) :: tail_io(:)
      integer :: i
      do i = 1 , N_nEq
         if ( abs(nEq_io(i)%b - mu - tail_io(1)%a) < mu_same ) then
            return
         end if
      end do
      i = 0
    end function fits_right

    subroutine my_setup(suffix,N_nEq,nEq_c,nEq_io)
      character(len=*), intent(in) :: suffix
      integer, intent(inout) :: N_nEq
      type(ts_nEq_c), pointer :: nEq_c(:)
      type(ts_c_io), pointer :: nEq_io(:)

      ! Local variables
      integer :: i
      character(len=C_N_NAME_LEN), allocatable :: tmp(:)

      N_nEq = fdf_nc_iotype('TS',suffix)
      if ( N_nEq < 1 ) call die('You must specify at least one non-equilbrium &
           &contour for the '//trim(suffix)//'.')
      allocate(tmp(N_nEq))

      tmp(1) = fdf_name_c_iotype('TS',suffix,1)
      do i = 2 , N_nEq
         tmp(i) = fdf_name_c_iotype('TS',suffix,i)
         if ( count(tmp(:i-1) == tmp(i)) /= 0 ) then
            call die('You cannot have two names from the bias-window &
                 &to be the same...')
         end if
      end do
      
      ! allocate all required objects
      nullify(nEq_io,nEq_c)
      allocate(nEq_io(N_nEq),nEq_c(N_nEq))
      
      do i = 1 , N_nEq

         ! assign pointer
         nEq_c(i)%c_io => nEq_io(i)
         ! read in the contour
         call ts_read_contour_block('TS',suffix,tmp(i),nEq_io(i), kT, Volt)
       
      end do
      deallocate(tmp)

      if ( N_nEq == 1 ) then
         call ts_fix_contour(nEq_io(1))
      end if
      do i = 1 , N_nEq - 1
         if ( 1 < i ) then
            call ts_fix_contour(nEq_io(i), &
                 prev=nEq_io(i-1), next=nEq_io(i+1))
         else
            call ts_fix_contour(nEq_io(i), next=nEq_io(i+1) )
         end if
      end do
      if ( N_nEq > 1 ) then
         call ts_fix_contour(nEq_io(N_neq), prev=nEq_io(N_nEq-1) )
      end if

      ! setup the contour
      do i = 1 , N_nEq
         ! allocate contour
         allocate(nEq_c(i)%c(nEq_c(i)%c_io%N),nEq_c(i)%w(nEq_c(i)%c_io%N))
         call setup_nEq_contour(nEq_c(i), kT, nEq_Eta)
      end do

      if ( nEq_c(1)%c_io%a > nEq_c(N_nEq)%c_io%b ) then
         call die('The non-equilibrium contours must be in increasing &
              energy. Even if your bias is negative. Please correct.')
      end if

    end subroutine my_setup

    subroutine allocate_and_assign_electrodes(nEq_seg,mus,left)
      type(ts_nEq_seg), intent(inout) :: nEq_seg
      logical, intent(in) :: mus(:), left
      integer :: i, j
      if ( left ) then
         allocate(nEq_seg%Els1(count(mus)))
      else
         allocate(nEq_seg%Els2(count(mus)))
      end if
      j = 0
      do i = 1 , size(mus)
         if ( mus(i) ) then
            j = j + 1
            if ( left ) then
               nEq_seg%Els1(j)%El => Elecs(i)
            else
               nEq_seg%Els2(j)%El => Elecs(i)
            end if
         end if
      end do
    end subroutine allocate_and_assign_electrodes

  end subroutine read_contour_neq_options

  ! This routine assures that we have setup all the 
  ! equilibrium contours for the passed electrode
  subroutine setup_nEq_contour(c, kT, Eta)
    type(ts_neq_c), intent(inout) :: c
    real(dp), intent(in) :: kT, Eta

    if ( c%c_io%part == 'line' ) then
       
       call contour_line(c,kT,Eta)
       
    else if ( c%c_io%part == 'tail' ) then
       
       call contour_tail(c,kT,Eta)

    else
       
       call die('Unrecognized contour type for the &
            &non-equilibrium part.')
       
    end if
    
! right-hand side of the center
    write(*,*)'TODO insert checks for contour'
    
  end subroutine setup_nEq_contour

  function segment_has_c(seg,i_c) result(has)
    type(ts_nEq_seg), intent(in) :: seg ! a segment that needs to be tested for 
                                        ! part
    integer, intent(in) :: i_c          ! the index of the part in the list of parts
    logical :: has
    has = allocated(seg%io)
    if ( has ) has = any(i_c == seg%io)
  end function segment_has_c

  function segment_has_tail_c(seg,i_c) result(has)
    type(ts_nEq_seg), intent(in) :: seg ! a segment that needs to be tested for 
                                        ! part
    integer, intent(in) :: i_c          ! the index of the part in the list of parts
    logical :: has
    has = any(i_c == seg%tail_io)
  end function segment_has_tail_c

  function segment_has_El(seg,El) result(has)
    type(ts_nEq_seg), intent(in) :: seg
    type(Elec), intent(in) :: El
    logical :: has
    integer :: i
    has = .false.
    do i = 1 , size(seg%Els1)
       if ( seg%Els1(i)%El .eq. El ) then
          has = .true.
          return
       end if
    end do

    do i = 1 , size(seg%Els2)
       if ( seg%Els2(i)%El .eq. El ) then
          has = .true.
          return
       end if
    end do
  end function segment_has_El
    

  subroutine weight_El_c(El,seg,i_c,i_e,kT, weight)
    use m_ts_aux, only : nf
    type(Elec), intent(in) :: El ! the electrode we calculate the contribution with respect to
    type(ts_nEq_seg), intent(in) :: seg ! the segment the energy-point is drawn from
    integer, intent(in) :: i_c ! the index of the contour (in the nEq_c array)
    integer, intent(in) :: i_e ! the index of the energy point in the contour
    real(dp), intent(in) :: kT ! the temperature
    real(dp), intent(out) :: weight ! the weight returned
    logical :: left
    ! local variables
    integer :: i

    if ( .not. (segment_has_c(seg,i_c) .or. &
         segment_has_El(seg,El)) ) then
       weight = 0._dp
       return
    end if
    ! check if the electrode is in the left segment
    left = .false.
    do i = 1 , size(seg%Els1)
       if ( seg%Els1(i)%El .eq. El ) then
          left = .true.
       end if
    end do
    
    weight = nEq_c(i_c)%w(i_e) * &
         nf(real(nEq_c(i_c)%c(i_e),dp),seg%mu1,seg%mu2,kT)
    if ( left ) weight = - weight

  end subroutine weight_El_c

  subroutine contour_line(c,kT,Eta)
    use m_integrate
    use m_gauss_quad
    type(ts_neq_c), intent(inout) :: c
    real(dp), intent(in) :: kT, Eta

    ! local variables
    real(dp) :: a,b, tmp
    real(dp), allocatable :: ce(:), cw(:)

    if ( c%c_io%part /= 'line' ) &
         call die('Contour is not a line')

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

    write(*,*) 'TODO correct TANH-sinh option gathering'

    select case ( method(c%c_io%method) )
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
       call die('Could not determine the line-integral')
    end select

    c%c = dcmplx(ce,Eta)
    c%w = dcmplx(cw,0._dp)

    deallocate(ce,cw)
    
  end subroutine contour_line

  subroutine contour_tail(c,kT,Eta)
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
    type(ts_neq_c), intent(inout) :: c
    real(dp), intent(in) :: kT, Eta

    ! local variables
    integer :: ioffset, infinity
    real(dp) :: a,b
    real(dp), allocatable :: ce(:), cw(:)

    if ( c%c_io%part /= 'tail' ) &
         call die('Contour is not a tail contour')

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b

    if ( b < 0._dp ) then
       call die('The non-equilbrium tail contours can only be &
            &defined with respect to the Fermi-level')
    end if

    write(*,*) 'TODO check the contours for the gaussian quadrature'

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))
    
    select case ( method(c%c_io%method) )
    case ( CC_G_NF_MIN:CC_G_NF_MAX )

       ! calculate the offset
       ioffset = nint(a/kT)
       if ( abs(ioffset * kT - a) > 1.e-7_dp ) then
          call die('The integer value of the kT offset for the &
               &Gauss-Fermi tail integral is not valid, please check input')
       end if
       if ( b > 30.5_dp * kT ) then
          infinity = huge(1)
       else
          infinity = nint(b/kT)
       end if
       if ( infinity > 30 ) infinity = huge(1)
       
       ! calculate the offset from the energy chemical potential tail
       select case ( infinity )
       case ( huge(1) )
          call GaussFermi_inf(ioffset,c%c_io%N,ce,cw)
       case ( 30 )
          call GaussFermi_30(ioffset,c%c_io%N,ce,cw)
       case ( 28 )
          call GaussFermi_28(ioffset,c%c_io%N,ce,cw)
       case ( 26 ) 
          call GaussFermi_26(ioffset,c%c_io%N,ce,cw)
       case ( 24 ) 
          call GaussFermi_24(ioffset,c%c_io%N,ce,cw)
       case ( 22 )
          call GaussFermi_22(ioffset,c%c_io%N,ce,cw)
       case ( 20 )
          call GaussFermi_20(ioffset,c%c_io%N,ce,cw)
       case ( 19 )
          call GaussFermi_19(ioffset,c%c_io%N,ce,cw)
       case ( 18 )
          call GaussFermi_18(ioffset,c%c_io%N,ce,cw)
       case ( 17 )
          call GaussFermi_17(ioffset,c%c_io%N,ce,cw)
       case default
          call die('Unknown tail integral ending')
       end select

       ! correct for the Gauss-Fermi unit-transformation
       ce = ce * kT
       cw = cw * kT

       ! move over the weights and the contour values
       c%c = dcmplx(ce,Eta)
       c%w = dcmplx(cw,0._dp)

    case default

       ! we revert so that we can actually use the line-integral
       c%c_io%part = 'line'

       call contour_line(c,kT,Eta)

    end select

    deallocate(ce,cw)

  end subroutine contour_tail


  function nEq_E(id,step) result(c)
    integer, intent(in) :: id
    integer, intent(in), optional :: step
    type(ts_c) :: c ! the configuration of the energy-segment
    integer :: lstep, i, PN
    lstep = 1
    if ( present(step) ) lstep = step
    PN = N_nEq_E()
    i = MOD(PN,lstep)
    if ( i /= 0 ) PN = PN + lstep - i
    do i = 1 , PN , lstep
       if ( i == id ) then
          c = get_c(id)
          return
       end if
    end do
    c = get_c(-1)
  end function nEq_E

  function get_c(id) result(c)
    integer, intent(in) :: id
    type(ts_c) :: c
    integer :: i,j,iE
    c%exist = .false.
    c%e     = dcmplx(0._dp,0._dp)
    c%idx   = 0
    if ( id < 1 ) return

    iE = 0
    do j = 1 , N_nEq ! number of contours
       do i = 1 , nEq_c(j)%c_io%N
          iE = iE + 1 
          if ( iE == id ) then
             c%exist = .true.
             c%e     = nEq_c(j)%c(i)
             c%idx(1) = 2 ! this designates the non-equilibrium contours
             c%idx(2) = j ! this designates the index of the non-equilibrium contour
             c%idx(3) = i ! this is the index of the non-equilibrium contour
             return
          end if
       end do
    end do

    do j = 1 , N_nEq_tail ! number of contours
       do i = 1 , nEq_tail_c(j)%c_io%N
          iE = iE + 1 
          if ( iE == id ) then
             c%exist = .true.
             c%e     = nEq_tail_c(j)%c(i)
             c%idx(1) = 3 ! this designates the tail non-equilibrium contours
             c%idx(2) = j ! this designates the index of the tail non-equilibrium contour
             c%idx(3) = i ! this is the index of the tail non-equilibrium contour
             return
          end if
       end do
    end do

  end function get_c

  function N_nEq_E() result(N)
    integer :: N
    N = N_nEq_window_E() + N_nEq_tail_E()
  end function N_nEq_E

  function N_nEq_window_E() result(N)
    integer :: N, i
    N = 0
    do i = 1 , N_nEq
       N = N + size(nEq_c(i)%c)
    end do
  end function N_nEq_window_E

  function N_nEq_tail_E() result(N)
    integer :: N, i
    N = 0
    do i = 1 , N_nEq_tail
       N = N + size(nEq_tail_c(i)%c)
    end do
  end function N_nEq_tail_E

  subroutine print_contour_neq_block(prefix)
    use parallel, only : IONode
    character(len=*), intent(in) :: prefix

    integer :: i

    if ( IONode ) then
       write(*,'(2a)') '%block ',trim(prefix)//'.Contours.Bias.Window'
       do i = 1 , N_nEq
          write(*,'(tr4,a)') trim(nEq_io(i)%name)
       end do
       write(*,'(2a,/)') '%endblock ',trim(prefix)//'.Contours.Bias.Window'
    end if

    do i = 1 , N_nEq
       call ts_print_contour_block(trim(prefix)//'.Contour.Bias.Window.',nEq_io(i))
    end do

    if ( IONode ) then
       write(*,'(/,2a)') '%block ',trim(prefix)//'.Contours.Bias.Tail'
       do i = 1 , N_tail
          write(*,'(tr4,a)') trim(tail_io(i)%name)
       end do
       write(*,'(2a,/)') '%endblock ',trim(prefix)//'.Contours.Bias.Tail'
    end if
    do i = 1 , N_tail
       call ts_print_contour_block(trim(prefix)//'.Contour.Bias.Tail.',tail_io(i))
    end do
  end subroutine print_contour_neq_block


  subroutine print_contour_neq_options(prefix)

    use parallel, only : IONode
    use units, only : eV

    use m_ts_io_contour

    character(len=*), intent(in) :: prefix
    character(len=200) :: chars
    integer :: i
    type(ts_c_opt_ll), pointer :: opt

    if ( .not. IONode ) return
    
    write(*,opt_n) '        >> non-Equilibrium contour << '
    write(*,opt_g_u) 'non-Equilibrium Greens function Eta',nEq_Eta/eV,'eV'
    do i = 1 , N_nEq
       chars = '  '//trim(nEq_io(i)%part)
       write(*,opt_c) 'Contour name',trim(prefix)//'.Contour.Bias.Window.'//trim(neq_io(i)%name)
       call write_e(trim(chars)//' contour E_min',neq_io(i)%a)
       call write_e(trim(chars)//' contour E_max',neq_io(i)%b)
       write(*,opt_int) trim(chars)//' contour points',neq_io(i)%N
       write(*,opt_c) trim(chars)//' contour method', &
            trim(longmethod2str(neq_io(i)))
       opt => neq_io(i)%opt
       do while ( associated(opt) )
          write(*,opt_c) '   Option for contour method', trim(opt%opt)
          opt => opt%next
       end do
    end do

    write(*,opt_n) '       > non-Equilibrium tail contour <'
    do i = 1 , N_neq_tail
       chars = '  '//trim(neq_tail_io(i)%part)
       !write(*,opt_c) 'Contour name',trim(prefix)//'.Contour.Bias.Tail.'//trim(neq_tail_io(i)%name)
       call write_e(trim(chars)//' contour E_min',neq_tail_io(i)%a)
       call write_e(trim(chars)//' contour E_max',neq_tail_io(i)%b)
       write(*,opt_int) trim(chars)//' contour points',neq_tail_io(i)%N
       write(*,opt_c) trim(chars)//' contour method', &
            trim(longmethod2str(neq_tail_io(i)))
       opt => neq_tail_io(i)%opt
       do while ( associated(opt) )
          write(*,opt_c) '   Option for contour method', trim(opt%opt)
          opt => opt%next
       end do
    end do

  end subroutine print_contour_neq_options

  subroutine io_contour_neq(slabel,suffix)
    use parallel, only : IONode
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

! *********************
! * LOCAL variables   *
! *********************
    character(len=25) :: tmp_suffix
    integer :: i

    if ( .not. IONode ) return

    do i = 1 , size(nEq_segs)
       
       if ( present(suffix) ) then
          write(tmp_suffix,'(a,i0)') trim(suffix)//'-',i
       else
          write(tmp_suffix,'(a,i0)') 'TSNEQ-',i
       end if
       call io_contour_neq_seg(nEq_segs(i),slabel,tmp_suffix)
       
    end do

  end subroutine io_contour_neq


  subroutine io_contour_neq_seg(seg,slabel,suffix)
    use parallel, only : IONode
    use units, only : eV
    type(ts_nEq_seg), intent(in) :: seg
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in) :: suffix

! *********************
! * LOCAL variables   *
! *********************
    character(len=200) :: fname
    integer :: i, unit
    
    if ( .not. IONode ) return
    
    fname = trim(slabel)//'.'//trim(suffix)

    call io_assign( unit )
    open( unit, file=fname, status='unknown' )
    write(unit,'(a)') '# Contour path for the non-equilibrium part'
    write(unit,'(a)') '# Segment between following chemical potentials in eV'
    write(unit,'(a,2(tr1,f10.5))') '#',seg%mu1/eV,seg%mu2/eV
    write(unit,'(a,a12,3(tr1,a13))') '#','Re(c)[eV]','Im(c)[eV]','Re(w)','Im(w)'

    if ( seg%tail_io(1) /= 0 ) then
       call io_contour_c(unit,nEq_tail_c(seg%tail_io(1)))
    end if
    do i = 1 , size(seg%io)
       call io_contour_c(unit,nEq_c(seg%io(i)),seg%mu1,seg%mu2)
    end do
    if ( seg%tail_io(2) /= 0 ) then
       call io_contour_c(unit,nEq_tail_c(seg%tail_io(2)))
    end if
    
    call io_close( unit )

  end subroutine io_contour_neq_seg

  subroutine io_contour_c(unit,c,mu1,mu2)
    use units, only : eV
    use m_ts_aux, only : nf
    integer, intent(in) :: unit
    type(ts_nEq_c), intent(in) :: c
    real(dp), intent(in), optional :: mu1, mu2
    integer :: i
    if ( present(mu1) .and. present(mu2) ) then
       write(*,*) 'TODO missing kT factor here'
       do i = 1 , size(c%c)
          write(unit,'(4(e13.6,tr1))') c%c(i)/eV,c%w(i) * &
               nf(real(c%c(i),dp),mu1,mu2,.025_dp*eV)
       end do
    else
       do i = 1 , size(c%c)
          write(unit,'(4(e13.6,tr1))') c%c(i)/eV,c%w(i)
       end do
    end if
  end subroutine io_contour_c
    
end module m_ts_contour_neq
