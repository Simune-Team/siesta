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
  use fdf, only : leqi

! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_electype

  use m_ts_chem_pot
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
  type(ts_cw), pointer, save, public :: nEq_c(:) => null(), nEq_tail_c(:) => null()
  ! this is the actual tail integral read in from the options.
  ! it is merely a placeholder before we revert to the nEq_tail array
  integer, save, private :: N_tail
  type(ts_c_io),  save, pointer, private :: tail_io(:) => null()
  type(ts_cw), save, pointer, private :: tail_c(:) => null()

  ! type to contain the information about each contour element.
  type :: ts_nEq_seg
     type(ts_mu), pointer :: mu1 => null(), mu2 => null()
     ! the indices for the nEq_io contours so that we don't look them up each time
     integer, allocatable :: io(:)
     ! the indices for the tail_io contours so that we don't look them up each time
     ! we always have two tails in any one segment
     integer :: tail_io(2) = 0
  end type ts_nEq_seg
  integer, save, public :: N_nEq_segs = 0
  type(ts_nEq_seg), save, pointer :: nEq_segs(:) => null()

  type :: ts_nEq_id
     integer :: ID = 0
     integer :: imu = 0
     integer :: iEl = 0
     type(ts_nEq_seg), pointer :: seg => null()
  end type ts_nEq_id
  integer, save, public :: N_nEq_id = 0
  type(ts_nEq_id), save, pointer :: nEq_ID(:) => null()

  ! The contour specific variables
  real(dp), save, public :: nEq_Eta

  ! this is heavily linked with the CONTOUR_EQ from m_ts_contour_eq
  integer, parameter, public :: CONTOUR_NEQ = 2
  integer, parameter, public :: CONTOUR_NEQ_TAIL = 3

  public :: read_contour_neq_options
  public :: print_contour_neq_options
  public :: print_contour_neq_block
  public :: io_contour_neq
  public :: N_nEq_E, N_nEq_window_E, N_nEq_tail_E
  public :: nEq_E
  public :: has_cE_neq
  public :: c2weight_neq
  public :: ID2mu

  private

contains

  subroutine read_contour_neq_options(N_Elec,Elecs,N_mu,mus, kT, Volt)

    use units, only : eV
    use fdf
    use m_ts_electype

    ! only to gain access to the chemical shifts
    integer, intent(in) :: N_Elec
    type(Elec), intent(in), target :: Elecs(N_Elec)
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in), target :: mus(N_mu)
    real(dp), intent(in) :: kT, Volt
    
    integer :: i, j, k, N, cur_mu, left, right
    integer, allocatable :: mus_tail(:)


    call fdf_obsolete('TS.biasContour.Eta')

    ! check that we in-fact have a bias calculation
    if ( N_mu < 2 ) then
       call die('Something has gone wrong. We can only find one chemical potential')
    end if

    ! broadening
    nEq_Eta = fdf_get('TS.Contours.nEq.Eta',0._dp*eV,'Ry')
    if ( nEq_Eta < 0._dp ) call die('ERROR: nEq_Eta < 0, we do not allow &
         &for using the advanced Greens function, please correct.')

    ! We only allow the user to either use the old input format, or the new
    ! per-electrode input

    ! Bias-window setup
    call my_setup('Bias.Window',N_nEq,nEq_c,nEq_io)
    if ( N_nEq < 1 ) then

       if ( N_Elec /= 2 .or. N_mu /= 2 ) then
          call die('When using other than 2 electrodes and &
               &chemical potentials you are forced to setup &
               &the contour input yourself.')
       end if

       ! We create the default version
       ! *** NOTE requires an even splitting of the bias
       N_nEq = 1
       nullify(nEq_io,nEq_c)
       allocate(nEq_io(N_nEq),nEq_c(N_nEq))
       nEq_c(1)%c_io => nEq_io(1)
       nEq_io(1)%part = 'line'
       nEq_io(1)%name = 'neq-1'
       nEq_io(1)%ca = '-|V|/2'
       nEq_io(1)%a  = -abs(Volt) * .5_dp
       nEq_io(1)%cb = '|V|/2'
       nEq_io(1)%b  =  abs(Volt) * .5_dp
       ! number of points
       nEq_io(1)%cd = '0.01 eV'
       nEq_io(1)%d = 0.01_dp * eV
       nEq_io(1)%method = 'simpson-mix'
       call ts_fix_contour(nEq_io(1))

       ! setup the contour
       allocate(nEq_c(1)%c(nEq_c(1)%c_io%N),nEq_c(1)%w(nEq_c(1)%c_io%N,1))
       call setup_nEq_contour(nEq_c(1), kT, nEq_Eta)

    end if

    ! Here we setup the tail integral
    call my_setup('Bias.Tail',N_tail,tail_c,tail_io)
    if ( N_tail > 1 ) &
         call die('You can currently only use one tail integral')
    if ( N_tail < 1 ) then

       if ( N_Elec /= 2 .or. N_mu /= 2 ) then
          call die('When using other than 2 electrodes and &
               &chemical potentials you are forced to setup &
               &the input yourself.')
       end if

       ! We create the default version
       ! *** NOTE requires an even splitting of the bias
       ! We create the default version
       ! *** NOTE requires an even splitting of the bias
       N_tail = 1
       nullify(tail_c,tail_io)
       allocate(tail_io(N_tail),tail_c(N_tail))
       tail_c(1)%c_io => tail_io(1)
       tail_io(1)%part = 'tail'
       tail_io(1)%name = 'neq-tail'
       tail_io(1)%ca = '0. kT'
       tail_io(1)%a  = 0._dp
       tail_io(1)%cb = '12 kT'
       tail_io(1)%b  =  12._dp * kT
       ! number of points
       tail_io(1)%cd = '0.01 eV'
       tail_io(1)%d =   0.01_dp * eV
       tail_io(1)%method = 'simpson-mix'
       call ts_fix_contour(tail_io(1))

       ! setup the contour
       allocate(tail_c(1)%c(tail_c(1)%c_io%N),tail_c(1)%w(tail_c(1)%c_io%N,1))
       call setup_nEq_contour(tail_c(1), kT, nEq_Eta)

    end if
    if ( tail_io(N_tail)%b < 10._dp * kT ) then
       call neq_die('The tail integrals used for the non-equilibrium tails &
            &are too close to the chemical potential. It must be at least &
            &10kT from mu.')
    end if

    ! Create all the different tail segments
    ! We have one tail in both ends and two tails for each middle segment
    ! We also check whether the tail integral fits in any of the contours
    N_nEq_tail = 2 + (N_mu - 2) * 2
    allocate(nEq_tail_io(N_nEq_tail))
    allocate(nEq_tail_c(N_nEq_tail))
    allocate(mus_tail(N_mu))
    mus_tail = 0
    j = 1
    do i = 1 , N_mu
       ! this should check that we are at an edge chemical potential
       if ( abs(minval(mus(:)%mu) - mus(i)%mu) < mu_same ) then
          mus_tail(i) = fits_left(mus(i),tail_io)

          call copy(tail_io(1),nEq_tail_io(j)) ! TODO N_tail > 1
          call assign_set_E(nEq_tail_c(j),nEq_tail_io(j), &
               tail_c(1),tail_io(:),mus(i), reverse = .true. ) ! TODO N_tail > 1
          j = j + 1

       else if ( abs(maxval(mus(:)%mu) - mus(i)%mu) < mu_same ) then
          mus_tail(i) = fits_right(mus(i),tail_io)

          call copy(tail_io(1),nEq_tail_io(j)) ! TODO N_tail > 1
          call assign_set_E(nEq_tail_c(j),nEq_tail_io(j), & 
               tail_c(1), tail_io(:), mus(i), reverse = .false. )
          j = j + 1

       else
          ! we are in a middle chemical potential
          mus_tail(i) = min(fits_left(mus(i),tail_io),fits_right(mus(i),tail_io))

          ! we need both the left and right tail
          call copy(tail_io(1),nEq_tail_io(j)) ! TODO N_tail > 1
          call assign_set_E(nEq_tail_c(j),nEq_tail_io(j), & 
               tail_c(1), tail_io(:), mus(i), reverse= .true.)
          j = j + 1

          call copy(tail_io(1),nEq_tail_io(j)) ! TODO N_tail > 1
          call assign_set_E(nEq_tail_c(j),nEq_tail_io(j), & 
               tail_c(1), tail_io(:), mus(i), reverse=.false.)
          j = j + 1

       end if
    end do
    ! check that a tail can be placed at all segments
    if ( any(mus_tail == 0) ) then
       call neq_die('No real axis contour tail fits with your chemical potentials')
    end if
    deallocate(mus_tail)


    ! Allocate all the segments, this comes from simple permutation rules
    ! 2 chemical potentials => 1 segment
    ! 3 chemical potentials => 3 segments
    ! 4 chemical potentials => 6 segments, etc.
    N_nEq_segs = 0 ! count
    do j = N_mu - 1 , 1 , -1
       N_nEq_segs = N_nEq_segs + j
    end do
    allocate(nEq_segs(N_nEq_segs))

    ! count the number of non-equilibrium segments concerning the
    ! electrodes density matrix update
    ! 2 electrodes, 2 mu => 2
    ! 3 electrodes, 2 mu => 3
    ! 4 electrodes, 2 mu => 4
    ! 5 electrodes, 2 mu => 5
    ! 3 electrodes, 3 mu => 6
    ! 6 electrodes, 2 mu => 6
    ! 7 electrodes, 2 mu => 7
    ! 4 electrodes, 3 mu => 8
    ! 5 electrodes, 3 mu => 10
    ! 4 electrodes, 4 mu => 12
    ! etc.
    ! For each chemical potential we need the contribution
    ! from each other electrode with different chemical potential
    N_nEq_id = 0
    do i = 1 , N_mu
       N_nEq_id = N_nEq_id + N_Elec - mus(i)%N_El
    end do
    if ( N_nEq_id < 1 ) then
       call neq_die('Could not find any non-equilibrium segments. &
            &Please correct.')
    end if
    allocate(nEq_id(N_nEq_id))
    

    ! populate the segments by their chemical potentials
    N_nEq_id = 0
    cur_mu = 1
    do i = 1 , N_mu - 1
       do j = i + 1 , N_mu
          ! the chemical potentials in the segments
          ! must be sorted! Otherwise the weights are incorrect
          if ( mus(i)%mu < mus(j)%mu ) then
             nEq_segs(cur_mu)%mu1 => mus(i) ! smallest
             nEq_segs(cur_mu)%mu2 => mus(j) ! biggest
          else
             nEq_segs(cur_mu)%mu1 => mus(j) ! smallest
             nEq_segs(cur_mu)%mu2 => mus(i) ! biggest
          end if

          ! create the ID's
          call add_ID(nEq_id,N_nEq_id,nEq_segs(cur_mu))
                   
          ! we need to find the tail contours that fit
          left  = fits_left(nEq_segs(cur_mu)%mu1,tail_io)
          right = fits_right(nEq_segs(cur_mu)%mu2,tail_io)
          !print *,nEq_segs(cur_mu)%mu1,left,right,nEq_segs(cur_mu)%mu2
          if ( left == 0 .or. right == 0 ) &
               call neq_die('Ensure that the tail integrals lower bound fits with &
               &a break at every chemical potential. This will ensure correct &
               &partitioning.')
          if ( left > right ) &
               call neq_die('The contours have not been sorted properly, please &
               &contact the developers')

          ! Allocate the pointers to the non-equilibrium contours
          N = abs(left-right) + 1
          allocate(nEq_segs(cur_mu)%io(N))
          do k = left , right
             nEq_segs(cur_mu)%io(k-left+1) = k
          end do
          nEq_segs(cur_mu)%tail_io(:) = 0
          do k = 1 , N_nEq_tail
             if ( abs(nEq_io(nEq_segs(cur_mu)%io(1))%a - nEq_tail_io(k)%b) < mu_same ) then
                nEq_segs(cur_mu)%tail_io(1) = k
             end if
             if ( abs(nEq_io(nEq_segs(cur_mu)%io(N))%b - nEq_tail_io(k)%a) < mu_same ) then
                nEq_segs(cur_mu)%tail_io(2) = k
             end if
          end do
          if ( any(nEq_segs(cur_mu)%tail_io == 0) ) then
             call neq_die('Could not find all tails')
          end if

          cur_mu = cur_mu + 1
       end do
    end do
    if ( N_nEq_id /= size(nEq_id) ) &
         call neq_die('Error in code')

    ! TODO correct empty cycles, i.e. if two line contours are neighbours 
    ! then we have overlying energy points...

  contains 

    subroutine assign_set_E(c,c_io,c_from,c_io_from,mu,reverse)
      type(ts_cw), intent(inout) :: c
      type(ts_c_io), target :: c_io
      type(ts_c_io) :: c_io_from(:)
      type(ts_cw), intent(in) :: c_from
      type(ts_mu), intent(in) :: mu
      logical, intent(in), optional :: reverse
      logical :: lreverse
      integer :: i, j
      ! assign
      c%c_io => c_io

      lreverse = .false.!mu%mu < 0._dp
      if ( present(reverse) ) lreverse = reverse

      ! correct the end points
      if ( lreverse ) then
         ! we are the lower tail (so reverse it)
         c_io%a = -c_io_from(size(c_io_from))%b + mu%mu
         c_io%b = -c_io_from(1)%a               + mu%mu
      else
         c_io%a =  c_io_from(1)%a               + mu%mu
         c_io%b =  c_io_from(size(c_io_from))%b + mu%mu
      end if

      ! create the contours
      allocate(c%c(c_io%N),c%w(c_io%N,1))
      if ( lreverse ) then
         do j = c_io%N , 1 , -1
            i = c_io%N - j + 1
            c%c(i) = dcmplx(-dreal(c_from%c(j)),dimag(c_from%c(j))) + mu%mu
            c%w(i,1) = c_from%w(j,1)
         end do
      else
         do i = 1 , c_io%N
            c%c(i)   = c_from%c(i) + mu%mu
            c%w(i,1) = c_from%w(i,1)
         end do
      end if
    end subroutine assign_set_E

    function fits_left(mu,tail_io) result(i)
      type(ts_mu), intent(in) :: mu
      type(ts_c_io), intent(in) :: tail_io(:)
      integer :: i
      do i = 1 , N_nEq
         if ( abs(mu%mu - nEq_io(i)%a - tail_io(1)%a) < mu_same ) then
            return
         end if
      end do
      i = 0
    end function fits_left

    function fits_right(mu,tail_io) result(i)
      type(ts_mu), intent(in) :: mu
      type(ts_c_io), intent(in) :: tail_io(:)
      integer :: i
      do i = 1 , N_nEq
         if ( abs(nEq_io(i)%b - mu%mu - tail_io(1)%a) < mu_same ) then
            return
         end if
      end do
      i = 0
    end function fits_right

    subroutine my_setup(suffix,N_nEq,nEq_c,nEq_io)
      character(len=*), intent(in) :: suffix
      integer, intent(inout) :: N_nEq
      type(ts_cw), pointer :: nEq_c(:)
      type(ts_c_io), pointer :: nEq_io(:)

      ! Local variables
      logical :: connected
      integer :: i, j
      character(len=C_N_NAME_LEN), allocatable :: tmp(:)

      N_nEq = fdf_nc_iotype('TS',suffix)
      if ( N_nEq < 1 ) return

      allocate(tmp(N_nEq))

      tmp(1) = fdf_name_c_iotype('TS',suffix,1)
      do i = 2 , N_nEq
         tmp(i) = fdf_name_c_iotype('TS',suffix,i)
         do j = 1 , i - 1
            if ( leqi(tmp(j),tmp(i)) ) then
               call neq_die('You cannot have two names from the bias-window &
                    &to be the same...')
            end if
         end do
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

      do i = 1 , N_nEq - 1
         if ( i == 1 ) then
            call ts_fix_contour(nEq_io(i), next=nEq_io(i+1), &
                 connected = connected )
         else
            call ts_fix_contour(nEq_io(i), &
                 prev=nEq_io(i-1), next=nEq_io(i+1), &
                 connected = connected )
         end if
         if ( .not. connected ) then
            call die('Contour: '//trim(nEq_io(i)%name)// &
                 ' and '//trim(nEq_io(i+1)%name)//' are not connected.')
         end if
      end do
      ! The fix_contour checks and corrects the neighbouring 
      ! contours.
      call ts_fix_contour(nEq_io(N_nEq))

      ! setup the contour
      do i = 1 , N_nEq
         if ( nEq_c(i)%c_io%N < 1 ) then
            write(*,*) 'Contour: '//trim(nEq_c(i)%c_io%Name)//' has 0 points.'
            write(*,*) 'Please ensure at least 1 point in each contour...'
            call die('Contour number of points, invalid')
         end if

         ! allocate contour
         allocate(nEq_c(i)%c(nEq_c(i)%c_io%N),nEq_c(i)%w(nEq_c(i)%c_io%N,1))
         call setup_nEq_contour(nEq_c(i), kT, nEq_Eta)
      end do

      if ( nEq_c(1)%c_io%a > nEq_c(N_nEq)%c_io%b ) then
         call neq_die('The non-equilibrium contours must be in increasing &
              energy. Even if your bias is negative. Please correct.')
      end if

    end subroutine my_setup

    subroutine add_ID(nEq_id,ID,nEq_seg)
      type(ts_nEq_id), intent(in out) :: nEq_id(:)
      integer, intent(inout) :: ID
      type(ts_nEq_seg), target :: nEq_seg
      integer :: i
      do i = 1 , nEq_seg%mu1%N_El
         ID = ID + 1
         if ( size(nEq_id) < ID ) then
            call die('Error in parsing number of non-equilibrium &
                 &electrode contributions')
         end if
         nEq_id(ID)%ID  = ID
         nEq_id(ID)%iEl = nEq_seg%mu1%el(i)
         nEq_id(ID)%imu = nEq_seg%mu2%ID
         nEq_id(ID)%seg => nEq_seg
      end do
      do i = 1 , nEq_seg%mu2%N_El
         ID = ID + 1
         if ( size(nEq_id) < ID ) then
            call die('Error in parsing number of non-equilibrium &
                 &electrode contributions')
         end if
         nEq_id(ID)%ID  = ID
         nEq_id(ID)%iEl = nEq_seg%mu2%el(i)
         nEq_id(ID)%imu = nEq_seg%mu1%ID
         nEq_id(ID)%seg => nEq_seg
      end do
    end subroutine add_ID

  end subroutine read_contour_neq_options

  ! This routine assures that we have setup all the 
  ! equilibrium contours for the passed electrode
  subroutine setup_nEq_contour(c, kT, Eta)
    type(ts_cw), intent(inout) :: c
    real(dp), intent(in) :: kT, Eta

    if ( leqi(c%c_io%part,'line') ) then
       
       call contour_line(c,kT,Eta)
       
    else if ( leqi(c%c_io%part,'tail') ) then
       
       call contour_tail(c,kT,Eta)

    else
       
       call neq_die('Unrecognized contour type for the &
            &non-equilibrium part.')
       
    end if
    
  end subroutine setup_nEq_contour

  function has_cE_neq(cE,iEl,ID) result(has)
    type(ts_c_idx), intent(in) :: cE
    integer, intent(in), optional :: iEl, ID
    logical :: has
    integer :: i
    has = .false.
    if ( cE%fake ) return
    if ( cE%idx(1) /= CONTOUR_NEQ .and. &
         cE%idx(1) /= CONTOUR_NEQ_TAIL ) return
    
    if ( present(ID) .and. present(iEl) ) then
       has = nEq_ID(ID)%iEl == iEl

    else if ( present(iEl) ) then
       
       select case ( cE%idx(1) ) 
       case ( CONTOUR_NEQ )
          do i = 1 , N_nEq_segs
             if ( seg_has_c(nEq_segs(i), cE%idx(2)) ) then
                has = has .or. seg_has_El(nEq_segs(i), iEl)
             end if
          end do
       case ( CONTOUR_NEQ_TAIL ) 
          do i = 1 , N_nEq_segs
             if ( seg_has_tail_c(nEq_segs(i), cE%idx(2)) ) then
                has = has .or. seg_has_El(nEq_segs(i), iEl)
             end if
          end do
       end select
       
    end if

  contains
    
    function seg_has_c(seg,i_c) result(has)
      type(ts_nEq_seg), intent(in) :: seg ! a segment that needs to be tested for 
                                          ! part
      integer, intent(in) :: i_c          ! the index of the part in the list of parts
      logical :: has
      has = allocated(seg%io)
      if ( has ) has = any(i_c == seg%io)
    end function seg_has_c

    function seg_has_tail_c(seg,i_c) result(has)
      type(ts_nEq_seg), intent(in) :: seg ! a segment that needs to be tested for 
                                          ! part
      integer, intent(in) :: i_c          ! the index of the part in the list of parts
      logical :: has
      has = any(i_c == seg%tail_io)
    end function seg_has_tail_c

    function seg_has_El(seg,iEl) result(has)
      type(ts_nEq_seg), intent(in) :: seg
      integer, intent(in) :: iEl
      logical :: has
      has = hasEl(seg%mu1,iEl)
      if ( has ) return
      has = hasEl(seg%mu2,iEl)
    end function seg_has_El

  end function has_cE_neq


  subroutine c2weight_neq(c,kT,iEl,ID, k,W,imu,ZW)
    use m_ts_aux, only : nf
    type(ts_c_idx), intent(in) :: c
    real(dp), intent(in) :: kT ! the temperature
    integer, intent(in) :: iEl, ID ! the electrode index (wrt. Elecs-array)
    real(dp), intent(in) :: k ! generic weight
    integer, intent(out) :: imu
    complex(dp), intent(out) :: W, ZW ! the weight returned
    ! local variables
    real(dp) :: E
    logical :: has_correct_weight, isLeft
    integer :: i
    type(ts_cw), pointer :: cw

    if ( .not. has_cE_neq(c,iEl,ID) ) then
       W  = 0._dp
       ZW = 0._dp
       return
    end if

    has_correct_weight = .false.
    if ( c%idx(1) == CONTOUR_NEQ ) then
       cw => nEq_c(c%idx(2))
    else if ( c%idx(1) == CONTOUR_NEQ_TAIL ) then
       cw => nEq_tail_c(c%idx(2))
       has_correct_weight = &
            any(method(cw%c_io) == (/(i,i=CC_G_NF_MIN,CC_G_NF_MAX)/)) 
    else
       print *,c%idx
       call die('c2weight_neq: Error in code')
    end if

    E = real(cw%c(c%idx(3)),dp)

    imu = nEq_ID(ID)%imu
    isLeft = hasEl(nEq_ID(ID)%seg%mu1,iEl)
    if ( .not. isLeft ) then
       if ( .not. hasEl(nEq_ID(ID)%seg%mu2,iEl) ) then
          call die('c2weight_neq: Error in code')
       end if
    end if

    ! Notice that we multiply with -i to move Gf.Gamma.Gf^\dagger
    ! to the negative imaginary part
    ! This means that we do not need to create different
    ! routines for the equilibrium density and the non-equilibrium
    ! density

    ! nf function is: nF(E-E1) - nF(E-E2) IMPORTANT
    if ( has_correct_weight ) then
       ! the gauss-fermi contour has the "correct" weight already...
       W = -k * cw%w(c%idx(3),1) * dcmplx(0._dp,-1._dp)
    else if ( isLeft ) then
       ! This weight will always be negative due to the sorting
       ! of the chemical potentials
       W =  k * cw%w(c%idx(3),1) * &
            nf(E, &
            nEq_ID(ID)%seg%mu1%mu, &
            nEq_ID(ID)%seg%mu2%mu, kT) * dcmplx(0._dp,-1._dp)
    else
       ! This weight will always be positive due to the sorting
       ! of the chemical potentials
       W =  k * cw%w(c%idx(3),1) * &
            nf(E, &
            nEq_ID(ID)%seg%mu2%mu, &
            nEq_ID(ID)%seg%mu1%mu, kT) * dcmplx(0._dp,-1._dp)
    end if

    if ( isLeft ) then
       ZW = E * W * nEq_ID(ID)%seg%mu2%N_el
    else
       ZW = E * W * nEq_ID(ID)%seg%mu1%N_el
    end if

  end subroutine c2weight_neq

  subroutine cseq2weight_neq(c,kT,seg,W)
    use m_ts_aux, only : nf
    type(ts_c_idx), intent(in) :: c
    real(dp), intent(in) :: kT ! the temperature
    type(ts_nEq_seg), intent(in) :: seg
    complex(dp), intent(out) :: W
    ! local variables
    real(dp) :: E
    type(ts_cw), pointer :: cw
    logical :: has_correct_weight
    integer :: i

    has_correct_weight = .false.
    if ( c%idx(1) == CONTOUR_NEQ ) then
       cw => nEq_c(c%idx(2))
    else if ( c%idx(1) == CONTOUR_NEQ_TAIL ) then
       cw => nEq_tail_c(c%idx(2))
       has_correct_weight = &
            any(method(cw%c_io) == (/(i,i=CC_G_NF_MIN,CC_G_NF_MAX)/)) 
    else
       call die('cseq2weight_neq: Error in code')
    end if

    E = real(cw%c(c%idx(3)),dp)

    if ( has_correct_weight ) then
       ! the gauss-fermi contour has the "correct" weight already...
       W =  cw%w(c%idx(3),1)
    else
       ! nf function is: nF(E-E1) - nF(E-E2) IMPORTANT
       ! We use this to get the positive weight (the mu's are sorted in descending order)
       W =  cw%w(c%idx(3),1) * &
            nf(E, &
            seg%mu2%mu, &
            seg%mu1%mu, kT)
    end if
    
  end subroutine cseq2weight_neq

  function ID2mu(ID) result(imu)
    integer, intent(in) :: ID
    integer :: imu
    imu = nEq_ID(ID)%imu
  end function ID2mu

  subroutine contour_line(c,kT,Eta)
    use m_integrate
    use m_gauss_quad
    type(ts_cw), intent(inout) :: c
    real(dp), intent(in) :: kT, Eta

    ! local variables
    character(len=c_N) :: tmpC
    real(dp) :: a,b, tmp
    real(dp), allocatable :: ce(:), cw(:)

    if ( .not. leqi(c%c_io%part,'line') ) &
         call die('Contour is not a line')

    if ( c%c_io%N < 1 ) then
       call die('Contour: '//trim(c%c_io%Name)//' has &
            an errorneous number of points.')
    end if

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

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

       ! we should also gain an option for this
       if ( c_io_has_opt(c%c_io,'precision') ) then
          tmpC = c_io_get_opt(c%c_io,'precision')
          read(tmpC,'(g20.10)') tmp
       else
          tmp = 2.e-2_dp * abs(b-a) / real(c%c_io%N,dp)
          write(tmpC,'(g20.10)') tmp
          call c_io_add_opt(c%c_io,'precision',tmpC)
       end if

       call TanhSinh_Exact(c%c_io%N,ce,cw,a,b, p=tmp)
       
    case default

       call die('Could not determine the line-integral')

    end select

    c%c = dcmplx(ce,Eta)
    c%w(:,1) = dcmplx(cw,0._dp)

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
    type(ts_cw), intent(inout) :: c
    real(dp), intent(in) :: kT, Eta

    ! local variables
    integer :: ioffset, infinity
    real(dp) :: a,b
    real(dp), allocatable :: ce(:), cw(:)

    if ( .not. leqi(c%c_io%part,'tail') ) &
         call die('Contour is not a tail contour')

    if ( c%c_io%N < 1 ) then
       call die('Contour: '//trim(c%c_io%Name)//' has &
            an errorneous number of points.')
    end if

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b

    if ( b < 0._dp ) then
       call die('The non-equilbrium tail contours can only be &
            &defined with respect to the Fermi-level')
    end if

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
       c%w(:,1) = dcmplx(cw,0._dp)

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
    type(ts_c_idx) :: c ! the configuration of the energy-segment
    integer :: lstep, i, PN
    lstep = 1
    if ( present(step) ) lstep = step
    PN = N_nEq_E()
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
  end function nEq_E

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
    do j = 1 , N_nEq ! number of contours
       if ( iE + nEq_c(j)%c_io%N < id ) then
          iE = iE + nEq_c(j)%c_io%N
          cycle
       end if
       i = id - iE
       if ( i <= nEq_c(j)%c_io%N ) then
          c%exist = .true.
          c%e     = nEq_c(j)%c(i)
          c%idx(1) = 2 ! designates the non-equilibrium contours
          c%idx(2) = j ! designates the index of the non-equilibrium contour
          c%idx(3) = i ! is the index of the non-equilibrium contour
          return
       end if
    end do

    do j = 1 , N_nEq_tail ! number of contours
       if ( iE + nEq_tail_c(j)%c_io%N < id ) then
          iE = iE + nEq_tail_c(j)%c_io%N
          cycle
       end if
       i = id - iE
       if ( i <= nEq_tail_c(j)%c_io%N ) then
          c%exist = .true.
          c%e     = nEq_tail_c(j)%c(i)
          c%idx(1) = 3 ! designates the tail non-equilibrium contours
          c%idx(2) = j ! designates the index of the tail non-equilibrium contour
          c%idx(3) = i ! is the index of the tail non-equilibrium contour
          return
       end if
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
       write(*,opt_c) 'Contour name',trim(prefix)// &
            '.Contour.Bias.Window.'//trim(neq_io(i)%name)
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
          if ( len_trim(opt%val) > 0 ) then
             write(*,opt_cc) '   Option for contour method',trim(opt%opt),trim(opt%val)
          else
             write(*,opt_c)  '   Option for contour method',trim(opt%opt)
          end if
          opt => opt%next
       end do
    end do

  end subroutine print_contour_neq_options

  subroutine io_contour_neq(slabel,kT,suffix)
    use parallel, only : IONode
    character(len=*), intent(in) :: slabel
    real(dp), intent(in) :: kT
    character(len=*), intent(in), optional :: suffix

! *********************
! * LOCAL variables   *
! *********************
    character(len=25) :: tmp_suffix
    integer :: i

    if ( .not. IONode ) return

    do i = 1 , N_nEq_segs
       
       if ( present(suffix) ) then
          write(tmp_suffix,'(a,i0)') trim(suffix)//'-',i
       else
          write(tmp_suffix,'(a,i0)') 'TSCCNEQ-',i
       end if
       call io_contour_neq_seg(kT,nEq_segs(i),slabel,tmp_suffix)
       
    end do

  end subroutine io_contour_neq


  subroutine io_contour_neq_seg(kT,seg,slabel,suffix)
    use parallel, only : IONode
    use units, only : eV
    real(dp), intent(in) :: kT
    type(ts_nEq_seg), intent(in) :: seg
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in) :: suffix

! *********************
! * LOCAL variables   *
! *********************
    character(len=200) :: fname
    integer        :: i, unit
    type(ts_c_idx) :: cidx
    
    if ( .not. IONode ) return
    
    fname = trim(slabel)//'.'//trim(suffix)

    call io_assign( unit )
    open( unit, file=fname, status='unknown' )
    write(unit,'(a)') '# Contour path for the non-equilibrium part'
    write(unit,'(a)') '# Segment between following chemical potentials:'
    write(unit,'(a,3(tr1,a))') '#',trim(seg%mu1%name),'and',trim(seg%mu2%name)
    write(unit,'(a,2(tr1,f10.5),tr1,a)') '#',seg%mu1%mu/eV,seg%mu2%mu/eV,'eV'

    write(unit,'(a,a12,2(tr1,a13))') '#','Re(c) [eV]','Im(c) [eV]','w'

    cidx%idx(1) = CONTOUR_NEQ_TAIL
    if ( seg%tail_io(1) /= 0 ) then
       cidx%idx(2) = seg%tail_io(1)
       call io_contour_c(unit,kT,seg,cidx)
    end if

    cidx%idx(1) = CONTOUR_NEQ
    do i = 1 , size(seg%io)
       cidx%idx(2) = seg%io(i)
       call io_contour_c(unit,kT,seg,cidx)
    end do

    cidx%idx(1) = CONTOUR_NEQ_TAIL
    if ( seg%tail_io(2) /= 0 ) then
       cidx%idx(2) = seg%tail_io(2)
       call io_contour_c(unit,kT,seg,cidx)
    end if
    
    call io_close( unit )

  end subroutine io_contour_neq_seg

  subroutine io_contour_c(unit,kT,seg,cidx)
    use units,    only : eV
    use m_ts_aux, only : nf
    integer, intent(in) :: unit
    real(dp), intent(in) :: kT
    type(ts_nEq_seg), intent(in) :: seg
    type(ts_c_idx), intent(inout) :: cidx
    type(ts_cw), pointer :: c
    integer :: i
    complex(dp) :: W

    if ( cidx%idx(1) == CONTOUR_NEQ ) then
       c => nEq_c(cidx%idx(2))
    else if ( cidx%idx(1) == CONTOUR_NEQ_TAIL ) then
       c => nEq_tail_c(cidx%idx(2))
    else
       call die('io_contour_c: Error in code')
    end if

    do i = 1 , size(c%c)
       cidx%e = c%c(i)
       cidx%idx(3) = i
       call cseq2weight_neq(cidx,kT,seg,W)
       if ( abs(aimag(W)) > 1.e-10_dp ) then
          call die('Error in the imaginary weight of &
               a non-equilibrium integration point.')
       end if
       write(unit,'(3(e13.6,tr1))') c%c(i) / eV, real(W,dp)
    end do

  end subroutine io_contour_c

  subroutine neq_die(msg)
    character(len=*), intent(in) :: msg
    
    write(*,*) 'Killing... printing out so-far gathered information'
    call print_contour_neq_options('TS')
    call print_contour_neq_block('TS')
    call die(msg)
  end subroutine neq_die
    
end module m_ts_contour_neq
