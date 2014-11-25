module m_tbt_contour

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

  implicit none

  ! Contour path
  integer, save, public :: N_tbt
  type(ts_c_io),  pointer, save, public :: tbt_io(:) => null()
  type(ts_cw), pointer, save, public :: tbt_c(:) => null()

  ! The contour specific variables
  real(dp), save, public :: tbt_Eta

  ! this is heavily linked with the CONTOUR_EQ from m_ts_contour_eq
  integer, parameter, public :: CONTOUR_TBT = 5

  public :: tbt_read_contour_options
  public :: print_contour_tbt_options
  public :: print_contour_tbt_block
  public :: io_contour_tbt
  public :: N_tbt_E
  public :: tbt_E

  private

contains

  subroutine tbt_read_contour_options(N_Elec,Elecs,N_mu,mus, kT)

    use parallel, only : Node
    use units, only : eV
    use fdf
    use m_ts_electype

    ! only to gain access to the chemical shifts
    integer, intent(in) :: N_Elec
    type(Elec), intent(in), target :: Elecs(N_Elec)
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in), target :: mus(N_mu)
    real(dp), intent(in) :: kT
    
    real(dp) :: Volt, tmp

    ! broadening
    tbt_Eta = fdf_get('TBT.Contours.Eta',0._dp,'Ry')
    if ( tbt_Eta < 0._dp ) call die('ERROR: Eta < 0, we do not allow &
         &for using the advanced Greens function, please correct.')

    ! We only allow the user to either use the old input format, or the new
    ! per-electrode input

    ! Bias-window setup
    call my_setup('Window',N_tbt,tbt_c,tbt_io)
    
    tmp = fdf_get('TS.Voltage',0._dp,'Ry')
    Volt = fdf_get('TBT.Voltage',tmp,'Ry')
    if ( abs(tmp-Volt) > 1.e-5_dp .and. Node == 0 ) then
       write(*,'(a)') '*** WARNING: transiesta and tbtrans bias &
            &not equivalent!'
       write(*,'(a)') '*** WARNING: Be sure to use an interpolation scheme!'
    end if

    if ( N_tbt < 1 ) then

       ! We should default to the bias window!

       if ( N_Elec /= 2 .or. N_mu /= 2 ) then
          call die('When using other than 2 electrodes and &
               &chemical potentials you are forced to setup &
               &the input yourself.')
       end if

       ! We create the default version
       ! *** NOTE requires an even splitting of the bias
       N_tbt = 1
       nullify(tbt_io,tbt_c)
       allocate(tbt_io(N_tbt),tbt_c(N_tbt))
       tbt_c(1)%c_io => tbt_io(1)
       tbt_io(1)%part = 'line'
       tbt_io(1)%name = 'neq'
       tbt_io(1)%ca = '-|V|/2 - 5. kT'
       tbt_io(1)%a  = -abs(Volt) * .5_dp - 5._dp * kT
       tbt_io(1)%cb = '|V|/2 + 5. kT'
       tbt_io(1)%b  =  abs(Volt) * .5_dp + 5._dp * kT
       ! number of points
       tbt_io(1)%cd = '0.01 eV'
       tbt_io(1)%d = 0.01_dp * eV
       tbt_io(1)%method = 'mid-rule'
       call ts_fix_contour(tbt_io(1))

       ! setup the contour
       allocate(tbt_c(1)%c(tbt_c(1)%c_io%N),tbt_c(1)%w(tbt_c(1)%c_io%N,1))
       call setup_tbt_contour(tbt_c(1), kT, tbt_Eta)

    end if

    ! TODO correct empty cycles, i.e. if two line contours are neighbours 
    ! then we have overlying energy points...

  contains 

    subroutine my_setup(suffix,N_tbt,tbt_c,tbt_io)
      character(len=*), intent(in) :: suffix
      integer, intent(inout) :: N_tbt
      type(ts_cw), pointer :: tbt_c(:)
      type(ts_c_io), pointer :: tbt_io(:)

      ! Local variables
      integer :: i, j
      character(len=C_N_NAME_LEN), allocatable :: tmp(:)
      logical :: connected

      N_tbt = fdf_nc_iotype('TBT',suffix)
      if ( N_tbt < 1 ) return

      allocate(tmp(N_tbt))

      tmp(1) = fdf_name_c_iotype('TBT',suffix,1)
      do i = 2 , N_tbt
         tmp(i) = fdf_name_c_iotype('TBT',suffix,i)
         do j = 1 , i - 1
            if ( leqi(tmp(j),tmp(i)) ) then
               call neq_die('You cannot have two names from the window &
                    &to be the same...')
            end if
         end do
      end do
      
      ! allocate all required objects
      nullify(tbt_io,tbt_c)
      allocate(tbt_io(N_tbt),tbt_c(N_tbt))
      
      do i = 1 , N_tbt

         ! assign pointer
         tbt_c(i)%c_io => tbt_io(i)
         ! read in the contour
         call ts_read_contour_block('TBT',suffix,tmp(i),tbt_io(i), kT, Volt)

      end do
      deallocate(tmp)

      do i = 1 , N_tbt - 1
         if ( i == 1 ) then
            call ts_fix_contour(tbt_io(i), next=tbt_io(i+1) , &
                 connected=connected )
         else
            call ts_fix_contour(tbt_io(i), &
                 prev=tbt_io(i-1), next=tbt_io(i+1), &
                 connected=connected )
         end if
         if ( Node == 0 .and. .not. connected ) then
            write(*,'(a)') 'tbtrans: *** Contours are not connected ***'
         end if
            
      end do
      call ts_fix_contour(tbt_io(N_tbt))

      ! setup the contour
      do i = 1 , N_tbt
         if ( tbt_c(i)%c_io%N < 1 ) then
            write(*,*) 'Contour: '//trim(tbt_c(i)%c_io%Name)//' has 0 points.'
            write(*,*) 'Please ensure at least 1 point in each contour...'
            call die('Contour number of points, invalid')
         end if

         ! allocate contour
         allocate(tbt_c(i)%c(tbt_c(i)%c_io%N),tbt_c(i)%w(tbt_c(i)%c_io%N,1))
         call setup_tbt_contour(tbt_c(i), kT, tbt_Eta)

      end do

    end subroutine my_setup

  end subroutine tbt_read_contour_options
                                         

  ! This routine assures that we have setup all the 
  ! equilibrium contours for the passed electrode
  subroutine setup_tbt_contour(c, kT, Eta)
    type(ts_cw), intent(inout) :: c
    real(dp), intent(in) :: kT, Eta

    if ( leqi(c%c_io%part,'line') ) then
       
       call contour_line(c,kT,Eta)
       
    else
       
       call neq_die('Unrecognized contour type for &
            &tbtrans, MUST be a line part.')
       
    end if
    
  end subroutine setup_tbt_contour

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

  function TBT_E(id,step) result(c)
    integer, intent(in) :: id
    integer, intent(in), optional :: step
    type(ts_c_idx) :: c ! the configuration of the energy-segment
    integer :: lstep, i, PN
    lstep = 1
    if ( present(step) ) lstep = step
    PN = N_tbt_E()
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
  end function TBT_E

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
    do j = 1 , N_tbt ! number of contours
       if ( iE + tbt_c(j)%c_io%N < id ) then
          iE = iE + tbt_c(j)%c_io%N
          cycle
       end if
       i = id - iE
       if ( i <= tbt_c(j)%c_io%N ) then
          c%exist = .true.
          c%e     = tbt_c(j)%c(i)
          c%idx(1) = 2 ! designates the non-equilibrium contours
          c%idx(2) = j ! designates the index of the non-equilibrium contour
          c%idx(3) = i ! is the index of the non-equilibrium contour
          return
       end if
    end do

  end function get_c

  function N_TBT_E() result(N)
    integer :: N, i
    N = 0
    do i = 1 , N_tbt
       N = N + size(tbt_c(i)%c)
    end do
  end function N_TBT_E

  subroutine print_contour_tbt_block(prefix)
    use parallel, only : IONode
    character(len=*), intent(in) :: prefix

    integer :: i

    if ( IONode ) then
       write(*,'(2a)') '%block ',trim(prefix)//'.Contours.Window'
       do i = 1 , N_tbt
          write(*,'(tr4,a)') trim(tbt_io(i)%name)
       end do
       write(*,'(2a,/)') '%endblock ',trim(prefix)//'.Contours.Window'
    end if

    do i = 1 , N_tbt
       call ts_print_contour_block(trim(prefix)//'.Contour.Window.',tbt_io(i))
    end do

  end subroutine print_contour_tbt_block


  subroutine print_contour_tbt_options(prefix)

    use parallel, only : IONode
    use units, only : eV

    use m_ts_io_contour

    character(len=*), intent(in) :: prefix
    character(len=200) :: chars
    integer :: i
    type(ts_c_opt_ll), pointer :: opt

    if ( .not. IONode ) return
    
    write(*,opt_n) '             >> TBtrans contour << '
    write(*,opt_g_u) 'non-Equilibrium Greens function Eta',tbt_Eta/eV,'eV'
    do i = 1 , N_tbt
       chars = '  '//trim(tbt_io(i)%part)
       write(*,opt_c) 'Contour name',trim(prefix)// &
            '.Contour.Window.'//trim(tbt_io(i)%name)
       call write_e(trim(chars)//' contour E_min',tbt_io(i)%a)
       call write_e(trim(chars)//' contour E_max',tbt_io(i)%b)
       write(*,opt_int) trim(chars)//' contour points',tbt_io(i)%N
       write(*,opt_c) trim(chars)//' contour method', &
            trim(longmethod2str(tbt_io(i)))
       opt => tbt_io(i)%opt
       do while ( associated(opt) )
          write(*,opt_c) '   Option for contour method', trim(opt%opt)
          opt => opt%next
       end do
    end do

  end subroutine print_contour_tbt_options

  subroutine io_contour_tbt(slabel,kT,suffix)
    use parallel, only : IONode
    character(len=*), intent(in) :: slabel
    real(dp), intent(in) :: kT
    character(len=*), intent(in), optional :: suffix

! *********************
! * LOCAL variables   *
! *********************
    integer :: iu, i
    type(ts_c_idx) :: cidx

    if ( .not. IONode ) return

    call io_assign( iu )
    open( iu, file=trim(slabel)//'.TBTCC', status='unknown' )
    write(iu,'(a)') '# Contour path for the transport part'
    write(iu,'(a,a12,2(tr1,a13))') '#','Re(c) [eV]','Im(c) [eV]'

    cidx%idx(1) = CONTOUR_TBT

    do i = 1 , N_tbt
       
       cidx%idx(2) = i
       call io_contour_c(iu,kT,cidx)
       
    end do

  end subroutine io_contour_tbt

  subroutine io_contour_c(iu,kT,cidx)
    use units,    only : eV
    use m_ts_aux, only : nf
    integer, intent(in) :: iu
    real(dp), intent(in) :: kT
    type(ts_c_idx), intent(inout) :: cidx
    type(ts_cw), pointer :: c
    integer :: i

    if ( cidx%idx(1) == CONTOUR_TBT ) then
       c => tbt_c(cidx%idx(2))
    else
       call die('io_contour_c: Error in code')
    end if

    do i = 1 , size(c%c)
       write(iu,'(2(e13.6,tr1))') c%c(i) / eV
    end do

  end subroutine io_contour_c

  subroutine neq_die(msg)
    character(len=*), intent(in) :: msg
    
    write(*,*) 'Killing... printing out so-far gathered information'
    call print_contour_tbt_options('TBT')
    call print_contour_tbt_block('TBT')
    call die(msg)

  end subroutine neq_die
    
end module m_tbt_contour
