! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_spin
  use precision, only: dp

  implicit none
  
  private

  logical, save, public :: NoMagn  = .false.
  logical, save, public :: SPpol   = .false.
  logical, save, public :: NonCol  = .false.
  logical, save, public :: SpOrb   = .false.
  logical, save, public :: TrSym  = .true.

  ! Number of spin components
  integer, save, public :: nspin

  integer, save, public :: spinor_dim   ! Spin dimension of electronic states
                                        ! Used for sizing eo and qo,
                                        ! efs and qfs
  integer, save, public :: h_spin_dim   ! Number of spin components in H and D
  integer, save, public :: e_spin_dim   ! Number of spin components in E_dm

  ! Different Fermi-levels for different fixed spin-components
  real(dp), pointer, save, public  :: efs(:)
  real(dp), pointer, save, public  :: qs(:)

  ! Whether we are performing spiral arrangement of spins
  logical, save, public :: Spiral  = .false.
  ! Pitch wave vector for spiral configuration
  real(dp), save, public :: qSpiral(3) = 0._dp

  public :: init_spin
  public :: print_spin
  public :: init_spiral

contains
  
  subroutine init_spin()
    
    use sys, only: die
    use fdf, only : fdf_get, leqi
    use alloc, only: re_alloc

    character(len=32) :: opt

    ! All default to false
    NoMagn = .false.
    SPpol  = .false.
    NonCol = .false.
    SpOrb  = .false.

    ! Time reversal symmetry
    TRSym  = fdf_get('TimeReversalSymmetry',.true.)

    ! Read in old flags:
    SPpol  = fdf_get('SpinPolarized',.false.)
    NonCol = fdf_get('NonCollinearSpin',.false.)
    SpOrb  = fdf_get('SpinOrbit',.false.)

    ! Set default option from "old" options
    if ( SpOrb ) then
       opt = 'spin-orbit'
    else if ( NonCol ) then
       opt = 'non-collinear'
    else if ( SPpol ) then
       opt = 'polarized'
    else
       opt = 'non-polarized'
    end if

    
    ! In order to enable text input (and obsolete the
    ! 4 different options we use a single value now)
    opt = fdf_get('Spin', opt)

    if ( leqi(opt, 'non-polarized') .or. &
         leqi(opt, 'non-polarised') .or. &
         leqi(opt, 'NP') .or. leqi(opt,'N-P') ) then
       NoMagn = .true.
    else if ( leqi(opt, 'polarized') .or. &
         leqi(opt, 'polarised') .or. leqi(opt, 'P') ) then
       SPpol = .true.
    else if ( leqi(opt, 'non-collinear') .or. &
         leqi(opt, 'NC') .or. leqi(opt, 'N-C') ) then
       NonCol = .true.
    else if ( leqi(opt, 'spin-orbit') .or. leqi(opt, 'S-O') .or. &
         leqi(opt, 'SOC') .or. leqi(opt, 'SO') ) then
       SpOrb = .true.
    else
       write(*,*) 'Unknown spin flag: ', trim(opt)
       call die('Spin: unknown flag, please assert the correct input.')
    end if

    ! Note that, in what follows,
    !   spinor_dim = min(h_spin_dim,2)
    !   e_spin_dim = min(h_spin_dim,4)
    !   nspin      = min(h_spin_dim,4)  ! Relevant for dhscf, local DM
    !      should probably be called nspin_grid
    !
    ! so everything can be determined if h_spin_dim is known.
    ! It is tempting to go back to the old 'nspin' overloading,
    ! making 'nspin' mean again 'h_spin_dim'.
    ! But this has to be done carefully, as some routines expect
    ! an argument 'nspin' which really means 'spinor_dim' (like diagon),
    ! and others (such as dhscf) expect 'nspin' to mean 'nspin_grid'.

    if ( SpOrb ) then
       NoMagn     = .false.
       SPpol      = .false.
       NonCol     = .false.
       SpOrb      = .true.
       TRSym      = .false.
       nspin      = 4
       h_spin_dim = 8
       e_spin_dim = 4
       spinor_dim = 2
       
    else if ( NonCol ) then
       NoMagn     = .false.
       SPpol      = .false.
       NonCol     = .true.
       SpOrb      = .false.
       TRSym      = .false.
       nspin      = 4
       h_spin_dim = 4
       e_spin_dim = 4
       spinor_dim = 2
       
    else if ( SPpol ) then
       NoMagn     = .false.
       SPpol      = .true.
       NonCol     = .false.
       SpOrb      = .false.
       TRSym      = .true.
       nspin      = 2
       h_spin_dim = 2
       e_spin_dim = 2
       spinor_dim = 2

    else 
       NoMagn     = .true.
       SPpol      = .false.
       NonCol     = .false.
       SpOrb      = .false.
       TRSym      = .true.
       nspin      = 1
       h_spin_dim = 1
       e_spin_dim = 1
       spinor_dim = 1

    end if

    nullify(efs,qs)
    call re_alloc(efs,1,spinor_dim,name="efs",routine="init_spin")
    call re_alloc(qs,1,spinor_dim,name="qs",routine="init_spin")

  end subroutine init_spin

  ! Print out spin-configuration options
  subroutine print_spin( )
    use parallel, only: IONode

    character(len=32) :: opt

    if ( .not. IONode ) return

    if ( SpOrb ) then
       opt        = 'spin-orbit'
    else if ( NonCol ) then
       opt        = 'non-collinear'
    else if ( SPpol ) then
       opt        = 'polarized'
    else 
       opt        = 'non-polarized'
    end if

    write(*,'(2a)')  'redata: Spin configuration               = ',trim(opt)
    write(*,'(a,i0)')'redata: Number of spin components        = ',h_spin_dim
    write(*,'(a,l1)')'redata: Time-Reversal Symmetry           = ',TRSym
    write(*,'(a,l1)')'redata: Spin-spiral                      = ',Spiral
    if ( Spiral .and. .not. NonCol ) then
       write(*,'(a)') 'redata: WARNING: spin-spiral requires non-collinear spin'
    end if

    if ( SpOrb ) then
       write(*,'(a)') repeat('#',60)
       write(*,'(a,t16,a,t60,a)') '#','Spin-orbit coupling is in beta','#'
       write(*,'(a,t13,a,t60,a)') '#','Several options may not be compatible','#'
       write(*,'(a)') repeat('#',60)
    else if ( NonCol ) then
       write(*,'(a)') repeat('#',60)
       write(*,'(a,t17,a,t60,a)') '#','Non-collinear spin is in beta','#'
       write(*,'(a,t13,a,t60,a)') '#','Several options may not be compatible','#'
       write(*,'(a)') repeat('#',60)
    end if

  end subroutine print_spin
  

  subroutine init_spiral( ucell )
    use fdf, only : fdf_get, leqi
    use fdf, only: block_fdf, parsed_line
    use fdf, only: fdf_block, fdf_bline
    use fdf, only: fdf_bnames, fdf_bvalues
    use units, only: Pi

    ! Unit cell lattice vectors
    real(dp), intent(in) :: ucell(3,3)
    
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    ! Reciprocal cell vectors
    real(dp) :: rcell(3,3), alat
    character(len=30) :: lattice

    ! read in lattice constant
    alat = fdf_get('LatticeConstant',0.0_dp,'Bohr')

    call reclat( ucell, rcell, 1 )

    Spiral = fdf_block('SpinSpiral', bfdf)

    if ( .not. Spiral ) return

    if (.not. fdf_bline(bfdf,pline)) &
         call die('init_spiral: ERROR in SpinSpiral block')

    ! Read lattice
    lattice = fdf_bnames(pline,1)
    ! Read pitch wave-vector
    qSpiral(1) = fdf_bvalues(pline,1)
    qSpiral(2) = fdf_bvalues(pline,2)
    qSpiral(3) = fdf_bvalues(pline,3)

    if ( leqi(lattice,'Cubic') ) then
       qSpiral(1) = Pi * qs(1) / alat
       qSpiral(2) = Pi * qs(2) / alat
       qSpiral(3) = Pi * qs(3) / alat
    else if ( leqi(lattice,'ReciprocalLatticeVectors') ) then
       qSpiral = matmul(rcell,qSpiral)
    else
       call die('init_spiral: ERROR: ReciprocalCoordinates must be' // &
            ' ''Cubic'' or ''ReciprocalLatticeVectors''')
    end if

  end subroutine init_spiral
  
end module m_spin
