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

  ! Number of spin components
  integer, save, public :: nspin

  integer, save, public :: spinor_dim   ! Spin dimension of electronic states
                                        ! Used for sizing eo and qo,
                                        ! efs and qfs
  integer, save, public :: h_spin_dim   ! Number of spin components in H and D
  integer, save, public :: e_spin_dim   ! Number of spin components in E_dm
  integer, save, public :: MColl        ! Spin size of Hamiltonian matrix that
                                        ! is to be diagonalized by cdiag

  ! Different Fermi-levels for different fixed spin-components
  real(dp), pointer, save, public  :: efs(:)
  real(dp), pointer, save, public  :: qs(:)

  ! Whether we are performing spiral arrangement of spins
  logical, save, public :: Spiral  = .false.
  ! Pitch wave vector for spiral configuration
  real(dp), save, public :: qSpiral(3) = 0._dp

  public :: init_spin

contains
  
  subroutine init_spin( ucell )
    
    use sys, only: die
    use fdf, only : fdf_get, leqi
    use alloc, only: re_alloc
    use parallel, only: IOnode, Nodes

    ! Unit cell lattice vectors
    real(dp), intent(in) :: ucell(3,3)

    character(len=32) :: opt

    logical :: TRSym

    ! All default to false
    NoMagn = .false.
    SPpol  = .false.
    NonCol = .false.
    SpOrb  = .false.

    ! Time reversal symmetry
    TRSym  = fdf_get('TimeReversalSymmetry',.false.)

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
    ! 4! different options we use a single value now)
    opt = fdf_get('Magnetization', opt)

    if ( leqi(opt, 'non-magnetic') .or. &
         leqi(opt, 'non-polarized') .or. &
         leqi(opt, 'non-polarised') .or. &
         leqi(opt, 'NP') .or. leqi(opt,'N-P') .or. &
         leqi(opt, 'NM') .or. leqi(opt, 'N-M') ) then
       NoMagn = .true.
    else if ( leqi(opt, 'magnetic') .or. &
         leqi(opt, 'M') .or. leqi(opt, 'polarized') .or. &
         leqi(opt, 'P') .or. leqi(opt, 'polarised') ) then
       SPpol = .true.
    else if ( leqi(opt, 'non-collinear') .or. &
         leqi(opt, 'NC') .or. leqi(opt, 'N-C') ) then
       NonCol = .true.
    else if ( leqi(opt, 'spin-orbit') .or. leqi(opt, 'S-O') .or. &
         leqi(opt, 'SOC') .or. leqi(opt, 'SO') ) then
       SpOrb = .true.
    end if
       

    if ( SpOrb ) then
       opt        = 'spin-orbit'
       NoMagn     = .false.
       SPpol      = .false.
       NonCol     = .false.
       SpOrb      = .true.
       TRSym      = .false.
       nspin      = 4
       h_spin_dim = 8
       e_spin_dim = 4
       spinor_dim = 2
       MColl      = 2
       
    else if ( NonCol ) then
       opt        = 'non-collinear'
       NoMagn     = .false.
       SPpol      = .false.
       NonCol     = .true.
       SpOrb      = .false.
       TRSym      = .false.
       nspin      = 4
       h_spin_dim = 4
       e_spin_dim = 4
       spinor_dim = 2
       MColl      = 2
       
    else if ( SPpol ) then
       opt        = 'polarized'
       NoMagn     = .false.
       SPpol      = .true.
       NonCol     = .false.
       SpOrb      = .false.
       TRSym      = .true.
       nspin      = 2
       h_spin_dim = 2
       e_spin_dim = 2
       spinor_dim = 2
       MColl      = 1

    else 
       opt        = 'non-polarized'
       NoMagn     = .true.
       SPpol      = .false.
       NonCol     = .false.
       SpOrb      = .false.
       TRSym      = .true.
       nspin      = 1
       h_spin_dim = 1
       e_spin_dim = 1
       spinor_dim = 1
       MColl      = 1

    end if

    ! Read spin-spiral settings
    call init_spiral()

    if ( IONode ) then
     write(*,'(2a)')  'redata: Magnetization                    = ',trim(opt)
     write(*,'(a,i0)')'redata: Number of spin components        = ',nspin
     write(*,'(a,l1)')'redata: Time-Reversal Symmetry           = ',TRSym
     write(*,'(a,l1)')'redata: Spin-spiral                      = ',Spiral
     if ( Spiral .and. .not. NonCol ) then
        write(*,'(a)') 'redata: WARNING: spin-spiral requires non-collinear spin'
     end if
    end if

    nullify(efs,qs)
    call re_alloc(efs,1,nspin,name="efs",routine="init_spin")
    call re_alloc(qs,1,nspin,name="qs",routine="init_spin")

  contains

    subroutine init_spiral
      use fdf, only: block_fdf, parsed_line
      use fdf, only: fdf_block, fdf_bline
      use fdf, only: fdf_bnames, fdf_bvalues
      use units, only: Pi
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
    
  end subroutine init_spin
  
end module m_spin
