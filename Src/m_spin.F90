module m_spin
  use precision, only: dp
  private

  ! CC RC Added
  logical, save, public           :: NoMagn    = .false.
  logical, save, public           :: SPpol     = .false.
  logical, save, public           :: NonCol    = .false.
  logical, save, public           :: SpOrb     = .false.
  logical, save, public           :: TRSym     = .false.
  ! CC RC Added

  ! Number of spin components
  integer, save, public           :: nspin

  ! CC RC  Added
  integer, save, public           :: spinor_dim     ! Spin dimension of electronic states
                                                    ! Used for sizing eo and qo,
                                                    ! efs and qfs
  integer, save, public           :: h_spin_dim     ! Number of spin components in H and D
  integer, save, public           :: e_spin_dim     ! Number of spin components in E_dm
                                                    ! Could possibly be made
                                                    ! equal to spinor_dim
  integer, save, public           :: MColl          ! Spin size of Hamiltonian matrix that  
                                                    ! is to be diagonalized by cdiag
  ! CC RC Added


  real(dp), pointer, save, public  :: efs(:)
  real(dp), pointer, save, public  :: qs(:)

  public :: init_spin

CONTAINS

  subroutine init_spin()
! CC RC Added: m_fdf_global and sys
!    use m_fdf_global
    use sys, only: die

    use fdf, only : fdf_get
    use alloc, only: re_alloc
    use parallel, only: IOnode, Nodes ! CC RC  Added: Nodes

    implicit none

!    logical  :: SPpol, NonCol

!    sppol  = fdf_get('SpinPolarized',.false.)
!    noncol = fdf_get('NonCollinearSpin',.false.)

! CC RC  name changed and added 
    SPpol  = fdf_get('SpinPolarized',.false.)
    NonCol = fdf_get('NonCollinearSpin',.false.)
    SpOrb  = fdf_get('SpinOrbit',.false.)
    TRSym  = fdf_get('TimeReversalSymmetry',.false.)    

! CC RC  Added SpOrb and related tags
    write(6,'(a,l)') 'm_spin: SpOrb = ', SpOrb
    if (SpOrb) then
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
    elseif (NonCol) then
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
    elseif (SPpol) then
#ifdef TRANSIESTA
       call die("TranSiesta does not implement spin-orbit yet")
#endif
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
    endif

! CC RC  Added
    if (ionode) then
       write(6,'(a,4x,l1)') 'redata: Non-magnetic run                = ',NoMagn
       write(6,'(a,4x,l1)') 'redata: SpinPolarized run                = ',SPpol
       write(6,'(a,4x,l1)') 'redata: Non-Collinear-spin run           = ',NonCol
       write(6,'(a,4x,l1)') 'redata: Spin-Orbit run                   = ',SpOrb
       write(6,'(a,4x,l1)') 'redata: Time-Reversal Symmetry           = ',TRSym
       write(6,'(a,4x,i1)') 'redata: Number of spin components        = ',nspin
    end if

!    if (ionode) then
!       write(6,'(a,4x,l1)') 'redata: Non-Collinear-spin run           = ',NonCol
!       write(6,'(a,4x,l1)') 'redata: SpinPolarized (Up/Down) run      = ',SPpol
!       write(6,'(a,4x,i1)') 'redata: Number of spin components        = ',nspin
!    end if

    nullify(efs,qs)
    call re_alloc(efs,1,nspin,name="efs",routine="init_spin")
    call re_alloc(qs,1,nspin,name="qs",routine="init_spin")

  end subroutine init_spin

end module m_spin
