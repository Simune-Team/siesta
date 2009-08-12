module m_spin
  use precision, only: dp
  private
  integer, save, public           :: nspin
  integer, save, public           :: MColl     ! used for matrix allocation
                                               ! if (NonCol) then MColl=2 
                                               !             else MColl=1
  logical, save, public           :: SPpol     = .false. 
                                        ! .true. if system is spin-polarized
  logical, save, public           :: NonCol    = .false. 
                                        ! .true. if system is non-collinear
  integer, save, public           :: NumSpin   = 1       
                                        ! formal number of spin components

  real(dp), pointer, save, public  :: efs(:)
  real(dp), pointer, save, public  :: qs(:)

  public :: init_spin

CONTAINS

  subroutine init_spin()
    use fdf, only : fdf_get
    use alloc, only: re_alloc
    use parallel, only: IOnode

    implicit none


    sppol  = fdf_get('SpinPolarized',.false.)
    noncol = fdf_get('NonCollinearSpin',.false.)

    if (NonCol) then

       nspin     = 4
       NumSpin   = 1
       MColl     = 2
       NonCol    = .true.
       SPpol     = .false.
    elseif (SPpol) then
       nspin     = 2
       NumSpin   = nspin
       MColl     = 1
       NonCol    = .false.
       SPpol     = .true.
    else 
       nspin     = 1
       NumSpin   = nspin
       MColl     = 1
       NonCol    = .false.
       SPpol     = .false.
    endif

    if (ionode) then
       write(6,'(a,4x,l1)') 'redata: SpinPolarized run                = ',SPpol
       write(6,'(a,4x,l1)') 'redata: Non-Collinear-spin run           = ',NonCol
       write(6,'(a,4x,i1)') 'redata: Number of spin components        = ',nspin
    end if

    nullify(efs,qs)
    call re_alloc(efs,1,nspin,name="efs",routine="init_spin")
    call re_alloc(qs,1,nspin,name="qs",routine="init_spin")

  end subroutine init_spin

end module m_spin
