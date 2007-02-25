module m_spin
  use precision, only: dp
  private
  integer, save, public            :: nspin
  real(dp), pointer, save, public  :: efs(:)
  real(dp), pointer, save, public  :: qs(:)

  public :: init_spin

CONTAINS

  subroutine init_spin()
    use m_fdf_global
    use alloc, only: re_alloc
    use parallel, only: IOnode

    implicit none

    logical  noncol, sppol

    call fdf_global_get(sppol,'SpinPolarized',.false.)
    call fdf_global_get(noncol,'NonCollinearSpin',.false.)

    if (noncol) then
       nspin = 4
    elseif (sppol) then
       nspin = 2
    else 
       nspin = 1
    endif

    if (IOnode) then
       write(6,'(a,4x,i1)')                                 &
            'redata: Number of spin components        = ',nspin
    endif

    nullify(efs,qs)
    call re_alloc(efs,1,nspin,name="efs",routine="init_spin")
    call re_alloc(qs,1,nspin,name="qs",routine="init_spin")

  end subroutine init_spin

end module m_spin
