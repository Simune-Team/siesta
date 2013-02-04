module m_pexsi
  public :: pexsi
CONTAINS
  subroutine pexsi(no_u, nspin,  &
       maxnh, numh, listhptr, listh, H, S)

    use precision, only: dp
    use m_hsx, only    : write_hs_formatted

    implicit          none

    integer           maxnh, no_u, nspin
    integer           listh(maxnh), numh(*), listhptr(*)
    real(dp)          H(maxnh,nspin), S(maxnh)

    external      :: timer

    call timer("pexsi", 1)

    call write_hs_formatted(no_u, nspin,  &
         maxnh, numh, listhptr, listh, H, S)

    call timer("pexsi", 2)

  end subroutine pexsi
end module m_pexsi
