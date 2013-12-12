! Auxilary functions for transiesta

module m_ts_aux
      
  use precision, only : dp

  implicit none

  interface nf
     module procedure nf_z
     module procedure nf_d
     module procedure nf2
  end interface nf

  private :: nf_z, nf_d

contains
  
  ! Fermi Function
  elemental function nf_d(d) result(nf)
    real(dp), intent(in) :: d
    real(dp) :: nf
    nf = 1._dp/(1._dp + exp(d))
  end function nf_d
  elemental function nf_z(z) result(nf)
    complex(dp), intent(in) :: z
    complex(dp) :: nf
    nf = dcmplx(1._dp,0._dp)/(dcmplx(1._dp,0._dp) + exp(z))
  end function nf_z

  ! Double fermi function
  ! calculates
  !   nF(E-E2) - nF(E-E1)
  elemental function nf2(E,E1,E2,kT) result(nf)
    real(dp), intent(in) :: E, E1,E2,kT
    real(dp) :: nf
    nf = nf_d((E-E2)/kT) - nf_d((E-E1)/kT)
  end function nf2

end module m_ts_aux
