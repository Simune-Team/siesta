! 
! Reprogrammed tbtrans
! This code has been fully created by;
! Nick Papior Andersen, nickpapior @ gmail.com
! The code has been constructed in 2014 and is
! meant to superseede any previous tbtrans versions.

! One major difference in tbtrans is that we never deal with
! xij arrays.
! We ONLY deal with unitcell offsets
! Hence just after reading in TSHS we convert it
! to the correct format (list_col contains the supercell
! index AND the column in the unitcell)
program tbtrans

  use class_Sparsity
  use class_dSpData1D
  use class_dSpData2D

  use m_tbt_options, only : kT
  use m_tbt_regions
  use m_tbt_hs, only : TSHS

  use m_tbtrans

  implicit none

  ! Initialize everything.
  call tbt_init()

  ! Call tbtrans
  call tbt(TSHS, kT)

  call tbt_end()

end program tbtrans

