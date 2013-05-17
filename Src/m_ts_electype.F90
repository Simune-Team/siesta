module m_ts_electype

  implicit none

  public :: Elec
  public :: HSfile
  public :: Atoms, UsedAtoms, TotUsedAtoms
  public :: Orbs, UsedOrbs, TotUsedOrbs
  public :: Rep
  public :: RepA1, RepA2, RepA3

  integer, parameter :: HSfile_len = 200

  type :: Elec
     sequence
     character(len=HSfile_len) :: HSfile
     integer :: Atoms, UsedAtoms
     integer :: Orbs, UsedOrbs
     integer :: RepA1 = 1, RepA2 = 1, RepA3 = 1
  end type Elec

contains
  
  elemental function HSfile(this)
    type(Elec), intent(in) :: this
    character(len=HSfile_len) :: HSfile
    HSfile = this%HSfile
  end function HSfile

  elemental function Atoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%Atoms
  end function Atoms

  elemental function UsedAtoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%UsedAtoms
  end function UsedAtoms
  elemental function TotUsedAtoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%UsedAtoms * Rep(this)
  end function TotUsedAtoms

  elemental function Rep(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = RepA1(this)*RepA2(this)*RepA3(this)
  end function Rep
  elemental function RepA1(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%RepA1
  end function RepA1
  elemental function RepA2(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%RepA2
  end function RepA2
  elemental function RepA3(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%RepA3
  end function RepA3

  elemental function Orbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%Orbs
  end function Orbs
  elemental function UsedOrbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%UsedOrbs
  end function UsedOrbs
  elemental function TotUsedOrbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%UsedOrbs * Rep(this)
  end function TotUsedOrbs

end module m_ts_electype
