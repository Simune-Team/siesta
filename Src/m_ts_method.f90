module m_ts_method
  
  implicit none

  public
  save
  
  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA.
  integer, parameter :: TS_SPARSITY = 1

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA as well 
  ! as the heavily optimized tri-diagonalization
  integer, parameter :: TS_SPARSITY_TRI = 2

  ! The default solution method (it will be reset
  ! after option reading)
  integer :: ts_method = TS_SPARSITY_TRI

  integer, pointer :: ts_a_type(:) => null()
  integer, pointer :: ts_o_type(:) => null()

  integer, parameter :: TYP_BUFFER = -1
  integer, parameter :: TYP_DEVICE = 0

contains

  subroutine ts_init_region_types(na_BufL,Elecs, &
       na_BufR, &
       na_u, lasto)

    use alloc
    use m_ts_electype

    type(Elec), intent(in) :: Elecs(:)
    integer,    intent(in) :: na_BufL, na_BufR
    integer,    intent(in) :: na_u, lasto(0:na_u)

    integer :: i, ia
    
    ! allocate regions
    call re_alloc(ts_a_type,1,na_u)
    call re_alloc(ts_o_type,1,lasto(na_u))

    ts_a_type(:) = TYP_DEVICE
    if ( na_BufL > 0 ) ts_a_type(1:na_BufL)       = TYP_BUFFER
    if ( na_BufR > 0 ) ts_a_type(na_u-na_BufR+1:) = TYP_BUFFER

    ! we have already set the type of the buffers and device
    do ia = 1 , na_u
       do i = 1 , size(Elecs)
          if ( Elecs(i)%idx_na <= ia .and. &
               ia < Elecs(i)%idx_na + TotUsedAtoms(Elecs(i)) ) then
             ts_a_type(ia) = i
             exit
          end if
       end do
       ts_o_type(lasto(ia-1)+1:lasto(ia)) = ts_a_type(ia)
    end do

  end subroutine ts_init_region_types

  elemental function get_orb_type(io) result(typ)
    use geom_helper, only : UCORB
    integer, intent(in) :: io
    integer :: typ
    typ = ts_o_type(ucorb(io,size(ts_o_type)))
  end function get_orb_type

  elemental function get_atom_type(ia) result(typ)
    integer, intent(in) :: ia
    integer :: typ
    typ = ts_a_type(ia)
  end function get_atom_type
    
end module m_ts_method
  
