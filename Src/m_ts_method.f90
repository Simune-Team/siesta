module m_ts_method

  use m_region
  
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

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA as well 
  ! as the MUMPS library
  integer, parameter :: TS_SPARSITY_MUMPS = 3

  ! The default solution method (it will be reset
  ! after option reading)
  integer :: ts_method = TS_SPARSITY_TRI

  ! The buffer atoms have type = -1, dev = 0, Electrodes = E_idx
  integer, parameter :: TYP_BUFFER = -1
  integer, parameter :: TYP_DEVICE = 0

  ! Containers for which atom/orbital is what
  integer, private :: no_u_TS = 0
  integer, allocatable, private :: a_type(:)
  integer, allocatable, private :: o_type(:)
  integer, allocatable, private :: a_offset(:)
  integer, allocatable, private :: o_offset(:)

  ! Containers for r_[ao]Buf%n for easier handling.
  integer :: na_Buf = 0
  integer :: no_Buf = 0
  type(tRgn) :: r_aBuf, r_oBuf ! the buffer region
  type(tRgn) :: r_aC,   r_oC   ! the entire calculation region

  ! We create a table for converting transiesta orbitals to 
  ! siesta orbitals:
  !   s_io = r_pvt%r(ts_io)
  ! This allows to easily convert from one to the other
  type(tRgn) :: r_pvt

contains

  subroutine ts_init_regions(prefix,N_Elec, Elecs, na_u, lasto)

    use alloc
    use fdf
    use fdf_extra, only : fdf_brange
    use m_ts_electype

    character(len=*), intent(in) :: prefix
    integer,    intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer,    intent(in) :: na_u, lasto(0:na_u)

    integer :: i, ia

    ! prepare to read in the data...
    character(len=len_trim(prefix)+13) :: bName
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=50) :: g
    type(tRgn) :: r_tmp

    no_u_TS = lasto(na_u)
    
    ! allocate regions
    if ( allocated(a_type) ) then

       deallocate(a_type,a_offset,o_type,o_offset)

       ! Clean regions
       call rgn_delete(r_aBuf,r_oBuf,r_aC,r_oC,r_pvt)

    end if

    allocate(a_type(na_u)   ,a_offset(na_u)    )
    allocate(o_type(no_u_TS),o_offset(no_u_TS) )
    a_type(:)   = TYP_DEVICE
    a_offset(:) = 0
    o_type(:)   = TYP_DEVICE
    o_offset(:) = 0
    
    ! we have already set the type of the buffers and device
    ! (as elecs might later on implement distributed
    !  positions we do it "stupidly")
    do i = 1 , N_Elec
       do ia = 1 , TotUsedAtoms(Elecs(i))
          call set_type(i,Elecs(i)%idx_a - 1 + ia,na_u,lasto)
       end do
    end do

    ! old options for buffer atoms
    call fdf_obsolete('TS.BufferAtomsLeft')
    call fdf_obsolete('TS.BufferAtomsRight')

    ! Read in TS.Atoms.Buffer
    bName = trim(prefix)//'.Atoms.Buffer'
    if ( fdf_block(bName,bfdf) ) then
    
       ! read by line and set them to be buffer atoms
       do while ( fdf_bline(bfdf,pline) ) 
          ! empty line
          if ( fdf_bnnames(pline) == 0 ) cycle
       
          g = fdf_bnames(pline,1)
          if ( leqi(g,'atom') ) then

             call fdf_brange(pline,r_tmp,1,na_u)
             if ( r_aBuf%n == 0 ) then
                call rgn_copy(r_tmp,r_aBuf)
             else if ( r_tmp%n > 0 ) then
                call rgn_union(r_aBuf,r_tmp,r_aBuf)
             end if
             
          end if
          
       end do

       if ( r_aBuf%n > 0 ) then
          do i = 1 , r_aBuf%n
             call set_type(TYP_BUFFER,r_aBuf%r(i),na_u,lasto)
          end do
       end if

    end if

    ! Create the calculation region
    call rgn_range(r_aC,1,na_u)
    call rgn_complement(r_aBuf,r_aC,r_tmp)
    call rgn_copy(r_tmp,r_aC)
    call rgn_delete(r_tmp)

    ! Sort the regions (faster look-ups)
    call rgn_sort(r_aBuf)
    call rgn_sort(r_aC)
    
    ! Convert atom regions to orbital regions
    call rgn_Atom2Orb(r_aBuf,na_u,lasto,r_oBuf)
    call rgn_Atom2Orb(r_aC,na_u,lasto,r_oC)
    ! Just tell them that they are sorted (they MUST be)
    call rgn_sort(r_oBuf)
    call rgn_sort(r_oC)

    ! Update counting buffers
    na_Buf = r_aBuf%n
    no_Buf = r_oBuf%n

    ! Create the "pivoting" array
    call rgn_range(r_pvt,1,lasto(na_u)-r_oBuf%n)
    ia = 0
    do i = 1 , lasto(na_u)
       if ( orb_type(i) == TYP_BUFFER ) cycle
       ia = ia + 1
       r_pvt%r(ia) = i
    end do

  contains

    subroutine set_type(typ,ia,na_u,lasto)
      integer, intent(in) :: typ, ia, na_u,lasto(0:na_u)
      integer :: i, no
      if ( a_type(ia) /= TYP_DEVICE ) then
         write(*,'(2(a,i0))') 'Trying to set atom ',ia,' to type: ',typ
         write(*,'(2(a,i0))') 'Atom ',ia,' is already: ',a_type(ia)

         call die('Error in setup. Atoms are having two types, check for &
              &electrode and buffer atom overlap...')
      end if
      a_type(ia) = typ
      o_type(lasto(ia-1)+1:lasto(ia)) = typ
      if ( typ == TYP_BUFFER ) then
         do i = ia , na_u
            a_offset(i) = a_offset(i) + 1
         end do
         no = lasto(ia) - lasto(ia-1)
         do i = lasto(ia-1) + 1 , lasto(na_u)
            o_offset(i) = o_offset(i) + no
         end do
      end if
    end subroutine set_type

  end subroutine ts_init_regions
  
  elemental function orb_type(io) result(typ)
    use geom_helper, only : UCORB
    integer, intent(in) :: io
    integer :: typ
    typ = o_type(ucorb(io,no_u_TS))
  end function orb_type

  elemental function orb_offset(io) result(off)
    use geom_helper, only : UCORB
    integer, intent(in) :: io
    integer :: off
    off = o_offset(ucorb(io,no_u_TS))
  end function orb_offset

  elemental function ts2s_orb(io) result(off)
    integer, intent(in) :: io
    integer :: off
    do off = io , no_u_TS
       if ( o_type(off) == TYP_BUFFER ) cycle ! the buffer atoms are NOT transiesta
       if ( off - o_offset(off) == io ) return
    end do
  end function ts2s_orb

  elemental function atom_offset(ia) result(off)
    use geom_helper, only : UCORB
    integer, intent(in) :: ia
    integer :: off
    off = a_offset(ia)
  end function atom_offset

  elemental function atom_type(ia) result(typ)
    integer, intent(in) :: ia
    integer :: typ
    typ = a_type(ia)
  end function atom_type

  elemental function a_isBuffer(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = a_type(ia) == TYP_BUFFER
  end function a_isBuffer

  elemental function a_isElec(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = a_type(ia) > 0
  end function a_isElec

  elemental function a_isDev(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = a_type(ia) == TYP_DEVICE
  end function a_isDev
    
end module m_ts_method
  
