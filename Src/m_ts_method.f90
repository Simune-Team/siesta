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

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA as well 
  ! as the MUMPS library
  integer, parameter :: TS_SPARSITY_MUMPS = 3

  ! The default solution method (it will be reset
  ! after option reading)
  integer :: ts_method = TS_SPARSITY_TRI

  integer, parameter :: TYP_BUFFER = -1
  integer, parameter :: TYP_DEVICE = 0

  integer :: na_Buf = 0
  integer :: no_Buf = 0

  integer, private :: ts_no_u = 0
  integer, pointer, private :: ts_a_type(:) => null()
  integer, pointer, private :: ts_o_type(:) => null()
  integer, pointer, private :: ts_a_offset(:) => null()
  integer, pointer, private :: ts_o_offset(:) => null()

contains

  subroutine ts_init_regions(prefix,N_Elec, Elecs, na_u, lasto)

    use alloc
    use fdf
    use m_ts_electype

    character(len=*), intent(in) :: prefix
    integer,    intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer,    intent(in) :: na_u, lasto(0:na_u)

    integer :: i, ia, ia1, ia2, ia3

    ! prepare to read in the data...
    character(len=len_trim(prefix)+13) :: bName
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=50) :: g
    
    ! allocate regions
    if ( .not. associated(ts_a_type) ) then
       call re_alloc(ts_a_type,1,na_u)
       call re_alloc(ts_a_offset,1,na_u)
       call re_alloc(ts_o_type,1,lasto(na_u))
       call re_alloc(ts_o_offset,1,lasto(na_u))
       ts_a_type(:) = TYP_DEVICE
       ts_a_offset(:) = 0
       ts_o_type(:) = TYP_DEVICE
       ts_o_offset(:) = 0
       ts_no_u = lasto(na_u)
    end if
    
    ! we have already set the type of the buffers and device
    do ia = 1 , na_u
       do i = 1 , N_Elec
          if ( AtomInElec(Elecs(i),ia) ) then
             call set_type(i,ia,na_u,lasto)
             exit
          end if
       end do
    end do

    ! Read in TS.Atoms.Buffer

    bName = trim(prefix)//'.Atoms.Buffer'
    if ( fdf_block(bName,bfdf) ) then
    
       ! read by line and set them to be buffer atoms
    do while ( fdf_bline(bfdf,pline) ) 
       ! empty line
       if ( fdf_bnnames(pline) == 0 ) cycle
       
       g = fdf_bnames(pline,1)
       if ( leqi(g,'position') ) then
          ! we have a position
          if ( fdf_bnnames(pline) == 1 ) then
             do i = 1 , fdf_bnintegers(pline)
                ia = fdf_bintegers(pline,i)
                if (ia < 0) ia = na_u + ia + 1
                call set_type(TYP_BUFFER,ia,na_u,lasto)
             end do
          else if ( fdf_bnnames(pline) == 3 ) then
             g = fdf_bnames(pline,2)
             if ( .not. leqi(g,'from') ) then
                call die('Error in block '//bName//': &
                     &position from <int> to <int> is ill formatted')
             end if
             g = fdf_bnames(pline,3)
             if ( .not. leqi(g,'to') ) then
                call die('Error in block '//bName//': &
                     &position from <int> to <int> is ill formatted')
             end if
             if ( fdf_bnintegers(pline) < 2 ) then
                call die('Error in block '//bName//': &
                     &position from <int> to <int> is ill formatted')
             end if
             ia1 = fdf_bintegers(pline,1)
             if (ia1 < 0) ia1 = na_u + ia1 + 1
             ia2 = fdf_bintegers(pline,2)
             if (ia2 < 0) ia2 = na_u + ia2 + 1
             ia3 = 1
             if ( fdf_bnintegers(pline) > 2 ) then
                ia3 = abs(fdf_bnintegers(pline,3))
             end if
             do ia = ia1,ia2,ia3
                call set_type(TYP_BUFFER,ia,na_u,lasto)
             end do
          end if
          
       end if
    end do

    end if
    ! old options for buffer atoms

    call fdf_obsolete('TS.BufferAtomsLeft')
    call fdf_obsolete('TS.BufferAtomsRight')

    ! Update counting buffers
    na_Buf = 0
    no_Buf = 0
    do ia = 1 , na_u
       if ( atom_type(ia) /= TYP_BUFFER ) cycle
       na_Buf = na_Buf + 1
       no_Buf = no_Buf + lasto(ia) - lasto(ia-1)
    end do

  contains

    subroutine set_type(typ,ia,na_u,lasto)
      integer, intent(in) :: typ, ia, na_u,lasto(0:na_u)
      integer :: i, no
      if ( ts_a_type(ia) /= TYP_DEVICE ) then
         write(*,'(2(a,i0))') 'Trying to set atom ',ia,' to type: ',typ
         write(*,'(2(a,i0))') 'Atom ',ia,' is already: ',ts_a_type(ia)

         call die('Error in setup. Atoms are having two types, check for &
              &electrode and buffer atom overlap...')
      end if
      ts_a_type(ia) = typ
      ts_o_type(lasto(ia-1)+1:lasto(ia)) = typ
      if ( typ == TYP_BUFFER ) then
         do i = ia , na_u
            ts_a_offset(i) = ts_a_offset(i) + 1
         end do
         no = lasto(ia) - lasto(ia-1)
         do i = lasto(ia-1) + 1 , lasto(na_u)
            ts_o_offset(i) = ts_o_offset(i) + no
         end do
      end if
    end subroutine set_type

  end subroutine ts_init_regions
  
  elemental function orb_type(io) result(typ)
    use geom_helper, only : UCORB
    integer, intent(in) :: io
    integer :: typ
    typ = ts_o_type(ucorb(io,ts_no_u))
  end function orb_type

  elemental function orb_offset(io) result(off)
    use geom_helper, only : UCORB
    integer, intent(in) :: io
    integer :: off
    off = ts_o_offset(ucorb(io,ts_no_u))
  end function orb_offset

  elemental function ts2s_orb(io) result(off)
    integer, intent(in) :: io
    integer :: off
    do off = io , ts_no_u
       if ( ts_o_type(off) == TYP_BUFFER ) cycle ! the buffer atoms are NOT transiesta
       if ( off - ts_o_offset(off) == io ) return
    end do
  end function ts2s_orb

  elemental function atom_offset(ia) result(off)
    use geom_helper, only : UCORB
    integer, intent(in) :: ia
    integer :: off
    off = ts_a_offset(ia)
  end function atom_offset

  elemental function atom_type(ia) result(typ)
    integer, intent(in) :: ia
    integer :: typ
    typ = ts_a_type(ia)
  end function atom_type

  elemental function a_isBuffer(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = ts_a_type(ia) == TYP_BUFFER
  end function a_isBuffer

  elemental function a_isElec(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = ts_a_type(ia) > 0
  end function a_isElec

  elemental function a_isDev(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = ts_a_type(ia) == TYP_DEVICE
  end function a_isDev
    
end module m_ts_method
  
