! Module for easy expansion/copying a DM file from another 
! geometry to the current one...

! This code has been fully created by:
!    Nick Papior Andersen, nickpapior@gmail.com
module m_handle_sparse

  use precision, only : dp
  use parallel, only : Node, Nodes
  use geom_helper, only : ucorb, iaorb

  implicit none

  private

  public :: bulk_expand
  public :: expand_spd2spd_2D
  public :: unfold_noauxiliary_supercell_SpD
  public :: copy_supercell_SpD
  public :: correct_supercell_SpD
  public :: reduce_spin_size

  interface copy_supercell_SpD
    module procedure copy_supercell_Sp1D
    module procedure copy_supercell_Sp2D
  end interface copy_supercell_SpD

  interface unfold_noauxiliary_supercell_SpD
    module procedure unfold_noauxiliary_supercell_Sp2D
  end interface unfold_noauxiliary_supercell_SpD

  interface correct_supercell_SpD
    module procedure correct_supercell_Sp1D
    module procedure correct_supercell_Sp2D
  end interface correct_supercell_SpD

contains

  subroutine bulk_expand(na_u,xa,lasto,cell,nsc,isc_off,DM_2D)

    use fdf
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use class_dSpData1D
    use m_os, only: file_exist
    use units, only: Ang
    use m_iodm, only: read_DM
    use m_ts_io, only: ts_read_TSHS

    ! The input parameters that govern the simulation
    integer, intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u), cell(3,3)
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))
    ! The two arrays that will be needed
    type(dSpData2D), intent(inout) :: DM_2D

    ! dummy arrays
    type(OrbitalDistribution) :: fake_dit
    type(Sparsity)  :: fsp
    type(Sparsity), pointer :: psp
    type(dSpData1D) :: tmp_1D
    type(dSpData2D) :: fDM_2D
    real(dp) :: fcell(3,3), fkdispl(3), fEf, fQtot, fTemp
    real(dp), pointer :: fxa(:,:) => null()
    integer :: fnsc(3), fna_u, fno_u, fnspin, fkcell(3,3), at, fn_s
    integer :: Tile(3), Reps(3), fnsc_DM(3)
    real(dp) :: def_xa_EPS, xa_EPS

    integer, pointer :: flasto(:) => null(), fisc_off(:,:) => null()

    integer :: allowed(na_u), itmp, s_na, c_na
    character(len=400) :: HSfile, DMfile, ln
    logical :: d_log1, d_log2, d_log3
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()

    ! If the block does not exist, we might as well return
    if ( .not. fdf_block('DM.Init.Bulk',bfdf) ) return

    if ( Node == 0 ) then
      write(*,'(/,a)') 'siesta: Initializing DM from bulk.'
      if ( product(nsc) == 1 ) then
        write(*,'(/,a)') 'siesta: *** WARNING *** Non-supercell calculation, &
            &will not be able to correctly handle cross-boundary connections.'
      end if
    end if

    ! Using this method we allow all interactions to be
    ! initialized (default)
    do at = 1 , na_u
      allowed(at) = at
    end do

    ! Default coordinate acceptance, 0.005 Ang
    def_xa_EPS = fdf_get('DM.Init.Bulk.Coord.Eps', 0.001_dp * Ang, 'Bohr')

    ! Read in each segment and copy data!
    do while ( fdf_bline(bfdf,pline) )
      if ( fdf_bnnames(pline) == 0 ) cycle ! skip empty lines

      ! Get the name of the segment that we will copy
      ln = ' '
      ln = fdf_bnames(pline,1)

      xa_EPS = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Coord.Eps',def_xa_EPS, 'Bohr')

      ! we search for the fdf-flags
      HSfile = ' '
      HSfile = fdf_get('DM.Init.Bulk.'//trim(ln),'NONE')
      ! Now we have all required information
      if ( .not. file_exist(HSfile, Bcast = .true. ) ) then
        write(*,*) trim(HSfile)
        call die('You at least need to supply the TSHS file for &
            &bulk segment '//trim(ln)//'.')
      end if

      ! Get the tiling/repetition
      if ( fdf_islist('DM.Init.Bulk.'//trim(ln)//'.Tile') ) then
        call fdf_list('DM.Init.Bulk.'//trim(ln)//'.Tile', 3, Tile)
      else
        Tile(:) = 1
      end if
      Tile(1) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Tile.A1',Tile(1))
      Tile(2) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Tile.A2',Tile(2))
      Tile(3) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Tile.A3',Tile(3))
      if ( fdf_islist('DM.Init.Bulk.'//trim(ln)//'.Repeat') ) then
        call fdf_list('DM.Init.Bulk.'//trim(ln)//'.Repeat', 3, Reps)
      else
        Reps(:) = 1
      end if
      Reps(1) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Repeat.A1',Reps(1))
      Reps(2) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Repeat.A2',Reps(2))
      Reps(3) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Repeat.A3',Reps(3))

      ! We first try and guess the DM file
      at = len_trim(HSfile)
      DMfile = fdf_get('DM.Init.Bulk.'//trim(ln)//'.DM', HSfile(1:at-4)//'DM')
      if ( .not. file_exist(DMfile, Bcast = .true.) ) then
        DMfile = fdf_get('DM.Init.Bulk.'//trim(ln)//'.DM', HSfile(1:at-4)//'TSDE')
      end if
      if ( .not. file_exist(DMfile, Bcast = .true. ) ) then
        call die('DM file could not be found, have you supplied an &
            &erroneous path?')
      end if

      ! Get the starting atom we need to copy into!
      at = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Atom.Insert',0)
      if ( at < 0 ) then
        at = na_u + at + 1
      end if
      if ( at <= 0 .or. na_u < at ) then
        print *,'Requested atom:',at
        call die('You need to supply the starting atom for the &
            &copy operation!')
      end if

      ! We have gathered all needed information!

      ! We require a TSHS file as that file contains
      ! all necessary information.
      call ts_read_TSHS(HSfile, &
          d_log1, d_log2, d_log3, &
          fcell, fnsc, fna_u, fno_u, fnspin,  &
          fkcell, fkdispl, &
          fxa, flasto, &
          fsp, fDM_2D, tmp_1D, fisc_off, &
          fEf, fQtot, fTemp, &
          itmp, itmp, &
          Bcast= .true. )
      fn_s = product(fnsc)
      ! 
      ! Clean-up, we do not need the Hamilton and overlap
      call delete(fDM_2D)
      call delete(tmp_1D)

      ! Get the options of how many atoms we actually
      ! copy
      s_na = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Atom.Start',1)
      if ( s_na < 0 ) then
        s_na = fna_u + s_na + 1
      end if
      c_na = fdf_get('DM.Init.Bulk.'//trim(ln)//'.Atom.Count',fna_u)
      if ( s_na <= 0 .or. s_na + c_na - 1 > fna_u ) then
        call die('You are requesting to copy more atoms than present &
            &in the file.')
      end if

      ! Luckily, the TSDE and the DM file
      ! exactly the same, except that TSDE has an extra EDM and Ef
      ! at the end, we do not care about that! :)
      ! read in DM file
      call read_DM( DMfile, fake_dit, fnsc_DM, fDM_2D, d_log1, Bcast=.true.)
      if ( fnsc_DM(1) == 0 ) fnsc_DM = fnsc
      if ( size(fDM_2D, 2) /= fnspin ) then
        call die('bulk_expand: DM and TSHS does not have the same spin')
      end if
      if ( nrows_g(fDM_2D) /= fno_u ) then
        call die('bulk_expand: DM and TSHS does not have the same no_u')
      end if
      psp => spar(fDM_2D)
      psp = fsp
      call delete(fsp)
      if ( .not. d_log1 ) then
        call die('Something went wrong, file not found?')
      end if

      if ( Node == 0 ) then
        ! Write out the settings
        itmp = product(Tile) * product(Reps) * c_na - 1
        write(*,'(a,i0,'' ,'',i0,2a)') &
            'siesta: Initializing bulk DM for atoms [ ',at,at+itmp,']', &
            ' using segment: '//trim(ln)
      end if

      call expand_spd2spd_2D(s_na,c_na,fna_u,flasto,fxa, &
          fDM_2D,&
          fcell, Tile, Reps, fn_s, fisc_off, &
          na_u,xa,lasto,DM_2D,cell,product(nsc),isc_off, at, xa_EPS, &
          print = .true., allowed_a = allowed)

      ! De-allocate before reading the next thing...
      call delete(fDM_2D)
      deallocate(fxa) ; nullify(fxa)
      deallocate(flasto) ; nullify(flasto)
      deallocate(fisc_off) ; nullify(fisc_off)

    end do
    if ( Node == 0 ) write(*,*) ! new-line

  end subroutine bulk_expand

  ! This will expand an inner sparse matrix to an outer sparse matrix
  ! by copying elements.
  subroutine expand_spd2spd_2D(s_a,na,na_i,lasto_i,xa_i,in,&
      cell_i, tile_i, repeat_i, n_s_i, sc_off_i, &
      na_o,xa_o,lasto_o,out,cell_o,n_s_o,sc_off_o, at, xa_EPS, &
      print, allowed_a)

    use units, only: Ang
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
#ifdef MPI
    use mpi_siesta
#endif
    use m_region

    ! Starting atom and number of atoms from the starting atom.
    integer, intent(in) :: s_a, na
    integer, intent(in) :: na_i, lasto_i(0:na_i), tile_i(3), repeat_i(3)
    real(dp), intent(in) :: xa_i(3,na_i), cell_i(3,3)
    integer, intent(in) :: n_s_i, sc_off_i(3,0:n_s_i-1)
    ! The density matrices that describes two
    ! different parts.
    type(dSpData2D), intent(inout) :: in, out
    integer, intent(in) :: na_o, lasto_o(0:na_o)
    real(dp), intent(in) :: xa_o(3,na_o), cell_o(3,3)
    integer, intent(in) :: n_s_o, sc_off_o(3,0:n_s_o-1)
    ! The atom that the in-sparsity pattern will start at
    integer, intent(in) :: at
    real(dp), intent(in) :: xa_EPS
    ! Whether we should print-out the information for the user
    logical, intent(in), optional :: print
    ! The updated elements can be chosen per atomic placement
    ! A list of allowed atoms can be supplied
    integer, intent(in), optional :: allowed_a(:)

    ! We are ready to check and copy the sparsity pattern...
    type(Sparsity), pointer :: sp_i, sp_o
    type(OrbitalDistribution), pointer :: dit_i, dit_o

    ! arrays for the sparsity patterns
    integer, pointer :: o_ptr(:), o_ncol(:), o_col(:)
    integer, pointer :: i_ptr(:), i_ncol(:), i_col(:)
    ! loop variables for the sparsity patterns
    integer :: ia_i, io_i, i_i, iat, io_o, lio_o, i_o
    integer :: i_s, i, ao, no_o, no_i, at_end, lat
    integer :: i1, i2, i3, o1, o2, o3
    integer :: orb_i, orb_o
    integer :: copy(3)

    real(dp) :: xc_i(3), xc_o(3), xj_i(3), xj_o(3)
    real(dp), pointer :: a_i(:,:), a_o(:,:)

    type(tRgn) :: rallow

#ifdef MPI
    integer :: tmp_copy(3)
    integer :: MPIerror
#endif

    o1 = product(tile_i)
    i1 = product(repeat_i)
    i = at - 1 + o1 * i1 * na
    if ( i > na_o ) then
      call die('Requested expansion region too large.')
    end if

    if ( present(allowed_a) ) then
      ! We copy over the allowed atoms
      call rgn_init(rallow, size(allowed_a))
      do i = 1 , rallow%n
        rallow%r(i) = allowed_a(i)
      end do
      call rgn_sort(rallow)
    else
      call rgn_range(rallow, at, at + o1 * i1 * na - 1)
    end if

    dit_i => dist(in)
    sp_i  => spar(in)
    a_i   => val(in)
    dit_o => dist(out)
    sp_o  => spar(out)
    a_o   => val(out)

    call attach(sp_i,n_col=i_ncol, list_ptr=i_ptr, list_col=i_col, &
        nrows_g=no_i)
    call attach(sp_o,n_col=o_ncol, list_ptr=o_ptr, list_col=o_col, &
        nrows_g=no_o)

    lat    = at
    at_end = at - 1 

    ! Loop on all equivalent atoms
    iat = at - 1
    ! We count number of copied data
    copy(:) = 0

    o3_loop: do o3 = 0 , tile_i(3) - 1
      o2_loop: do o2 = 0 , tile_i(2) - 1
        o1_loop: do o1 = 0 , tile_i(1) - 1

          at_end = at_end + product(repeat_i) * na

          ! We loop over the input SP which we will copy
          do ia_i = s_a , s_a + na - 1

            i3_loop: do i3 = 0 , repeat_i(3) - 1
              i2_loop: do i2 = 0 , repeat_i(2) - 1
                i1_loop: do i1 = 0 , repeat_i(1) - 1

                  ! Step atom index in out-put array
                  iat = iat + 1

                  ! The expanded atomic position
                  if ( na == 1 ) then

                    ! As we only compare one atom, we force the alignment.
                    ! Hence we only compare orbital distances.
                    ! Hence, we need not expand the cell!

                  else

                    ! Calculate the actual position of this expanded atom
                    xc_i(:) = xa_i(:,ia_i) - xa_i(:,s_a) + &
                        cell_i(:,1)*i1+cell_i(:,2)*i2+cell_i(:,3)*i3

                    xc_o(:) = xa_o(:,iat) - xa_o(:,lat)

                    if ( maxval(abs(xc_o - xc_i)) > xa_EPS ) then
                      write(*,'(2(tr1,a,tr1,i0),3(tr1,f10.5),3(tr1,i2))') &
                          'ia',ia_i,'matching',iat, (xc_o-xc_i)/Ang, i1, i2, i3
                      call die('Atomic coordinates do not coincide, &
                          &have you employed correct ordering for expansion A1->A2->A3?')
                    end if

                  end if

                  ! Check that the number of orbitals is the same
                  if ( lasto_i(ia_i) - lasto_i(ia_i-1) /= &
                      lasto_o(iat) - lasto_o(iat-1) ) then
                    write(*,*) 'ia',ia_i,' has',lasto_i(ia_i) - lasto_i(ia_i-1), &
                        'ca',iat,' has',lasto_o(iat) - lasto_o(iat-1), &
                        ' orbitals'
                    call die('Not the same atom! Please correct')
                  end if

                  ! loop over the orbitals of this atom
                  do io_i = lasto_i(ia_i-1) + 1 , lasto_i(ia_i)

                    ! The equivalent orbital in the out-put sparsity
                    ! pattern is this:
                    io_o = lasto_o(iat-1) + io_i - lasto_i(ia_i-1)
                    lio_o = index_global_to_local(dit_o,io_o,Node)
                    ! If the orbital does not exist on the current node
                    ! we simply skip it...
                    if ( lio_o <= 0 ) cycle
                    !print '(a,i0,2(a,i4))','Node: ',Node, &
                    !     ' orbital: ',io_o,' local: ',lio_o

                    ! Now we can actually do something!

                    ! Capture number of possible orbitals that we have
                    copy(3) = copy(3) + o_ncol(lio_o)

                    ! Loop over indices in this row
                    ! the large sparsity pattern must be the largest
                    ! hence we have this as the outer loop
                    do i_o = o_ptr(lio_o) + 1 , o_ptr(lio_o) + o_ncol(lio_o)

                      ! First we figure out which atomic position this
                      ! corresponds to:
                      ao = iaorb(o_col(i_o),lasto_o)
                      ! Do not allow overwriting DM outside of region.
                      if ( .not. in_rgn(rallow, ao) ) cycle

                      i_s   = (o_col(i_o)-1) / no_o
                      orb_o = ucorb(o_col(i_o),no_o) - lasto_o(ao-1)

                      ! We always compare distances between 
                      ! two orbitals, if they match, then we
                      ! must be having the same orbital connection
                      xj_o(:) = xa_o(:,ao) - xa_o(:,iat) + &
                          cell_o(:,1) * sc_off_o(1,i_s) + &
                          cell_o(:,2) * sc_off_o(2,i_s) + &
                          cell_o(:,3) * sc_off_o(3,i_s) 

                      ! Now we need to figure out all the orbitals
                      ! that has the same meaning in both sparsity patterns
                      ! Get the cell-offset

                      do i_i = i_ptr(io_i) + 1 , i_ptr(io_i) + i_ncol(io_i)

                        i     = iaorb(i_col(i_i),lasto_i)
                        orb_i = ucorb(i_col(i_i),no_i) - lasto_i(i-1)

                        ! Different orbitals can have the same
                        ! orbital center, so we check the orbital
                        ! index to be the same as well...
                        if ( orb_i /= orb_o ) cycle

                        i_s   = (i_col(i_i)-1) / no_i

                        xj_i(:) = xa_i(:,i) - xa_i(:,ia_i) + &
                            cell_i(:,1) * sc_off_i(1,i_s) + &
                            cell_i(:,2) * sc_off_i(2,i_s) + &
                            cell_i(:,3) * sc_off_i(3,i_s)

                        ! If they are not equivalent we will not do anything
                        if ( maxval(abs(xj_o - xj_i)) > xa_EPS ) cycle

                        ! WUHUU, we have the equivalent atom and equivalent
                        ! orbital connection. We copy data now!

                        if ( lat <= ao .and. ao <= at_end ) then
                          copy(1) = copy(1) + 1 ! diagonal contribution
                        else
                          copy(2) = copy(2) + 1 ! off-diagonal contribution
                        end if

                        a_o(i_o,:) = a_i(i_i,:)

                      end do
                    end do
                  end do

                end do i1_loop
              end do i2_loop
            end do i3_loop

          end do

          lat = lat + product(repeat_i) * na

        end do o1_loop
      end do o2_loop
    end do o3_loop

    call rgn_delete(rallow)

    if ( .not. present(print) ) return
    if ( .not. print ) return

#ifdef MPI
    tmp_copy = copy
    call MPI_Reduce(tmp_copy,copy,3,MPI_Integer, MPI_Sum, &
        0, MPI_Comm_World, MPIerror)
#endif

    if ( Node == 0 ) then
      write(*,'(a,''[ '',i0,'', '',i0,''] out of '',i0,a)') &
          'Expanded in total [UC,SC] ',copy,' possible elements.'
    end if

  end subroutine expand_spd2spd_2D


  ! Copy one super-cell sparse pattern to another super-cell sparse
  ! pattern. This will copy [in] to [out].
  ! The data D2_in contains one sparsity pattern with nsc_in supercells.
  ! See sparse_matrices.F90 for details regarding the sparsity pattern layout.
  ! The data D2_out contains the "new" sparsity pattern we want to retain.
  ! I.e. if the number of supercells change, we only retain the equivalent supercell
  ! elements. And if the number of non-zero elements change we retain only the
  ! new ones (discarding the older ones).
  subroutine copy_supercell_Sp2D(d2_out, nsc_out, d2_in, nsc_in)
    
    use class_Sparsity
    use class_dSpData2D
#ifdef MPI
    use mpi_siesta
#endif

    ! output D2 sparse pattern
    type(dSpData2D), intent(inout) :: d2_out
    ! output D2 number of supercells
    integer, intent(in) :: nsc_out(3)
    ! input D2 sparse pattern
    type(dSpData2D), intent(inout) :: d2_in
    ! input D2 number of supercells
    integer, intent(in) :: nsc_in(3)

    ! We are ready to check and copy the sparsity pattern...
    type(Sparsity), pointer :: sp_i, sp_o

    ! Local variables
    integer :: no_u, no_l
    integer :: io, i, i_ind, o_ind, isc(3)
    integer :: i_is, o_is, o_hsc(3), new_col
    integer :: discarded(2), dim_min

    ! arrays for the sparsity patterns
    integer, pointer :: o_ptr(:), o_ncol(:), o_col(:)
    integer, pointer :: i_ptr(:), i_ncol(:), i_col(:)
    real(dp), pointer :: a_o(:,:), a_i(:,:)
    integer, allocatable :: o_isc(:,:,:), i_isc(:,:)
    integer, allocatable :: out_index(:)

    if ( all(nsc_out == nsc_in) ) return

    sp_i => spar(d2_in)
    a_i => val(d2_in)
    sp_o => spar(d2_out)
    a_o => val(d2_out)
    
    call attach(sp_i,n_col=i_ncol, list_ptr=i_ptr, list_col=i_col, &
        nrows=no_l, nrows_g=no_u)
    call attach(sp_o,n_col=o_ncol, list_ptr=o_ptr, list_col=o_col, &
        nrows=io, nrows_g=i)

    dim_min = min(size(a_i,2), size(a_o,2))
    if ( no_u /= i ) &
        call die('copy_supercell_sp_d2: error in number of global orbitals.')
    if ( no_l /= io ) &
        call die('copy_supercell_sp_d2: error in number of local orbitals.')

    ! Now create the conversion tables
    call generate_isc(nsc_out, o_isc)
    call generate_linear_isc(nsc_in, i_isc)

    ! Allocate look-up table
    allocate(out_index(product(nsc_out)*no_u))
    ! Set all elements to zero
    out_index(:) = 0

    ! Count the number of discarded non-zero elements
    !   (1) is the supercell discarded (for missing supercells)
    !   (2) is because the orbital interaction does not exist
    discarded = 0

    ! We need to check whether in-put SC is too large
    o_hsc = nsc_out / 2

    ! Since we are going to copy, then we have to set it to zero...
    a_o(:,:) = 0._dp

    do io = 1, no_l
      
      ! copy output lookup table to not search every element
      do i = 1, o_ncol(io)
        out_index(o_col(o_ptr(io)+i)) = o_ptr(io) + i
      end do

      ! Now we can do the copy...
      inner_columns: do i = 1, i_ncol(io)
        i_ind = i_ptr(io) + i
        i_is = (i_col(i_ind)-1) / no_u

        ! Get isc
        isc = i_isc(:, i_is)

        ! Check that the out supercell exists
        if ( any(abs(isc) > o_hsc) ) then
          discarded(1) = discarded(1) + 1
          cycle inner_columns
        end if

        ! We know the supercell exists, lets see if the orbital connection
        ! exists.
        o_is = o_isc(isc(1),isc(2),isc(3))
        
        ! Transfer the orbital index to the correct supercell
        o_ind = out_index(o_is * no_u + ucorb(i_col(i_ind), no_u))
        
        if ( o_ind == 0 ) then
          ! The orbital interaction does not exist
          discarded(2) = discarded(2) + 1
        else
          a_o(o_ind,1:dim_min) = a_i(i_ind,1:dim_min)
        end if
        
      end do inner_columns

      ! restore lookup table
      do i = 1, o_ncol(io)
        out_index(o_col(o_ptr(io)+i)) = 0
      end do
      
    end do

    deallocate(o_isc, i_isc, out_index)

  end subroutine copy_supercell_Sp2D

  subroutine copy_supercell_Sp1D(d1_out, nsc_out, d1_in, nsc_in)
    
    use class_Sparsity
    use class_dSpData1D
#ifdef MPI
    use mpi_siesta
#endif

    ! output D1 sparse pattern
    type(dSpData1D), intent(inout) :: d1_out
    ! output D1 number of supercells
    integer, intent(in) :: nsc_out(3)
    ! input D1 sparse pattern
    type(dSpData1D), intent(inout) :: d1_in
    ! input D1 number of supercells
    integer, intent(in) :: nsc_in(3)

    ! We are ready to check and copy the sparsity pattern...
    type(Sparsity), pointer :: sp_i, sp_o

    ! Local variables
    integer :: no_u, no_l
    integer :: io, i, i_ind, o_ind, isc(3)
    integer :: i_is, o_is, o_hsc(3), new_col
    integer :: discarded(2), dim_min

    ! arrays for the sparsity patterns
    integer, pointer :: o_ptr(:), o_ncol(:), o_col(:)
    integer, pointer :: i_ptr(:), i_ncol(:), i_col(:)
    real(dp), pointer :: a_o(:), a_i(:)
    integer, allocatable :: o_isc(:,:,:), i_isc(:,:)
    integer, allocatable :: out_index(:)

    if ( all(nsc_out == nsc_in) ) return

    sp_i => spar(d1_in)
    a_i => val(d1_in)
    sp_o => spar(d1_out)
    a_o => val(d1_out)
    
    call attach(sp_i,n_col=i_ncol, list_ptr=i_ptr, list_col=i_col, &
        nrows=no_l, nrows_g=no_u)
    call attach(sp_o,n_col=o_ncol, list_ptr=o_ptr, list_col=o_col, &
        nrows=io, nrows_g=i)

    if ( no_u /= i ) &
        call die('copy_supercell_sp_d2: error in number of global orbitals.')
    if ( no_l /= io ) &
        call die('copy_supercell_sp_d2: error in number of local orbitals.')

    ! Now create the conversion tables
    call generate_isc(nsc_out, o_isc)
    call generate_linear_isc(nsc_in, i_isc)

    ! Allocate look-up table
    allocate(out_index(product(nsc_out)*no_u))
    ! Set all elements to zero
    out_index(:) = 0

    ! Count the number of discarded non-zero elements
    !   (1) is the supercell discarded (for missing supercells)
    !   (2) is because the orbital interaction does not exist
    discarded = 0

    ! We need to check whether in-put SC is too large
    o_hsc = nsc_out / 2

    ! Since we are going to copy, then we have to set it to zero...
    a_o(:) = 0._dp

    do io = 1, no_l
      
      ! copy output lookup table to not search every element
      do i = 1, o_ncol(io)
        out_index(o_col(o_ptr(io)+i)) = o_ptr(io) + i
      end do

      ! Now we can do the copy...
      inner_columns: do i = 1, i_ncol(io)
        i_ind = i_ptr(io) + i
        i_is = (i_col(i_ind)-1) / no_u

        ! Get isc
        isc = i_isc(:, i_is)

        ! Check that the out supercell exists
        if ( any(abs(isc) > o_hsc) ) then
          discarded(1) = discarded(1) + 1
          cycle inner_columns
        end if

        ! We know the supercell exists, lets see if the orbital connection
        ! exists.
        o_is = o_isc(isc(1),isc(2),isc(3))
        
        ! Transfer the orbital index to the correct supercell
        o_ind = out_index(o_is * no_u + ucorb(i_col(i_ind), no_u))
        
        if ( o_ind == 0 ) then
          ! The orbital interaction does not exist
          discarded(2) = discarded(2) + 1
        else
          a_o(o_ind) = a_i(i_ind)
        end if
        
      end do inner_columns

      ! restore lookup table
      do i = 1, o_ncol(io)
        out_index(o_col(o_ptr(io)+i)) = 0
      end do
      
    end do

    deallocate(o_isc, i_isc, out_index)

  end subroutine copy_supercell_Sp1D

  !> Expand an nsc == 1 DM to an nsc_ == * supercell.
  !>
  !> In cases where the DM is constructed from all(nsc == 1) we know that
  !> all elements in the supercell has the same elements as in folded DM.
  !>
  !> However, IFF one has a DM with folded elements (S(io,io) > 1.)
  !> we have a problem because the DM elements are made of sums of
  !> supercell and primary unit-cell contributians.
  subroutine unfold_noauxiliary_supercell_Sp2D(sp_sc, D2)

    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData2D
    use class_dData2D
#ifdef MPI
    use mpi_siesta
#endif

    !> Input supercell sparsity pattern (this will contain periodic connections)
    type(Sparsity), intent(inout) :: sp_sc
    !> Input/Output 2D data with associated non-supercell sparse pattern, upon exit
    !> this contains as many non-zero elements as in sp_sc with expanded elements.
    type(dSpData2D), intent(inout) :: D2

    ! Local variables
    type(OrbitalDistribution) :: dit
    type(Sparsity), pointer :: sp
    type(dData2D) :: A2D

    real(dp), pointer :: A2(:,:), A2_sc(:,:)
    integer :: no_u, no_l
    integer :: io, jo, ind, scind

    ! arrays for the sparsity patterns
    integer, pointer :: ptr(:), ncol(:), col(:)
    integer, pointer :: scptr(:), scncol(:), sccol(:)

    ! This will silently assume same sizes in the sparse patterns

    sp => spar(D2)
    A2 => val(D2)
    call attach(sp,n_col=ncol, list_ptr=ptr, list_col=col, &
        nrows=no_l, nrows_g=no_u)

    ! Supercell sparsity pattern
    call attach(sp_sc, n_col=scncol, list_ptr=scptr, list_col=sccol, &
        nnzs=io)

    ! Create array that hosts the new data
    ind = size(A2, dim=2)
    call newdData2D(A2D, io, ind,"(unfold 2D)")
    A2_sc => val(A2D)
    A2_sc(:,:) = 0._dp

    do io = 1, no_l

      ! Loop supercell sparse pattern
      do scind = scptr(io) + 1, scptr(io) + scncol(io)

        jo = ucorb(sccol(scind), no_u)
        
        do ind = ptr(io) + 1, ptr(io) + ncol(io)
          if ( col(ind) == jo ) then
            A2_sc(scind, :) = A2(ind, :)
            exit
          end if
        end do
        
      end do
    end do

    dit = dist(D2)
    call newdSpData2D(sp_sc, A2D, dit, D2, name="Unfolded Sp2D")

    call delete(dit)
    call delete(A2D)

  end subroutine unfold_noauxiliary_supercell_Sp2D

  ! Correct a sparse pattern (in-place) from an old NSC to a new NSC
  ! A supercell sparse matrix layout is described in sparse_matrices.F90
  ! Each column value also holds the supercell index by an offset equal to
  ! the supercell index X no_u such that all columns and supercells are
  ! unique.
  ! Here we take one sparse data pattern and change all supercell indices
  ! from the nsc_old indices to nsc_new.
  ! If nsc_old == nsc_new, nothing will happen.
  subroutine correct_supercell_Sp2D(nsc_old, D2, nsc_new)

    use class_Sparsity
    use class_dSpData2D
#ifdef MPI
    use mpi_siesta
#endif

    type(dSpData2D), intent(inout) :: D2
    integer, intent(in) :: nsc_old(3)
    integer, intent(in) :: nsc_new(3)

    ! Local variables
    type(Sparsity), pointer :: sp
    integer :: no_u, no_l
    integer :: io, ind, isc(3)
    integer :: old_n_s, new_outside
    integer :: old_is, new_is, new_hsc(3)

    ! arrays for the sparsity patterns
    integer, pointer :: ptr(:), ncol(:), col(:)
    integer, allocatable :: old_isc(:,:), new_isc(:,:,:)
    real(dp), pointer :: a2(:,:)

    if ( all(nsc_old == nsc_new) ) return

    sp => spar(D2)
    call attach(sp,n_col=ncol, list_ptr=ptr, list_col=col, &
        nrows=no_l, nrows_g=no_u)

    ! Required to set removed elements to 0
    A2 => val(D2)

    ! Now create the conversion tables
    !> The linear isc is a list of index offsets (equivalent to isc_off)
    !> The order of indices may be found in [[atomlist:superx]] (or in the ORB_INDX output file)
    call generate_linear_isc(nsc_old, old_isc)
    !> Create the isc_off in supercell format such that:
    !>   new_isc(1, 1, 1) yields the index for the [1, 1, 1] supercell.
    call generate_isc(nsc_new, new_isc)

    ! Calculate total number of supercells
    old_n_s = product(nsc_old)
    new_outside = product(nsc_new) * no_u

    ! We need to check whether in-put SC is too large
    new_hsc = nsc_new / 2
    do io = 1, no_l

      ! Now we can do the copy...
      inner_columns: do ind = ptr(io) + 1, ptr(io) + ncol(io)
        !> Calculate the old supercell index (supercell offset is old_isc(:, old_is))
        old_is = (col(ind)-1) / no_u
        isc = old_isc(:, old_is)

        if ( old_is >= old_n_s ) then ! it should be removed
          ! Set it to zero
          A2(ind, :) = 0._dp
          ! Also set the column index to a position outside the new sparse pattern
          col(ind) = ucorb(col(ind), no_u) + new_outside
          cycle inner_columns
        end if

        ! Check that the out supercell exists
        if ( any(abs(isc) > new_hsc) ) then
          ! Set it to zero
          A2(ind, :) = 0._dp
          ! Also set the column index to a position higher
          col(ind) = ucorb(col(ind), no_u) + new_outside
          cycle inner_columns
        end if

        ! We know the supercell exists, so convert column
        new_is = new_isc(isc(1),isc(2),isc(3))
        
        ! Transfer the orbital index to the correct supercell
        col(ind) = ucorb(col(ind), no_u) + new_is * no_u
        
      end do inner_columns
      
    end do

    deallocate(old_isc, new_isc)
    
  end subroutine correct_supercell_Sp2D

  subroutine correct_supercell_Sp1D(nsc_old, D1, nsc_new)

    use class_Sparsity
    use class_dSpData1D
#ifdef MPI
    use mpi_siesta
#endif

    type(dSpData1D), intent(inout) :: D1
    integer, intent(in) :: nsc_old(3)
    integer, intent(in) :: nsc_new(3)

    ! Local variables
    type(Sparsity), pointer :: sp
    integer :: no_u, no_l
    integer :: io, ind, isc(3)
    integer :: old_n_s, new_outside
    integer :: old_is, new_is, new_hsc(3)

    ! arrays for the sparsity patterns
    integer, pointer :: ptr(:), ncol(:), col(:)
    integer, allocatable :: old_isc(:,:), new_isc(:,:,:)
    real(dp), pointer :: A1(:)

    if ( all(nsc_old == nsc_new) ) return

    sp => spar(D1)
    call attach(sp,n_col=ncol, list_ptr=ptr, list_col=col, &
        nrows=no_l, nrows_g=no_u)

    ! Required to set removed elements to 0
    A1 => val(D1)

    ! Now create the conversion tables
    !> The linear isc is a list of index offsets (equivalent to isc_off)
    !> The order of indices may be found in [[atomlist:superx]] (or in the ORB_INDX output file)
    call generate_linear_isc(nsc_old, old_isc)
    !> Create the isc_off in supercell format such that:
    !>   new_isc(1, 1, 1) yields the index for the [1, 1, 1] supercell.
    call generate_isc(nsc_new, new_isc)

    ! Calculate total number of supercells
    old_n_s = product(nsc_old)
    new_outside = product(nsc_new) * no_u

    ! We need to check whether in-put SC is too large
    new_hsc = nsc_new / 2
    do io = 1, no_l

      ! Now we can do the copy...
      inner_columns: do ind = ptr(io) + 1, ptr(io) + ncol(io)
        !> Calculate the old supercell index (supercell offset is old_isc(:, old_is))
        old_is = (col(ind)-1) / no_u
        isc = old_isc(:, old_is)

        ! Simple error handling in case a user has supplied a too large column index
        if ( old_is >= old_n_s ) then ! it should be removed
          ! Set it to zero
          A1(ind) = 0._dp
          ! Also set the column index to a position outside the new sparse pattern
          ! since it isn't defined in the old sparse pattern, it can't be defined in
          ! the new one.
          col(ind) = ucorb(col(ind), no_u) + new_outside
          cycle inner_columns
        end if

        ! Check that the out supercell exists
        if ( any(abs(isc) > new_hsc) ) then
          ! Set it to zero
          A1(ind) = 0._dp
          ! Also set the column index to a position higher
          col(ind) = ucorb(col(ind), no_u) + new_outside
          cycle inner_columns
        end if

        ! We know the supercell exists, so convert column
        new_is = new_isc(isc(1),isc(2),isc(3))
        
        ! Transfer the orbital index to the correct supercell
        col(ind) = ucorb(col(ind), no_u) + new_is * no_u
        
      end do inner_columns
      
    end do

    deallocate(old_isc, new_isc)
    
  end subroutine correct_supercell_Sp1D

  !> Generate supercells in order of the indices 2 dimensions of size `3`, `0:product(nsc)-1`
  !>
  !> The given supercell offsets for a given supercell is:
  !>
  !> ```fortran
  !>   isc(:, is)
  !> ```
  !>
  !> where `is` is the supercell index `(list_col(ind) - 1)/no_u`.
  !> The order is equivalent to those generated in [[atomlist::superx]].
  subroutine generate_linear_isc(nsc, isc)
    integer, intent(in) :: nsc(3)
    integer, allocatable :: isc(:,:)

    integer :: x, y, z, i
    integer :: nx, ny, nz

    allocate(isc(3,0:product(nsc)-1))

    i = 0
    do z = 0, nsc(3) - 1
      nz = linear2pm(z, nsc(3))
      do y = 0, nsc(2) - 1
        ny = linear2pm(y, nsc(2))
        do x = 0, nsc(1) - 1
          nx = linear2pm(x, nsc(1))
          isc(1,i) = nx
          isc(2,i) = ny
          isc(3,i) = nz
          i = i + 1
        end do
      end do
    end do

  end subroutine generate_linear_isc

  !> Generate a supercell index array comprising of 3 dimensions of size `nx`, `ny` and `nz`
  !>
  !> To obtain the index of a given supercell simply do:
  !>
  !> ```fortran
  !>   isc(ix, iy, iz)
  !> ```
  !>
  !> where `ix`, `iy` and `iz` are the supercells for the individiual
  !> lattice vector directions. The order of supercells is equivalent to
  !> those generated in [[atomlist::superx]].
  subroutine generate_isc(nsc, isc)
    integer, intent(in) :: nsc(3)
    integer, allocatable :: isc(:,:,:)

    integer :: x, y, z, i
    integer :: nx, ny, nz
    integer :: hsc(3)

    ! nsc % 2 == 1, nsc / 2 % 2 == 0
    hsc = nsc / 2

    allocate(isc(-hsc(1):hsc(1),-hsc(2):hsc(2),-hsc(3):hsc(3)))

    i = 0
    do z = 0, nsc(3) - 1
      nz = linear2pm(z, nsc(3))
      do y = 0, nsc(2) - 1
        ny = linear2pm(y, nsc(2))
        do x = 0, nsc(1) - 1
          nx = linear2pm(x, nsc(1))
          isc(nx,ny,nz) = i
          i = i + 1
        end do
      end do
    end do

  end subroutine generate_isc

  pure function linear2pm(i,n) result(j)
    integer, intent(in) :: i, n
    integer :: j
    if ( i > n / 2 ) then
      j = -n + i
    else
      j = i
    end if
  end function linear2pm

  subroutine reduce_spin_size(ispin,H_2D,S_1D,Ef)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D
    integer, intent(in) :: ispin
    type(dSpData2D), intent(inout) :: H_2D
    type(dSpData1D), intent(inout), optional :: S_1D
    real(dp), intent(in), optional :: Ef
    type(dSpData2D) :: tmp

    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp

    integer :: dim_spin
    real(dp), pointer :: H_orig(:,:), H_new(:,:), S_orig(:)

    ! we also shift to the Fermi level
    H_orig => val(H_2D)
    if ( present(S_1D) ) then
      S_orig => val(S_1D)
    end if

    ! In case there only is one spin channel
    if ( spar_dim(H_2D) == 1 ) then
      dim_spin = 2
    else
      dim_spin = 1
    end if
    if ( size(H_orig,dim=dim_spin) == 1 ) then
      if ( present(Ef) .and. present(S_1D) ) then
!$OMP parallel workshare default(shared)
        H_orig(:,1) = H_orig(:,1) - Ef * S_orig(:)
!$OMP end parallel workshare
      end if
    else

      ! The sparsity pattern associated with the Hamiltonian
      dit => dist(H_2D)
      sp => spar(H_2D)
      call newdSpData2D(sp,1,dit,tmp)
      H_new => val(tmp)
      if ( present(Ef) .and. present(S_1D) ) then
!$OMP parallel workshare default(shared)
        H_new(:,1) = H_orig(:,ispin) - Ef * S_orig(:)
!$OMP end parallel workshare
      else
!$OMP parallel workshare default(shared)
        H_new(:,1) = H_orig(:,ispin)
!$OMP end parallel workshare
      end if

      ! Copy/delete/clean to the old array
      H_2D = tmp
      call delete(tmp)

    end if

  end subroutine reduce_spin_size

end module m_handle_sparse
