! Module for easy expansion/copying a DM file from another 
! geometry to the current one...

! This code has been fully created by:
!    Nick Papior Andersen, nickpapior@gmail.com
module m_handle_sparse

  use precision, only : dp
  use parallel, only : Node, Nodes

  use geom_helper, only : ucorb, iaorb
  use class_OrbitalDistribution
  use class_Sparsity
  use class_dSpData2D

  implicit none

contains

  subroutine bulk_expand(na_u,xa,lasto,cell,nsc,isc_off,DM_2D)
    use fdf
    use class_dSpData1D
    use m_os, only : file_exist
    use m_iodm, only : read_DM
    use m_ts_io, only : ts_read_TSHS

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
    integer :: fnsc(3), fna_u, fno_u, fnspin, fkscell(3,3), at, fn_s
    integer :: ORep(3), IRep(3)

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

    ! Read in each segment and copy data!
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle ! skip empty lines

       ! Get the name of the segment that we will copy
       ln = ' '
       ln = fdf_bnames(pline,1)
       
       ! we search for the fdf-flags
       HSfile = ' '
       HSfile = fdf_get('DM.Init.Bulk.'//trim(ln),'NONE')
       ! Now we have all required information
       if ( .not. file_exist(HSfile, Bcast = .true. ) ) then
          write(*,*) trim(HSfile)
          call die('You at least need to supply the TSHS file for &
               &bulk segment '//trim(ln)//'.')
       end if

       ! Get the repetitions
       ORep(1) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.ORep.A1',1)
       ORep(2) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.ORep.A2',1)
       ORep(3) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.ORep.A3',1)
       IRep(1) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.IRep.A1',1)
       IRep(2) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.IRep.A2',1)
       IRep(3) = fdf_get('DM.Init.Bulk.'//trim(ln)//'.IRep.A3',1)

       ! We first try and guess the DM file
       at = len_trim(HSfile)
       DMfile = fdf_get('DM.Init.Bulk.'//trim(ln)//'.DM', &
            HSfile(1:at-4)//'DM')
       if ( .not. file_exist(DMfile, Bcast = .true.) ) then
          DMfile = fdf_get('DM.Init.Bulk.'//trim(ln)//'.DM', &
               HSfile(1:at-4)//'TSDE')
       end if
       if ( .not. file_exist(DMfile, Bcast = .true. ) ) then
          call die('DM file could not be found, have you supplied an &
               &erroneous path?')
       end if

       ! Get the starting atom we need to copy into!
       at = fdf_get('DM.Init.Bulk.'//trim(ln)//'.InsertAt',0)
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
            fkscell, fkdispl, &
            fxa, flasto, &
            fsp, fDM_2D, tmp_1D, fisc_off, &
            fEf, fQtot, fTemp, &
            itmp, itmp, &
            Bcast= .true. )
       fn_s = product(fnsc)

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
       call read_DM( DMfile, fake_dit, fDM_2D, d_log1, Bcast=.true.)
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
          itmp = product(ORep) * product(IRep) * c_na - 1
          write(*,'(a,i0,'' ,'',i0,2a)') &
               'siesta: Initializing bulk DM for atoms [ ',at,at+itmp,']', &
               ' using segment: '//trim(ln)
       end if

       call expand_spd2spd_2D(s_na,c_na,fna_u,flasto,fxa, &
            fDM_2D,&
            fcell, ORep, IRep, fn_s, fisc_off, &
            na_u,xa,lasto,DM_2D,cell,product(nsc),isc_off, at, &
            print = .true., allowed_a = allowed)

       ! De-allocate before reading the next thing...
       call delete(fDM_2D)
       deallocate(fxa) ; nullify(fxa)
       deallocate(flasto) ; nullify(flasto)
       deallocate(fisc_off) ; nullify(fisc_off)
       
    end do
    if ( Node == 0 ) write(*,*) ! new-line
    
  end subroutine bulk_expand

  subroutine expand_spd2spd_2D(s_a,na,na_i,lasto_i,xa_i,in,&
       cell_i, orep_i, irep_i, n_s_i, sc_off_i, &
       na_o,xa_o,lasto_o,out,cell_o,n_s_o,sc_off_o, at, &
       print, allowed_a)

#ifdef MPI
    use mpi_siesta
#endif

    ! Starting atom and number of atoms from the starting atom.
    integer, intent(in) :: s_a, na
    integer, intent(in) :: na_i, lasto_i(0:na_i), orep_i(3), irep_i(3)
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
    real(dp), parameter :: xa_EPS = 1.e-3_dp

    integer, allocatable :: lallow(:)

#ifdef MPI
    integer :: tmp_copy(3)
    integer :: MPIerror
#endif

    o1 = product(orep_i)
    i1 = product(irep_i)
    i = at - 1 + o1 * i1 * na
    if ( i > na_o ) then
       call die('Requested expansion region too large.')
    end if

    if ( present(allowed_a) ) then
       ! We copy over the allowed atoms
       allocate(lallow(size(allowed_a)))
       lallow(:) = allowed_a(:)
    else
       ! We only allow copying the diagonal entries
       allocate(lallow(o1 * i1 * na))
       do i = 1 , o1 * i1 * na
          lallow(i) = at + i - 1
       end do
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

    o3_loop: do o3 = 0 , orep_i(3) - 1
    o2_loop: do o2 = 0 , orep_i(2) - 1
    o1_loop: do o1 = 0 , orep_i(1) - 1

    at_end = at_end + product(irep_i) * na
    
    ! We loop over the input SP which we will copy
    do ia_i = s_a , s_a + na - 1

     i3_loop: do i3 = 0 , irep_i(3) - 1
     i2_loop: do i2 = 0 , irep_i(2) - 1
     i1_loop: do i1 = 0 , irep_i(1) - 1

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
            print *,'ia',ia_i,'matching',iat,xc_o-xc_i,i1,i2,i3
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
        if ( .not. any(ao == lallow) ) cycle

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

    lat = lat + product(irep_i) * na

    end do o1_loop
    end do o2_loop
    end do o3_loop

    deallocate(lallow)

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


  subroutine reduce_spin_size(ispin,H_2D,S_1D,Ef)
    use class_dSpData1D
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
