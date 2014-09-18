!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
!
module m_sparse

  use precision, only : dp
  use geom_helper, only : ucorb, iaorb
  use parallel, only : Node, Nodes

  implicit none

  private

  public :: nsc_to_offset
  public :: offset2idx

  public :: calc_nsc
  public :: list_col_correct

  interface xij_offset
     module procedure xij_offset_sp, xij_offset_direct
  end interface xij_offset
  public :: xij_offset

  interface offset_xij
     module procedure offset_xij_sp, offset_xij_direct
  end interface offset_xij
  public :: offset_xij

  integer, parameter, public :: TRANSFER_ALL = -999999

contains

  subroutine nsc_to_offset(nsc,sc_off)
    ! Number of supercells in each direction
    integer, intent(in) :: nsc(3)
    ! index of offsets
    integer, pointer :: sc_off(:,:)

    ! ** local variables
    integer :: n_s
    integer :: ia, ib, ic, is

    nullify(sc_off)

    n_s = product(nsc)
    allocate(sc_off(3,n_s))

    ! Create offsets
    do ic = -nsc(3)/2 , nsc(3)/2
    do ib = -nsc(2)/2 , nsc(2)/2
    do ia = -nsc(1)/2 , nsc(1)/2
       is = offset2idx(nsc,(/ia,ib,ic/))
       sc_off(1,is) = ia
       sc_off(2,is) = ib
       sc_off(3,is) = ic
    end do
    end do
    end do
    
  end subroutine nsc_to_offset

  ! From lnsc we calculate the supercell index
  function offset2idx(nsc,tm) result(idx)
    integer, intent(in) :: nsc(3), tm(3)
    integer :: idx
    integer :: ia, ib, ic

    idx = 0
    do ic = -nsc(3)/2 , nsc(3)/2
    do ib = -nsc(2)/2 , nsc(2)/2
    do ia = -nsc(1)/2 , nsc(1)/2
       idx = idx + 1
       if ( tm(1) == ia .and. &
            tm(2) == ib .and. &
            tm(3) == ic ) return
    end do
    end do
    end do
    idx = 0

  end function offset2idx


  subroutine calc_nsc(ucell,na_u,xa,lasto,no_l,no_u,n_nzs, &
       ncol,l_ptr,l_col,xij,nsc,Bcast)
    use cellSubs, only : reclat
#ifdef MPI
    use parallelsubs, only : LocalToGlobalOrb
    use mpi_siesta
#endif

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: ucell(3,3) ! The unit cell of system
    integer, intent(in)  :: na_u ! Unit cell atoms
    real(dp), intent(in) :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in)  :: lasto(0:na_u) ! last orbital number of equivalent atom
    integer, intent(in)  :: no_l,no_u ! Unit cell orbitals
    integer, intent(in)  :: n_nzs ! Hamiltonian size
    integer, intent(in)  :: ncol(no_l), l_ptr(no_l)
    integer, intent(in)  :: l_col(n_nzs)
    real(dp), intent(in) :: xij(3,n_nzs) ! differences with unitcell, differences with unitcell
    logical, intent(in), optional :: Bcast
! ***********************
! * OUTPUT variables    *
! ***********************
    integer, intent(out) :: nsc(3)

! ***********************
! * LOCAL variables     *
! ***********************
    logical :: lBcast
    real(dp) :: recell(3,3)
    real(dp) :: xijo(3), xc
    integer :: ia, ja
    integer :: lio, io, jo, ind
#ifdef MPI
    integer :: MPIerror, tnsc(3)
#endif

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    ! Initialize the transfer cell to:
    nsc(:) = 0

    ! Prepare the cell to calculate the index of the atom
    call reclat(ucell,recell,0) ! Without 2*Pi
    
    do lio = 1 , no_l
#ifdef MPI
       if ( lBcast ) then
          call LocalToGlobalOrb(lio,Node,Nodes,io)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)
          jo = ucorb(l_col(ind),no_u)
          ja = iaorb(jo,lasto)
          xijo(:) = xij(:,ind) - ( xa(:,ja)-xa(:,ia) )

          ! Loop over directions
          ! recell is already without 2*Pi
          xc = sum(xijo(:) * recell(:,1))
          nsc(1) = max(nsc(1),abs(nint(xc)))
          xc = sum(xijo(:) * recell(:,2))
          nsc(2) = max(nsc(2),abs(nint(xc)))
          xc = sum(xijo(:) * recell(:,3))
          nsc(3) = max(nsc(3),abs(nint(xc)))
       end do
    end do
    
#ifdef MPI
    if ( lBcast ) then
       tnsc = nsc
       call MPI_AllReduce(tnsc,nsc,3,MPI_Integer, &
            MPI_Max,MPI_Comm_World,MPIerror)
    end if
#endif
    
    ! Calculate the actual number of super-cells
    nsc(:) = 2 * nsc(:) + 1
    
  end subroutine calc_nsc

  ! Create an index array containing the unit-cell expansions
  ! This means we can do this:
  !   xij(:,ind) = ucell(:,1) * tm(1) + ucell(:,2) * tm(2) + ucell(:,3) * tm(3) + xa(:,iaorb(jo))-xa(:,iaorb(io))
  ! to create the full xij array...
  subroutine list_col_correct(ucell,na_u,no_l,no_u,n_nzs, &
       lasto, xa, ncol, l_ptr, l_col, xij, nsc, Bcast)
#ifdef MPI
    use parallelsubs, only : LocalToGlobalOrb
#endif
    use cellSubs, only : reclat

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: ucell(3,3) ! The unit cell of system
    integer, intent(in)  :: na_u ! Unit cell atoms
    integer, intent(in)  :: no_l,no_u ! Unit cell orbitals
    integer, intent(in)  :: n_nzs ! Hamiltonian size
    integer, intent(in)  :: lasto(0:na_u) ! Last orbital of atom
    real(dp), intent(in) :: xa(3,na_u) ! atom positions
    integer, intent(in)  :: ncol(no_l), l_ptr(no_l)
    ! We correct the indices here to conform with the corrected offsets
    integer, intent(inout)  :: l_col(n_nzs)
    real(dp), intent(in) :: xij(3,n_nzs) ! differences with unitcell, differences with unitcell
    ! Number of supercells in each direction ** MUST be corrected using nsc_to_offests
    integer, intent(in)  :: nsc(3) ! Number of supercells in each direction
    logical, intent(in), optional :: Bcast
! ***********************
! * LOCAL variables     *
! ***********************
    logical :: lBcast
    integer :: lio, io, jo, ind, ia, ja, is, tm(3)
    real(dp) :: xijo(3), rcell(3,3)

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    
    ! Prepare the cell to calculate the index of the atom
    call reclat(ucell,rcell,0) ! Without 2*Pi

    !err = 0._dp
    do lio = 1 , no_u
#ifdef MPI
       if ( lBcast ) then
          call LocalToGlobalOrb(lio,Node,Nodes,io)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)

          jo = ucorb(l_col(ind),no_u)
          ja = iaorb(jo,lasto)

          xijo(:) = xij(:,ind) - ( xa(:,ja) - xa(:,ia) )
          tm(:) = nint( matmul(xijo,rcell) )

          ! get supercell index
          is = offset2idx(nsc,tm) - 1

          if ( is < 0 ) then
             ! Index not found
             call die('Error in index, possible shifting')
          end if

          ! Correct index
          l_col(ind) = is * no_u + jo

          ! The error on this conversion
          ! tends to be on the order of 1e-14
          ! xijo = ucell(:,1) * tm(1) &
          !      + ucell(:,2) * tm(2) &
          !      + ucell(:,3) * tm(3) &
          !      + xa(:,ja) - xa(:,ia)
          ! err = max(maxval(abs(xijo - xij(:,ind))),err)
          ! This error could be "important" when
          ! calculating self-energies

       end do
    end do

  end subroutine list_col_correct


  subroutine xij_offset_sp(ucell,nsc, na_u,xa,lasto, &
       block_dist,sp,n_nzs,xij,isc_off,Bcast)
    
    use cellSubs, only : reclat
    use intrinsic_missing, only : VNORM
    use class_OrbitalDistribution
    use class_Sparsity
    use alloc
#ifdef MPI
    use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    ! Unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of super-cells in each direction
    integer, intent(in) :: nsc(3)
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    real(dp), intent(in) :: xa(3,na_u)
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    ! The distribution for the sparsity-pattern
    type(OrbitalDistribution), intent(inout) :: block_dist
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! number of non-zero elements in H
    integer, intent(in) :: n_nzs
    ! vectors from i-J
    real(dp), intent(in) :: xij(3,n_nzs)
    logical, intent(in), optional :: Bcast
! **********************
! * OUTPUT variables   *
! **********************
    integer, pointer :: isc_off(:,:)

! **********************
! * LOCAL variables    *
! **********************
    logical :: lBcast
#ifdef MPI
    integer, allocatable :: ioff_2(:,:)
    integer :: MPIerror, iNode
#endif
    integer :: n_s, tm(3), is
    integer :: no_l, no_u
    integer :: lio, io, ind, ia, ja
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    real(dp) :: xijo(3), rcell(3,3)

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    ! Attach to the sparsity pattern
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    ! Number of super cells
    n_s = product(nsc)
    ! Calculate the offsets (instead of using xij)
    call re_alloc(isc_off,1,3,1,n_s)

    ! Initialize the integer offset
    isc_off(:,:) = 0
    
    ! Prepare the cell to calculate the index of the atom
    call reclat(ucell,rcell,0) ! Without 2*Pi

    do lio = 1 , no_l

       if ( ncol(lio) == 0 ) cycle

#ifdef MPI
       if ( lBcast ) then
          ! Transfer to global index
          io = index_local_to_global(block_dist,lio,Node)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)

          ja = iaorb(ucorb(l_col(ind),no_u),lasto)

          ! the supercell index (counting from one)
          is = (l_col(ind) - 1)/no_u + 1

          xijo(:) = xij(:,ind) - ( xa(:,ja) - xa(:,ia) )

          tm(:) = nint( matmul(xijo,rcell) )

          if ( any(tm(:) /= isc_off(:,is)) .and. any(isc_off(:,is)/=0) ) then
             write(*,'(2(a10,2(tr1,i6)))')'r,C',io,l_col(ind),'ia,ja',ia,ja
             write(*,'(a,i3,2(tr2,3(tr1,i3)))') 'is, tm, old_tm: ',is, tm, isc_off(:,is)
             write(*,'(a,4(tr1,f10.5))') 'xij - ucell*tm, |V|: ',xijo(:) - &
                  ucell(:,1) * tm(1) - ucell(:,2) * tm(2) - ucell(:,3) * tm(3), vnorm(xijo)
             write(*,'(2(a10,3(tr1,f10.5)))') 'xijo: ',xijo(:), &
                  'ucell: ',ucell(:,1) * tm(1) + ucell(:,2) * tm(2) + ucell(:,3) * tm(3)
             call die('Error on transfer matrix indices...')
          else
             isc_off(:,is) = tm(:)
          end if

       end do
    end do

#ifdef MPI
    if ( lBcast ) then
    ! Reduce all sc_off
    allocate(ioff_2(3,n_s))

    do iNode = 0 , Nodes - 1
       if ( iNode == Node ) ioff_2(:,:) = isc_off(:,:)
       call MPI_Bcast(ioff_2(1,1),3*n_s, &
            MPI_Integer, iNode, MPI_Comm_World, MPIerror)

       ! Search for differing sc_off
       do is = 1 , n_s
          if ( all(isc_off(:,is) == 0) ) then
             ! it hasn't been found locally, just copy
             isc_off(:,is) = ioff_2(:,is)
          else if ( all(ioff_2(:,is) == 0) ) then
             ! do nothing, this is not set on other node
          else if ( any(isc_off(:,is) /= ioff_2(:,is)) ) then
             print '(2(tr1,a,3(tr1,i2)))','Local:',isc_off(:,is),'Found:',ioff_2(:,is)
             ! We do not have the same supercell indices
             call die('Supercell reduction is erroneous. &
                  &Please consult the developers.')
          end if
          
       end do

    end do

    deallocate(ioff_2)
    end if
#endif
    
  end subroutine xij_offset_sp

  subroutine xij_offset_direct(ucell,nsc,na_u,xa,lasto, &
       no_l,no_u,n_nzs,ncol,l_ptr,l_col,xij,isc_off,Bcast)
    
    use cellSubs, only : reclat
    use intrinsic_missing, only : VNORM
    use alloc
#ifdef MPI
    use parallelsubs, only : LocalToGlobalOrb
    use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    ! Unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of super-cells in each direction
    integer, intent(in) :: nsc(3)
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    real(dp), intent(in) :: xa(3,na_u)
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    ! Sizes
    integer, intent(in) :: no_l, no_u, n_nzs
    ! sparsity pattern
    integer, intent(in) :: ncol(no_l), l_ptr(no_l), l_col(n_nzs)
    ! vectors from i-J
    real(dp), intent(in) :: xij(3,n_nzs)
    logical, intent(in), optional :: Bcast

! **********************
! * OUTPUT variables   *
! **********************
    integer, pointer :: isc_off(:,:)

! **********************
! * LOCAL variables    *
! **********************
#ifdef MPI
    integer, allocatable :: ioff_2(:,:)
    integer :: MPIerror, iNode
#endif
    logical :: lBcast
    integer :: n_s, is, tm(3)
    integer :: lio, io, ind, ia, ja
    real(dp) :: xijo(3), rcell(3,3)

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    ! Number of super cells
    n_s = product(nsc)
    ! Calculate the offsets (instead of using xij)
    call re_alloc(isc_off,1,3,1,n_s)

    ! Initialize the integer offset
    isc_off(:,:) = 0
    
    ! Prepare the cell to calculate the index of the atom
    call reclat(ucell,rcell,0) ! Without 2*Pi

    do lio = 1 , no_l

       if ( ncol(lio) == 0 ) cycle

#ifdef MPI
       if ( lBcast ) then
          call LocalToGlobalOrb(lio,Node,Nodes,io)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)

          ja = iaorb(ucorb(l_col(ind),no_u),lasto)

          ! the supercell index (counting from 1)
          is = (l_col(ind) - 1)/no_u + 1

          xijo(:) = xij(:,ind) - ( xa(:,ja) - xa(:,ia) )

          tm(:) = nint( matmul(xijo,rcell) )

          if ( any(tm(:) /= isc_off(:,is)) .and. any(isc_off(:,is)/=0) ) then
             write(*,'(2(a10,2(tr1,i6)))')'r,C',io,l_col(ind),'ia,ja',ia,ja
             write(*,'(a,i3,2(tr2,3(tr1,i3)))') 'is, tm, old_tm: ',is, tm, isc_off(:,is)
             write(*,'(a,4(tr1,f10.5))') 'xij - ucell*tm, |V|: ',xijo(:) - &
                  ucell(:,1) * tm(1) - ucell(:,2) * tm(2) - ucell(:,3) * tm(3), vnorm(xijo)
             write(*,'(2(a10,3(tr1,f10.5)))') 'xijo: ',xijo(:), &
                  'ucell: ',ucell(:,1) * tm(1) + ucell(:,2) * tm(2) + ucell(:,3) * tm(3)
             call die('Error on transfer matrix indices...')
          else
             isc_off(:,is) = tm(:)
          end if

       end do
    end do

#ifdef MPI
    if ( lBcast ) then
    ! Reduce all sc_off
    allocate(ioff_2(3,n_s))

    do iNode = 0 , Nodes - 1
       if ( iNode == Node ) ioff_2(:,:) = isc_off(:,:)
       call MPI_Bcast(ioff_2(1,1),3*n_s, &
            MPI_Integer, iNode, MPI_Comm_World, MPIerror)

       ! Search for differing sc_off
       do is = 1 , n_s
          if ( all(isc_off(:,is) == 0) ) then
             ! it hasn't been found locally, just copy
             isc_off(:,is) = ioff_2(:,is)
          else if ( all(ioff_2(:,is) == 0) ) then
             ! do nothing, this is not set on other node
          else if ( any(isc_off(:,is) /= ioff_2(:,is)) ) then
             print '(2(tr1,a,3(tr1,i2)))','Local:',isc_off(:,is),'Found:',ioff_2(:,is)
             ! We do not have the same supercell indices
             call die('Supercell reduction is erroneous. &
                  &Please consult the developers.')
          end if
          
       end do

    end do

    deallocate(ioff_2)

    end if

#endif

  end subroutine xij_offset_direct

  subroutine offset_xij_sp(ucell,n_s,isc_off, na_u,xa,lasto, &
       block_dist,sp,n_nzs,xij,Bcast)
    
    use class_OrbitalDistribution
    use class_Sparsity
    use alloc

! **********************
! * INPUT variables    *
! **********************
    ! Unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of super-cells
    integer, intent(in) :: n_s
    ! Integer supercells
    integer, intent(in) :: isc_off(3,0:n_s-1)
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    real(dp), intent(in) :: xa(3,na_u)
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    ! The distribution for the sparsity-pattern
    type(OrbitalDistribution), intent(inout) :: block_dist
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    integer, intent(in) :: n_nzs
    logical, intent(in), optional :: Bcast
! **********************
! * OUTPUT variables   *
! **********************
    real(dp), intent(out) :: xij(3,n_nzs)

! **********************
! * LOCAL variables    *
! **********************
    logical :: lBcast
    integer :: tm(3), is
    integer :: no_l, no_u
    integer :: lio, io, ind, ia, ja
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    ! Attach to the sparsity pattern
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    do lio = 1 , no_l

       if ( ncol(lio) == 0 ) cycle

#ifdef MPI
       if ( lBcast ) then
          ! Transfer to global index
          io = index_local_to_global(block_dist,lio,Node)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)

          ja = iaorb(ucorb(l_col(ind),no_u),lasto)

          ! the supercell index (counting from zero)
          is = (l_col(ind) - 1)/no_u

          tm(:) = isc_off(:,is)
          xij(:,ind) = tm(1)*ucell(:,1)+tm(2)*ucell(:,2)+tm(3)*ucell(:,3)
          xij(:,ind) = xij(:,ind) + xa(:,ja) - xa(:,ia)

       end do
    end do

  end subroutine offset_xij_sp

  subroutine offset_xij_direct(ucell,n_s,isc_off, na_u,xa,lasto, &
       no_l,no_u,n_nzs,ncol,l_ptr,l_col,xij,Bcast)
    
    use alloc
#ifdef MPI
    use parallelsubs, only : LocalToGlobalOrb
#endif

! **********************
! * INPUT variables    *
! **********************
    ! Unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of super-cells
    integer, intent(in) :: n_s
    ! Integer supercells
    integer, intent(in) :: isc_off(3,0:n_s-1)
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    real(dp), intent(in) :: xa(3,na_u)
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    ! The sparsity pattern
    integer, intent(in) :: no_l, no_u, n_nzs, ncol(no_l), l_ptr(no_l), l_col(n_nzs)
    ! vectors from i-J
    logical, intent(in), optional :: Bcast
! **********************
! * OUTPUT variables   *
! **********************
    real(dp), intent(out) :: xij(3,n_nzs)

! **********************
! * LOCAL variables    *
! **********************
    logical :: lBcast
    integer :: tm(3), is
    integer :: lio, io, ind, ia, ja

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    do lio = 1 , no_l

       if ( ncol(lio) == 0 ) cycle

#ifdef MPI
       if ( lBcast ) then
          call LocalToGlobalOrb(lio,Node,Nodes,io)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)

          ja = iaorb(ucorb(l_col(ind),no_u),lasto)

          ! the supercell index (counting from zero)
          is = (l_col(ind) - 1)/no_u

          tm(:) = isc_off(:,is)
          xij(:,ind) = tm(1)*ucell(:,1)+tm(2)*ucell(:,2)+tm(3)*ucell(:,3)
          xij(:,ind) = xij(:,ind) + xa(:,ja) - xa(:,ia)

       end do
    end do

  end subroutine offset_xij_direct

end module m_sparse
  
