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

  subroutine expand_spd2spd_2D(na_i,lasto_i,xa_i,in,&
       cell_i, rep_i, n_s_i, sc_off_i, &
       na_o,xa_o,lasto_o,out,cell_o,n_s_o,sc_off_o, at, &
       print, allowed_a)

#ifdef MPI
    use mpi_siesta
#endif

    integer, intent(in) :: na_i, lasto_i(0:na_i), rep_i(3)
    real(dp), intent(in) :: xa_i(3,na_i), cell_i(3,3)
    integer, intent(in) :: n_s_i, sc_off_i(3,n_s_i)
    ! The density matrices that describes two
    ! different parts.
    type(dSpData2D), intent(inout) :: in, out
    integer, intent(in) :: na_o, lasto_o(0:na_o)
    real(dp), intent(in) :: xa_o(3,na_o), cell_o(3,3)
    integer, intent(in) :: n_s_o, sc_off_o(3,n_s_o)
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
    integer :: i_s, i, ao, no_o, no_i, i1, i2, i3, at_end
    integer :: orb_i, orb_o
    integer :: copy(2)

    real(dp) :: xc_i(3), xc_o(3), xj_i(3), xj_o(3)
    real(dp), pointer :: a_i(:,:), a_o(:,:)
    real(dp), parameter :: xa_EPS = 1.e-3_dp

    integer, allocatable :: lallow(:)

#ifdef MPI
    integer :: MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif

    i = at - 1 + product(rep_i) * na_i
    if ( i > na_o ) then
       call die('Requested expansion region too large.')
    end if

    if ( present(allowed_a) ) then
       ! We copy over the allowed atoms
       allocate(lallow(size(allowed_a)))
       lallow(:) = allowed_a(:)
    else
       ! We only allow copying the diagonal entries
       allocate(lallow(product(rep_i)*na_i))
       do i = 1 , product(rep_i) * na_i
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

    at_end = at + na_i * product(rep_i) - 1

    ! Loop on all equivalent atoms
    iat = at - 1
    ! We count number of copied data
    copy(:) = 0
    
    ! We loop over the input SP which we will copy
    do ia_i = 1 , na_i

     i3_loop: do i3 = 0 , rep_i(3) - 1
     i2_loop: do i2 = 0 , rep_i(2) - 1
     i1_loop: do i1 = 0 , rep_i(1) - 1

      ! The expanded atomic position
      xc_i(:) = xa_i(:,ia_i) - xa_i(:,1) + &
           cell_i(:,1)*i1+cell_i(:,2)*i2+cell_i(:,3)*i3

      ! Step atom index in out-put array
      iat = iat + 1

      xc_o(:) = xa_o(:,iat) - xa_o(:,at)
      
      if ( maxval(abs(xc_o - xc_i)) > xa_EPS ) then
         print *,'ia',ia_i,'matching',iat,xc_o-xc_i
         call die('Atomic coordinates does not coincide, &
              &have you employed correct ordering for expansion A1->A2->A3?')
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

       ! Loop over indices in this row
       ! the large sparsity pattern must be the largest
       ! hence we have this as the outer loop
       do i_o = o_ptr(lio_o) + 1 , o_ptr(lio_o) + o_ncol(lio_o)

        ! First we figure out which atomic position this
        ! corresponds to:
        i_s = (o_col(i_o)-1)/no_o + 1
        ao = iaorb(o_col(i_o),lasto_o)
        orb_o = ucorb(o_col(i_o),no_o) - lasto_o(ao-1)
        ! Do not allow overwriting DM outside of region.
        if ( .not. any(ao==lallow) ) cycle
        xj_o(:) = xa_o(:,ao) - xa_o(:,at) + &
             cell_o(:,1) * sc_off_o(1,i_s) + &
             cell_o(:,2) * sc_off_o(2,i_s) + &
             cell_o(:,3) * sc_off_o(3,i_s)        
 
        ! Now we need to figure out all the orbitals
        ! that has the same meaning in both sparsity patterns
        ! Get the cell-offset

        do i_i = i_ptr(io_i) + 1 , i_ptr(io_i) + i_ncol(io_i)

           i_s = (i_col(i_i)-1)/no_i + 1
           i = iaorb(i_col(i_i),lasto_i)
           orb_i = ucorb(i_col(i_i),no_i) - lasto_i(i-1)

           ! Different orbitals can have the same
           ! orbital center, so we check the orbital
           ! index to be the same as well...
           if ( orb_i /= orb_o ) cycle
           
           xj_i(:) = xa_i(:,i) - xa_i(:,1) + &
                cell_i(:,1) * sc_off_i(1,i_s) + &
                cell_i(:,2) * sc_off_i(2,i_s) + &
                cell_i(:,3) * sc_off_i(3,i_s) + &
                cell_i(:,1)*i1+cell_i(:,2)*i2+cell_i(:,3)*i3
           
           ! If they are not equivalent we will not do anything
           if ( maxval(abs(xj_o - xj_i)) > xa_EPS ) cycle

           ! WUHUU, we have the equivalent atom and equivalent
           ! orbital connection. We copy data now!

           if ( at <= ao .and. ao <= at_end ) then
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

    deallocate(lallow)

    if ( .not. present(print) ) return
    if ( .not. print ) return

    io_i = copy(1)
    io_o = copy(2)
    if ( Node == 0 ) then
       write(*,'(a,''[ '',i0,'', '',i0,'']'',a,i0)') &
            'Expanding ',copy,' elements on node ',0
    end if
    
#ifdef MPI
    if ( Node == 0 ) then
       do iat = 1 , Nodes - 1
          call MPI_Recv(copy,2,MPI_Integer, &
               iat, 0, MPI_Comm_World, MPIstatus, MPIerror)
          write(*,'(a,''[ '',i0,'', '',i0,'']'',a,i0)') &
               'Expanding ',copy,' elements on node ',iat
          io_i = io_i + copy(1)
          io_o = io_o + copy(2)
       end do
    else
       call MPI_Send( copy , 2, MPI_Integer, &
            0, 0, MPI_Comm_World, MPIerror)
    end if
#endif
    
    if ( Node == 0 ) then
       copy(1) = io_i
       copy(2) = io_o
       write(*,'(a,''[ '',i0,'', '',i0,'']'',a)') &
            'Expanded in total ',copy,' elements.'
    end if
    
  end subroutine expand_spd2spd_2D

end module m_handle_sparse
