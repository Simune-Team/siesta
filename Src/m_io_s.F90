! Module for easy reading of sparsity patterns and 
! data distributed in sparsity patterns.

! Ideally all IO routines should utilize these routines

! Fully implemented by Nick Papior Andersen

! SIESTA io-module
module m_io_s

  use class_OrbitalDistribution
  use class_Sparsity
  use precision, only : sp, dp
  use parallel, only : IONode, Node, Nodes
#ifdef MPI
  use mpi_siesta, only : MPI_Bcast, MPI_AllReduce, MPI_Sum, MPI_Max
  use mpi_siesta, only : MPI_Comm_World, MPI_Comm_Self
  use mpi_siesta, only : MPI_Integer, MPI_Double_Precision
  use mpi_siesta, only : MPI_Success, MPI_Status_Size
#endif

  implicit none 

  private

  public :: io_read_Sp
  public :: io_read_d1D, io_read_d2D

contains

  ! Reads in a sparsity pattern at the
  ! current position in the file (iu)
  ! The sparsity pattern "sp" will be returned
  ! as populated.
  ! If dist is supplied it will distribute
  ! the sparsity pattern as supplied (this implies Bcast = .true.)
  ! Else if Bcast is true it will b-cast the sparsity 
  ! pattern fully.
  subroutine io_read_Sp(iu, no, sp, tag, dit, Bcast)

    ! File handle
    integer, intent(in) :: iu
    ! Number of orbitals readed
    integer, intent(in) :: no
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(inout), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast

    ! Local variables for reading the values
    integer, pointer :: ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    integer :: io, n_nzs, ind, nl
    logical :: ldit, lBcast
    integer, pointer :: buf(:) => null()
#ifdef MPI
    integer, pointer :: buf2(:) => null()
    integer :: gio, max_n
    integer :: MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq, BNode
#endif

#ifdef MPI
    ldit = present(dit)
#else
    ldit = .false.
#endif
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    if ( IONode ) then

       allocate(buf(no))

       ! First read in number of non-zero 
       ! entries per orbital
       read(iu) buf

#ifdef MPI
       max_n = maxval(buf)
#endif

    end if

    nl = no

    ! If a distribution is present, then do something
#ifdef MPI
    if ( ldit ) then

       nl = 0
       ! First count number of local entries
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          if ( BNode == Node ) nl = nl + 1
       end do

       ! allocate room for the number of columns in
       ! each row
       allocate(ncol(nl))
       MPIerror = MPI_Success - 1

       ! distribute it
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)

          if ( BNode == Node ) then
             ! We have the node, get the local orbital
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                ncol(io) = buf(gio)
             else
                call MPI_Recv( ncol(io) , 1, MPI_Integer, &
                     0, gio, MPI_Comm_World, MPIstatus, MPIerror )
             end if
          else if ( Node == 0 ) then
             call MPI_ISSend( buf(gio) , 1, MPI_Integer, &
                  BNode, gio, MPI_Comm_World, MPIreq, MPIerror)
          end if
          
       end do

       if ( IONode .and. MPIerror /= MPI_Success - 1 ) then
          call MPI_Wait(MPIreq,MPIstatus,MPIerror)
       end if

    else if ( lBcast ) then

       ! Everything should be b-casted
       if ( IONode ) then
          ncol => buf
       else
          allocate(ncol(nl))
       end if
       
       ! Bcast everything
       call MPI_Bcast(ncol(1),nl,MPI_Integer,0,MPI_Comm_World, &
            MPIError)
       
    else if ( .not. IONode ) then
       ! no distribution, no b-cast.
       ! The sparsity pattern will only exist on the IONode
       return ! the sparsity pattern will not be created then...
    end if
#else
    ! Point to the buffer
    ncol => buf
#endif

    ! Allocate pointer
    allocate(l_ptr(nl))
    
    l_ptr(1) = 0
    do io = 2 , nl
       l_ptr(io) = l_ptr(io-1) + ncol(io-1)
    end do
    
    ! Number of local non-zero elements
    ! (also works for any bcast methods)
    n_nzs = l_ptr(nl) + ncol(nl)

    ! Allocate space
    allocate(l_col(n_nzs))

#ifdef MPI
    if ( ldit ) then
       
       ! We have a distributed read
       if ( IONode ) then
          allocate(buf2(max_n))
       end if

       MPIerror = MPI_Success - 1

       ! read in columns
       ind = 0
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          
          if ( BNode == Node ) then

             ! We have the node, get the local orbital
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                read(iu) l_col(ind+1:ind+ncol(io))
             else
                ! No need to block the receive...
                call MPI_IRecv( l_col(ind+1) , ncol(io), MPI_Integer, &
                     0, gio, MPI_Comm_World, MPIreq, MPIerror )
             end if

             ind = ind + ncol(io)

          else if ( Node == 0 ) then

             read(iu) buf2(1:buf(gio))
             call MPI_Send( buf2(1) , buf(gio), MPI_Integer, &
                  BNode, gio, MPI_Comm_World, MPIerror)

          end if
          
       end do

       if ( IONode ) then
          deallocate(buf,buf2)
       else if ( MPIerror /= MPI_Success - 1 ) then
          call MPI_Wait(MPIreq,MPIstatus,MPIerror)
       end if

    else if ( lBcast ) then

       if ( IONode ) then
          
          ind = 0
          do gio = 1 , no
             read(iu) l_col(ind+1:ind+ncol(gio))
             ind = ind + ncol(gio)
          end do

       end if

       ! Bcast
       call MPI_Bcast(l_col(1),n_nzs,MPI_Integer,0,MPI_Comm_World, &
            MPIError)

    end if

#else
    
    ind = 0
    do io = 1 , no
       read(iu) l_col(ind+1:ind+ncol(io))
       ind = ind + ncol(io)
    end do

#endif

    ! Create the sparsity pattern
    call newSparsity(sp,nl,no, &
         n_nzs, ncol, l_ptr, l_col, trim(tag))
    
    deallocate(ncol,l_ptr,l_col)
    
  end subroutine io_read_Sp

  ! Writes a sparsity pattern at the
  ! current position in the file (iu)
  ! If dist is supplied it will write a distributed sparsity pattern
  subroutine io_write_Sp(iu, sp, dit)

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! distribution
    type(OrbitalDistribution), intent(inout), optional :: dit

    ! Local variables for reading the values
    integer, pointer :: buf(:) => null()
    integer, pointer :: ncol(:) => null(), l_col(:) => null()

    integer :: lno, no, io, max_n, ind
    logical :: ldit
#ifdef MPI
    integer :: gio
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif

    ! Get the sparsity sizes
    call attach(sp,n_col=ncol, list_col=l_col, nrows=lno,nrows_g=no)

    ldit = present(dit)
    if ( ldit ) ldit = lno /= no

    if ( ldit ) then

#ifdef MPI
       ! Allocate the full ncol
       allocate(buf(no))
       MPIerror = MPI_Success - 1
          
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                buf(gio) = ncol(io)
             else
                call MPI_ISSend( ncol(io) , 1, MPI_Integer, &
                     0, gio, MPI_Comm_World, MPIreq, MPIerror)
             end if
          else if ( Node == 0 ) then
             call MPI_Recv( buf(gio) , 1, MPI_Integer, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
          end if
       end do

       if ( .not. Node == 0 .and. MPIerror /= MPI_Success - 1) then
          ! Wait for the last one to not send
          ! to messages with the same tag...
          call MPI_Wait(MPIreq,MPIstatus,MPIerror)
       end if
#else
       buf => ncol
#endif

    else
       
       buf => ncol
       
    end if
    
    if ( Node == 0 ) then
       
       write(iu) buf

       ! Retrive the maximum number of non-zero
       ! elements in each row
       max_n = maxval(buf)

    end if

#ifdef MPI
    if ( ldit ) then
       deallocate(buf)
    end if
#endif

#ifdef MPI
    ! Write the list_col array
    if ( ldit ) then

       ! The ionode now has the maximum retrieved array
       if ( Node == 0 ) then
          nullify(buf)
          allocate(buf(max_n))
       end if
       MPIerror = MPI_Success - 1

       ! Loop size
       ind = 0
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                write(iu) l_col(ind+1:ind+ncol(io))
             else
                call MPI_ISSend( l_col(ind+1) , ncol(io), &
                     MPI_Integer, 0, gio, MPI_Comm_World, MPIreq, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Integer, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: io_write_Sp')
             call MPI_Get_Count(MPIstatus, MPI_Integer, io, MPIerror)
             write(iu) buf(1:io)
          end if
       end do
       
       if ( Node == 0 ) then
          deallocate(buf)
       else if ( MPIerror /= MPI_Success - 1 ) then
          ! Wait for the last one to not send
          ! to messages with the same tag...
          call MPI_Wait(MPIreq,MPIstatus,MPIerror)
       end if
    else

       ind = 0
       do io = 1 , no
          write(iu) l_col(ind+1:ind+ncol(io))
          ind = ind + ncol(io)
       end do
       
    end if

#else
    
    ind = 0
    do io = 1 , no
       write(iu) l_col(ind+1:ind+ncol(io))
       ind = ind + ncol(io)
    end do
#endif
    
  end subroutine io_write_Sp

  subroutine io_read_d1D(iu, sp, dSp1D, tag, dit, Bcast)

    use class_dSpData1D

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Data array with sparsity pattern
    type(dSpData1D), intent(inout) :: dSp1D
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(inout), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast

    ! Local variables for reading the values
    type(OrbitalDistribution) :: fdit
    real(dp), pointer :: a(:) => null()
    integer, pointer :: ncol(:) => null()

    integer :: io, lno, no, ind, n_nzs
    logical :: ldit, lBcast
    integer, pointer :: ibuf(:) => null()
#ifdef MPI
    real(dp), pointer :: buf(:) => null()
    integer :: max_n, gio
    integer :: MPIerror, BNode, MPIreq, MPIstatus(MPI_STATUS_SIZE)
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)
    ldit = .not. lno == no

    if ( ldit ) then
       call newdSpData1D(sp,dit,dSp1D,name=trim(tag))
       
       if ( lno == no ) then
          ibuf => ncol
       else
#ifdef MPI
          allocate(ibuf(no))

          do gio = 1 , no
             BNode = node_handling_element(dit,gio)
             
             if ( Node == BNode ) then
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
                   ibuf(gio) = ncol(io)
                else
                   call MPI_ISSend( ncol(io) , 1, MPI_Integer, &
                        0, gio, MPI_Comm_World, MPIreq, MPIerror)
                end if
             else if ( Node == 0 ) then
                call MPI_Recv( ibuf(gio) , 1, MPI_Integer, &
                     BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             end if
          end do
#else
          call die('Error in distribution, io_read_d1D')
#endif
       end if

    else
       ! Create the Fake distribution
#ifdef MPI
       call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
       call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
       call newdSpData1D(sp,fdit,dSp1D,name=trim(tag))
    end if

    a => val(dSp1D)

    if ( ldit ) then
#ifdef MPI
       
       ! Allocate the maximum number of entries
       if ( IONode ) then
          max_n = maxval(ibuf)
          allocate(buf(max_n))
       end if
       MPIerror = MPI_Success - 1

       ! Loop size
       ind = 0
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                read(iu) a(ind+1:ind+ncol(io))
             else
                call MPI_IRecv( a(ind+1) , ncol(io), MPI_Double_Precision, &
                     0, gio, MPI_Comm_World, MPIreq, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             read(iu) buf(1:ibuf(gio))
             call MPI_Send( buf(1) , ibuf(gio) , MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
          end if
       end do

       if ( Node /= 0 .and. MPIerror /= MPI_Success - 1 ) then
          call MPI_Wait(MPIreq,MPIstatus,MPIerror)
       end if
#else
       call die('Error in distribution for, io_read_d1D')
#endif
    else

       if ( IONode ) then
          ind = 0
          do io = 1 , no
             read(iu) a(ind+1:ind+ncol(io))
             ind = ind + ncol(io)
          end do
       end if

    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then
       call MPI_Bcast(a(1),n_nzs,MPI_Double_Precision, &
            0, MPI_Comm_World, MPIError)
    
    else if ( .not. IONode ) then
       return ! the sparsity pattern will not be created then...
    end if
#endif

  end subroutine io_read_d1D

  subroutine io_write_d1D(iu, dSp1D)

    use class_dSpData1D

    ! File handle
    integer, intent(in) :: iu
    ! Data array with sparsity pattern (and an attached distribution)
    type(dSpData1D), intent(inout) :: dSp1D

    ! Local variables for reading the values
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: a(:) => null()
    integer, pointer :: ncol(:) => null()

    integer :: io, lno, no, ind
    logical :: ldit
#ifdef MPI
    real(dp), pointer :: buf(:) => null()
    integer :: gio, max_n
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol)

    ldit = .not. lno == no

    ! Retrieve data
    a => val(dSp1D)

    if ( ldit ) then
       
#ifdef MPI
       ! Get the maximum size of the array
       io = maxval(ncol)
       call MPI_Reduce(io,max_n,1,MPI_Integer, &
            MPI_Max, 0, MPI_Comm_World, MPIerror)

       ! The ionode now has the maximum retrieved array
       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if
       MPIerror = MPI_Success - 1

       ! Loop size
       ind = 0
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                write(iu) a(ind+1:ind+ncol(io))
             else
                call MPI_ISSend( a(ind+1) , ncol(io), &
                     MPI_Double_Precision, 0, gio, MPI_Comm_World, MPIreq, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: io_write_dSp1D')
             call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
             write(iu) buf(1:io)
          end if
       end do
       
       if ( Node == 0 ) then
          deallocate(buf)
       else if ( MPIerror /= MPI_Success - 1 ) then
          ! Wait for the last one to not send
          ! to messages with the same tag...
          call MPI_Wait(MPIreq,MPIstatus,MPIerror)
       end if
#else
       call die('Error in io_write_d1D')
#endif
    else

       ind = 0
       do io = 1 , no
          write(iu) a(ind+1:ind+ncol(io))
          ind = ind + ncol(io)
       end do
       
    end if
    
  end subroutine io_write_d1D

  subroutine io_read_d2D(iu, sp, dSp2D, dim2, tag, &
       sparsity_dim, dit, Bcast)
    
    use class_dSpData2D

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Data array with sparsity pattern
    type(dSpData2D), intent(inout) :: dSp2D
    ! The non-sparse dimension
    integer, intent(in) :: dim2
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! This denotes the sparsity dimension (either 1 or 2)
    integer, intent(in), optional :: sparsity_dim
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(inout), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast

    ! Local variables for reading the values
    type(OrbitalDistribution) :: fdit
    real(dp), pointer :: a(:,:) => null()
    integer, pointer :: ncol(:) => null(), ibuf(:) => null()

    integer :: io, lno, no, s, ind, n_nzs
    integer :: sp_dim
    logical :: ldit, lBcast
#ifdef MPI
    real(dp), pointer :: buf(:) => null()
    integer :: gio, max_n
    integer :: MPIerror, BNode, MPIreq, MPIstatus(MPI_STATUS_SIZE)
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    sp_dim = 1
    if ( present(sparsity_dim) ) sp_dim = sparsity_dim

    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,nnzs=n_nzs)
    ldit = .not. lno == no

    if ( ldit ) then
       call newdSpData2D(sp,dim2,dit,dSp2D,name=trim(tag), &
            sparsity_dim=sp_dim)

       if ( lno == no ) then
          ibuf => ncol
       else
#ifdef MPI
          allocate(ibuf(no))

          do gio = 1 , no
             BNode = node_handling_element(dit,gio)
             
             if ( Node == BNode ) then
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
                   ibuf(gio) = ncol(io)
                else
                   call MPI_ISSend( ncol(io) , 1, MPI_Integer, &
                        0, gio, MPI_Comm_World, MPIreq, MPIerror)
                end if
             else if ( Node == 0 ) then
                call MPI_Recv( ibuf(gio) , 1, MPI_Integer, &
                     BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             end if
          end do
#else
          call die('Error in distribution, io_read_d2D')
#endif
       end if

    else
       ! Create the Fake distribution
#ifdef MPI
       call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
       call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
       call newdSpData2D(sp,dim2,fdit,dSp2D,name=trim(tag), &
            sparsity_dim=sp_dim)
    end if

    a => val(dSp2D)

    if ( ldit ) then

#ifdef MPI
       ! Allocate maximum number of entries
       if ( IONode ) then
          max_n = maxval(ibuf)
          allocate(buf(max_n*dim2))
       end if
       MPIerror = MPI_Success - 1

    if ( sp_dim == 2 ) then ! collapsed IO
       ! Loop size
       ind = 0
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                read(iu) a(1:dim2,ind+1:ind+ncol(io))
             else
                call MPI_IRecv( a(1,ind+1) , dim2*ncol(io), MPI_Double_Precision, &
                     0, gio, MPI_Comm_World, MPIreq, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             read(iu) buf(1:ibuf(gio)*dim2)
             call MPI_Send( buf(1) , ibuf(gio)*dim2 , MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
          end if
       end do

       if ( Node /= 0 .and. MPIerror /= MPI_Success - 1 ) then
          call MPI_Wait(MPIreq,MPIstatus,MPIerror)
       end if

    else ! non-collapsed IO

       do s = 1 , dim2
          ind = 0
          do gio = 1 , no
             BNode = node_handling_element(dit,gio)
             
             if ( Node == BNode ) then
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
                   read(iu) a(ind+1:ind+ncol(io),s)
                else
                   call MPI_IRecv( a(ind+1,s) , ncol(io), MPI_Double_Precision, &
                        0, gio, MPI_Comm_World, MPIreq, MPIerror)
                end if
                ind = ind + ncol(io)
             else if ( Node == 0 ) then
                read(iu) buf(1:ibuf(gio))
                call MPI_Send( buf(1) , ibuf(gio) , MPI_Double_Precision, &
                     BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             end if
          end do

          if ( Node /= 0 .and. MPIerror /= MPI_Success - 1 ) then
             call MPI_Wait(MPIreq,MPIstatus,MPIerror)
          end if

       end do

    end if

#else
       call die('Error in distribution for, io_read_d2D')
#endif
    else
       
       if ( IONode ) then
          if ( sp_dim == 2 ) then ! collapsed IO
             ind = 0
             do io = 1 , no
                read(iu) (a(1:dim2,ind+s),s=1,ncol(io))
                ind = ind + ncol(io)
             end do
          else ! non-collapsed IO
             do s = 1 , dim2 
                ind = 0
                do io = 1 , no
                   read(iu) a(ind+1:ind+ncol(io),s)
                   ind = ind + ncol(io)
                end do
             end do
          end if
       end if
    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then
       call MPI_Bcast(a(1,1),dim2*n_nzs,MPI_Double_Precision, &
            0, MPI_Comm_World, MPIError)
    else if ( .not. IONode ) then
       return ! the sparsity pattern will not be created then...
    end if
#endif

  end subroutine io_read_d2D

  subroutine io_write_d2D(iu, dSp2D)

    use class_dSpData2D

    ! File handle
    integer, intent(in) :: iu
    ! Data array with sparsity pattern (and an attached distribution)
    type(dSpData2D), intent(inout) :: dSp2D

    ! Local variables for reading the values
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: a(:,:) => null()
    integer, pointer :: ncol(:) => null()

    integer :: io, lno, no, ind
    integer :: id2, dim2, sp_dim, n_nzs
    logical :: ldit
#ifdef MPI
    real(dp), pointer :: buf(:) => null()
    integer :: gio, max_n
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)
    
    ldit = .not. lno == no
    
    ! Retrieve data
    a => val(dSp2D)
    if ( size(a,dim=2) == n_nzs ) then
       sp_dim = 2
       dim2 = size(a,dim=1)
    else
       sp_dim = 1
       dim2 = size(a,dim=2)
    end if

    if ( ldit ) then
       
#ifdef MPI
       ! Get the maximum size of the array
       io = maxval(ncol)
       call MPI_Reduce(io,max_n,1,MPI_Integer, &
            MPI_Max, 0, MPI_Comm_World, MPIerror)

       ! The ionode now has the maximum retrieved array
       if ( Node == 0 ) then
          allocate(buf(max_n*dim2))
       end if
       MPIerror = MPI_Success - 1

       ! Loop size
    if ( sp_dim == 1 ) then

       do id2 = 1 , dim2
          ind = 0
          do gio = 1 , no
             BNode = node_handling_element(dit,gio)
             
             if ( Node == BNode ) then
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
                   write(iu) a(ind+1:ind+ncol(io),id2)
                else
                   call MPI_ISSend( a(ind+1,id2) , ncol(io), &
                        MPI_Double_Precision, 0, gio, MPI_Comm_World, MPIreq, MPIerror)
                end if
                ind = ind + ncol(io)
             else if ( Node == 0 ) then
                call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                     BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
                if ( MPIerror /= MPI_Success ) &
                     call die('Error in code (1): io_write_dSp2D')
                call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
                write(iu) buf(1:io)
             end if
          end do ! gio
       end do ! id2

    else
       
       ind = 0
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                write(iu) a(1:dim2,ind+1:ind+ncol(io))
             else
                call MPI_ISSend( a(1,ind+1) , dim2*ncol(io), &
                     MPI_Double_Precision, 0, gio, MPI_Comm_World, MPIreq, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code (2): io_write_dSp2D')
             call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
             write(iu) buf(1:io)
          end if
       end do

    end if
    
       if ( Node == 0 ) then
          deallocate(buf)
       else if ( MPIerror /= MPI_Success - 1 ) then
          ! Wait for the last one to not send
          ! to messages with the same tag...
          call MPI_Wait(MPIreq,MPIstatus,MPIerror)
       end if
#else
       call die('Error in io_write_d2D')
#endif
    else
       
       if ( sp_dim == 1 ) then
          do id2 = 1 , dim2
             ind = 0
             do io = 1 , no
                write(iu) a(ind+1:ind+ncol(io),id2)
                ind = ind + ncol(io)
             end do
          end do
       else
          ind = 0
          do io = 1 , no
             write(iu) a(1:dim2,ind+1:ind+ncol(io))
             ind = ind + ncol(io)
          end do
       end if
       
    end if
    
  end subroutine io_write_d2D

end module m_io_s
  
