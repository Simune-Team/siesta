! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module m_ncdf_io

  use precision, only : sp, dp
  use parallel, only : Node, Nodes
#ifdef NCDF_4
  use variable
  use dictionary
  use nf_ncdf, ncdf_parallel => parallel
#endif
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  private

#ifdef NCDF_4

  public :: cdf_w_sp
  public :: cdf_w_d1D, cdf_w_d2D

#endif

contains

#ifdef NCDF_4

  subroutine cdf_w_Sp(ncdf,dit,sp)
    use m_io_s,only : Node_Sp_gncol
    use class_OrbitalDistribution
    use class_Sparsity
  
    type(hNCDF), intent(inout) :: ncdf
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp

    ! Local variables
    integer, pointer :: ncol(:), l_col(:)
    integer, allocatable :: buf(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n
#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif
    integer :: n_nnzs, g_nzs

    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

#ifdef MPI
    if ( no_l /= no_u ) then
       call MPI_Reduce(n_nzs,g_nzs,1,MPI_Integer, &
            MPI_Sum, 0, MPI_Comm_World, MPIerror)
    else
       g_nzs = n_nzs
    end if
#else
    g_nzs = nnzs(sp)
#endif

    ! Read dimension nnzs
    call ncdf_inq_dim(ncdf,'nnzs',len=n_nnzs)
    if ( Node == 0 .and. n_nnzs /= g_nzs ) then
       call die('Number of non-zero elements is not equivalent.')
    end if
    
    if ( parallel_io(ncdf) ) then

       print *,'Here:',node,' this isformatted erroneously'

       ! Calculate the offset for the global
       ! partitioning
       gind = global_offset(dit,no_l)

       ! we write it using MPI
       call ncdf_par_access(ncdf,name='n_col',access=NF90_COLLECTIVE)
       call ncdf_put_var(ncdf,'n_col',ncol,start=(/gind+1/), &
            count=(/no_l/))

       gind = global_offset(dit,n_nnzs)

       ! we write it using MPI
       call ncdf_par_access(ncdf,name='list_col',access=NF90_COLLECTIVE)
       call ncdf_put_var(ncdf,'list_col',l_col,start=(/gind+1/), &
            count=(/n_nzs/))

    else if ( no_l == no_u ) then

       call ncdf_put_var(ncdf,'n_col',ncol)
       call ncdf_put_var(ncdf,'list_col',l_col)

    else

#ifdef MPI
       allocate(buf(no_u))
       call Node_Sp_gncol(0,sp,dit,no_u,buf)
       call ncdf_put_var(ncdf,'n_col',buf)
#else
       call ncdf_put_var(ncdf,'n_col',ncol)
#endif

       if ( Node == 0 ) then
#ifdef MPI
          max_n = maxval(buf)
          deallocate(buf)
#else
          max_n = maxval(ncol)
#endif
          allocate(buf(max_n))
       end if

       ! Write list_col
#ifdef MPI

       MPIreq = MPI_REQUEST_NULL

       ! Loop size
       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)

          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                call ncdf_put_var(ncdf,'list_col',l_col(ind+1:ind+ncol(io)), &
                     count=(/ncol(io)/),start=(/gind/))
                gind = gind + ncol(io)
             else
                call MPI_Send( l_col(ind+1) , ncol(io), MPI_Integer, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Integer, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: ncdf_write_Sp')
             call MPI_Get_Count(MPIstatus, MPI_Integer, io, MPIerror)
             call ncdf_put_var(ncdf,'list_col',buf(1:io), &
                  count=(/io/),start=(/gind/))
             gind = gind + io
          end if
       end do

       if ( Node == 0 ) then
          deallocate(buf)
       end if

#else
       call ncdf_put_var(ncdf,'list_col',l_col)
#endif

    end if

  end subroutine cdf_w_Sp

  subroutine cdf_w_d1D(ncdf,vname,dSp1D)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData1D), intent(inout) :: dSp1D

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n
    real(dp), pointer :: a(:)
    real(dp), allocatable :: buf(:)

#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

    a => val(dSp1D)
    
    if ( parallel_io(ncdf) ) then

       call ncdf_par_access(ncdf,name=vname,access=NF90_COLLECTIVE)

       ! Calculate the processors offset
       gind = global_offset(dit,n_nzs)

       ! we write it using MPI
       call ncdf_put_var(ncdf,vname,a,start=(/gind+1/), &
            count=(/n_nzs/))

    else

#ifdef MPI
       io = maxval(ncol)
       call MPI_Reduce(io,max_n,1,MPI_Integer, &
            MPI_Max, 0, MPI_Comm_World, MPIerror)
       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if

       MPIreq = MPI_REQUEST_NULL

       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                call ncdf_put_var(ncdf,vname,a(ind+1:ind+ncol(io)), &
                     count=(/ncol(io)/),start=(/gind/))
                gind = gind + ncol(io)
             else
                call MPI_Send( a(ind+1) , ncol(io), MPI_Double_Precision, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: ncdf_write_d1D')
             call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
             call ncdf_put_var(ncdf,vname,buf(1:io), &
                  count=(/io/),start=(/gind/))
             gind = gind + io
          end if
       end do

#else
       call ncdf_put_var(ncdf,trim(vname),a)
#endif

    end if
    
  end subroutine cdf_w_d1D

  subroutine cdf_w_d2D(ncdf,vname,dSp2D)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData2D), intent(inout) :: dSp2D

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n
    integer :: id2, dim2, sp_dim
    real(dp), pointer :: a(:,:)
    real(dp), allocatable :: buf(:)

#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

    a => val(dSp2D)

    if ( size(a,dim=2) == n_nzs ) then
       sp_dim = 2
       dim2 = size(a,dim=1)
    else
       sp_dim = 1
       dim2 = size(a,dim=2)
    end if

    if ( parallel_io(ncdf) ) then

       call ncdf_par_access(ncdf,name=vname,access=NF90_COLLECTIVE)

       ! we write it using MPI
       gind = global_offset(dit,n_nzs)

       if ( sp_dim == 1 ) then
          do io = 1 , dim2
             call ncdf_put_var(ncdf,trim(vname),a(:,io), &
                  start=(/gind+1,io/), count=(/n_nzs/))
          end do
       else
          call ncdf_put_var(ncdf,trim(vname),a,start=(/1,gind+1/), &
               count=(/dim2,n_nzs/))
       end if

    else

#ifdef MPI
       io = maxval(ncol)
       call MPI_Reduce(io,max_n,1,MPI_Integer, &
            MPI_Max, 0, MPI_Comm_World, MPIerror)
       max_n = max_n * dim2
       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if

       MPIreq = MPI_REQUEST_NULL

    if ( sp_dim == 1 ) then

       do id2 = 1 , dim2
          ind = 0
          gind = 1
          do gio = 1 , no_u
             BNode = node_handling_element(dit,gio)
             
             if ( Node == BNode ) then
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
                   call ncdf_put_var(ncdf,trim(vname),a(ind+1:ind+ncol(io),id2), &
                        count=(/ncol(io)/),start=(/gind,id2/))
                   gind = gind + ncol(io)
                else
                   call MPI_Send( a(ind+1,id2) , ncol(io), MPI_Double_Precision, &
                        0, gio, MPI_Comm_World, MPIerror)
                   
                end if
                ind = ind + ncol(io)
             else if ( Node == 0 ) then
                call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                     BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
                if ( MPIerror /= MPI_Success ) &
                     call die('Error in code (1): ncdf_write_d2D')
                call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
                call ncdf_put_var(ncdf,trim(vname),buf(1:io), &
                     count=(/io/),start=(/gind,id2/))
                gind = gind + io
             end if
          end do ! gio
       end do ! id2

    else

       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                call ncdf_put_var(ncdf,trim(vname),a(1:dim2,ind+1:ind+ncol(io)), &
                     count=(/dim2,ncol(io)/),start=(/1,gind/))
                gind = gind + ncol(io)
             else
                call MPI_Send( a(1,ind+1) , dim2*ncol(io), MPI_Double_Precision, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: ncdf_write_d2D')
             call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
             io = io / dim2
             call ncdf_put_var(ncdf,trim(vname),reshape(buf(1:io*dim2),(/dim2,io/)), &
                  count=(/dim2,io/),start=(/1,gind/))
             gind = gind + io
          end if
       end do

    end if
#else
       call ncdf_put_var(ncdf,trim(vname),a)
#endif

    end if
    
  end subroutine cdf_w_d2D
  
#else
  subroutine dummy()
    !pass
  end subroutine dummy
#endif

end module m_ncdf_io
