! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module m_ncdf_io

  use precision, only : sp, dp, grid_p
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
  public :: cdf_init_mesh
  public :: cdf_w_grid

  public :: cdf_w_sp
  public :: cdf_w_d1D, cdf_w_d2D

#endif

  ! Create the type to hold the mesh distribution data
  type :: tMeshDist
     private
     ! number of mesh divisions in each axis
     integer :: m(3)
     ! mesh box bounds of each node in each direction
            ! box(1,iAxis,iNode)=lower bounds
            ! box(2,iAxis,iNode)=upper bounds
     integer, allocatable :: box(:,:,:)
  end type tMeshDist
  type(tMeshDist), save :: distr

contains

#ifdef NCDF_4

  subroutine cdf_init_mesh(mesh,nsm)
    
    use parallel, ProcY => ProcessorY

    ! mesh divisions (fine points), number of fine points per big points
    integer, intent(in) :: mesh(3), nsm

    ! Local quantities
    integer :: nm(3), nmesh, ntot
    integer :: ProcZ, blocY, blocZ, nremY, nremZ
    integer :: dimX, dimY, dimZ
    integer :: PP, iniY, iniZ, PY, PZ

    if ( .not. allocated(distr%box) ) then
       allocate(distr%box(2,3,0:Nodes-1))
    end if
    
    nm(1:3) = mesh(1:3) / nsm
    nmesh = product(mesh)
    distr%m(1:3) = mesh(1:3)
    
    ProcZ = Nodes/ProcY
    
    blocY = nm(2)/ProcY
    nremY = nm(2) - blocY*ProcY
    blocZ = nm(3)/ProcZ
    nremZ = nm(3) - blocZ*ProcZ

    dimX = nm(1) * nsm
    
    ntot = 0
    
    PP   = 0
    iniY = 1
    do PY = 1, ProcY
       
       dimY = blocY
       if ( PY <= nremY ) dimY = dimY + 1  ! Add extra points starting from the first nodes
       dimY = dimY * nsm                 ! For fine points
       
        iniZ = 1
        do PZ = 1, ProcZ
           dimZ = blocZ
           if ( PZ <= nremZ ) dimZ = dimZ + 1
           dimZ = dimZ*nsm                 ! For fine points
           
           distr%box(1,1,PP) = 1
           distr%box(2,1,PP) = dimX
           distr%box(1,2,PP) = iniY
           distr%box(2,2,PP) = iniY + dimY - 1
           distr%box(1,3,PP) = iniZ
           distr%box(2,3,PP) = iniZ + dimZ - 1

           ntot = ntot + dimX * dimY * dimZ
           
           iniZ = iniZ + dimZ
           PP   = PP + 1
           
        end do

        iniY = iniY + dimY

     end do
     
     if (ntot /= nmesh) then
        if (Node == 0) then
           write(6,*) "Nominal npt: ", nmesh, " /= assigned npt:", ntot
        end if
        call die()
     end if
     
  end subroutine cdf_init_mesh

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
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
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
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
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
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
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
  
  subroutine cdf_w_grid(ncdf,name,mesh,lnpt,grid,idx)

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: name
    integer, intent(in) :: mesh(3), lnpt
    real(grid_p), intent(in) :: grid(lnpt)
    integer, intent(in), optional :: idx

#ifdef MPI
    integer :: MPIstat(MPI_STATUS_SIZE)
    integer :: MPIerror, mnpt
    integer :: lb(3), nel(3), iN, inpt
    real(grid_p), allocatable :: gb(:)
#endif

#ifdef MPI

    if ( parallel_io(ncdf) ) then

       ! Ensure collective writing
       call ncdf_par_access(ncdf,name=name,access=NF90_COLLECTIVE)

       lb(:)  = distr%box(1,:,Node)
       nel(:) = distr%box(2,:,Node) - lb(:) + 1

       if ( present(idx) ) then
          call ncdf_put_var(ncdf,name,grid, &
               start=(/lb(1),lb(2),lb(3),idx/), &
               count=(/nel(1),nel(2),nel(3),1/) )          
       else
          call ncdf_put_var(ncdf,name,grid, start=lb, count=nel )
       end if
   
    else

       mnpt = 0
       do iN = 0 , Nodes - 1
          nel(:) = distr%box(2,:,iN) - distr%box(1,:,iN) + 1
          mnpt = max(mnpt,product(nel))
       end do
    
       ! The main node can safely write the data...
       if ( Node == 0 ) then
          
          allocate(gb(mnpt))
          
          ! First save it's own data
          lb(:) = distr%box(1,:,0) 
          nel(:) = distr%box(2,:,0) - distr%box(1,:,0) + 1
          
          if ( present(idx) ) then
             call ncdf_put_var(ncdf,name,grid, &
                  start=(/lb(1),lb(2),lb(3),idx/), &
                  count=(/nel(1),nel(2),nel(3),1/) )
          else
             call ncdf_put_var(ncdf,name,grid, &
                  start=lb, count=nel )
          end if

          ! Loop on the remaining nodes
          do iN = 1 , Nodes - 1

             ! we retrieve data from the iN'th node.
             call MPI_Recv(gb,mnpt,MPI_grid_real,iN,iN, &
                  MPI_Comm_World,MPIstat, MPIerror)
             
             ! Just make sure we only pass the correct size
             call MPI_Get_Count(MPIstat, MPI_Grid_Real, inpt, MPIerror)
             
             lb(:) = distr%box(1,:,iN) 
             nel(:) = distr%box(2,:,iN) - lb(:) + 1
             if ( inpt /= product(nel) ) then
                call die('Error when receiving the distributed &
                     &grid data for writing to the NetCDF file.')
             end if
             
             if ( present(idx) ) then
                call ncdf_put_var(ncdf,name,gb(1:inpt), &
                     start=(/lb(1),lb(2),lb(3),idx/), &
                     count=(/nel(1),nel(2),nel(3),1/) )
             else
                call ncdf_put_var(ncdf,name,gb(1:inpt), &
                     start=lb, count=nel )
             end if
          
          end do
          
          deallocate(gb)
       else
          call MPI_Send(grid,lnpt,MPI_grid_real,0, &
               Node,MPI_Comm_World,MPIerror)
       end if
       
    end if

#else
    if ( present(idx) ) then
       call ncdf_put_var(ncdf,name,grid,start=(/1,1,1,idx/))
    else
       call ncdf_put_var(ncdf,name,grid)
    end if
#endif
  
  end subroutine cdf_w_grid

#else
  subroutine dummy()
    !pass
  end subroutine dummy
#endif

end module m_ncdf_io
