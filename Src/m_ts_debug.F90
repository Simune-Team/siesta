module m_ts_debug

  use precision, only :dp

  implicit none

  private :: dp

#ifdef TRANSIESTA_DEBUG
contains

  subroutine write_TriMat(iu,tri)
    use class_zTriMat
    integer, intent(inout) :: iu
    type(zTriMat), intent(inout) :: tri
    complex(dp), pointer :: z(:)
    integer :: i,j, p, np, n,idx, ntmp

    iu = iu + 1
    ! Number of parts
    np = parts(tri)
    n = 0
    
    do p = 1 , np - 1

       z => val(tri,p,p)
       do j = 1 , nrows_g(tri,p)
          idx = (j-1)*nrows_g(tri,p)
          do i = 1 , nrows_g(tri,p)
             write(iu,'(2(tr1,i5),2(tr1,g10.6))') n+i, n+j,z(idx+i)
          end do
       end do

       z => val(tri,p,p+1)
       ntmp = n+nrows_g(tri,p)
       do j = 1 , nrows_g(tri,p)
          idx = (j-1)*nrows_g(tri,p+1)
          do i = 1 , nrows_g(tri,p+1)
             write(iu,'(2(tr1,i5),2(tr1,g10.6))') n+ntmp+i, n+j,z(idx+i)
          end do
       end do

       z => val(tri,p+1,p)
       ntmp = n+nrows_g(tri,p)
       do j = 1 , nrows_g(tri,p+1)
          idx = (j-1)*nrows_g(tri,p)
          do i = 1 , nrows_g(tri,p)
             write(iu,'(2(tr1,i5),2(tr1,g10.6))') n+i, n+ntmp+j,z(idx+i)
          end do
       end do
       
       n = n + nrows_g(tri,p)
    end do

    z => val(tri,np,np)
    do j = 1 , nrows_g(tri,np)
       idx = (j-1)*nrows_g(tri,np)
       do i = 1 , nrows_g(tri,np)
          write(iu,'(2(tr1,i5),2(tr1,g10.6))') n+i, n+j,z(idx+i)
       end do
    end do
  end subroutine write_TriMat

  subroutine write_Full(iu,no,GF)
    integer, intent(inout) :: iu
    integer, intent(in) :: no
    complex(dp), intent(inout) :: GF(no,no)
    integer :: i,j, p, np, n,idx
    iu = iu + 1
    
    do j = 1 , no
       do i = 1 , no
          write(iu,'(2(tr1,i5),2(tr1,g10.6))') i, j,GF(i,j)
       end do
    end do

  end subroutine write_Full

  
  subroutine sp_to_file(u,sp)
    use class_Sparsity
    use geom_helper, only : UCORB
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
#endif
    integer, intent(in) :: u
    type(Sparsity), intent(inout) :: sp
    integer :: io,jo,j,ind
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)

    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    write(u,'(i5)') nrows(sp)

    do io = 1 , nrows(sp)
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)
          ind = l_ptr(io) + j
          jo = UCORB(l_col(ind),nrows(sp))
          write(u,'(3(tr1,i5))') io,jo,1
       end do
    end do
    if ( ind /= nnzs(sp) ) then
       call die('Have not looped through all things')
    end if

#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,io)
#endif

  end subroutine sp_to_file

#endif

end module m_ts_debug
