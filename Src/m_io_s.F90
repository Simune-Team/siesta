! Module for easy reading of sparsity patterns and 
! data distributed in sparsity patterns.

! Ideally all IO routines should utilize these routines

! Fully implemented by Nick Papior Andersen

! SIESTA io-module
module m_io_s

  use alloc
  use class_OrbitalDistribution
  use class_Sparsity
  use precision, only : sp, dp
  use parallel, only : IONode, Node, Nodes
#ifdef MPI
  use mpi_siesta, only : MPI_Bcast, MPI_AllReduce, MPI_Sum
  use mpi_siesta, only : MPI_Comm_World, MPI_Comm_Self
  use mpi_siesta, only : MPI_Integer, MPI_Double_Precision
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
  subroutine io_read_Sp(iu, no, sp, name, dit, Bcast)

    ! File handle
    integer, intent(in) :: iu
    ! Number of orbitals readed
    integer, intent(in) :: no
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The name of the sparsity pattern
    character(len=*), intent(in) :: name
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(inout), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast

    ! Local variables for reading the values
    integer, pointer :: ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    integer :: io, n_nzs, ind
    logical :: ldit, lBcast
#ifdef MPI
    integer :: MPIerror
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    ! Allocate arrays
    call re_alloc(ncol,1,no)
    call re_alloc(l_ptr,1,no)

    if ( IONode ) then

       read(iu) ncol
       n_nzs = sum(ncol)
       call re_alloc(l_col,1,n_nzs)
       ind = 0
       do io = 1 , no
          read(iu) l_col(ind+1:ind+ncol(io))
          ind = ind + ncol(io)
       end do

    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( ldit ) then
       call die('Not implemented yet!')
    else if ( lBcast ) then
       ! Bcast everything
       call MPI_Bcast(ncol(1),no,MPI_Integer,0,MPI_Comm_World, &
            MPIError)
       l_ptr(1) = 0
       do io = 2 , no
          l_ptr(io) = l_ptr(io-1) + ncol(io-1)
       end do
       n_nzs = l_ptr(no) + ncol(no)
       if ( .not. IONode ) call re_alloc(l_col,1,n_nzs)
       call MPI_Bcast(l_col(1),n_nzs,MPI_Integer,0,MPI_Comm_World, &
            MPIError)
    else if ( .not. IONode ) then
       return ! the sparsity pattern will not be created then...
    end if
#else
    l_ptr(1) = 0
    do io = 2 , no
       l_ptr(io) = l_ptr(io-1) + ncol(io-1)
    end do
#endif

    ! Create the sparsity pattern
    call newSparsity(sp,no,no, &
         n_nzs, ncol, l_ptr, l_col, trim(name))

    call de_alloc(ncol)
    call de_alloc(l_ptr)
    call de_alloc(l_col)
    
  end subroutine io_read_Sp

  subroutine io_read_d1D(iu, sp, dSp1D, name, dit, Bcast)

    use class_dSpData1D

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Data array with sparsity pattern
    type(dSpData1D), intent(inout) :: dSp1D
    ! The name of the sparsity pattern
    character(len=*), intent(in) :: name
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
#ifdef MPI
    integer :: MPIerror
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)

    if ( ldit ) then
       call newdSpData1D(sp,dit,dSp1D,name=trim(name))
       if ( lno /= no ) then
          ! We have a problem... We haven't made
          ! this routine distribute the data yet...
          call die('Not implemented yet!')
       end if
    else
       ! Create the Fake distribution
#ifdef MPI
       call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
       call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
       call newdSpData1D(sp,fdit,dSp1D,name=trim(name))
    end if

    a => val(dSp1D)

    if ( IONode ) then
       ind = 0
       do io = 1 , no
          read(iu) a(ind+1:ind+ncol(io))
          ind = ind + ncol(io)
       end do
    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( ldit .or. lBcast ) then
       call MPI_Bcast(a(1),n_nzs,MPI_Double_Precision,0,MPI_Comm_World, &
            MPIError)
    
    else if ( .not. IONode ) then
       return ! the sparsity pattern will not be created then...
    end if
#endif

  end subroutine io_read_d1D

  subroutine io_read_d2D(iu, sp, dSp2D, dim2, name, &
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
    ! The name of the sparsity pattern
    character(len=*), intent(in) :: name
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
    integer, pointer :: ncol(:) => null()

    integer :: io, lno, no, s, ind, n_nzs
    integer :: sp_dim
    logical :: ldit, lBcast
#ifdef MPI
    integer :: MPIerror
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    sp_dim = 1
    if ( present(sparsity_dim) ) sp_dim = sparsity_dim

    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,nnzs=n_nzs)

    if ( ldit ) then
       call newdSpData2D(sp,dim2,dit,dSp2D,name=trim(name), &
            sparsity_dim=sp_dim)
       if ( lno /= no ) then
          ! We have a problem... We haven't made
          ! this routine distribute the data yet...
          call die('Not implemented yet!')
       end if
    else
       ! Create the Fake distribution
#ifdef MPI
       call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
       call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
       call newdSpData2D(sp,dim2,fdit,dSp2D,name=trim(name), &
            sparsity_dim=sp_dim)
    end if

    a => val(dSp2D)

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

    ! If a distribution is present, then do something
#ifdef MPI
    if ( ldit .or. lBcast ) then
       call MPI_Bcast(a(1,1),dim2*n_nzs,MPI_Double_Precision,0,MPI_Comm_World, &
            MPIError)
    else if ( .not. IONode ) then
       return ! the sparsity pattern will not be created then...
    end if
#endif

  end subroutine io_read_d2D

end module m_io_s
  
