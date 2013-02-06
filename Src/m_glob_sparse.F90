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
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
module m_glob_sparse
!
! Routines that are used for the globalization and distribution of the Hamiltonian 
! and scattering matrices into a full size sparse matrix instead of distributed sparse
!  matrices. It has several options for creating different kinds of matrices.
!
! call glob_sparse_arrays(no_l,no_u,no_s,maxnh, &
!       numh ,listhptr ,listh ,xij , Gamma, &
!       maxnhg, &
!       numhg,listhptrg,listhg,xijg)
! is used to generate the fully globalized versions of those with suffix "g".
! This glob_sparse_arrays is needed before using glob_sparse_matrix.
!
! call glob_sparse_matrix(no_l,no_u,no_s, &
!       maxnh,  numh , listhptr , H(:,ispin) , &
!       maxnhg, numhg, listhptrg, Hf)
! is used to create distribute the sparse matrix H to their full non-distributed
! sparse matrix representation.
! 
! Together with m_hs_matrix this module provides an easy interface to generate 
! any k-point matrix in full form during MPI runs.

! This module also creates generic routines for "sparse array" allocation.
! Meaning that one can create the globalized numh array, then the listhptr, etc.
! This will leverage some memory problems when large systems are used.
! Especially it is justified when performing IO-operations.
!
!
  implicit none

  private

#ifdef MPI
  public :: glob_sparse_numh
  public :: glob_sparse_listh
  public :: glob_sparse_listhptr
  public :: glob_sparse_xij

  public :: glob_sparse_arrays
  public :: glob_sparse_arrays_dealloc
  public :: glob_sparse_matrix
  public :: glob_sparse_matrix_dealloc
  public :: glob_sparse_matrices
  public :: glob_sparse_matrices_dealloc
#endif

contains

#ifdef MPI

! Globalization of 'numh' array
  subroutine glob_sparse_numh(no_l,no_u,numh,numhg)

    use parallel,     only : Node,Nodes
    use mpi_siesta,   only : MPI_Comm_World,MPI_Integer
    use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_l              ! no. orbs. in unit cell (local)
    integer, intent(in) :: no_u              ! no. orbs. in unit cell (global)
    integer, intent(in) :: numh(no_l)        ! Number of nonzero elements of each row
! ******************************
! * OUTPUT variables           *
! * The globalized equivalents *
! ******************************
    integer, allocatable,  intent(out) :: numhg(:)
        
! ************************
! * LOCAL variables      *
! ************************
    integer :: io, iio
    integer :: BNode
    integer :: MPIerror

    allocate(numhg(no_u))
    call memory('A','I',no_u,'globArrays')
    do io = 1,no_u
       call WhichNodeOrb(io,Nodes,BNode)
       if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
       end if
       call MPI_Bcast(numhg(io),1,MPI_Integer,BNode,MPI_Comm_World, &
            MPIError)
    end do
  end subroutine glob_sparse_numh

! Globalization of 'listhptr' array
  subroutine glob_sparse_listhptr(no_u,numhg,listhptrg)

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_u              ! no. orbs. in unit cell (global)
    integer, intent(in) :: numhg(no_u)       ! Number of nonzero elements of each row, ALREADY GLOBALIZED
! ******************************
! * OUTPUT variables           *
! * The globalized equivalents *
! ******************************
    integer, allocatable,  intent(out) :: listhptrg(:)
        
! ************************
! * LOCAL variables      *
! ************************
    integer :: io

    allocate(listhptrg(no_u))
    call memory('A','I',no_u,'globArrays')
    listhptrg(1) = 0
    do io = 2,no_u
       listhptrg(io) = listhptrg(io-1) + numhg(io-1)
    end do

  end subroutine glob_sparse_listhptr

! Globalize the 'listh' array 
  subroutine glob_sparse_listh(no_l,no_u, maxnh, &
       numh ,listhptr , listh, &
       numhg,listhptrg, maxnhg, listhg)

    use parallel,     only : Node,Nodes
    use mpi_siesta,   only : MPI_Comm_World,MPI_Integer
    use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_l              ! no. orbs. in unit cell (local)
    integer, intent(in) :: no_u              ! no. orbs. in unit cell (global)
    integer, intent(in) :: maxnh             ! Maximum number of nonzero elements of 
!                                              each row of hamiltonian matrix
    integer, intent(in) :: numh(no_l)        ! Number of nonzero elements of each row
    integer, intent(in) :: listhptr(no_l)    ! Pointer to each row (-1) of the
!                                              hamiltonian matrix
!                                              column indexes for each matrix row
    integer, intent(in) :: listh(maxnh)      ! Nonzero hamiltonian-matrix element
!                                              column indexes for each matrix row
!                                              of hamiltonian matrix
    integer, intent(in) :: numhg(no_u)       ! The globalized equivalent...
    integer, intent(in) :: listhptrg(no_u)   ! The globalized equivalent...
! ******************************
! * OUTPUT variables           *
! * The globalized equivalents *
! ******************************
    integer, intent(out)               :: maxnhg
    integer, allocatable,  intent(out) :: listhg(:)
        
! ************************
! * LOCAL variables      *
! ************************
    integer :: io, jo, iio
    integer :: BNode
    integer :: MPIerror

! Globalize maxnh
    maxnhg = listhptrg(no_u) + numhg(no_u)
    
    allocate(listhg(maxnhg))
    call memory('A','I',maxnhg,'globArrays')
    do io = 1 , no_u
       call WhichNodeOrb(io,Nodes,BNode)
       if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
               listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
       end if
       call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_Integer, &
            BNode,MPI_Comm_World,MPIError)
    end do

  end subroutine glob_sparse_listh

! Globalize the 'xij' array 
  subroutine glob_sparse_xij(no_l,no_u, maxnh, &
       numh ,listhptr, xij, Gamma, &
       numhg,listhptrg,maxnhg, xijg)

    use precision,    only : dp
    use parallel,     only : Node,Nodes
    use mpi_siesta,   only : MPI_Comm_World,DAT_Double
    use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_l              ! no. orbs. in unit cell (local)
    integer, intent(in) :: no_u              ! no. orbs. in unit cell (global)
    integer, intent(in) :: maxnh             ! Maximum number of nonzero elements of 
!                                              each row of hamiltonian matrix
    integer, intent(in) :: numh(no_l)        ! Number of nonzero elements of each row
    integer, intent(in) :: listhptr(no_l)    ! Pointer to each row (-1) of the
!                                              hamiltonian matrix
!                                              column indexes for each matrix row
    real(dp), intent(in):: xij(3,maxnh)      ! Vectors between orbital centers (sparse)
    logical, intent(in) :: Gamma             ! A gamma calculation? Does xij exist?
    integer, intent(in) :: numhg(no_u)       ! The globalized equivalent...
    integer, intent(in) :: listhptrg(no_u)   ! The globalized equivalent...
! ******************************
! * OUTPUT variables           *
! * The globalized equivalents *
! ******************************
    integer, intent(out)               :: maxnhg
    real(dp), allocatable, intent(out) :: xijg(:,:)
        
! ************************
! * LOCAL variables      *
! ************************
    integer :: io, jo, iio
    integer :: BNode
    integer :: MPIerror

! Globalize maxnh
    maxnhg = listhptrg(no_u) + numhg(no_u)

! Globalize xij
    if ( .not. Gamma ) then
       allocate(xijg(3,maxnhg))
       call memory('A','D',3*maxnhg,'globArrays')
       do io = 1 , no_u
          call WhichNodeOrb(io,Nodes,BNode)
          if ( Node .eq. BNode ) then
             call GlobalToLocalOrb(io,Node,Nodes,iio)
             do jo = 1 , numh(iio)
                xijg(1:3,listhptrg(io)+jo) = xij(1:3,listhptr(iio)+jo)
             end do
          end if
          call MPI_Bcast(xijg(1,listhptrg(io)+1),3*numhg(io) &
               ,DAT_double,BNode,MPI_Comm_World,MPIerror)
       end do
    end if

  end subroutine glob_sparse_xij


  ! Routine for globalizing lists within the unit cell.
  ! These are the arrays needed when using the glob_HS.
  ! This routine will globalize every array needed for constructing the full
  ! sparse matrices
  subroutine glob_sparse_arrays(no_l,no_u,no_s,maxnh, &
       numh ,listhptr ,listh ,xij , Gamma,&
       numhg,listhptrg,maxnhg,listhg,xijg)

    use precision,    only : dp  

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_l              ! no. orbs. in unit cell (local)
    integer, intent(in) :: no_u              ! no. orbs. in unit cell (global)
    integer, intent(in) :: no_s              ! no. orbs. in supercell
    integer, intent(in) :: maxnh             ! Maximum number of nonzero elements of 
!                                              each row of hamiltonian matrix
    integer, intent(in) :: numh(no_l)        ! Number of nonzero elements of each row
    integer, intent(in) :: listhptr(no_l)    ! Pointer to each row (-1) of the
!                                              hamiltonian matrix
!                                              column indexes for each matrix row
    integer, intent(in) :: listh(maxnh)      ! Nonzero hamiltonian-matrix element
!                                              column indexes for each matrix row
!                                              of hamiltonian matrix
    real(dp), intent(in) :: xij(3,maxnh)     ! Vectors between orbital centers (sparse)
    logical, intent(in)  :: Gamma            ! Is is a Gamma calculation
! ******************************
! * OUTPUT variables           *
! * The globalized equivalents *
! ******************************
    integer, allocatable,  intent(out) :: numhg(:)
    integer, allocatable,  intent(out) :: listhptrg(:)
    integer, intent(out)               :: maxnhg
    integer, allocatable,  intent(out) :: listhg(:)
    real(dp), allocatable, intent(out) :: xijg(:,:)
        
    call glob_sparse_numh(no_l,no_u,numh,numhg)

    call glob_sparse_listhptr(no_u,numhg,listhptrg)

    call glob_sparse_listh(no_l,no_u, maxnh, &
         numh ,listhptr , listh, &
         numhg,listhptrg, maxnhg, listhg)

    call glob_sparse_xij(no_l,no_u, maxnh, &
         numh ,listhptr , xij, Gamma, &
         numhg,listhptrg, maxnhg, xijg)
    
  end subroutine glob_sparse_arrays


  ! Routine for deallocating all auxillary sparse lists
  ! You can call this routine even with already deallocated arrays,
  ! if they are allocated they will be de-allocated.
  subroutine glob_sparse_arrays_dealloc(no_u, Gamma, &
       maxnhg, &
       numhg,listhptrg,listhg,xijg)

    use precision,    only : dp  

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_u   ! no. orbs. in unit cell (global)
    logical, intent(in) :: Gamma  ! Is is a Gamma calculation
    integer, intent(in) :: maxnhg ! Number of elements in the full hamiltonian
! *******************************
! * OUTPUT variables            *
! * The deallocated equivalents *
! *******************************
    integer, allocatable,  intent(in out) :: numhg(:)
    integer, allocatable,  intent(in out) :: listhptrg(:)
    integer, allocatable,  intent(in out) :: listhg(:)
    real(dp), allocatable, intent(in out) :: xijg(:,:)
        
    if ( allocated(numhg) ) then
       call memory('D','I',no_u,'globArrays')
       deallocate(numhg)
    end if
    
    if ( allocated(listhptrg) ) then
       call memory('D','I',no_u,'globArrays')
       deallocate(listhptrg)
    end if
    
    if ( allocated(listhg) ) then
       call memory('D','I',maxnhg,'globArrays')
       deallocate(listhg)
    end if
    
    if ( .not. Gamma .and. allocated(xijg) ) then
       call memory('D','D',3*maxnhg,'globArrays')
       deallocate(xijg)
    end if
    
  end subroutine glob_sparse_arrays_dealloc


  subroutine glob_sparse_matrices(no_l,no_u,no_s, &
       maxnh,  numh , listhptr , H ,S , &
       maxnhg, numhg, listhptrg, Hg,Sg)

    use parallel,     only : Node,Nodes
    use precision,    only : dp  
    use mpi_siesta,   only : MPI_Comm_World,DAT_double
    use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_l              ! no. orbs. in unit cell (local)
    integer, intent(in) :: no_u              ! no. orbs. in unit cell (global)
    integer, intent(in) :: no_s              ! no. orbs. in supercell
    integer, intent(in) :: maxnh             ! Maximum number of nonzero elements of 
!                                              each row of hamiltonian matrix
    integer, intent(in) :: numh(no_l)        ! Number of nonzero elements of each row
    integer, intent(in) :: listhptr(no_l)    ! Pointer to each row (-1) of the
    real(dp), intent(in):: H(maxnh)          ! The local Hamiltonian in sparse format
    real(dp), intent(in):: S(maxnh)          ! The local overlap in sparse format
    integer, intent(in) :: maxnhg            ! Sum reduction of maxnh
    integer, intent(in) :: numhg(no_u)       ! Number of nonzero elements of each row
    integer, intent(in) :: listhptrg(no_u)   ! Pointer to each row (-1) of the

! ******************************
! * OUTPUT variables           *
! ******************************
    real(dp), intent(out) :: Hg(maxnhg), Sg(maxnhg)
        
! ************************
! * LOCAL variables      *
! ************************
    integer :: io, jo, iio
    integer :: BNode
    integer :: MPIerror
    
    do io = 1 , no_u
       call WhichNodeOrb(io,Nodes,BNode)
       if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numh(iio)
             Hg(listhptrg(io)+jo) = H(listhptr(iio)+jo)
             Sg(listhptrg(io)+jo) = S(listhptr(iio)+jo)
          end do
       endif
       call MPI_Bcast(Hg(listhptrg(io)+1),numhg(io),DAT_double,BNode,MPI_Comm_World,MPIerror)
       call MPI_Bcast(Sg(listhptrg(io)+1),numhg(io),DAT_double,BNode,MPI_Comm_World,MPIerror)
    end do
    
  end subroutine glob_sparse_matrices

  ! Routine for deallocating the full sparse matrices H and S
  subroutine glob_sparse_matrices_dealloc(maxnhg,Hg,Sg)

    use precision,    only : dp  

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: maxnhg ! Number of elements in the full hamiltonian

! *******************************
! * OUTPUT variables            *
! * The deallocated equivalents *
! *******************************
    real(dp), allocatable, intent(in out) :: Hg(:), Sg(:)
        
    if ( allocated(Hg) ) then
       call memory('D','D',maxnhg,'globArrays')
       deallocate(Hg)
    end if

    if ( allocated(Sg) ) then
       call memory('D','D',maxnhg,'globArrays')
       deallocate(Sg)
    end if
    
  end subroutine glob_sparse_matrices_dealloc


  subroutine glob_sparse_matrix(no_l,no_u,no_s, &
       maxnh,  numh , listhptr , H  , &
       maxnhg, numhg, listhptrg, Hg)

    use parallel,     only : Node,Nodes
    use precision,    only : dp  
    use mpi_siesta,   only : MPI_Comm_World,DAT_double
    use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_l              ! no. orbs. in unit cell (local)
    integer, intent(in) :: no_u              ! no. orbs. in unit cell (global)
    integer, intent(in) :: no_s              ! no. orbs. in supercell
    integer, intent(in) :: maxnh             ! Maximum number of nonzero elements of 
!                                              each row of hamiltonian matrix
    integer, intent(in) :: numh(no_l)        ! Number of nonzero elements of each row
    integer, intent(in) :: listhptr(no_l)    ! Pointer to each row (-1) of the
    real(dp), intent(in):: H(maxnh)          ! The local Hamiltonian in sparse format
    integer, intent(in) :: maxnhg            ! Sum reduction of maxnh
    integer, intent(in) :: numhg(no_u)       ! Number of nonzero elements of each row
    integer, intent(in) :: listhptrg(no_u)   ! Pointer to each row (-1) of the

! ******************************
! * OUTPUT variables           *
! ******************************
    real(dp), intent(out) :: Hg(maxnhg)
        
! ************************
! * LOCAL variables      *
! ************************
    integer :: io, jo, iio
    integer :: BNode
    integer :: MPIerror
    
    do io = 1 , no_u
       call WhichNodeOrb(io,Nodes,BNode)
       if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numh(iio)
             Hg(listhptrg(io)+jo) = H(listhptr(iio)+jo)
          end do
       endif
       call MPI_Bcast(Hg(listhptrg(io)+1),numhg(io),DAT_double,BNode,MPI_Comm_World,MPIerror)
    end do
    
  end subroutine glob_sparse_matrix

! Routine for deallocating one full sparse matrices 
  subroutine glob_sparse_matrix_dealloc(maxnhg,Hg)

    use precision,    only : dp  

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: maxnhg ! Number of elements in the full hamiltonian

! *******************************
! * OUTPUT variables            *
! * The deallocated equivalents *
! *******************************
    real(dp), allocatable, intent(in out) :: Hg(:)
        
    if ( allocated(Hg) ) then
       call memory('D','D',maxnhg,'globArrays')
       deallocate(Hg)
    end if
    
  end subroutine glob_sparse_matrix_dealloc

#else
  subroutine dummy()
  end subroutine dummy
#endif


end module m_glob_sparse
  
