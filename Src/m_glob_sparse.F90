!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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
!       numh ,listh ,listhptr ,xij , Gamma, &
!       maxnhg, &
!       numhg,listhg,listhptrg,xijg)
! is used to generate the fully globalized versions of those with suffix "g".
! This glob_sparse_arrays is needed before using glob_sparse_matrix.
!
! call glob_sparse_matrix(no_l,no_u,no_s, &
!       maxnh,  numh , listhptr , H ,S , &
!       maxnhg, numhg, listhptrg, Hf,Sf)
! is used to create distribute the sparse matrix H and S to their full non-distributed
! sparse matrix representation.
! 
! Together with m_hs_matrix this module provides an easy interface to generate 
! any k-point matrix in full form during MPI runs.
!
!
  implicit none

  private

#ifdef MPI
  public :: glob_sparse_arrays
  public :: glob_sparse_arrays_dealloc
  public :: glob_sparse_matrix
#endif

contains


#ifdef MPI
  ! Routine for globalizing lists within the unit cell.
  ! These are the arrays needed when using the glob_HS.
  subroutine glob_sparse_arrays(no_l,no_u,no_s,maxnh, &
       numh ,listh ,listhptr ,xij , Gamma,&
       maxnhg, &
       numhg,listhg,listhptrg,xijg)

    use parallel,     only : Node,Nodes
    use precision,    only : dp  
    use mpi_siesta,   only : MPI_Comm_World,MPI_Integer
    use mpi_siesta,   only : DAT_double => MPI_double_precision
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
    integer, intent(in) :: listh(maxnh)      ! Nonzero hamiltonian-matrix element
!                                              column indexes for each matrix row
!                                              of hamiltonian matrix
    integer, intent(in) :: listhptr(no_l)    ! Pointer to each row (-1) of the
!                                              hamiltonian matrix
!                                              column indexes for each matrix row
    real(dp), intent(in) :: xij(3,maxnh)     ! Vectors between orbital centers (sparse)
    logical, intent(in)  :: Gamma            ! Is is a Gamma calculation
! ******************************
! * OUTPUT variables           *
! * The globalized equivalents *
! ******************************
    integer, intent(out)               :: maxnhg
    integer, allocatable,  intent(out) :: numhg(:)
    integer, allocatable,  intent(out) :: listhg(:)
    integer, allocatable,  intent(out) :: listhptrg(:)
    real(dp), allocatable, intent(out) :: xijg(:,:)
        
! ************************
! * LOCAL variables      *
! ************************
    integer :: io, jo, iio
    integer :: BNode
    integer :: MPIerror

! Globalize numh
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
         
! Globalize listhptr
    allocate(listhptrg(no_u))
    call memory('A','I',no_u,'globArrays')
    listhptrg(1) = 0
    do io = 2 , no_u
       listhptrg(io) = listhptrg(io-1) + numhg(io-1)
    end do
    
! Globalize maxnh
    maxnhg = listhptrg(no_u) + numhg(no_u)
    
! Globalize listh
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

  end subroutine glob_sparse_arrays


  ! Routine for globalizing lists within the unit cell.
  ! These are the arrays needed when using the glob_HS.
  subroutine glob_sparse_arrays_dealloc(no_u, Gamma, &
       maxnhg, &
       numhg,listhg,listhptrg,xijg)

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
    integer, allocatable,  intent(in out) :: listhg(:)
    integer, allocatable,  intent(in out) :: listhptrg(:)
    real(dp), allocatable, intent(in out) :: xijg(:,:)
        
    call memory('D','I',no_u,'globArrays')
    deallocate(numhg)
    
    call memory('D','I',no_u,'globArrays')
    deallocate(listhptrg)
    
    call memory('D','I',maxnhg,'globArrays')
    deallocate(listhg)
    
    if ( .not. Gamma ) then
       call memory('D','D',3*maxnhg,'globArrays')
       deallocate(xijg)
    end if
    
  end subroutine glob_sparse_arrays_dealloc

  subroutine glob_sparse_matrix(no_l,no_u,no_s, &
       maxnh,  numh , listhptr , H ,S , &
       maxnhg, numhg, listhptrg, Hf,Sf)

    use parallel,     only : Node,Nodes
    use precision,    only : dp  
    use mpi_siesta,   only : MPI_Comm_World
    use mpi_siesta,   only : DAT_double => MPI_double_precision
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
    real(dp), intent(out) :: Hf(maxnhg), Sf(maxnhg)
        
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
             Hf(listhptrg(io)+jo) = H(listhptr(iio)+jo)
             Sf(listhptrg(io)+jo) = S(listhptr(iio)+jo)
          end do
       endif
       call MPI_Bcast(Hf(listhptrg(io)+1),numhg(io),DAT_double,BNode,MPI_Comm_World,MPIerror)
       call MPI_Bcast(Sf(listhptrg(io)+1),numhg(io),DAT_double,BNode,MPI_Comm_World,MPIerror)
    end do
    
  end subroutine glob_sparse_matrix

#else
  subroutine dummy()
  end subroutine dummy
#endif


end module m_glob_sparse
  
