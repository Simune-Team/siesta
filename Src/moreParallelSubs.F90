! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2008.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
!******************************************************************************
! MODULE moreParallelSubs
! Provides some utility routines for parallel execution
! Written by J.M.Soler. Feb.2008
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! miscAllReduce : Reduces a miscellaneous set of variables and arrays
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use alloc,      only: de_alloc        ! De-allocation routine
! use alloc,      only: re_alloc        ! Re-allocation routine
! use sys,        only: die             ! Termination routine
!
!   USED module parameters:
! use precision,  only: dp              ! Real double precision type
!
!   USED MPI procedures, interfaces and types:
! use mpi_siesta, only: MPI_AllReduce
! use mpi_siesta, only: MPI_COMM_WORLD
! use mpi_siesta, only: MPI_Integer
! use mpi_siesta, only: MPI_double_precision
!
!   EXTERNAL procedures used:
! none
!
!******************************************************************************
!
! SUBROUTINE miscAllReduce( op, a0, b0, c0, d0, e0, f0, &
!                               a1, b1, c1, &
!                               a2, b2, &
!                               a3 )
!
! Reduces a miscellaneous set of variables and arrays
!------------------------------ INPUT -----------------------------------------
! character(len=*) op  ! Operation to be applied ('sum'|'prod'|'max'|'min')
!-------------------- OPTIONAL INPUT and OUTPUT -------------------------------
! type                 :: a0, b0, c0, d0, e0, f0  ! Scalar arguments
! type,dimension(:)    :: a1, b1, c1              ! Vector arguments
! type,dimension(:,:)  :: a2, b2                  ! Rank-2 array arguments
! type,dimension(:,:,:):: a3                      ! Rank-3 array arguments
!   Where type=(integer|double precision)
!   All arguments must be of the same type
!----------------------------- UNITS ------------------------------------------
! Units of all arguments are arbitrary
!----------------------------- USAGE ------------------------------------------
! use precision,        only: dp
! use moreParallelSubs, only: miscAllReduce
! real(dp):: Ek, Eh, polarization(3), stress(3,3)
! real(dp),allocatable:: force(:,:)
! ... Find Ek, Eh, polarization, and stress. Allocate and find force.
! call miscAllReduce( 'sum', Ek, Eh, a1=polarization, a2=force, b2=stress )
!---------------------------- BEHAVIOUR ---------------------------------------
! In serial execution, it does nothing.
! If op/=('sum'|'prod'|'max'|'min'), it stops with an error message.
!--------------------------- ALGORITHMS ---------------------------------------
! All present variables and arrays are packed consecutively in a buffer vector
! The buffer is reduced with MPI_AllReduce and unpacked back into the present
!   variables and arrays.
!
!******************************************************************************

MODULE moreParallelSubs

! Used module procedures
  use alloc,     only: de_alloc  ! De-allocation routine
  use alloc,     only: re_alloc  ! Re-allocation routine
  use sys,       only: die       ! Termination routine

! Used module parameters
  use precision, only: dp        ! Real double precision type

! MPI interfaces and types
#ifdef MPI
  use mpi_siesta
#endif

! All public procedures (there are no public types, parameters, or variables):
PUBLIC:: &
  miscAllReduce  ! Reduces a miscellaneous set of variables and arrays

PRIVATE ! Nothing is declared public beyond this point

! Overload procedure name
! Note: a compiler warning may occur because calling miscAllReduce without any
!       optional arguments would not allow to resolve the right module proced.
!       But such a call would make no sense and it will never happen.
  interface miscAllReduce
    module procedure      &
!      miscAllReduceInt,   &! Integer version
      miscAllReduceDouble  ! Double precision real version
  end interface miscAllReduce

CONTAINS

!******************************************************************************

SUBROUTINE miscAllReduceInt( op, a0, b0, c0, d0, e0, f0, &
                                 a1, b1, c1, &
                                 a2, b2, &
                                 a3 )

  implicit none
  character(len=*),intent(in)    :: op  ! ('sum'|'prod'|'max'|'min')
  integer,optional,intent(inout)                 :: a0, b0, c0, d0, e0, f0
  integer,optional,intent(inout),dimension(:)    :: a1, b1, c1
  integer,optional,intent(inout),dimension(:,:)  :: a2, b2
  integer,optional,intent(inout),dimension(:,:,:):: a3

  character(len=*),parameter:: myName  = 'miscAllReduceInt '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: MPIerror, m, n
  integer,pointer:: recvBuff(:)=>null(), sendBuff(:)=>null()

! Do nothing unless execution is parallel
#ifdef MPI

! Find total size of variables and arrays to be gathered
  n = 0
  if (present(a0)) n = n + 1
  if (present(b0)) n = n + 1
  if (present(c0)) n = n + 1
  if (present(d0)) n = n + 1
  if (present(e0)) n = n + 1
  if (present(f0)) n = n + 1
  if (present(a1)) n = n + size(a1)
  if (present(b1)) n = n + size(b1)
  if (present(c1)) n = n + size(c1)
  if (present(a2)) n = n + size(a2)
  if (present(b2)) n = n + size(b2)
  if (present(a3)) n = n + size(a3)

! Allocate buffers to send and receive
  call re_alloc( sendBuff, 1,n, name=myName//'sendBuff' )
  call re_alloc( recvBuff, 1,n, name=myName//'recvBuff' )

! Pack all variables and arrays into sendBuff
  n = 0
  if (present(a0)) then
    sendBuff(n+1) = a0
    n = n + 1
  end if
  if (present(b0)) then
    sendBuff(n+1) = b0
    n = n + 1
  end if
  if (present(c0)) then
    sendBuff(n+1) = c0
    n = n + 1
  end if
  if (present(d0)) then
    sendBuff(n+1) = d0
    n = n + 1
  end if
  if (present(e0)) then
    sendBuff(n+1) = e0
    n = n + 1
  end if
  if (present(f0)) then
    sendBuff(n+1) = f0
    n = n + 1
  end if
  if (present(a1)) then
    m = size(a1)
    sendBuff(n+1:n+m) = a1
    n = n + m
  end if
  if (present(b1)) then
    m = size(b1)
    sendBuff(n+1:n+m) = b1
    n = n + m
  end if
  if (present(c1)) then
    m = size(c1)
    sendBuff(n+1:n+m) = c1
    n = n + m
  end if
  if (present(a2)) then
    m = size(a2)
    sendBuff(n+1:n+m) = reshape( a2, (/m/) )
    n = n + m
  end if
  if (present(b2)) then
    m = size(b2)
    sendBuff(n+1:n+m) = reshape( b2, (/m/) )
    n = n + m
  end if
  if (present(a3)) then
    m = size(a3)
    sendBuff(n+1:n+m) = reshape( a3, (/m/) )
    n = n + m
  end if

! Reduce sendBuff into recvBuff
  if (op=='sum' .or. op=='Sum' .or. op=='SUM') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_Integer, &
                        MPI_Sum, MPI_Comm_World, MPIerror )
  else if (op=='prod' .or. op=='Prod' .or. op=='PROD') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_Integer, &
                        MPI_Prod, MPI_Comm_World, MPIerror )
  else if (op=='max' .or. op=='Max' .or. op=='MAX') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_Integer, &
                        MPI_Max, MPI_Comm_World, MPIerror )
  else if (op=='min' .or. op=='Min' .or. op=='MIN') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_Integer, &
                        MPI_Min, MPI_Comm_World, MPIerror )
  else
    call die(errHead//'unknown operator')
  end if

! Unpack recvBuff
  n = 0
  if (present(a0)) then
    a0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(b0)) then
    b0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(c0)) then
    c0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(d0)) then
    d0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(e0)) then
    e0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(f0)) then
    f0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(a1)) then
    m = size(a1)
    a1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(b1)) then
    m = size(b1)
    b1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(c1)) then
    m = size(c1)
    c1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(a2)) then
    m = size(a2)
    a2 = reshape( recvBuff(n+1:n+m), shape(a2) )
    n = n + m
  end if
  if (present(b2)) then
    m = size(b2)
    b2 = reshape( recvBuff(n+1:n+m), shape(b2) )
    n = n + m
  end if
  if (present(a3)) then
    m = size(a3)
    a3 = reshape( recvBuff(n+1:n+m), shape(a3) )
    n = n + m
  end if

  call de_alloc( recvBuff )
  call de_alloc( sendBuff )

#endif

END SUBROUTINE miscAllReduceInt

!******************************************************************************

SUBROUTINE miscAllReduceDouble( op, a0, b0, c0, d0, e0, f0, &
                                    a1, b1, c1, &
                                    a2, b2, &
                                    a3 )

  implicit none
  character(len=*),intent(in)    :: op  ! ('sum'|'prod'|'max'|'min')
  real(dp),optional,intent(inout)                 :: a0, b0, c0, d0, e0, f0
  real(dp),optional,intent(inout),dimension(:)    :: a1, b1, c1
  real(dp),optional,intent(inout),dimension(:,:)  :: a2, b2
  real(dp),optional,intent(inout),dimension(:,:,:):: a3

  character(len=*),parameter:: myName  = 'miscAllReduceInt '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: MPIerror, m, n
  real(dp),pointer:: recvBuff(:)=>null(), sendBuff(:)=>null()

! Do nothing unless execution is parallel
#ifdef MPI

! Find total size of variables and arrays to be reduced
  n = 0
  if (present(a0)) n = n + 1
  if (present(b0)) n = n + 1
  if (present(c0)) n = n + 1
  if (present(d0)) n = n + 1
  if (present(e0)) n = n + 1
  if (present(f0)) n = n + 1
  if (present(a1)) n = n + size(a1)
  if (present(b1)) n = n + size(b1)
  if (present(c1)) n = n + size(c1)
  if (present(a2)) n = n + size(a2)
  if (present(b2)) n = n + size(b2)
  if (present(a3)) n = n + size(a3)

! Allocate buffers to send and receive
  call re_alloc( sendBuff, 1,n, name=myName//'sendBuff' )
  call re_alloc( recvBuff, 1,n, name=myName//'recvBuff' )

! Pack all variables and arrays into sendBuff
  n = 0
  if (present(a0)) then
    sendBuff(n+1) = a0
    n = n + 1
  end if
  if (present(b0)) then
    sendBuff(n+1) = b0
    n = n + 1
  end if
  if (present(c0)) then
    sendBuff(n+1) = c0
    n = n + 1
  end if
  if (present(d0)) then
    sendBuff(n+1) = d0
    n = n + 1
  end if
  if (present(e0)) then
    sendBuff(n+1) = e0
    n = n + 1
  end if
  if (present(f0)) then
    sendBuff(n+1) = f0
    n = n + 1
  end if
  if (present(a1)) then
    m = size(a1)
    sendBuff(n+1:n+m) = a1
    n = n + m
  end if
  if (present(b1)) then
    m = size(b1)
    sendBuff(n+1:n+m) = b1
    n = n + m
  end if
  if (present(c1)) then
    m = size(c1)
    sendBuff(n+1:n+m) = c1
    n = n + m
  end if
  if (present(a2)) then
    m = size(a2)
    sendBuff(n+1:n+m) = reshape( a2, (/m/) )
    n = n + m
  end if
  if (present(b2)) then
    m = size(b2)
    sendBuff(n+1:n+m) = reshape( b2, (/m/) )
    n = n + m
  end if
  if (present(a3)) then
    m = size(a3)
    sendBuff(n+1:n+m) = reshape( a3, (/m/) )
    n = n + m
  end if

! Reduce sendBuff into recvBuff
  if (op=='sum' .or. op=='Sum' .or. op=='SUM') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_double_precision, &
                        MPI_Sum, MPI_Comm_World, MPIerror )
  else if (op=='prod' .or. op=='Prod' .or. op=='PROD') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_double_precision, &
                        MPI_Prod, MPI_Comm_World, MPIerror )
  else if (op=='max' .or. op=='Max' .or. op=='MAX') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_double_precision, &
                        MPI_Max, MPI_Comm_World, MPIerror )
  else if (op=='min' .or. op=='Min' .or. op=='MIN') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_double_precision, &
                        MPI_Min, MPI_Comm_World, MPIerror )
  else
    call die(errHead//'unknown operator')
  end if

! Unpack recvBuff
  n = 0
  if (present(a0)) then
    a0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(b0)) then
    b0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(c0)) then
    c0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(d0)) then
    d0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(e0)) then
    e0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(f0)) then
    f0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(a1)) then
    m = size(a1)
    a1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(b1)) then
    m = size(b1)
    b1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(c1)) then
    m = size(c1)
    c1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(a2)) then
    m = size(a2)
    a2 = reshape( recvBuff(n+1:n+m), shape(a2) )
    n = n + m
  end if
  if (present(b2)) then
    m = size(b2)
    b2 = reshape( recvBuff(n+1:n+m), shape(b2) )
    n = n + m
  end if
  if (present(a3)) then
    m = size(a3)
    a3 = reshape( recvBuff(n+1:n+m), shape(a3) )
    n = n + m
  end if

  call de_alloc( recvBuff )
  call de_alloc( sendBuff )

#endif

END SUBROUTINE miscAllReduceDouble

!******************************************************************************

END MODULE moreParallelSubs
