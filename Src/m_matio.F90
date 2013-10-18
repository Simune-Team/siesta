!     
!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996- .
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     
      module m_matio

      implicit none

      private

      PUBLIC :: write_mat

CONTAINS

!------------------------------------------------------
  subroutine write_mat (maxnd, no_l, nspin, &
       numd, listdptr, listd, mat, userfile)

    use mpi
    use parallel, only: mpi_comm => SIESTA_comm, blocksize
    use alloc

    integer, parameter :: dp = selected_real_kind(10,100)


    integer, intent(in) :: maxnd
    integer, intent(in) :: no_l
    integer, intent(in) :: nspin
    integer, intent(in) :: numd(1:no_l)
    integer, intent(in) :: listdptr(1:no_l)
    integer, intent(in) :: listd(maxnd)
    real(dp), intent(in) :: mat(maxnd,nspin)

    character(len=*), intent(in) :: userfile

    integer :: no_u, m, ml, im, ndmaxg, unit1, is
    integer :: n_l, nsize_l, n_g, nsize_g, node, myrank, nprocs
#ifdef MPI
    integer  :: MPIerror, stat(MPI_STATUS_SIZE)
    real(dp), dimension(:), pointer :: buffer => null()
    integer,  dimension(:), pointer :: ibuffer => null()
#endif
    integer, dimension(:), pointer  :: numdg => null()

    call timer("WriteMat",1)

#ifdef MPI
    call MPI_Comm_Size( MPI_Comm, nprocs, MPIerror )
    call MPI_Comm_Rank( MPI_Comm, myrank, MPIerror )
#else
    nprocs = 1
    myrank = 0
#endif

!     Find total number of orbitals over all Nodes
!     *** Do we want *all*reduce?
#ifdef MPI
    call MPI_AllReduce(no_l,no_u,1,MPI_integer,MPI_sum,MPI_Comm,MPIerror)
#else
    no_u = no_l
#endif

    if (myrank.eq.0) then
       print *, "Blocksize, no_u: ", blocksize, no_u

       call io_assign(unit1)
       open( unit1, file=trim(userfile), form="formatted", status='unknown' )
       rewind(unit1)
       write(unit1,*) no_u, nspin, blocksize
       
       call re_alloc( numdg, 1, no_u, 'numdg', 'write_mat' )
       print *, "Size of numdg: ", size(numdg)
    endif

!     Get info about numd
    n_g = 0
    n_l = 0
    node = -1
    DO
!!       call mpi_barrier(mpi_comm, mpierror)
       node = node + 1
       if (node == nprocs) node = 0

       print *, " node: ", node, " myrank: ", myrank

       if (myrank == node) then
          nsize_l = min(blocksize,no_l-n_l)
          print *, "myrank: ", myrank, " processing size:", nsize_l
          if (node==0) then
             nsize_g = min(blocksize,no_u-n_g)
             print *, "myrank: ", myrank, " will just copy: (l,g)", &
                      nsize_l, nsize_g
             numdg(n_g+1:n_g+nsize_g) = numd(n_l+1:n_l+nsize_l)
             n_g = n_g + nsize_g
             print *, "root has received so far: ", n_g
          else
             print *, "myrank: ", myrank, " will send to 0: ", nsize_l
             call MPI_Send(numd(n_l+1),nsize_l,MPI_integer, &
                  0,1,MPI_Comm,MPIerror)
             print *, "myrank: ", myrank, " completed send ", nsize_l
          endif
          n_l = n_l + nsize_l

       else if (myrank == 0) then
          nsize_g = min(blocksize,no_u-n_g)
          print *, "root will receive from ", node, " size: ", nsize_g
          print *, "will put it starting at: ", n_g+1
!          call sub(numdg(n_g+1:),nsize_g)
          call MPI_Recv(numdg(n_g+1:),nsize_g,MPI_integer, &
                node,1,MPI_Comm,stat,MPIerror)
          n_g = n_g + nsize_g
          print *, "root has received so far: ", n_g
       endif

       if (myrank == 0) then
          if (n_g == no_u) EXIT
       else
          if (n_l == no_l) then
             print *, "rank ", myrank, " exiting loop"
             EXIT
          endif
       endif
          
    enddo
          
!     Write out numd array
      if (myrank.eq.0) then
!         ndmaxg = maxval(numdg(1:no_u))
         write(unit1,*) (numdg(m),m=1,no_u)
         call io_close(unit1)
!!$#ifdef MPI
!!$         nullify(buffer,ibuffer)
!!$         call re_alloc( buffer,  1, ndmaxg, 'buffer',  'write_dm' )
!!$         call re_alloc( ibuffer, 1, ndmaxg, 'ibuffer', 'write_dm' )
!!$#endif
      endif

    CONTAINS
    subroutine sub(a,n)
      integer, intent(in) :: a(*)
      integer :: n
      print *, a(1:n)
    end subroutine sub

    end subroutine write_mat

  end module m_matio
