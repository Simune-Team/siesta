subroutine setup_ordern_indexes(no_l,no_u,Nodes)

  use spatial, only:  nL2G, nG2L, nNode, nOrbPerNode
  use domain_decom, only :  use_dd_perm, dd_nuo, dd_perm, &
                            ulimit, llimit, dd_invp, dd_nnode
  use alloc, only: re_alloc

#ifdef MPI
  use mpi_siesta, only: MPI_Comm_World, MPI_Integer
  use mpi_siesta, only: MPI_AllGather, MPI_AllGatherV
#endif

  implicit none

  integer, intent(in) :: no_l, no_u
  integer, intent(in) :: Nodes

#ifdef MPI
  integer, dimension(:), allocatable :: nl2gtmp, counts, displs
  integer :: i, MPIerr
#endif
  integer :: LOrb, GOrb

  call re_alloc( nL2G, 1, no_u, 1, Nodes, 'nL2G', 'setup_ordern_indexes' )
  call re_alloc( nG2L, 1, no_u, 'nG2L', 'setup_ordern_indexes' )
  call re_alloc( nNode, 1, no_u, 'nNode', 'setup_ordern_indexes' )
  call re_alloc( nOrbPerNode, 1, Nodes, 'nOrbPerNode','setup_ordern_indexes' )

! Global to local index and node mapping
  do Gorb = 1, no_u
     if (use_dd_perm) then
        Lorb = dd_perm(GOrb)
     else
        if (GOrb.ge.llimit .and. GOrb.lt.ulimit) then
           LOrb = GOrb - llimit + 1
        else
           LOrb = 0
        endif
     endif
     nG2L(Gorb) = LOrb
     nNode(Gorb) = dd_nnode(Gorb)
  enddo

! Local to global index
! (Local node only)

#ifdef MPI
  allocate(nl2gtmp(1:no_u))
  nl2gtmp(1:no_u) = 0
#endif

  do LOrb = 1, no_l
     if (use_dd_perm) then
        GOrb = dd_invp(LOrb)
     else
        GOrb = LOrb + llimit - 1
     endif
     !     nL2G(LOrb,Node+1) = Gorb
#ifdef MPI
     nl2gtmp(LOrb) = Gorb
#else
     nL2G(LOrb,1) = Gorb
#endif
  enddo

#ifdef MPI
  ! Now all_gatherv
  allocate(counts(1:Nodes),displs(1:Nodes))
  do i = 1, Nodes
     counts(i) = no_u
     displs(i) = (i-1)*no_u
  enddo
  call mpi_allgatherv(nl2gtmp,no_u,MPI_Integer,nl2G(1,1),counts,displs, &
                      MPI_integer,MPI_Comm_World, MPIerr)
  deallocate(nl2gtmp,counts,displs)
#endif

! Number of orbitals per node

  ! nOrbPerNode(Node+1) = dd_nuo
#ifdef MPI
  ! Now all_gather...
  call mpi_allgather(dd_nuo,1,MPI_Integer,nOrbPerNode,1, &
                      MPI_integer,MPI_Comm_World, MPIerr)
#else
  nOrbPerNode(1) = dd_nuo
#endif

end subroutine setup_ordern_indexes
