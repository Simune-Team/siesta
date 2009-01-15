      module m_norm
      CONTAINS
      function compute_norm(DM,nnz,nspin) result(norm)
      use precision, only: dp
      use sparse_matrices, only: S
      use parallel, only: node
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      integer, intent(in)   ::  nnz, nspin
      real(dp), intent(in)  ::  DM(:,:)
      real(dp)  ::       norm

#ifdef MPI
      real(dp):: buffer1   ! Temporary for the result of a globalize operation
#endif
      integer   ::  io, ispin

      print *, "*node, nnz, nspin:", node, nnz, nspin
      print *, "*node, sizeDM_1, sizeDM_2", node,
     $     size(DM,dim=1), size(DM,dim=2)
      print *, "*node, sizeS", node, size(S)
      norm = 0.0_dp
      do ispin = 1,min(nspin,2)
        do io = 1,nnz
          norm = norm + DM(io,ispin) * s(io)
        enddo
      enddo
#ifdef MPI
        call globalize_sum(norm,buffer1)
        norm = buffer1
#endif
      end function compute_norm
      end module m_norm
