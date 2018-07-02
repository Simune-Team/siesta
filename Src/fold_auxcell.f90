module m_fold_auxcell
  public :: fold_sparse_arrays
  private
  CONTAINS
    subroutine fold_sparse_arrays(no_l,no_u,numh,listhptr,nnz,listh, &
               numh_u,listhptr_u,nnz_u,listh_u,ind2ind_u)

! Fold-in a sparse matrix from the auxiliary supercell into the
! unit cell (see picture)

      ! Number of orbitals handled by this node
      integer, intent(in) :: no_l
      ! Number of orbitals in unit cell
      integer, intent(in) :: no_u
      
      ! Number of interactions per orbital; pointer to beginning of sparse array data
      integer, intent(in) :: numh(no_l), listhptr(no_l)
      ! Total number of interactions
      integer, intent(in) :: nnz
      ! Columns of the (rectangular) supercell-based matrix
      integer, intent(in) :: listh(nnz)

      ! Output: All interactions are collapsed to the unit cell orbitals

      ! Number of interactions per orbital; pointer to beginning of sparse array data
      integer, intent(out) :: numh_u(no_l), listhptr_u(no_l)
      ! Total number of interactions
      integer, intent(out) :: nnz_u
      ! Columns of the (square) unit cell-based matrix
      integer, allocatable, intent(out) :: listh_u(:)
      ! Mapper of indexes from supercell-sparse-array to unit-cell sparse array
      integer, allocatable, intent(out) :: ind2ind_u(:)


      ! Local variables
      integer :: iuo, ind, j, jo, j_u, ind_u, juo
      integer, allocatable :: mask(:)

      allocate(mask(no_u))

      ! Find out how many "folded" interactions there are
      nnz_u = 0
      do iuo = 1,no_l
         mask = 0
         do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = mod(jo-1,no_u) + 1
            mask(juo) = 1    ! Found one
         enddo
         numh_u(iuo) = sum(mask)
         nnz_u = nnz_u + numh_u(iuo)
      enddo

      allocate(listh_u(nnz_u))
      allocate(ind2ind_u(nnz))

      ! Generate folded pointer array
      listhptr_u(1) = 0
      do iuo = 2, no_l
         listhptr_u(iuo) = listhptr_u(iuo-1) + numh_u(iuo-1)
      enddo

      ! Complete the mapping
      do iuo = 1,no_l
         mask(:) = 0
         ind_u = listhptr_u(iuo)
         ind = listhptr(iuo)
         do j = 1,numh(iuo)
            ind = ind + 1
            jo = listh(ind)
            juo = mod(jo-1,no_u) + 1
            if (mask(juo) > 0) then
               ! juo already seen and its place recorded
               ind2ind_u(ind) = mask(juo)
            else
               ind_u = ind_u + 1
               listh_u(ind_u) = juo
               ind2ind_u(ind) = ind_u   ! map
               mask(juo) = ind_u  ! record 
            endif
         enddo
      enddo
    end subroutine fold_sparse_arrays
end module m_fold_auxcell
