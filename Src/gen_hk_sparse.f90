!
  ! A couple of routines to showcase the main issue in the implementation of
  ! k-points in the ELSI solver for Siesta
!
subroutine fold_sparse_arrays(no_l,no_u,no_s,numh,listhptr,nnz,listh, &
    indxuo,numh_u,listhptr_u,nnz_u,listh_u,ind2ind_u)

! Fold-in a sparse matrix from the auxiliary supercell into the
! unit cell (see picture)

      ! Number of orbitals handled by this node
      integer, intent(in) :: no_l
      ! Number of orbitals in unit cell
      integer, intent(in) :: no_u
      ! Number of orbitals in auxiliary supercell
      integer, intent(in) :: no_s
      
      ! Number of interactions per orbital; pointer to beginning of sparse array data
      integer, intent(in) :: numh(no_l), listhptr(no_l)
      ! Total number of interactions
      integer, intent(in) :: nnz
      ! Columns of the (rectangular) supercell-based matrix
      integer, intent(in) :: listh(nnz)

      ! Mapper of supercell-to-unit cell orbital indexes
      integer, intent(in) :: indxuo(no_s)


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
            juo = indxuo(jo)
            mask(juo) = 1    ! Found one
         enddo
         numh_u(iuo) = count(mask>0)
         nnz_u = nnz_u + numh_u(iuo)
      enddo
      
      allocate(listh_u(nnz_u))
      allocate(ind2ind_u(nnz))

      ! Generate folded pointer array
      listhptr_u(1) = 0
      do iuo = 2, no_l
         listhptr_u(iuo) = listhptr_u(iuo-1) + numh_u(iuo)
      enddo

      ! Complete the mapping
      do iuo = 1,no_l
         mask(:) = 0
         ind_u = listhptr_u(iuo)
         ind = listhptr(iuo)
         do j = 1,numh(iuo)
            ind = ind + 1
            jo = listh(ind)
            juo = indxuo(jo)
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
      
subroutine compute_sparse_hk_sk(no_l,no_u,no_s,numh,listhptr,nnz, &
                                listh,indxuo,nspin,H,S,xij,nk,kpoint)

      ! Number of orbitals handled by this node
      integer, intent(in) :: no_l
      ! Number of orbitals in unit cell
      integer, intent(in) :: no_u
      ! Number of orbitals in auxiliary supercell
      integer, intent(in) :: no_s

      ! Number of interactions per orbital; pointer to beginning of sparse array data
      integer, intent(in) :: numh(no_l), listhptr(no_l)
      ! Total number of interactions
      integer, intent(in) :: nnz
      ! Columns of the (rectangular) supercell-based matrix
      integer, intent(in) :: listh(nnz)
      ! Mapper of supercell-to-unit cell orbital indexes
      integer, intent(in) :: indxuo(no_s)

      integer, intent(in)  :: nspin         ! Number of spins (can be 1 or 2 in this simple version)
      real(dp), intent(in) :: H(nnz,nspin)  ! Hamiltonian sparse matrix
      real(dp), intent(in) :: S(nnz)        ! Overlap sparse matrix
      real(dp), intent(in) :: xij(3,nnz)    ! Vectors between orbital centers (sparse)

      integer, intent(in)  :: nk            ! Number of k-points
      real(dp), intent(in) :: kpoint(3,nk)  ! k-point vectors
      !  units: xij and kpoint must be in reciprocal coordinates of each other, so that xij*kpoint is a phase
      !         (see usage below)

      

      integer, allocatable :: numh_u(:), listhptr_u(:), listh_u(:)
      integer, allocatable :: ind2ind_u(:)

      complex(dp), allocatable :: Hk(:,:)
      complex(dp), allocatable :: Sk(:)
      
      integer ik, iuo, jo, juo, ind, ind_u, j, nnz_u
      real(dp) :: kxij
      complex(dp) :: kphs  ! phase factor
      
      call fold_sparse_arrays(no_l,no_u,no_s,numh,listhptr,nnz,listh, &
                             indxuo,numh_u,listhptr_u,nnz_u,listh_u,ind2ind_u)

      allocate(Hk(nnz_u,nspin), Sk(nnz_u,nspin))

      do ispin = 1, nspin  ! Serial over spins for now...
         
         do ik = 1, nk
     
            Sk = 0
            Hk = 0
            do iuo = 1,no_l
               do j = 1,numh(iuo)
                  ind = listhptr(iuo) + j
                  jo = listh(ind)
                  ind_u = ind2ind_u(ind)
                  juo = indxuo(jo)

                  kxij = kpoint(1,ik) * xij(1,ind) +    &
                         kpoint(2,ik) * xij(2,ind) +    &
                         kpoint(3,ik) * xij(3,ind) 
                  kphs = cdexp(dcmplx(0.0_dp, -1.0_dp)*kxij)

                  Sk(ind_u) = Sk(ind_u) + S(ind)*kphs
                  Hk(ind_u,ispin) = Hk(ind_u,ispin) + H(ind,ispin)*kphs
                  
               enddo 
            enddo
            ! Use Hk, Sk ...
            ! In ELSI:
            !   nnz_g = mpi_all_reduce(nnz_u, MPI_SUM....)
            !
            !   call elsi_set_csc(elsi_h, nnz_g, nnz_u, no_l, listh_u, listhptr_u)

                       ! Accumulate DM contrib for this k-point into full DM...
         enddo
      enddo

           
    end subroutine compute_sparse_hk_sk
