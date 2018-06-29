!
! Alberto Garcia, Feb. 2018
!
! ELSI DM-based interface to Siesta. It uses the sparse matrices from Siesta,
! and obtains the DM and EDM matrices in sparse form.
!
! The elsi_solver routine is in principle able to perform (spin-polarized)
! calculations for real matrices (i.e., at the Gamma point), and periodic
! calculations for complex matrices (i.e., multiple k-points).
!
! The structure of the solver routine is such that it will detect when it is
! called for the first scf step, so it can perform any needed initialization.
! The same idea can be used for the diagonalization mode, with (more extensive)
! appropriate changes.
!
! The module also exports the "elsi_finalize_scfloop" routine, to be called from
! the appropriate place. Some variables are kept at the module level for this.
!
! Usage: Compile Siesta with -DSIESTA__ELSI
!        Define
!           SolutionMethod ELSI
!        in the fdf file
! TODO:
!        -  Do not get EDM in every SCF step
!        -  MPI.Nprocs.SIESTA is not working -- maybe we will remove that feature
!        -  Spin works
!        -  k-points to be implemented soon
!        -  Add documentation
!
module m_elsi_interface

#ifdef SIESTA__ELSI

  use precision, only: dp
  use elsi

  implicit none

  private

  integer, parameter :: ELPA_SOLVER       = 1 ! solver
  integer, parameter :: OMM_SOLVER        = 2 ! solver
  integer, parameter :: PEXSI_SOLVER      = 3 ! solver
  integer, parameter :: SIPS_SOLVER       = 5 ! solver
  integer, parameter :: MULTI_PROC        = 1 ! parallel_mode
  integer, parameter :: PEXSI_CSC         = 1 ! distribution
  integer, parameter :: SIESTA_CSC        = 2 ! distribution
  integer, parameter :: GAUSSIAN          = 0 ! broadening
  integer, parameter :: FERMI             = 1 ! broadening
  integer, parameter :: METHFESSEL_PAXTON = 2 ! broadening
  integer, parameter :: CUBIC             = 3 ! broadening
  integer, parameter :: COLD              = 4 ! broadening
  integer, parameter :: ELSI_NOT_SET      = -910910

  type(elsi_handle) :: elsi_h

  integer :: elsi_global_comm
  integer :: which_solver

  real(dp), allocatable :: v_old(:,:)

  !  public :: elsi_solver
  public :: getdm_elsi
  public :: elsi_finalize_scfloop
  public :: elsi_save_potential

CONTAINS

subroutine getdm_elsi(iscf, no_s, nspin, no_l, maxnh, no_u,  &
     numh, listhptr, listh, H, S, qtot, temp, &
     xijo, indxuo, nkpnt, kpoint, kweight,    &
     eo, qo, Dscf, Escf, ef, Entropy, occtol, neigwanted)

  !
  ! Analogous to 'diagon', it dispatches ELSI solver routines as needed
  !

  use m_fold_auxcell, only: fold_sparse_arrays ! Could be called in state_init

! Missing for now
!  eo(no_u,nspin,nk)  : Eigenvalues 
!  qo(no_u,nspin,nk)  : Occupations of eigenstates


     real(dp), intent(inout) :: H(:,:), S(:)    ! Note!
     integer, intent(in) ::  iscf, maxnh, no_u, no_l, no_s, nkpnt
     integer, intent(in) ::  neigwanted, nspin
     integer, intent(in) ::  indxuo(no_s), listh(maxnh), numh(no_l), listhptr(no_l)

      real(dp), intent(out) ::  Dscf(maxnh,nspin), ef, Escf(maxnh,nspin), Entropy
      real(dp), intent(out) ::  eo(no_u,nspin,nkpnt), qo(no_u,nspin,nkpnt)
      real(dp), intent(in)  ::  kpoint(3,nkpnt), qtot, temp, kweight(nkpnt), occtol,xijo(3,maxnh)

      logical :: gamma, using_aux_cell
      integer, allocatable :: numh_u(:), listhptr_u(:), listh_u(:)
      integer, allocatable :: ind2ind_u(:)

      ! Real or complex, as appropriate
      real(dp), allocatable, dimension(:,:)  :: Dscf_k, Escf_k, Hk
      real(dp), allocatable, dimension(:)    :: Sk
      
      integer :: iuo, ispin, j, ind, ind_u, nnz_u
      !complex(dp) :: kphs

      external die

      gamma = ((nkpnt == 1) .and. (sum(abs(kpoint(:,1))) == 0.0_dp))
      using_aux_cell =  (no_s /= no_u)

      if (gamma) then
         if  (.not. using_aux_cell) then

            call elsi_solver(iscf, no_u, no_l, nspin, &
                     maxnh, listhptr, listh, qtot, temp, &
                     H, S, Dscf, Escf, ef, Entropy)

         else

            allocate(numh_u(no_l), listhptr_u(no_l))
            call fold_sparse_arrays(no_l,no_u,no_s,numh,listhptr,maxnh,listh, &
                             indxuo,numh_u,listhptr_u,nnz_u,listh_u,ind2ind_u)

            allocate(Hk(nnz_u,nspin), Sk(nnz_u))
            Sk = 0
            Hk = 0
            do ispin = 1, nspin  ! Serial over spins for now...
               do iuo = 1,no_l
                  do j = 1,numh(iuo)
                     ind = listhptr(iuo) + j
                     ind_u = ind2ind_u(ind)

                  !kxij = kpoint(1,ik) * xij(1,ind) +    &
                  !       kpoint(2,ik) * xij(2,ind) +    &
                  !       kpoint(3,ik) * xij(3,ind) 
                  !kphs = cdexp(dcmplx(0.0_dp, -1.0_dp)*kxij)

                     Sk(ind_u) = Sk(ind_u) + S(ind)     !*kphs
                     Hk(ind_u,ispin) = Hk(ind_u,ispin) + H(ind,ispin)   !*kphs
                  
                  enddo
               enddo
            enddo

            allocate(Dscf_k(nnz_u,nspin), Escf_k(nnz_u,nspin))
            call elsi_solver(iscf, no_u, no_l, nspin, &
                     nnz_u, listhptr_u, listh_u, qtot, temp, &
                     Hk, Sk, Dscf_k, Escf_k, ef, Entropy)
            deallocate(Hk, Sk)


            ! Unfold
            do ispin = 1, nspin
               do iuo = 1,no_l
                  do j = 1,numh(iuo)
                     ind = listhptr(iuo) + j
                     ind_u = ind2ind_u(ind)

                     !kxij = kpoint(1,ik) * xij(1,ind) +    &
                     !kpoint(2,ik) * xij(2,ind) +    &
                     !kpoint(3,ik) * xij(3,ind) 
                     !ckxij = cos(kxij)
                     !skxij = sin(kxij)
!                     Dscf(ind,ispin)=Dscf_k(ind_u,ispin)*ckxij - &
!                                     Dscf(2,juo,iuo)*skxij
                     Dscf(ind,ispin) = Dscf_k(ind_u,ispin)
                     Escf(ind,ispin) = Escf_k(ind_u,ispin)
                  enddo
               enddo
            enddo
            deallocate(Dscf_k, Escf_k)
            deallocate(numh_u, listhptr_u)
            deallocate(listh_u, ind2ind_u)
         endif
      else
         ! Always using aux cell...
         call die("Cannot do k-points with ELSI yet")
         ! Fold arrays
         !call elsi_solver_complex( )
         ! UNFOLD !!!
      endif

end subroutine getdm_elsi

! This version assumes *the same* distributions for Siesta (setup_H et al) and ELSI
! operations.  ELSI tasks 1D block-cyclic distributed CSC/CSR matrices as its
! input/output.
!
subroutine elsi_solver(iscf, n_basis, n_basis_l, n_spin, nnz_l, row_ptr, &
  col_idx, qtot, temp, ham, ovlp, dm, edm, ef, ets)

  use fdf,         only: fdf_get
  use m_mpi_utils, only: globalize_sum
  use parallel,    only: BlockSize
#ifdef MPI
  use mpi_siesta
#endif
  use class_Distribution
  use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
  use alloc
  use m_mpi_utils, only: globalize_sum
  
  implicit none

  integer,  intent(in)    :: iscf      ! SCF step counter
  integer,  intent(in)    :: n_basis   ! Global basis
  integer,  intent(in)    :: n_basis_l ! Local basis
  integer,  intent(in)    :: n_spin
  integer,  intent(in)    :: nnz_l     ! Local nonzero
  integer,  intent(in)    :: row_ptr(n_basis_l)
  integer,  intent(in), target    :: col_idx(nnz_l)
  real(dp), intent(in)    :: qtot
  real(dp), intent(in)    :: temp
  real(dp), intent(inout), target :: ham(nnz_l,n_spin)
  real(dp), intent(inout), target :: ovlp(nnz_l)
  real(dp), intent(out)   :: dm(nnz_l,n_spin)
  real(dp), intent(out)   :: edm(nnz_l,n_spin)
  real(dp), intent(out)   :: ef        ! Fermi energy
  real(dp), intent(out)   :: ets       ! Entropy/k, dimensionless

  integer :: mpirank
  integer :: ierr
  integer :: n_state
  integer :: nnz_g
  integer :: which_broad
  integer :: mp_order
  integer :: out_level
  integer :: out_json
  integer :: elpa_flavor
  integer :: omm_flavor
  integer :: omm_n_elpa
  integer :: pexsi_tasks_per_pole
  integer :: pexsi_tasks_symbolic
  integer :: sips_n_slice
  integer :: sips_n_elpa

  real(dp) :: omm_tol
  real(dp) :: energy

  character(len=5) :: solver_string
  character(len=5) :: broad_string

  integer, allocatable, dimension(:) :: row_ptr2
  integer, allocatable, dimension(:), target :: numh

  integer :: World_group
  integer :: elsi_Spatial_comm, elsi_Spatial_group
  integer :: elsi_Spin_comm

  type(distribution) :: dist1
  type(distribution), target  :: dist2_spin(2)
  type(distribution), pointer :: dist2
  
  type(aux_matrix), allocatable, target :: m1_spin(:)
  type(aux_matrix) :: m2
  type(aux_matrix), pointer :: m1

  integer :: color, spatial_rank, spin_rank
  integer :: npTotal, npPerSpin
  integer, allocatable, dimension(:) :: global_ranks_in_world
  integer, allocatable, dimension(:) :: ranks_in_world  ! should be 'spatial_ranks...'
  integer, allocatable, dimension(:,:) :: ranks_in_World_Spin

  integer :: elsi_spin

  integer :: numColLocal   ! change name later ....
  integer :: nnzLocal   ! change name later ....
  integer :: nnz
  integer, pointer  :: colptrLocal(:) => null()
  integer  :: i, ih, ispin
  real(dp) :: ets_spin

  integer, pointer  :: rowindLocal(:)
  real(dp), pointer :: SnzvalLocal(:)
  real(dp), pointer :: HnzvalLocal(:)
  real(dp), pointer :: DMnzvalLocal(:) => null()
  real(dp), pointer :: EDMnzvalLocal(:) => null()
  
  external :: timer

#ifndef MPI
  call die("This ELSI solver interface needs MPI")
#endif

  ! Global communicator is a duplicate of passed communicator
  call MPI_Comm_Dup(MPI_Comm_DFT, elsi_global_comm, ierr)
  call MPI_Comm_Rank(elsi_global_comm, mpirank, ierr)
  call MPI_Comm_Size(elsi_global_comm, npTotal, ierr)

  call timer("elsi", 1)

  ! Initialization
  if (iscf == 1) then

    ! Get ELSI options
    solver_string        = fdf_get("ELSISolver", "elpa")
    out_level            = fdf_get("ELSIOutputLevel", 0)
    out_json             = fdf_get("ELSIOutputJson", 1)
    broad_string         = fdf_get("ELSIBroadeningMethod", "fermi")
    mp_order             = fdf_get("ELSIBroadeningMPOrder", 1)
    elpa_flavor          = fdf_get("ELSIELPAFlavor", 2)
    omm_flavor           = fdf_get("ELSIOMMFlavor", 0)
    omm_n_elpa           = fdf_get("ELSIOMMELPASteps", 3)
    omm_tol              = fdf_get("ELSIOMMTolerance", 1.0e-9_dp)
    pexsi_tasks_per_pole = fdf_get("ELSIPEXSITasksPerPole", ELSI_NOT_SET)
    pexsi_tasks_symbolic = fdf_get("ELSIPEXSITasksSymbolic", 1)
    sips_n_slice         = fdf_get("ELSISIPSSlices", ELSI_NOT_SET)
    sips_n_elpa          = fdf_get("ELSISIPSELPASteps", 2)

    select case (solver_string)
    case ("elpa", "ELPA")
      which_solver = ELPA_SOLVER
    case ("omm", "OMM")
      which_solver = OMM_SOLVER
    case ("pexsi", "PEXSI")
      which_solver = PEXSI_SOLVER
    case ("sips", "SIPS")
      which_solver = SIPS_SOLVER
    case default
      which_solver = ELPA_SOLVER
    end select

    select case (broad_string)
    case ("gauss")
      which_broad = GAUSSIAN
    case ("fermi")
      which_broad = FERMI
    case ("mp")
      which_broad = METHFESSEL_PAXTON
    case ("cubic")
      which_broad = CUBIC
    case ("cold")
      which_broad = COLD
    case default
      which_broad = FERMI
    end select

    ! Number of states to solve when calling an eigensolver
    n_state = min(n_basis, n_basis/2+5)

    ! Now we have all ingredients to initialize ELSI
    call elsi_init(elsi_h, which_solver, MULTI_PROC, SIESTA_CSC, n_basis, &
      qtot, n_state)

    ! Output
    call elsi_set_output(elsi_h, out_level)
    call elsi_set_output_log(elsi_h, out_json)
    call elsi_set_write_unit(elsi_h, 6)

    ! Broadening
    call elsi_set_mu_broaden_scheme(elsi_h, which_broad)
    call elsi_set_mu_broaden_width(elsi_h, temp)
    call elsi_set_mu_mp_order(elsi_h, mp_order)

    ! Solver settings
    call elsi_set_elpa_solver(elsi_h, elpa_flavor)
    call elsi_set_omm_flavor(elsi_h, omm_flavor)
    call elsi_set_omm_n_elpa(elsi_h, omm_n_elpa)
    call elsi_set_omm_tol(elsi_h, omm_tol)

    if (pexsi_tasks_per_pole /= ELSI_NOT_SET) then
      call elsi_set_pexsi_np_per_pole(elsi_h, pexsi_tasks_per_pole)
    end if

    call elsi_set_pexsi_np_symbo(elsi_h, pexsi_tasks_symbolic)
    call elsi_set_pexsi_temp(elsi_h, temp)
!    call elsi_set_pexsi_gap(elsi_h, gap)
!    call elsi_set_pexsi_delta_e(elsi_h, delta_e)
    call elsi_set_pexsi_mu_min(elsi_h, -10.0_dp)
    call elsi_set_pexsi_mu_max(elsi_h, 10.0_dp)

    if (sips_n_slice /= ELSI_NOT_SET) then
      call elsi_set_sips_n_slice(elsi_h, sips_n_slice)
    end if

    call elsi_set_sips_n_elpa(elsi_h, sips_n_elpa)
    call elsi_set_sips_interval(elsi_h, -10.0_dp, 10.0_dp)

 endif   ! iscf == 1
  
    if (n_spin == 1) then
       
       ! Sparsity pattern
       call globalize_sum(nnz_l, nnz_g, comm=elsi_global_comm)

       allocate(row_ptr2(n_basis_l+1))
       row_ptr2(1:n_basis_l) = row_ptr(1:n_basis_l)+1
       row_ptr2(n_basis_l+1) = nnz_l+1

       call elsi_set_csc(elsi_h, nnz_g, nnz_l, n_basis_l, col_idx, row_ptr2)
       deallocate(row_ptr2)

       call elsi_set_csc_blk(elsi_h, BlockSize)
       call elsi_set_mpi(elsi_h, elsi_global_comm)

    else
       
       ! MPI logic for spin polarization

       ! Re-create numh, as we use it in the transfer
       allocate(numh(n_basis_l))
       numh(1) = row_ptr(2)
       do i = 2, n_basis_l-1
          numh(i) = row_ptr(i+1)-row_ptr(i)
       enddo
       numh(n_basis_l) = nnz_l - row_ptr(n_basis_l)

       ! Define the original distribution (over all the available nodes)

       allocate(global_ranks_in_world(npTotal))
       global_ranks_in_world = (/ (i, i=0, npTotal-1) /)
       call newDistribution(dist1,elsi_global_Comm,global_ranks_in_world, &
            TYPE_BLOCK_CYCLIC,BlockSize,"global dist")
       deallocate(global_ranks_in_world)
       call MPI_Barrier(elsi_global_Comm,ierr)


       ! "Row" communicator for independent PEXSI operations on each spin
       ! The name refers to "spatial" degrees of freedom.
       color = mod(mpirank,n_spin)    ! {0,1} for n_spin = 2, or {0} for n_spin = 1
       call MPI_Comm_Split(elsi_global_comm, color, mpirank, elsi_Spatial_Comm, ierr)
       ! "Column" communicator for spin reductions
       color = mpirank/n_spin
       call MPI_Comm_Split(elsi_global_comm, color, mpirank, elsi_Spin_Comm, ierr)

       call mpi_comm_rank( elsi_Spatial_Comm, spatial_rank, ierr )
       call mpi_comm_rank( elsi_Spin_Comm, spin_rank, ierr )
    
       ! Include the actual world ranks in the distribution objects

       npPerSpin = npTotal/n_spin

       call MPI_Comm_Group(elsi_Spatial_Comm, elsi_Spatial_Group, Ierr)
       allocate(ranks_in_world(npPerSpin))
       call MPI_Comm_Group(elsi_global_comm, World_Group, Ierr)
       call MPI_Group_translate_ranks( elsi_Spatial_Group, npPerSpin, &
            (/ (i,i=0,npPerSpin-1) /), &
            World_Group, ranks_in_world, ierr )

       call MPI_Group_Free(elsi_Spatial_Group, ierr)
       call MPI_Group_Free(World_Group, ierr)
       
       allocate (ranks_in_World_Spin(npPerSpin,n_spin))
       call MPI_AllGather(ranks_in_world,npPerSpin,MPI_integer,&
            Ranks_in_World_Spin(1,1),npPerSpin, &
            MPI_integer,elsi_Spin_Comm,ierr)

       ! Create distributions known to all nodes (we used allgather...)
       do ispin = 1, n_spin
          call newDistribution(dist2_spin(ispin), elsi_global_Comm, &
               Ranks_in_World_Spin(:,ispin),  &
               TYPE_BLOCK_CYCLIC, blockSize, "SPATIAL dist")
       enddo
       deallocate(ranks_in_world,Ranks_in_World_Spin)
       call MPI_Barrier(elsi_global_Comm,ierr)

       elsi_spin = spin_rank+1  ! {1,2}

       
       ! This is done serially, each time filling one spin set
       ! Note that **all processes** need to have the same m1
       ! but we probably do not need two copies of m1

       ! Just do allocate(m1%vals(2)) outside the loop
       
       allocate(m1_spin(n_spin))
       do ispin = 1, n_spin

          m1 => m1_spin(ispin)

          m1%norbs = n_basis
          m1%no_l  = n_basis_l
          m1%nnzl  = nnz_l
          m1%numcols => numh  ! numh
          m1%cols    => col_idx
          allocate(m1%vals(2))
          m1%vals(1)%data => ovlp(:)
          m1%vals(2)%data => ham(:,ispin)

          call timer("redist_orbs_fwd", 1)
 
          dist2 => dist2_spin(ispin)
          call redistribute_spmatrix(n_basis,m1,dist1,m2,dist2,elsi_global_Comm)
   
          call timer("redist_orbs_fwd", 2)

          if (elsi_spin == ispin) then  ! Each team gets their own data

             !nrows = m2%norbs          ! or simply 'norbs'
             numColLocal = m2%no_l
             nnzLocal    = m2%nnzl
             call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,elsi_Spatial_Comm,ierr)
             call re_alloc(colptrLocal,1,numColLocal+1,"colptrLocal","elsi_solver")
             colptrLocal(1) = 1
             do ih = 1,numColLocal
                colptrLocal(ih+1) = colptrLocal(ih) + m2%numcols(ih)
             enddo

             rowindLocal => m2%cols
             SnzvalLocal => m2%vals(1)%data
             HnzvalLocal => m2%vals(2)%data

             call re_alloc(DMnzvalLocal,1,nnzLocal,"DMnzvalLocal","elsi_solver")
             call re_alloc(EDMnzvalLocal,1,nnzLocal,"EDMnzvalLocal","elsi_solver")
          endif
       enddo

       call elsi_set_csc(elsi_h, nnz, nnzLocal, numColLocal, rowindLocal, colPtrLocal)
       call de_alloc(colPtrLocal,"colPtrLocal","elsi_solver")
       
       call elsi_set_csc_blk(elsi_h, BlockSize)
       call elsi_set_spin(elsi_h, n_spin, elsi_spin)
       call elsi_set_mpi(elsi_h, elsi_Spatial_comm)
       call elsi_set_mpi_global(elsi_h, elsi_global_comm)

    endif  ! n_spin
    
  call timer("elsi-solver", 1)

  if (n_spin == 1) then
     call elsi_dm_real_sparse(elsi_h, ham, ovlp, DM, energy)
     call elsi_get_edm_real_sparse(elsi_h, EDM)
     call elsi_get_entropy(elsi_h, ets)
  else
     ! Solve DM, and get (at every step for now) EDM, Fermi energy, and entropy
     ! Presumably energy is already summed over spins
     call elsi_dm_real_sparse(elsi_h, HnzValLocal, SnzValLocal, DMnzvallocal, energy)
     call elsi_get_edm_real_sparse(elsi_h, EDMnzvallocal)
     !... we need to sum this over spins
     call elsi_get_entropy(elsi_h, ets_spin)
     call globalize_sum(ets_spin, ets, comm=elsi_Spin_comm)
  endif

  ! We assume also that ef and ets are known to all nodes
  call elsi_get_mu(elsi_h, ef)


  ets = ets/temp

  call timer("elsi-solver", 2)

  if ( n_spin == 2) then
     ! Now we need to redistribute back

     do ispin = 1, n_spin

        m1 => m1_spin(ispin)

        if (elsi_spin == ispin) then
           ! Prepare m2 to transfer
           ! The other fields are the same
           call de_alloc(m2%vals(1)%data,"m2%vals(1)%data","elsi_solver")
           call de_alloc(m2%vals(2)%data,"m2%vals(2)%data","elsi_solver")

           m2%vals(1)%data => DMnzvalLocal(1:nnzLocal)
           m2%vals(2)%data => EDMnzvalLocal(1:nnzLocal)

        endif

        ! Prepare m1 to receive the results

        nullify(m1%vals(1)%data)    ! formerly pointing to S
        nullify(m1%vals(2)%data)    ! formerly pointing to H
        deallocate(m1%vals)
        nullify(m1%numcols)         ! formerly pointing to numH
        nullify(m1%cols)            ! formerly pointing to listH

        call timer("redist_orbs_bck", 1)
        dist2 => dist2_spin(ispin)
        call redistribute_spmatrix(n_basis,m2,dist2,m1,dist1,elsi_global_Comm)
        call timer("redist_orbs_bck", 2)

        if (elsi_spin == ispin) then
           ! Each team deallocates during "its" spin cycle
           call de_alloc(DMnzvalLocal, "DMnzvalLocal", "elsi_solver")
           call de_alloc(EDMnzvalLocal,"EDMnzvalLocal","elsi_solver")

           nullify(m2%vals(1)%data)    ! formerly pointing to DM
           nullify(m2%vals(2)%data)    ! formerly pointing to EDM
           deallocate(m2%vals)
           ! allocated in the direct transfer
           call de_alloc(m2%numcols,"m2%numcols","elsi_solver")
           call de_alloc(m2%cols,   "m2%cols",   "elsi_solver")
        endif


        ! In future, m1%vals(1,2) could be pointing to DM and EDM,
        ! and the 'redistribute' routine check whether the vals arrays are
        ! associated, to use them instead of allocating them.
        DM(:,ispin)  = m1%vals(1)%data(:)    
        EDM(:,ispin) = m1%vals(2)%data(:)    
        ! Check no_l
        if (n_basis_l /= m1%no_l) then
           call die("Mismatch in no_l")
        endif
        ! Check listH
        if (any(col_idx(:) /= m1%cols(:))) then
           call die("Mismatch in listH")
        endif

        ! Do this here if using two copies of m1. Otherwise, put it
        ! outside the loop
        call de_alloc(m1%vals(1)%data,"m1%vals(1)%data","elsi_solver")
        call de_alloc(m1%vals(2)%data,"m1%vals(2)%data","elsi_solver")
        deallocate(m1%vals)
        ! allocated in the direct transfer
        call de_alloc(m1%numcols,"m1%numcols","elsi_solver") 
        call de_alloc(m1%cols,   "m1%cols",   "elsi_solver")

     enddo

     call MPI_Comm_Free(elsi_Spatial_comm, ierr)
     call MPI_Comm_Free(elsi_Spin_comm, ierr)
     
  endif

  call timer("elsi", 2)

  
end subroutine elsi_solver

! Clean up:  Finalize ELSI instance and free MPI communicator.
!
subroutine elsi_finalize_scfloop()

  integer :: ierr

  ! Make which_solver invalid
  which_solver = ELSI_NOT_SET

  if (allocated(v_old)) then
    deallocate(v_old)
  end if

  call elsi_finalize(elsi_h)

  call MPI_Comm_Free(elsi_global_comm, ierr)

end subroutine elsi_finalize_scfloop

! Save Hartree + XC potential and find minimum and maximum change of it between
! two SCF iterations.
!
subroutine elsi_save_potential(n_pts, n_spin, v_scf)

  use m_mpi_utils, only: globalize_min, globalize_max

  integer,  intent(in) :: n_pts
  integer,  intent(in) :: n_spin
  real(dp), intent(in) :: v_scf(n_pts,n_spin)

  real(dp) :: mu_min
  real(dp) :: mu_max
  real(dp) :: dv_min
  real(dp) :: dv_max
  real(dp) :: tmp

  if (which_solver == PEXSI_SOLVER) then
    if (.not. allocated(v_old)) then
      allocate(v_old(n_pts,n_spin))

      v_old = v_scf

      mu_min = -10.0_dp
      mu_max = 10.0_dp
    else
      call elsi_get_pexsi_mu_min(elsi_h, mu_min)
      call elsi_get_pexsi_mu_max(elsi_h, mu_max)

      v_old = v_scf-v_old

      ! Get minimum and maximum of change of total potential
      tmp = minval(v_old)

      call globalize_min(tmp, dv_min, comm=elsi_global_comm)

      tmp = maxval(v_old)

      call globalize_max(tmp, dv_max, comm=elsi_global_comm)

      mu_min = mu_min+dv_min
      mu_max = mu_max+dv_max
    end if

    ! Adjust chemical potential range for PEXSI
    call elsi_set_pexsi_mu_min(elsi_h, mu_min)
    call elsi_set_pexsi_mu_max(elsi_h, mu_max)
  end if

end subroutine elsi_save_potential

#endif  /* SIESTA__ELSI */

end module m_elsi_interface
