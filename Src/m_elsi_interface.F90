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
     xijo, nkpnt, kpoint, kweight,    &
     eo, qo, Dscf, Escf, ef, Entropy, occtol, neigwanted)

  !
  ! Analogous to 'diagon', it dispatches ELSI solver routines as needed
  !

  use m_fold_auxcell, only: fold_sparse_arrays ! Could be called in state_init

! Missing for now
!  eo(no_u,nspin,nk)  : Eigenvalues 
!  qo(no_u,nspin,nk)  : Occupations of eigenstates


     real(dp), intent(inout) :: H(:,:), S(:)    ! Note: we might overwrite these
     integer, intent(in) ::  iscf, maxnh, no_u, no_l, no_s, nkpnt
     integer, intent(in) ::  neigwanted, nspin
     integer, intent(in) ::  listh(maxnh), numh(no_l), listhptr(no_l)
     real(dp), intent(in)  ::  kpoint(3,nkpnt), qtot, temp, kweight(nkpnt), occtol,xijo(3,maxnh)

      real(dp), intent(out) ::  Dscf(maxnh,nspin), ef, Escf(maxnh,nspin), Entropy
      real(dp), intent(out) ::  eo(no_u,nspin,nkpnt), qo(no_u,nspin,nkpnt)

      logical :: gamma, using_aux_cell
      integer, allocatable :: numh_u(:), listhptr_u(:), listh_u(:)
      integer, allocatable :: ind2ind_u(:)

      ! Real or complex, as appropriate
      real(dp), allocatable, dimension(:,:)  :: Dscf_k, Escf_k, Hk
      real(dp), allocatable, dimension(:)    :: Sk
      
      integer :: iuo, ispin, j, ind, ind_u, nnz_u

      external die

      gamma = ((nkpnt == 1) .and. (sum(abs(kpoint(:,1))) == 0.0_dp))
      using_aux_cell =  (no_s /= no_u)

      if (gamma) then
         if  (.not. using_aux_cell) then

            call elsi_real_solver(iscf, no_u, no_l, nspin, &
                     maxnh, listhptr, listh, qtot, temp, &
                     H, S, Dscf, Escf, ef, Entropy)

         else

            allocate(numh_u(no_l), listhptr_u(no_l))
            call fold_sparse_arrays(no_l,no_u,  &
                                    numh,listhptr,maxnh,listh, &
                                    numh_u,listhptr_u,nnz_u,listh_u,ind2ind_u)

            allocate(Hk(nnz_u,nspin), Sk(nnz_u))
            Sk = 0
            Hk = 0
               do iuo = 1,no_l
                  do j = 1,numh(iuo)
                     ind = listhptr(iuo) + j
                     ind_u = ind2ind_u(ind)
                     Sk(ind_u) = Sk(ind_u) + S(ind) 
                     do ispin = 1, nspin 
                        Hk(ind_u,ispin) = Hk(ind_u,ispin) + H(ind,ispin) 
                     enddo
                  enddo
               enddo

            allocate(Dscf_k(nnz_u,nspin), Escf_k(nnz_u,nspin))
            call elsi_real_solver(iscf, no_u, no_l, nspin, &
                     nnz_u, listhptr_u, listh_u, qtot, temp, &
                     Hk, Sk, Dscf_k, Escf_k, ef, Entropy)
            deallocate(Hk, Sk)
            deallocate(numh_u, listhptr_u, listh_u)


            ! Unfold 
            do iuo = 1,no_l
               do j = 1,numh(iuo)
                  ind = listhptr(iuo) + j
                  ind_u = ind2ind_u(ind)
                  do ispin = 1, nspin
                     Dscf(ind,ispin) = Dscf_k(ind_u,ispin)
                     Escf(ind,ispin) = Escf_k(ind_u,ispin)
                  enddo
               enddo
            enddo
            deallocate(Dscf_k, Escf_k)
            deallocate(ind2ind_u)

         endif  ! using auxiliary supercell with Gamma sampling
         
      else
         ! We need more preparation
         call elsi_kpoints_dispatcher(iscf, no_s, nspin, no_l, maxnh, no_u,  &
              numh, listhptr, listh, H, S, qtot, temp, &
              xijo, nkpnt, kpoint, kweight,    &
              eo, qo, Dscf, Escf, ef, Entropy, occtol, neigwanted)
         
      endif

end subroutine getdm_elsi

!======================================================================================
! ELSI takes 1D block-cyclic distributed CSC/CSR matrices as its
! input/output.
!
! But note taht this version assumes *the same* distributions for Siesta (setup_H et al) and ELSI
! operations.  
!
subroutine elsi_real_solver(iscf, n_basis, n_basis_l, n_spin, nnz_l, row_ptr, &
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

  integer :: elsi_Spatial_comm, elsi_Spin_comm

  type(distribution) :: dist_global
  type(distribution) :: dist_spin(2)
  
  type(aux_matrix) :: pkg_global, pkg_spin   ! Packages for transfer

  integer :: my_spin

  integer :: my_no_l  
  integer :: my_nnz_l 
  integer :: my_nnz
  integer, pointer  :: my_row_ptr2(:) => null()
  integer  :: i, ih, ispin, spin_rank
  real(dp) :: ets_spin

  integer, pointer  :: my_col_idx(:)
  real(dp), pointer :: my_S(:)
  real(dp), pointer :: my_H(:)
  real(dp), pointer :: my_DM(:) => null()
  real(dp), pointer :: my_EDM(:) => null()
  
  external :: timer

#ifndef MPI
  call die("This ELSI solver interface needs MPI")
#endif

  ! Global communicator is a duplicate of passed communicator
  call MPI_Comm_Dup(MPI_Comm_DFT, elsi_global_comm, ierr)

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

       ! Split the communicator in spins and get distribution objects
       ! for the data redistribution needed
       ! Note that dist_spin is an array
       call get_spin_comms_and_dists(elsi_global_comm,elsi_global_comm, &
            blocksize, n_spin, &
            dist_global,dist_spin, elsi_spatial_comm, elsi_spin_comm)

       ! Find out which spin team we are in, and tag the spin we work on
       call mpi_comm_rank( elsi_Spin_Comm, spin_rank, ierr )
       my_spin = spin_rank+1  ! {1,2}

       
       ! This is done serially, each time filling one spin set
       ! Note that **all processes** need to have the same pkg_global
       
       do ispin = 1, n_spin

          ! Load pkg_global data package
          pkg_global%norbs = n_basis
          pkg_global%no_l  = n_basis_l
          pkg_global%nnzl  = nnz_l
          pkg_global%numcols => numh  
          pkg_global%cols    => col_idx

          allocate(pkg_global%vals(2))
          ! Link the vals items to the appropriate arrays (no extra memory here)
          pkg_global%vals(1)%data => ovlp(:)
          ! Note that we *cannot* say  => ham(:,my_spin)
          ! and avoid the sequential loop, as then half the processors will send
          ! the information for 'spin up' and the other half the information for 'spin down',
          ! which is *not* what we want.
          pkg_global%vals(2)%data => ham(:,ispin)

          call timer("redist_orbs_fwd", 1)

          ! We are doing the transfers sequentially. One spin team is
          ! 'idle' (in the receiving side) in each pass, as the dist_spin(ispin) distribution
          ! does not involve them.
          
          call redistribute_spmatrix(n_basis,pkg_global,dist_global, &
                                             pkg_spin,dist_spin(ispin),elsi_global_Comm)
   
          call timer("redist_orbs_fwd", 2)

          if (my_spin == ispin) then  ! Each team gets their own data

             !nrows = pkg_spin%norbs          ! or simply 'norbs'
             my_no_l = pkg_spin%no_l
             my_nnz_l    = pkg_spin%nnzl
             call MPI_AllReduce(my_nnz_l,my_nnz,1,MPI_integer,MPI_sum,elsi_Spatial_Comm,ierr)
             ! generate off-by-one row pointer
             call re_alloc(my_row_ptr2,1,my_no_l+1,"my_row_ptr2","elsi_solver")
             my_row_ptr2(1) = 1
             do ih = 1,my_no_l
                my_row_ptr2(ih+1) = my_row_ptr2(ih) + pkg_spin%numcols(ih)
             enddo

             my_col_idx => pkg_spin%cols
             my_S => pkg_spin%vals(1)%data
             my_H => pkg_spin%vals(2)%data

             call re_alloc(my_DM,1,my_nnz_l,"my_DM","elsi_solver")
             call re_alloc(my_EDM,1,my_nnz_l,"my_EDM","elsi_solver")
          endif

          ! Clean pkg_global
          nullify(pkg_global%vals(1)%data)  
          nullify(pkg_global%vals(2)%data)  
          deallocate(pkg_global%vals)
          nullify(pkg_global%numcols)   
          nullify(pkg_global%cols)      

       enddo

       call elsi_set_csc(elsi_h, my_nnz, my_nnz_l, my_no_l, my_col_idx, my_row_ptr2)
       call de_alloc(my_row_ptr2,"my_row_ptr2","elsi_solver")
       
       call elsi_set_csc_blk(elsi_h, BlockSize)
       call elsi_set_spin(elsi_h, n_spin, my_spin)
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
     ! Energy is already summed over spins
     call elsi_dm_real_sparse(elsi_h, my_H, my_S, my_DM, energy)
     call elsi_get_edm_real_sparse(elsi_h, my_EDM)
     !... but we still need to sum the entropy over spins
     call elsi_get_entropy(elsi_h, ets_spin)
     call globalize_sum(ets_spin, ets, comm=elsi_Spin_comm)
  endif

  call elsi_get_mu(elsi_h, ef)
  ets = ets/temp
  
  ! Ef, energy, and ets are known to all nodes

  call timer("elsi-solver", 2)

  if ( n_spin == 2) then
     ! Now we need to redistribute back

     do ispin = 1, n_spin

        if (my_spin == ispin) then
           ! Prepare pkg_spin to transfer the right spin information
           ! The other fields (numcols, cols) are the same and are still there
           ! Deallocate my_S and my_H
           call de_alloc(pkg_spin%vals(1)%data,"pkg_spin%vals(1)%data","elsi_solver")
           call de_alloc(pkg_spin%vals(2)%data,"pkg_spin%vals(2)%data","elsi_solver")

           pkg_spin%vals(1)%data => my_DM(1:my_nnz_l)
           pkg_spin%vals(2)%data => my_EDM(1:my_nnz_l)

        endif

        ! pkg_global is clean now
        call timer("redist_orbs_bck", 1)
        call redistribute_spmatrix(n_basis,pkg_spin,dist_spin(ispin) &
                                          ,pkg_global,dist_global,elsi_global_Comm)
        call timer("redist_orbs_bck", 2)

        ! Clean pkg_spin
        if (my_spin == ispin) then
           ! Each team deallocates during "its" spin cycle
           call de_alloc(my_DM, "my_DM", "elsi_solver")
           call de_alloc(my_EDM,"my_EDM","elsi_solver")

           nullify(pkg_spin%vals(1)%data)    ! formerly pointing to DM
           nullify(pkg_spin%vals(2)%data)    ! formerly pointing to EDM
           deallocate(pkg_spin%vals)
           ! allocated in the direct transfer
           call de_alloc(pkg_spin%numcols,"pkg_spin%numcols","elsi_solver")
           call de_alloc(pkg_spin%cols,   "pkg_spin%cols",   "elsi_solver")
        endif


        ! In future, pkg_global%vals(1,2) could be pointing to DM and EDM,
        ! and the 'redistribute' routine check whether the vals arrays are
        ! associated, to use them instead of allocating them.
        DM(:,ispin)  = pkg_global%vals(1)%data(:)    
        EDM(:,ispin) = pkg_global%vals(2)%data(:)    
        ! Check no_l
        if (n_basis_l /= pkg_global%no_l) then
           call die("Mismatch in no_l")
        endif
        ! Check listH
        if (any(col_idx(:) /= pkg_global%cols(:))) then
           call die("Mismatch in listH")
        endif

        ! Clean pkg_global
        ! allocated by the transfer routine, but we did not actually
        ! look at them
        call de_alloc(pkg_global%numcols,"pkg_global%numcols","elsi_solver") 
        call de_alloc(pkg_global%cols,   "pkg_global%cols",   "elsi_solver")

        call de_alloc(pkg_global%vals(1)%data,"pkg_global%vals(1)%data","elsi_solver")
        call de_alloc(pkg_global%vals(2)%data,"pkg_global%vals(2)%data","elsi_solver")
        deallocate(pkg_global%vals)

     enddo

     call MPI_Comm_Free(elsi_Spatial_comm, ierr)
     call MPI_Comm_Free(elsi_Spin_comm, ierr)
     
  endif

  call timer("elsi", 2)

  
end subroutine elsi_real_solver

subroutine elsi_kpoints_dispatcher(iscf, no_s, nspin, no_l, maxnh, no_u,  &
     numh, listhptr, listh, H, S, qtot, temp, &
     xijo, nkpnt, kpoint, kweight,    &
     eo, qo, Dscf, Escf, ef, Entropy, occtol, neigwanted)

  use mpi_siesta, only: mpi_comm_dft
  use mpi
  use parallel, only: blocksize
  use class_Distribution
    use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
    use alloc

  !
  ! K-point redistribution, Hk and Sk building, and call to complex ELSI solver
  !

  use m_fold_auxcell, only: fold_sparse_arrays ! Could be called in state_init

! Missing for now
!  eo(no_u,nspin,nk)  : Eigenvalues 
!  qo(no_u,nspin,nk)  : Occupations of eigenstates


     real(dp), intent(in), target :: H(:,:), S(:)    ! Note that now we do not change them
     integer, intent(in) ::  iscf, maxnh, no_u, no_l, no_s, nkpnt
     integer, intent(in) ::  neigwanted, nspin
     integer, intent(in) ::  listhptr(no_l)
     integer, intent(in), target ::  listh(maxnh), numh(no_l)

      real(dp), intent(out) ::  Dscf(maxnh,nspin), ef, Escf(maxnh,nspin), Entropy
      real(dp), intent(out) ::  eo(no_u,nspin,nkpnt), qo(no_u,nspin,nkpnt)
      real(dp), intent(in)  ::  kpoint(3,nkpnt), qtot, temp, kweight(nkpnt), occtol,xijo(3,maxnh)

      
      integer :: mpirank, kcolrank, npGlobal
      integer :: npPerK, color, my_kpt_n
      
      integer :: kpt_comm, kpt_col_comm
      integer :: elsi_global_comm

      integer :: Global_Group, kpt_Group

      integer, allocatable :: ranks_in_world(:), ranks_in_world_AllK(:,:)
      type(distribution) :: dist_global
      type(distribution), allocatable :: dist_k(:)

      type(aux_matrix) :: pkg_global, pkg_k
      integer :: nvals
      real(dp), allocatable, target :: xijo_transp(:,:)
      real(dp), pointer :: my_xijo_transp(:,:) => null()
      real(dp), allocatable :: my_xij(:,:)

      integer :: my_no_l, my_nnz_l
      integer, pointer :: my_numh(:)
      integer, pointer :: my_listh(:)
      integer, pointer :: my_listhptr(:) => null()
      real(dp), pointer :: my_S(:)
      real(dp), pointer :: my_H(:,:)

      real(dp), pointer :: my_Escf(:,:) => null()
      real(dp), pointer :: my_Dscf(:,:) => null()
      real(dp), pointer :: my_Escf_reduced(:,:) => null()
      real(dp), pointer :: my_Dscf_reduced(:,:) => null()

      complex(dp), pointer :: DM_k(:,:)
      complex(dp), pointer :: EDM_k(:,:)
      complex(dp), allocatable :: Hk(:,:)
      complex(dp), allocatable :: Sk(:)

      integer :: iuo, j, ind, ind_u, ispin 
      real(dp) :: kxij, my_kpoint(3)
      complex(dp) :: kphs
      
      integer :: i, ih, ik, ierr

      integer, allocatable :: numh_u(:), listhptr_u(:), listh_u(:), ind2ind_u(:)
      integer :: nnz_u
      
      external die

      ! Split the global communicator
      ! Re-distribute H and S to the k-point (and spin) teams
      ! Generate Hk, Sk
      ! Call elsi_complex_solver
      ! Construct and re-distribute global DM and EDM

        ! Global communicator is a duplicate of passed communicator
      call MPI_Comm_Dup(MPI_Comm_DFT, elsi_global_comm, ierr)
      call MPI_Comm_Rank(elsi_global_comm, mpirank, ierr)
      call MPI_Comm_Size(elsi_global_comm, npGlobal, ierr)

      ! Split elsi_global_comm

      npPerK    = npGlobal/nkpnt

      ! Options for split: As many groups as nkpnt, so numbering is trivial
      color = mpirank/npPerK !  :  color 0: 0,1,2,3  ; color 1: 4,5,6,7
      call MPI_Comm_Split(elsi_global_comm, color, mpirank, kpt_Comm, ierr)
      my_kpt_n = 1 + color   !       1 + mpirank/npPerK
      ! Column communicator
      color = mod(mpirank, npPerK) ! :  color 0: 0,4  1: 1,5  2: 2,6  3: 3,7
      call MPI_Comm_Split(elsi_global_comm, color, mpirank, kpt_col_Comm, ierr) ! OK to use mpirank: just ordering
      call MPI_Comm_Rank(kpt_col_comm, kcolrank, ierr)

      print *, mpirank, "| ", "k-point ", my_kpt_n, " rank in col:", kcolrank
      
      ! Create distribution for k-point group
      call MPI_Comm_Group(elsi_global_comm, Global_Group, Ierr)
      call MPI_Comm_Group(kpt_Comm, kpt_Group, Ierr)
      allocate(ranks_in_world(npPerK))
      call MPI_Group_translate_ranks( kpt_Group, npPerK, &
           (/ (i,i=0,npPerK-1) /), &
           Global_Group, ranks_in_world, ierr )

      call MPI_Group_Free(kpt_Group, ierr)

      allocate (ranks_in_World_AllK(npPerK,nkpnt))
       ! We need everybody to have this information, as all nodes
       ! are part of the global distribution side of the communication
       ! (global_comm is the common-address space)

       ! Use the k-column communicator
       call MPI_AllGather(ranks_in_world,npPerK,MPI_integer,&
            Ranks_in_World_AllK(1,1),npPerK, &
            MPI_integer,kpt_col_Comm,ierr)

       allocate(dist_k(nkpnt))
       ! Create distributions known to all nodes  (again, those in base_comm)
       do ik = 1, nkpnt
          call newDistribution(dist_k(ik), elsi_global_Comm, &
               Ranks_in_World_AllK(:,ik),  &
               TYPE_BLOCK_CYCLIC, blockSize, "kpt dist")
       enddo
       deallocate(ranks_in_world,Ranks_in_World_AllK)
       call MPI_Barrier(elsi_global_Comm,ierr)

      call newDistribution(dist_global,elsi_global_Comm, (/ (i, i=0, npGlobal-1) /), &
           TYPE_BLOCK_CYCLIC,BlockSize,"global dist")
      call MPI_Barrier(elsi_global_Comm,ierr)

      ! Redistribute arrays

      ! No need to go sequentially over k-points
      ! Load pkg_global data package
          pkg_global%norbs = no_u
          pkg_global%no_l  = no_l
          pkg_global%nnzl  = maxnh
          pkg_global%numcols => numh  
          pkg_global%cols    => listh

          nvals = 1 + 3 + nspin   ! S, xijo (tranposed), H(:,nspin)
          allocate(pkg_global%vals(nvals))
          ! Link the vals items to the appropriate arrays (no extra memory here)
          pkg_global%vals(1)%data => S(:)
          call transpose(xijo,xijo_transp)
          do i = 1, 3
             pkg_global%vals(1+i)%data => xijo_transp(:,i)
          enddo
          do ispin = 1, nspin
             pkg_global%vals(4+ispin)%data => H(:,ispin)
          enddo

          call timer("redist_orbs_fwd", 1)
          call redistribute_spmatrix(no_u,pkg_global,dist_global, &
                                             pkg_k,dist_k(my_kpt_n),elsi_global_Comm)
          call timer("redist_orbs_fwd", 2)

          !------------------------------------------
          ! Unpack info: real S and H (and index arrays) distributed over each kpt_comm
          my_no_l   = pkg_k%no_l
          my_nnz_l  = pkg_k%nnzl
          my_numh => pkg_k%numcols

          my_listh => pkg_k%cols
          my_S => pkg_k%vals(1)%data

          do i = 1, 3
             my_xijo_transp(1:my_nnz_l,i:i) => pkg_k%vals(1+i)%data
          enddo

          do ispin = 1, nspin
             my_H(1:my_nnz_l,ispin:ispin) => pkg_k%vals(4+ispin)%data
          enddo
          !------------------------------------------
          
          ! Now we could clear pkg_k--
           !  nullify(pkg_k%numcols)
           !  nullify(pkg_k%cols)
           !  deallocate(pkg_k%vals)
             !
             ! but we need to remember to deallocate the actual arrays after use
             ! it is probably safer to keep the pkg references
          !---------------------------

          ! Clean pkg_global -- This is safe, as we use pointers to data only
          do i = 1, size(pkg_global%vals)
             nullify(pkg_global%vals(i)%data)
          enddo
          deallocate(pkg_global%vals)
          nullify(pkg_global%numcols)   
          nullify(pkg_global%cols)      

          deallocate(xijo_transp)  ! Auxiliary array used for sending
          
          ! generate listhptr for folding/unfolding operations
          call re_alloc(my_listhptr,1,my_no_l,"my_listhptr","elsi_solver")
          my_listhptr(1) = 0
          do ih = 2,my_no_l
             my_listhptr(ih) = my_listhptr(ih-1) + my_numh(ih-1)
          enddo

          ! The folded variables (*_u) are local to each team
                
          allocate(numh_u(my_no_l), listhptr_u(my_no_l))
          call fold_sparse_arrays(my_no_l,no_u, &
                              my_numh,my_listhptr,my_nnz_l,my_listh, &
                              numh_u,listhptr_u,nnz_u,listh_u,ind2ind_u)

      ! 
      call transpose(my_xijo_transp,my_xij)
      my_kpoint(:) = kpoint(:,my_kpt_n)
      print *, mpirank, "| ", "Doing k-point ", my_kpt_n
      
      allocate(Hk(nnz_u,nspin), Sk(nnz_u))
      Sk = 0
      Hk = 0
      do iuo = 1,my_no_l
         do j = 1,my_numh(iuo)
            ind = my_listhptr(iuo) + j
            ind_u = ind2ind_u(ind)

            kxij = my_kpoint(1) * my_xij(1,ind) +    &
                   my_kpoint(2) * my_xij(2,ind) +    &
                   my_kpoint(3) * my_xij(3,ind) 
            kphs = cdexp(dcmplx(0.0_dp, -1.0_dp)*kxij)

            Sk(ind_u) = Sk(ind_u) + my_S(ind) *kphs
            do ispin = 1, nspin
               Hk(ind_u,ispin) = Hk(ind_u,ispin) + my_H(ind,ispin) *kphs
            enddo
         enddo
      enddo
            
      ! Deallocate un-needed bits of pkg_k  (my_S, my_H, my_xijo_transp)
      do i = 1, size(pkg_k%vals)
         call de_alloc(pkg_k%vals(i)%data,"pkg_k%vals%data","kpoints_dispatcher")
      enddo
      deallocate(pkg_k%vals)
      !       

      ! Prepare arrays for holding results
      call re_alloc(DM_k,1,nnz_u,1,nspin,"DM_k","elsi_solver")
      call re_alloc(EDM_k,1,nnz_u,1,nspin,"EDM_k","elsi_solver")

      call elsi_complex_solver(iscf, no_u, my_no_l, nspin, nnz_u, numh_u, listhptr_u, &
                               listh_u, qtot, temp, Hk, Sk, DM_k, EDM_k, Ef, Entropy,  &
                               nkpnt, my_kpt_n, kpoint(:,my_kpt_n), kweight(my_kpt_n),    &
                               elsi_global_comm, kpt_Comm )
      
      deallocate(listhptr_u, numh_u, listh_u)
      deallocate(Hk,Sk)

      ! Re-create DM and EDM
      ! Possibly reduce entropy
      ! Unfold within a given k

      ! Prepare arrays for holding results: Note sizes: these are folded out
      call re_alloc(my_Dscf,1,my_nnz_l,1,nspin,"my_Dscf","elsi_solver")
      call re_alloc(my_Escf,1,my_nnz_l,1,nspin,"my_Escf","elsi_solver")

      do iuo = 1, my_no_l
         do j = 1, my_numh(iuo)
            ind = my_listhptr(iuo) + j
            ind_u = ind2ind_u(ind)

            kxij = my_kpoint(1) * my_xij(1,ind) +    &
                 my_kpoint(2) * my_xij(2,ind) +    &
                 my_kpoint(3) * my_xij(3,ind) 
            kphs = cdexp(dcmplx(0.0_dp, +1.0_dp)*kxij)  ! ** CHECK phase sign

            do ispin = 1, nspin
               my_Dscf(ind,ispin) = real( DM_k(ind_u,ispin) * kphs )
               my_Escf(ind,ispin) = real( EDM_k(ind_u,ispin) * kphs )
            enddo
         enddo
      enddo

      deallocate(my_xij)
      deallocate(ind2ind_u)
      call de_alloc(DM_k,"DM_k","elsi_solver")
      call de_alloc(EDM_k,"EDM_k","elsi_solver")

      if (my_kpt_n == 1 ) then
         ! Prepare arrays for holding reduced data
         call re_alloc(my_Dscf_reduced,1,my_nnz_l,1,nspin,"my_Dscf_reduced","elsi_solver")
         call re_alloc(my_Escf_reduced,1,my_nnz_l,1,nspin,"my_Escf_reduced","elsi_solver")
         if (kcolrank /= 0) call die("Rank 0 in kpt_comm not doing kpt 1") 
      endif

      ! Use k-point column communicator, and reduce to rank 0,
      ! which *should* correspond to kpt=1... (checked above)
      call MPI_Reduce( my_Dscf, my_Dscf_reduced, nspin*my_nnz_l, MPI_Double_Precision, &
           MPI_Sum, 0, kpt_col_comm, ierr )
      call MPI_Reduce( my_Escf, my_Escf_reduced, nspin*my_nnz_l, MPI_Double_Precision, &
           MPI_Sum, 0, kpt_col_comm, ierr )

      call de_alloc(my_Dscf,"my_Dscf")
      call de_alloc(my_Escf,"my_Escf")
      
      ! redistribute to global distribution, only from the first k-point

      if (my_kpt_n == 1) then

         ! This bits are still there *** update if we ever clean after use
         ! pkg_k%norbs = no_u
         ! pkg_k%no_l  = my_no_l
         ! pkg_k%nnzl  = my_nnz_l
         ! pkg_k%numcols => my_numh  
         ! pkg_k%cols    => my_listh

         nvals = 2*nspin   ! DM, EDM
         allocate(pkg_k%vals(nvals))
         do ispin = 1, nspin
            pkg_k%vals(ispin)%data => my_Dscf_reduced(:,ispin)
            pkg_k%vals(nspin+ispin)%data => my_Escf_reduced(:,ispin)
         enddo

      endif

      ! Everybody participates in the transfer in the receiving side, but
      ! only kpt=1 in the sending side, because we are using dist_k(1)
      call timer("redist_dm-edm_bck", 1)
      call redistribute_spmatrix(no_u,pkg_k,dist_k(1), &
           pkg_global,dist_global,elsi_global_Comm)
      call timer("redist_dm-edm_bck", 2)

      call de_alloc(my_Dscf_reduced,"my_Dscf_reduced")
      call de_alloc(my_Escf_reduced,"my_Escf_reduced")
      !Clean pkg_k in sender
      if (my_kpt_n == 1) then
         deallocate(pkg_k%vals)   ! just pointers
      endif

      !Clean all pkg_k: These were actually allocated in the forward transfer
      call de_alloc(pkg_k%numcols,"pkg_k%numcols")
      call de_alloc(pkg_k%cols,"pkg_k%cols")

      ! Unpack data (all processors, original global distribution)
        ! In future, pkg_global%vals(:) could be pointing to DM and EDM,
        ! and the 'redistribute' routine check whether the vals arrays are
        ! associated, to use them instead of allocating them.
         do ispin = 1, nspin
            Dscf(:,ispin)  = pkg_global%vals(ispin)%data(:)    
            Escf(:,ispin) = pkg_global%vals(nspin+ispin)%data(:)    
         enddo
        ! Check no_l
        if (no_l /= pkg_global%no_l) then
           call die("Mismatch in no_l at end of kpoints-dispatcher")
        endif
        ! Check listH
        if (any(listh(:) /= pkg_global%cols(:))) then
           call die("Mismatch in listH at end of kpoints-dispatcher")
        endif

        ! Clean pkg_global
        call de_alloc(pkg_global%numcols,"pkg_global%numcols","elsi_solver") 
        call de_alloc(pkg_global%cols,   "pkg_global%cols",   "elsi_solver")
        do i = 1, size(pkg_global%vals)
           call de_alloc(pkg_global%vals(i)%data,"pkg_global%vals%data","kpoints_dispatcher")
        enddo
        deallocate(pkg_global%vals)

        ! Possible reduction of entropy?
        
     call MPI_Comm_Free(kpt_comm, ierr)
     call MPI_Comm_Free(kpt_col_comm, ierr)
                  
end subroutine elsi_kpoints_dispatcher

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

! To-Do: Check the implications of distinguishing between 'global_comm'
! (common space for communications) and 'base_comm' (to be split)
! For the real case they can be the same
!
! This routine should be called by everybody in global_comm

subroutine get_spin_comms_and_dists(global_comm, base_comm, blocksize, n_spin,&
     dist_global, dist_spin, spatial_comm, spin_comm)

  use mpi_siesta
  use class_Distribution    ! distribution, newDistribution, types
  
  integer, intent(in)  :: global_comm    ! Communicator for global distribution
                                         ! Needed to include the actual global ranks
                                         ! in the distribution objects

  integer, intent(in)  :: base_comm      ! Communicator to be split
  integer, intent(out) :: spatial_comm   ! "Orbital" comm (one for each spin team)
  integer, intent(out) :: spin_comm      ! "Inner" spin comm (for reductions)

  integer, intent(in)  :: blocksize
  integer, intent(in)  :: n_spin              ! should be 2 in this version

  ! Note inout for bud objects
  type(distribution), intent(inout)  :: dist_global         ! global distribution object
  type(distribution), intent(inout)  :: dist_spin(2) ! per-spin objects

  ! Local variables
  integer :: color, base_rank
  integer :: npGlobal, npPerSpin
  integer, allocatable, dimension(:) :: global_ranks_in_world
  integer, allocatable, dimension(:) :: ranks_in_world  ! should be 'spatial_ranks...'
  integer, allocatable, dimension(:,:) :: ranks_in_World_Spin
  integer :: global_group, Spatial_group
  integer :: i, ispin, ierr
  
       call MPI_Comm_Size(global_comm, npGlobal, ierr)
         ! This distribution could be created by the caller...
       allocate(global_ranks_in_world(npGlobal))
       global_ranks_in_world = (/ (i, i=0, npGlobal-1) /)
       call newDistribution(dist_global,global_Comm,global_ranks_in_world, &
            TYPE_BLOCK_CYCLIC,BlockSize,"global dist")
       deallocate(global_ranks_in_world)
       call MPI_Barrier(global_Comm,ierr)

       ! Split the base communicator
       call MPI_Comm_Rank(base_comm, base_rank, ierr)
       ! "Row" communicator for independent operations on each spin
       ! The name refers to "spatial" degrees of freedom.
       color = mod(base_rank,n_spin)    ! {0,1} for n_spin = 2, or {0} for n_spin = 1
       call MPI_Comm_Split(base_comm, color, base_rank, Spatial_Comm, ierr)
       ! "Column" communicator for spin reductions
       color = base_rank/n_spin
       call MPI_Comm_Split(base_comm, color, base_rank, Spin_Comm, ierr)

       call MPI_Comm_Size(spatial_comm, npPerSpin, ierr)   ! number of procs in each spin team

       call MPI_Comm_Group(Spatial_Comm, Spatial_Group, Ierr)
       allocate(ranks_in_world(npPerSpin))
       call MPI_Comm_Group(global_comm, Global_Group, Ierr)
       call MPI_Group_translate_ranks( Spatial_Group, npPerSpin, &
            (/ (i,i=0,npPerSpin-1) /), &
            Global_Group, ranks_in_world, ierr )

       call MPI_Group_Free(Spatial_Group, ierr)
       call MPI_Group_Free(Global_Group, ierr)
       
       allocate (ranks_in_World_Spin(npPerSpin,n_spin))
       ! We need everybody to have this information, as all nodes
       ! are part of the global distribution side of the communication
       ! (global_comm is the common-address space)

       ! But note that this will tell only those in base_comm...
       ! This routine might not be the end of the story for cascading splits
       ! (k-points and spin)
       ! We might need an extra 'broadcast' over the k-point 'column' comm.
       call MPI_AllGather(ranks_in_world,npPerSpin,MPI_integer,&
            Ranks_in_World_Spin(1,1),npPerSpin, &
            MPI_integer,Spin_Comm,ierr)

       ! Create distributions known to all nodes  (again, those in base_comm)
       do ispin = 1, n_spin
          call newDistribution(dist_spin(ispin), global_Comm, &
               Ranks_in_World_Spin(:,ispin),  &
               TYPE_BLOCK_CYCLIC, blockSize, "SPATIAL dist")
       enddo
       deallocate(ranks_in_world,Ranks_in_World_Spin)
       call MPI_Barrier(global_Comm,ierr)
       
end subroutine get_spin_comms_and_dists

#endif  /* SIESTA__ELSI */

subroutine elsi_complex_solver(iscf, n_basis, n_basis_l, n_spin, nnz_l, numh, row_ptr, &
     col_idx, qtot, temp, ham, ovlp, dm, edm, ef, ets, &
     nkpnt, kpt_n, kpt, weight, elsi_global_comm, kpt_comm)

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
  integer,  intent(in), target    :: numh(n_basis_l)
  integer,  intent(in)    :: row_ptr(n_basis_l)
  integer,  intent(in), target    :: col_idx(nnz_l)
  real(dp), intent(in)    :: qtot
  real(dp), intent(in)    :: temp
  complex(dp), intent(inout), target :: ham(nnz_l,n_spin)
  complex(dp), intent(inout), target :: ovlp(nnz_l)
  complex(dp), intent(out)   :: dm(nnz_l,n_spin)
  complex(dp), intent(out)   :: edm(nnz_l,n_spin)
  real(dp), intent(out)   :: ef        ! Fermi energy
  real(dp), intent(out)   :: ets       ! Entropy/k, dimensionless
  integer,  intent(in)    :: nkpnt     ! number of k-points
  integer,  intent(in)    :: kpt_n
  real(dp), intent(in)    :: kpt(3:)
  real(dp), intent(in)    :: weight
  integer,  intent(in)    :: elsi_global_comm
  integer,  intent(in)    :: kpt_comm
  
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

  integer :: elsi_Spatial_comm, elsi_Spin_comm

  type(distribution) :: dist_global
  type(distribution) :: dist_spin(2)
  
  type(aux_matrix) :: pkg_global, pkg_spin   ! Packages for transfer

  integer :: my_spin

  integer :: my_no_l  
  integer :: my_nnz_l 
  integer :: my_nnz
  integer, pointer  :: my_row_ptr2(:) => null()
  integer  :: i, ih, ispin, spin_rank
  real(dp) :: ets_spin

  integer, pointer  :: my_col_idx(:)
  complex(dp), pointer :: my_S(:)
  complex(dp), pointer :: my_H(:)
  complex(dp), pointer :: my_DM(:) => null()
  complex(dp), pointer :: my_EDM(:) => null()
  
  external :: timer

#ifndef MPI
  call die("This ELSI solver interface needs MPI")
#endif

  call timer("elsi-complex-solver", 1)

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
       call globalize_sum(nnz_l, nnz_g, comm=kpt_comm)

       allocate(row_ptr2(n_basis_l+1))
       row_ptr2(1:n_basis_l) = row_ptr(1:n_basis_l)+1
       row_ptr2(n_basis_l+1) = nnz_l+1

       call elsi_set_csc(elsi_h, nnz_g, nnz_l, n_basis_l, col_idx, row_ptr2)
       deallocate(row_ptr2)

       call elsi_set_csc_blk(elsi_h, BlockSize)
       call elsi_set_kpoint(elsi_h, nkpnt, kpt_n, weight)
       call elsi_set_mpi(elsi_h, kpt_comm)
       call elsi_set_mpi_global(elsi_h, elsi_global_comm)

    else
       
       ! MPI logic for spin polarization

       ! Split the communicator in spins and get distribution objects
       ! for the data redistribution needed
       ! Note that dist_spin is an array
       call get_spin_comms_and_dists(elsi_global_comm,kpt_comm, &  !! **** kpt_comm as global?
            blocksize, n_spin, &
            dist_global,dist_spin, elsi_spatial_comm, elsi_spin_comm)

       ! Find out which spin team we are in, and tag the spin we work on
       call mpi_comm_rank( elsi_Spin_Comm, spin_rank, ierr )
       my_spin = spin_rank+1  ! {1,2}

       
       ! This is done serially, each time filling one spin set
       ! Note that **all processes** need to have the same pkg_global
       
       do ispin = 1, n_spin

          ! Load pkg_global data package
          pkg_global%norbs = n_basis
          pkg_global%no_l  = n_basis_l
          pkg_global%nnzl  = nnz_l
          pkg_global%numcols => numh  
          pkg_global%cols    => col_idx

          allocate(pkg_global%complex_vals(2))
          ! Link the vals items to the appropriate arrays (no extra memory here)
          pkg_global%complex_vals(1)%data => ovlp(:)
          ! Note that we *cannot* say  => ham(:,my_spin)
          ! and avoid the sequential loop, as then half the processors will send
          ! the information for 'spin up' and the other half the information for 'spin down',
          ! which is *not* what we want.
          pkg_global%complex_vals(2)%data => ham(:,ispin)

          call timer("redist_orbs_fwd", 1)

          ! We are doing the transfers sequentially. One spin team is
          ! 'idle' (in the receiving side) in each pass, as the dist_spin(ispin) distribution
          ! does not involve them.
          
          call redistribute_spmatrix(n_basis,pkg_global,dist_global, &
                                             pkg_spin,dist_spin(ispin),elsi_global_Comm)
   
          call timer("redist_orbs_fwd", 2)

          if (my_spin == ispin) then  ! Each team gets their own data

             !nrows = pkg_spin%norbs          ! or simply 'norbs'
             my_no_l = pkg_spin%no_l
             my_nnz_l    = pkg_spin%nnzl
             call MPI_AllReduce(my_nnz_l,my_nnz,1,MPI_integer,MPI_sum,elsi_Spatial_Comm,ierr)
             ! generate off-by-one row pointer
             call re_alloc(my_row_ptr2,1,my_no_l+1,"my_row_ptr2","elsi_solver")
             my_row_ptr2(1) = 1
             do ih = 1,my_no_l
                my_row_ptr2(ih+1) = my_row_ptr2(ih) + pkg_spin%numcols(ih)
             enddo

             my_col_idx => pkg_spin%cols
             my_S => pkg_spin%complex_vals(1)%data
             my_H => pkg_spin%complex_vals(2)%data

             call re_alloc(my_DM,1,my_nnz_l,"my_DM","elsi_solver")
             call re_alloc(my_EDM,1,my_nnz_l,"my_EDM","elsi_solver")
          endif

          ! Clean pkg_global
          nullify(pkg_global%complex_vals(1)%data)  
          nullify(pkg_global%complex_vals(2)%data)  
          deallocate(pkg_global%complex_vals)
          nullify(pkg_global%numcols)   
          nullify(pkg_global%cols)      

       enddo

       call elsi_set_csc(elsi_h, my_nnz, my_nnz_l, my_no_l, my_col_idx, my_row_ptr2)
       call de_alloc(my_row_ptr2,"my_row_ptr2","elsi_solver")
       
       call elsi_set_csc_blk(elsi_h, BlockSize)
       call elsi_set_spin(elsi_h, n_spin, my_spin)
       call elsi_set_kpoint(elsi_h, nkpnt, kpt_n, weight)
       call elsi_set_mpi(elsi_h, elsi_Spatial_comm)
       call elsi_set_mpi_global(elsi_h, elsi_global_comm)

    endif  ! n_spin
    
  call timer("elsi-solver", 1)

  if (n_spin == 1) then
     call elsi_dm_complex_sparse(elsi_h, ham, ovlp, DM, energy)
     call elsi_get_edm_complex_sparse(elsi_h, EDM)
     call elsi_get_entropy(elsi_h, ets)
  else
     ! Solve DM, and get (at every step for now) EDM, Fermi energy, and entropy
     ! Energy is already summed over spins
     call elsi_dm_complex_sparse(elsi_h, my_H, my_S, my_DM, energy)
     call elsi_get_edm_complex_sparse(elsi_h, my_EDM)
     !... but we still need to sum the entropy over spins
     call elsi_get_entropy(elsi_h, ets_spin)
     call globalize_sum(ets_spin, ets, comm=elsi_Spin_comm)
     ! And over kpoints??
  endif

  call elsi_get_mu(elsi_h, ef)
  ets = ets/temp
  
  ! Ef, energy, and ets are known to all nodes

  call timer("elsi-solver", 2)

  if ( n_spin == 2) then
     ! Now we need to redistribute back

     do ispin = 1, n_spin

        if (my_spin == ispin) then
           ! Prepare pkg_spin to transfer the right spin information
           ! The other fields (numcols, cols) are the same and are still there
           ! Deallocate my_S and my_H
           call de_alloc(pkg_spin%complex_vals(1)%data,"pkg_spin%vals(1)%data","elsi_solver")
           call de_alloc(pkg_spin%complex_vals(2)%data,"pkg_spin%vals(2)%data","elsi_solver")

           pkg_spin%complex_vals(1)%data => my_DM(1:my_nnz_l)
           pkg_spin%complex_vals(2)%data => my_EDM(1:my_nnz_l)

        endif

        ! pkg_global is clean now
        call timer("redist_orbs_bck", 1)
        call redistribute_spmatrix(n_basis,pkg_spin,dist_spin(ispin) &
                                          ,pkg_global,dist_global,elsi_global_Comm)
        call timer("redist_orbs_bck", 2)

        ! Clean pkg_spin
        if (my_spin == ispin) then
           ! Each team deallocates during "its" spin cycle
           call de_alloc(my_DM, "my_DM", "elsi_solver")
           call de_alloc(my_EDM,"my_EDM","elsi_solver")

           nullify(pkg_spin%complex_vals(1)%data)    ! formerly pointing to DM
           nullify(pkg_spin%complex_vals(2)%data)    ! formerly pointing to EDM
           deallocate(pkg_spin%complex_vals)
           ! allocated in the direct transfer
           call de_alloc(pkg_spin%numcols,"pkg_spin%numcols","elsi_solver")
           call de_alloc(pkg_spin%cols,   "pkg_spin%cols",   "elsi_solver")
        endif


        ! In future, pkg_global%vals(1,2) could be pointing to DM and EDM,
        ! and the 'redistribute' routine check whether the vals arrays are
        ! associated, to use them instead of allocating them.
        DM(:,ispin)  = pkg_global%complex_vals(1)%data(:)    
        EDM(:,ispin) = pkg_global%complex_vals(2)%data(:)    
        ! Check no_l
        if (n_basis_l /= pkg_global%no_l) then
           call die("Mismatch in no_l")
        endif
        ! Check listH
        if (any(col_idx(:) /= pkg_global%cols(:))) then
           call die("Mismatch in listH")
        endif

        ! Clean pkg_global
        ! allocated by the transfer routine, but we did not actually
        ! look at them
        call de_alloc(pkg_global%numcols,"pkg_global%numcols","elsi_solver") 
        call de_alloc(pkg_global%cols,   "pkg_global%cols",   "elsi_solver")

        call de_alloc(pkg_global%complex_vals(1)%data,"pkg_global%vals(1)%data","elsi_solver")
        call de_alloc(pkg_global%complex_vals(2)%data,"pkg_global%vals(2)%data","elsi_solver")
        deallocate(pkg_global%complex_vals)

     enddo

     call MPI_Comm_Free(elsi_Spatial_comm, ierr)
     call MPI_Comm_Free(elsi_Spin_comm, ierr)
     
  endif

  call timer("elsi-complex-solver", 2)

  
end subroutine elsi_complex_solver

subroutine transpose(a,b)
  real(dp), intent(in) :: a(:,:)
  real(dp), allocatable, intent(out) :: b(:,:)

  integer n, m, i, j
  
  n = size(a,1)
  m = size(a,2)
  allocate(b(m,n))
  do i = 1, n
     do j = 1, m
        b(j, i) = a(i, j)
     enddo
  enddo
end subroutine transpose

end module m_elsi_interface



