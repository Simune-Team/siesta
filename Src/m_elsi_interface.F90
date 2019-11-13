!
! Alberto Garcia, 2018-
! with help by Victor Yu and Victor M. Garcia-Suarez
!
! ELSI DM-based interface to Siesta. It uses the sparse matrices from Siesta,
! and obtains the DM (and optionally the EDM) matrices in sparse form.
!
! This interface does not generate eigenvalues nor eigenvectors, even if
! a diagonalization-based ELSI solver (i.e., ELPA) is used.
!
! The elsi_getdm routine is in principle able to perform (spin-polarized)
! calculations for real matrices (i.e., at the Gamma point), and periodic
! calculations for complex matrices (i.e., multiple k-points), including spin.
!
! The structure of the solver routine is such that it will detect when it is
! called for the first scf step, so it can perform any needed initialization.
!
! The module also exports the "elsi_finalize_scfloop" routine, to be called from
! the appropriate place. Some variables are kept at the module level for this.
!
! This interface has been tested with ELSI-v2.0.2 --> 2.3.1 
!
! Usage: Compile Siesta with -DSIESTA__ELSI
!        Define
!           SolutionMethod ELSI
!        in the fdf file

module m_elsi_interface

#if SIESTA__ELSI

  use precision, only: dp
  use parallel, only: ionode
  use units, only:    eV
  use elsi

  implicit none

  private

  integer, parameter :: ELPA_SOLVER       = 1 ! solver
  integer, parameter :: OMM_SOLVER        = 2 ! solver
  integer, parameter :: PEXSI_SOLVER      = 3 ! solver
  integer, parameter :: SIPS_SOLVER       = 5 ! solver
  integer, parameter :: NTPOLY_SOLVER     = 6 ! solver

  integer, parameter :: MULTI_PROC        = 1 ! parallel_mode
  integer, parameter :: SIESTA_CSC        = 2 ! distribution

  integer, parameter :: GAUSSIAN          = 0 ! broadening
  integer, parameter :: FERMI             = 1 ! broadening
  integer, parameter :: METHFESSEL_PAXTON = 2 ! broadening
  integer, parameter :: CUBIC             = 3 ! broadening
  integer, parameter :: COLD              = 4 ! broadening

  integer, parameter :: ELSI_NOT_SET      = -910910

  type(elsi_handle) :: elsi_h

  integer :: elsi_global_comm    ! Used by all routines. Freed at end of scf loop

  integer :: which_solver
  integer :: which_broad
  integer :: out_level
  integer :: out_json
  integer :: mp_order
  integer :: illcond_check
  real(dp) :: illcond_tol
  
  integer :: elpa_flavor
  integer :: elpa_gpu
  integer :: elpa_n_single
  integer :: elpa_autotune

  integer  :: omm_flavor
  integer  :: omm_n_elpa
  real(dp) :: omm_tol
  
  integer :: pexsi_tasks_per_pole
  integer :: pexsi_n_pole
  integer :: pexsi_n_mu
  integer :: pexsi_tasks_symbolic
  real(dp) :: pexsi_inertia_tol
  real(dp) :: pexsi_initial_mu_min
  real(dp) :: pexsi_initial_mu_max

  integer :: sips_n_slice
  integer :: sips_n_elpa

  integer  :: ntpoly_method
  real(dp) :: ntpoly_filter
  real(dp) :: ntpoly_tol

  character(len=6) :: solver_string
  character(len=5) :: broad_string

  real(dp), allocatable :: v_old(:,:)    ! For mu update in PEXSI solver
  real(dp), allocatable :: delta_v(:,:)   ! For mu update in PEXSI solver
  real(dp)  :: dv_min, dv_max             ! For mu update in PEXSI solver
  real(dp)  :: mu_min, mu_max

  public :: elsi_getdm
  public :: elsi_finalize_scfloop
  public :: elsi_save_potential

CONTAINS

subroutine elsi_getdm(iscf, no_s, nspin, no_l, maxnh, no_u,  &
     numh, listhptr, listh, H, S, qtot, temp, &
     xijo, nkpnt, kpoint, kweight,    &
     Dscf, ef, Entropy, occtol, neigwanted, Get_EDM_Only)

  !
  ! Analogous to 'diagon', it dispatches ELSI solver routines as needed
  !

  use m_fold_auxcell, only: fold_sparse_arrays ! Could be called in state_init


     real(dp), intent(inout) :: H(:,:), S(:)    ! Note: we might overwrite these
     integer, intent(in) ::  iscf, maxnh, no_u, no_l, no_s, nkpnt
     integer, intent(in) ::  neigwanted, nspin
     integer, intent(in) ::  listh(maxnh), numh(no_l), listhptr(no_l)
     real(dp), intent(in)  ::  kpoint(3,nkpnt), qtot, temp, kweight(nkpnt), occtol,xijo(3,maxnh)

      real(dp), intent(out) ::  Dscf(maxnh,nspin), ef, Entropy

      ! Interim flag to just get the EDM
      ! Note that, if .true., the DM is NOT obtained
      ! This needs to be refactored
      
      logical, intent(in) :: Get_EDM_Only 


      logical :: gamma, using_aux_cell
      integer, allocatable :: numh_u(:), listhptr_u(:), listh_u(:)
      integer, allocatable :: ind2ind_u(:)

      !
      real(dp), allocatable, dimension(:,:)  :: Dscf_u, H_u
      real(dp), allocatable, dimension(:)    :: S_u

      integer :: iuo, ispin, j, ind, ind_u, nnz_u

      external die

      gamma = ((nkpnt == 1) .and. (sum(abs(kpoint(:,1))) == 0.0_dp))
      using_aux_cell =  (no_s /= no_u)

      if (gamma) then

         if  (.not. using_aux_cell) then

            call elsi_real_solver(iscf, no_u, no_l, nspin, &
                     maxnh, listhptr, listh, qtot, temp, &
                     H, S, Dscf, ef, Entropy, Get_EDM_Only)

         else

            ! We fold-in the sparse arrays so that they add up all the
            ! matrix elements whose 'column' supercell index maps to the same
            ! unit-cell column. Note that they are still sparse.

            ! First, determine the appropriate index arrays (ending in _u), and
            ! a mapper ind2ind_u

            allocate(numh_u(no_l), listhptr_u(no_l))
            call fold_sparse_arrays(no_l,no_u,  &
                                    numh,listhptr,maxnh,listh, &
                                    numh_u,listhptr_u,nnz_u,listh_u,ind2ind_u)

            allocate(H_u(nnz_u,nspin), S_u(nnz_u))
            S_u = 0
            H_u = 0
               do iuo = 1,no_l
                  do j = 1,numh(iuo)
                     ind = listhptr(iuo) + j
                     ind_u = ind2ind_u(ind)
                     S_u(ind_u) = S_u(ind_u) + S(ind)
                     do ispin = 1, nspin
                        H_u(ind_u,ispin) = H_u(ind_u,ispin) + H(ind,ispin)
                     enddo
                  enddo
               enddo

            ! We can now call the standard real solver routine
            allocate(Dscf_u(nnz_u,nspin))

            call elsi_real_solver(iscf, no_u, no_l, nspin, &
                     nnz_u, listhptr_u, listh_u, qtot, temp, &
                     H_u, S_u, Dscf_u, ef, Entropy, Get_EDM_Only)

            deallocate(H_u, S_u)
            deallocate(numh_u, listhptr_u, listh_u)


            ! Unfold. We put the same '_u' DM entry into all image slots.
            do iuo = 1,no_l
               do j = 1,numh(iuo)
                  ind = listhptr(iuo) + j
                  ind_u = ind2ind_u(ind)
                  do ispin = 1, nspin
                     Dscf(ind,ispin) = Dscf_u(ind_u,ispin)
                  enddo
               enddo
            enddo
            deallocate(Dscf_u)
            deallocate(ind2ind_u)

         endif  ! using auxiliary supercell with Gamma sampling

      else

         ! We need more preparation
         call elsi_kpoints_dispatcher(iscf, no_s, nspin, no_l, maxnh, no_u,  &
              numh, listhptr, listh, H, S, qtot, temp, &
              xijo, nkpnt, kpoint, kweight,    &
              Dscf, ef, Entropy, occtol, neigwanted, Get_EDM_Only)

      endif

end subroutine elsi_getdm

subroutine elsi_get_opts()

  use fdf,         only: fdf_get

  solver_string        = fdf_get("ELSI-Solver", "elpa")
  out_level            = fdf_get("ELSI-Output-Level", 0)
  out_json             = fdf_get("ELSI-Output-Json", 1)
  broad_string         = fdf_get("ELSI-Broadening-Method", "fermi")
  mp_order             = fdf_get("ELSI-Broadening-MPOrder", 1)
  illcond_check        = fdf_get("ELSI-Ill-Condition-Check", 0 )
  illcond_tol          = fdf_get("ELSI-Ill-Condition-Tolerance", 1.0e-5_dp )
  
  elpa_flavor          = fdf_get("ELSI-ELPA-Flavor", 2)
  elpa_gpu             = fdf_get("ELSI-ELPA-GPU", 0)
  elpa_n_single        = fdf_get("ELSI-ELPA-N-single-precision", 0)
  elpa_autotune        = fdf_get("ELSI-ELPA-Autotune", 0)

  omm_flavor           = fdf_get("ELSI-OMM-Flavor", 0)
  omm_n_elpa           = fdf_get("ELSI-OMM-ELPA-Steps", 3)
  omm_tol              = fdf_get("ELSI-OMM-Tolerance", 1.0e-9_dp)

  pexsi_tasks_per_pole = fdf_get("ELSI-PEXSI-Tasks-Per-Pole", ELSI_NOT_SET)
  pexsi_tasks_symbolic = fdf_get("ELSI-PEXSI-Tasks-Symbolic", 1)
  pexsi_n_pole         = fdf_get("ELSI-PEXSI-Number-Of-Poles", 20)
  pexsi_n_mu           = fdf_get("ELSI-PEXSI-Number-Of-Mu-Points", 2)
  pexsi_inertia_tol    = fdf_get("ELSI-PEXSI-Inertia-Tolerance", 0.05_dp)
  pexsi_initial_mu_min = fdf_get("ELSI-PEXSI-Initial-Mu-Min", -1.0_dp, 'Ry')
  pexsi_initial_mu_max = fdf_get("ELSI-PEXSI-Initial-Mu-Max", 0.0_dp, 'Ry')


  sips_n_slice         = fdf_get("ELSI-SIPS-Slices", ELSI_NOT_SET)
  sips_n_elpa          = fdf_get("ELSI-SIPS-ELPA-Steps", 2)

  ntpoly_method        = fdf_get("ELSI-NTPOLY-Method", 2)
  ntpoly_filter        = fdf_get("ELSI-NTPOLY-Filter", 1.0e-9_dp)
  ntpoly_tol           = fdf_get("ELSI-NTPOLY-Tolerance", 1.0e-6_dp)

  select case (solver_string)
  case ("elpa", "ELPA")
    which_solver = ELPA_SOLVER
  case ("omm", "OMM")
    which_solver = OMM_SOLVER
  case ("pexsi", "PEXSI")
    which_solver = PEXSI_SOLVER
  case ("sips", "SIPS", "SIPs")
    which_solver = SIPS_SOLVER
  case ("ntpoly", "NTPOLY", "NTPoly")
    which_solver = NTPOLY_SOLVER
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

end subroutine elsi_get_opts

!======================================================================================
! ELSI takes 1D block-cyclic distributed CSC/CSR matrices as its
! input/output.
!
! But note taht this version assumes *the same* distributions for Siesta (setup_H et al) and ELSI
! operations.
!
subroutine elsi_real_solver(iscf, n_basis, n_basis_l, n_spin, nnz_l, row_ptr, &
  col_idx, qtot, temp, ham, ovlp, DM, ef, ets, Get_EDM_Only)

  use fdf,         only: fdf_get
  use m_mpi_utils, only: globalize_sum
  use parallel,    only: BlockSize
#ifdef MPI
  use mpi_siesta
#endif
  use class_Distribution
  use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
  use alloc, only: de_alloc  ! To deallocate some pointers in transfer matrices

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
  real(dp), intent(out)   :: DM(nnz_l,n_spin)    ! It can be the DM or the EDM
  real(dp), intent(out)   :: ef        ! Fermi energy
  real(dp), intent(out)   :: ets       ! Entropy/k, dimensionless

  integer :: ierr
  integer :: n_state
  integer :: nnz_g

  real(dp) :: energy

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
  integer, allocatable :: my_row_ptr2(:)
  integer  :: i, ih, ispin, spin_rank

  integer, pointer  :: my_col_idx(:)
  real(dp), pointer :: my_S(:)
  real(dp), pointer :: my_H(:)
  real(dp), allocatable, target :: my_DM(:) 

  integer :: date_stamp

  logical :: Get_EDM_Only

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
    call elsi_get_opts()
      
    ! Number of states to solve when calling an eigensolver
    n_state = min(n_basis, n_basis/2+5)

    ! Now we have all ingredients to initialize ELSI
    call elsi_init(elsi_h, which_solver, MULTI_PROC, SIESTA_CSC, n_basis, &
      qtot, n_state)

    ! Output
    if (ionode) then
      call elsi_set_output(elsi_h, out_level)
      call elsi_set_output_log(elsi_h, out_json)
      call elsi_set_write_unit(elsi_h, 6)
    endif

    ! Possible ill-conditioning of S
    call elsi_set_illcond_check(elsi_h, illcond_check)
    call elsi_set_illcond_tol(elsi_h, illcond_tol)

    ! Broadening
    call elsi_set_mu_broaden_scheme(elsi_h, which_broad)
    call elsi_set_mu_broaden_width(elsi_h, temp)
    call elsi_set_mu_mp_order(elsi_h, mp_order)

    ! Solver settings
    call elsi_set_elpa_solver(elsi_h, elpa_flavor)
    call elsi_set_elpa_n_single(elsi_h, elpa_n_single)
    call elsi_set_elpa_autotune(elsi_h, elpa_autotune)
    call elsi_set_elpa_gpu(elsi_h, elpa_gpu)

    call elsi_set_omm_flavor(elsi_h, omm_flavor)
    call elsi_set_omm_n_elpa(elsi_h, omm_n_elpa)
    call elsi_set_omm_tol(elsi_h, omm_tol)

! --- PEXSI
    
    if (pexsi_tasks_per_pole /= ELSI_NOT_SET) then
      call elsi_set_pexsi_np_per_pole(elsi_h, pexsi_tasks_per_pole)
    end if

    call elsi_set_pexsi_n_mu(elsi_h, pexsi_n_mu)
    call elsi_set_pexsi_n_pole(elsi_h, pexsi_n_pole)
    call elsi_set_pexsi_inertia_tol(elsi_h, pexsi_inertia_tol)
    
    call elsi_set_pexsi_np_symbo(elsi_h, pexsi_tasks_symbolic)
    call elsi_set_pexsi_temp(elsi_h, temp)

! --- SIPs
    
    if (sips_n_slice /= ELSI_NOT_SET) then
      call elsi_set_sips_n_slice(elsi_h, sips_n_slice)
    end if

    call elsi_set_sips_n_elpa(elsi_h, sips_n_elpa)
    call elsi_set_sips_ev_min(elsi_h, -10.0_dp)
    call elsi_set_sips_ev_max(elsi_h,  10.0_dp)

    call elsi_set_ntpoly_method(elsi_h, ntpoly_method)
    call elsi_set_ntpoly_filter(elsi_h, ntpoly_filter)
    call elsi_set_ntpoly_tol(elsi_h, ntpoly_tol)

 endif

 if ( (which_solver == PEXSI_SOLVER) .and. &
      .not. Get_EDM_only) then
    ! Set the proper bounds for the chemical potential
    if (iscf == 1) then
       call elsi_set_pexsi_mu_min(elsi_h, pexsi_initial_mu_min)
       call elsi_set_pexsi_mu_max(elsi_h, pexsi_initial_mu_max)
    else
       call elsi_get_pexsi_mu_min(elsi_h, mu_min)
       call elsi_get_pexsi_mu_max(elsi_h, mu_max)
       if (ionode) then
          print *, "*-- current mu_min, mu_max:", mu_min/eV, mu_max/eV
       endif
       mu_min = mu_min+dv_min
       mu_max = mu_max+dv_max
       ! Adjust chemical potential range for PEXSI
       call elsi_set_pexsi_mu_min(elsi_h, mu_min)
       call elsi_set_pexsi_mu_max(elsi_h, mu_max)
       if (ionode) then
          print *, "*-- updated mu_min, mu_max:", mu_min/eV, mu_max/eV
       endif

    endif   ! iscf == 1
 end if

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
             allocate(my_row_ptr2(my_no_l+1))
             my_row_ptr2(1) = 1
             do ih = 1,my_no_l
                my_row_ptr2(ih+1) = my_row_ptr2(ih) + pkg_spin%numcols(ih)
             enddo

             my_col_idx => pkg_spin%cols
             my_S => pkg_spin%vals(1)%data
             my_H => pkg_spin%vals(2)%data

             allocate(my_DM(my_nnz_l))
          endif

          ! Clean pkg_global
          nullify(pkg_global%vals(1)%data)   ! They were just pointing
          nullify(pkg_global%vals(2)%data)
          deallocate(pkg_global%vals)
          nullify(pkg_global%numcols)
          nullify(pkg_global%cols)

       enddo

       call elsi_set_csc(elsi_h, my_nnz, my_nnz_l, my_no_l, my_col_idx, my_row_ptr2)
       deallocate(my_row_ptr2)

       call elsi_set_csc_blk(elsi_h, BlockSize)
       call elsi_set_spin(elsi_h, n_spin, my_spin)
       call elsi_set_mpi(elsi_h, elsi_Spatial_comm)
       call elsi_set_mpi_global(elsi_h, elsi_global_comm)

    endif  ! n_spin

  call timer("elsi-solver", 1)

  if (n_spin == 1) then
     if (.not.Get_EDM_Only) then
       call elsi_dm_real_sparse(elsi_h, ham, ovlp, DM, energy)
       call elsi_get_entropy(elsi_h, ets)
     else
       call elsi_get_edm_real_sparse(elsi_h, DM)
     endif
  else
     ! Solve DM, and get (at every step for now) EDM, Fermi energy, and entropy
     ! Energy and entropy are already summed over spins
     if (.not.Get_EDM_Only) then
       call elsi_dm_real_sparse(elsi_h, my_H, my_S, my_DM, energy)
       call elsi_get_entropy(elsi_h, ets)
     else
       call elsi_get_edm_real_sparse(elsi_h, my_DM)
     endif
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
           deallocate(pkg_spin%vals)
           allocate(pkg_spin%vals(1))

           pkg_spin%vals(1)%data => my_DM(1:my_nnz_l)

        endif

        ! pkg_global is clean now
        call timer("redist_orbs_bck", 1)
        call redistribute_spmatrix(n_basis,pkg_spin,dist_spin(ispin) &
                                          ,pkg_global,dist_global,elsi_global_Comm)
        call timer("redist_orbs_bck", 2)

        ! Clean pkg_spin
        if (my_spin == ispin) then
           ! Each team deallocates during "its" spin cycle
           deallocate(my_DM)

           nullify(pkg_spin%vals(1)%data)    ! formerly pointing to DM
           deallocate(pkg_spin%vals)
           ! allocated in the direct transfer
           call de_alloc(pkg_spin%numcols,"pkg_spin%numcols","elsi_solver")
           call de_alloc(pkg_spin%cols,   "pkg_spin%cols",   "elsi_solver")
        endif


        ! In future, pkg_global%vals could be pointing to DM (or EDM)
        ! and the 'redistribute' routine check whether the vals arrays are
        ! associated, to use them instead of allocating them.
        DM(:,ispin)  = pkg_global%vals(1)%data(:)
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
        if (size(pkg_global%vals) /= 1) call die("pkg_global has two vals fields...")
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
     Dscf, ef, Entropy, occtol, neigwanted, Get_EDM_Only)

  use mpi_siesta, only: mpi_comm_dft
  use mpi
  use parallel, only: blocksize
  use class_Distribution
    use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
    use alloc, only: de_alloc  ! To deallocate some pointers in transfer matrices

  !
  ! K-point redistribution, Hk and Sk building, and call to complex ELSI solver
  !

  use m_fold_auxcell, only: fold_sparse_arrays ! Could be called in state_init

     real(dp), intent(in), target :: H(:,:), S(:)    ! Note that now we do not change them
     integer, intent(in) ::  iscf, maxnh, no_u, no_l, no_s, nkpnt
     integer, intent(in) ::  neigwanted, nspin
     integer, intent(in) ::  listhptr(no_l)
     integer, intent(in), target ::  listh(maxnh), numh(no_l)

      real(dp), intent(out) ::  Dscf(maxnh,nspin), ef, Entropy
      real(dp), intent(in)  ::  kpoint(3,nkpnt), qtot, temp, kweight(nkpnt), occtol,xijo(3,maxnh)


      integer :: mpirank, kcolrank, npGlobal
      integer :: npPerK, color, my_kpt_n

      integer :: kpt_comm, kpt_col_comm
      integer :: Global_Group, kpt_Group

      integer, allocatable :: ranks_in_world(:), ranks_in_world_AllK(:,:)
      type(distribution) :: dist_global
      type(distribution), allocatable :: dist_k(:)

      type(aux_matrix) :: pkg_global, pkg_k
      integer :: nvals
      real(dp), allocatable, target :: xijo_transp(:,:)
      real(dp), allocatable :: my_xijo_transp(:,:)
      real(dp), allocatable :: my_xij(:,:)

      integer :: my_no_l, my_nnz_l
      integer, pointer :: my_numh(:)
      integer, pointer :: my_listh(:)
      integer, allocatable :: my_listhptr(:)
      real(dp), allocatable :: my_S(:)
      real(dp), allocatable :: my_H(:,:)
      real(dp), pointer :: buffer(:)  ! for unpacking help

      real(dp), allocatable :: my_Dscf(:,:)
      real(dp), allocatable, target :: my_Dscf_reduced(:,:)

      complex(dp), allocatable :: DM_k(:,:)
      complex(dp), allocatable :: Hk(:,:)
      complex(dp), allocatable :: Sk(:)

      integer :: iuo, j, ind, ind_u, ispin
      real(dp) :: kxij, my_kpoint(3)
      complex(dp) :: kphs

      integer :: i, ih, ik, ierr

      integer, allocatable :: numh_u(:), listhptr_u(:), listh_u(:), ind2ind_u(:)
      integer :: nnz_u

      logical :: Get_EDM_Only

      external die

      ! Split the global communicator
      ! Re-distribute H and S to the k-point (and spin) teams
      ! Generate Hk, Sk
      ! Call elsi_complex_solver
      ! Construct and re-distribute global DM (or EDM)

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

      !print *, mpirank, "| ", "k-point ", my_kpt_n, " rank in col:", kcolrank

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

        ! Do sequential
        do ik=1, nkpnt

           call timer("redist_orbs_fwd", 1)
           call redistribute_spmatrix(no_u,pkg_global,dist_global, &
                                             pkg_k,dist_k(ik),elsi_global_Comm)
           call timer("redist_orbs_fwd", 2)

           if (my_kpt_n == ik) then
              !------------------------------------------
              ! Unpack info: real S and H (and index arrays) distributed over each kpt_comm
              my_no_l   = pkg_k%no_l
              my_nnz_l  = pkg_k%nnzl
              my_numh => pkg_k%numcols

              my_listh => pkg_k%cols

              allocate(my_S(my_nnz_l))
              my_S(:) = pkg_k%vals(1)%data(:)
              call de_alloc(pkg_k%vals(1)%data)

              allocate(my_xijo_transp(my_nnz_l,3))
              do i = 1, 3
                 buffer => pkg_k%vals(1+i)%data
                 my_xijo_transp(:,i) = buffer(:)
                 call de_alloc(pkg_k%vals(1+i)%data)
              enddo

              allocate(my_H(my_nnz_l,nspin))
              do ispin = 1, nspin
                 buffer => pkg_k%vals(4+ispin)%data
                 my_H(:,ispin) = buffer(:)
                 call de_alloc(pkg_k%vals(4+ispin)%data)
              enddo
              deallocate(pkg_k%vals)

              ! Now we could clear the rest of pkg_k--
              !  nullify(pkg_k%numcols)
              !  nullify(pkg_k%cols)
              !
              ! but we need to remember to deallocate the actual arrays after use
              ! it is probably safer to keep the pkg references
              !---------------------------

           endif
        enddo   !ik

          ! Clean pkg_global -- This is safe, as we use pointers to data only
          do i = 1, size(pkg_global%vals)
             nullify(pkg_global%vals(i)%data)
          enddo
          deallocate(pkg_global%vals)
          nullify(pkg_global%numcols)
          nullify(pkg_global%cols)

          deallocate(xijo_transp)  ! Auxiliary array used for sending

          ! generate listhptr for folding/unfolding operations
          allocate(my_listhptr(my_no_l))
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
          deallocate(my_xijo_transp)

      my_kpoint(:) = kpoint(:,my_kpt_n)

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

      deallocate(my_S,my_H)
      !

      !print *, mpirank, "| ", "k-point ", my_kpt_n, " Done folding"
      ! Prepare arrays for holding results
      allocate(DM_k(nnz_u,nspin))

      call elsi_complex_solver(iscf, no_u, my_no_l, nspin, nnz_u, numh_u, listhptr_u, &
                               listh_u, qtot, temp, Hk, Sk, DM_k, Ef, Entropy,  &
                               nkpnt, my_kpt_n, kpoint(:,my_kpt_n), kweight(my_kpt_n),    &
                               kpt_Comm, Get_EDM_Only )

      !print *, mpirank, "| ", "k-point ", my_kpt_n, " Done elsi_complex_solver"
      deallocate(listhptr_u, numh_u, listh_u)
      deallocate(Hk,Sk)

      ! Re-create DM (or EDM):
      ! Unfold within a given k
      ! Add up all the k contributions with the appropriate phases

      ! Prepare arrays for holding results: Note sizes: these are folded out
      allocate(my_Dscf(my_nnz_l,nspin))

      do iuo = 1, my_no_l
         do j = 1, my_numh(iuo)
            ind = my_listhptr(iuo) + j
            ind_u = ind2ind_u(ind)

            kxij = my_kpoint(1) * my_xij(1,ind) +    &
                 my_kpoint(2) * my_xij(2,ind) +    &
                 my_kpoint(3) * my_xij(3,ind)
            kphs = cdexp(dcmplx(0.0_dp, +1.0_dp)*kxij)

            do ispin = 1, nspin
               my_Dscf(ind,ispin) = real( DM_k(ind_u,ispin) * kphs, kind=dp )
            enddo
         enddo
      enddo

      deallocate(my_listhptr)

      ! Apply the k-point weight to the (apparently normalized) DM and EDM that
      ! come out of ELSI

      my_Dscf = kweight(my_kpt_n) * my_Dscf

      !print *, mpirank, "| ", "k-point ", my_kpt_n, " Done generating my_Dscf"
      deallocate(my_xij)
      deallocate(ind2ind_u)
      deallocate(DM_k)

      if (my_kpt_n == 1 ) then
         ! Prepare arrays for holding reduced data
         allocate(my_Dscf_reduced(my_nnz_l,nspin))
         if (kcolrank /= 0) call die("Rank 0 in kpt_comm not doing kpt 1")
      else
         ! These should not be referenced
         allocate(my_Dscf_reduced(1,1))
      endif

      ! Use k-point column communicator, and reduce to rank 0,
      ! which *should* correspond to kpt=1... (checked above)
      call MPI_Reduce( my_Dscf, my_Dscf_reduced, nspin*my_nnz_l, MPI_Double_Precision, &
           MPI_Sum, 0, kpt_col_comm, ierr )

      !print *, mpirank, "| ", "k-point ", my_kpt_n, " Done reducing my_Dscf"
      deallocate(my_Dscf)

      ! redistribute to global distribution, only from the first k-point

      if (my_kpt_n == 1) then

         ! These bits are still there *** update if we ever clean pkgs after use
         ! pkg_k%norbs = no_u
         ! pkg_k%no_l  = my_no_l
         ! pkg_k%nnzl  = my_nnz_l
         ! pkg_k%numcols => my_numh
         ! pkg_k%cols    => my_listh

         nvals = nspin   ! DM (or EDM)
         allocate(pkg_k%vals(nvals))
         do ispin = 1, nspin
            pkg_k%vals(ispin)%data => my_Dscf_reduced(:,ispin)
         enddo

      endif

      ! Everybody participates in the transfer in the receiving side, but
      ! only kpt=1 in the sending side, because we are using dist_k(1)
      call timer("redist_dm-edm_bck", 1)
      call redistribute_spmatrix(no_u,pkg_k,dist_k(1), &
           pkg_global,dist_global,elsi_global_Comm)
      call timer("redist_dm-edm_bck", 2)

      ! Deallocate aux arrays
      deallocate(my_Dscf_reduced)
      if (my_kpt_n == 1) then
         !Clean pkg_k in sender
         deallocate(pkg_k%vals)   ! just pointers
      endif

      !Clean all pkg_k: These were actually allocated in the forward transfer
      call de_alloc(pkg_k%numcols,"pkg_k%numcols")
      call de_alloc(pkg_k%cols,"pkg_k%cols")

      !print *, mpirank, "| ", "k-point ", my_kpt_n, " Done cleaning pkg_k"

      ! Unpack data (all processors, original global distribution)
        ! In future, pkg_global%vals(:) could be pointing to DM and EDM,
        ! and the 'redistribute' routine check whether the vals arrays are
        ! associated, to use them instead of allocating them.
         do ispin = 1, nspin
            Dscf(:,ispin)  = pkg_global%vals(ispin)%data(:)
         enddo
        ! Check no_l
        if (no_l /= pkg_global%no_l) then
           call die("Mismatch in no_l at end of kpoints-dispatcher")
        endif
        ! Check listH
        if (any(listh(:) /= pkg_global%cols(:))) then
           call die("Mismatch in listH at end of kpoints-dispatcher")
        endif
        !print *, mpirank, "| ", "k-point ", my_kpt_n, " Done Dscf"

        ! Clean pkg_global
        call de_alloc(pkg_global%numcols,"pkg_global%numcols","elsi_solver")
        call de_alloc(pkg_global%cols,   "pkg_global%cols",   "elsi_solver")
        do i = 1, size(pkg_global%vals)
           call de_alloc(pkg_global%vals(i)%data,"pkg_global%vals%data","kpoints_dispatcher")
        enddo
        deallocate(pkg_global%vals)
        !print *, mpirank, "| ", "k-point ", my_kpt_n, " Done cleaning pkg_global"

        ! Reduction of entropy over kpt_col_comm is not necessary

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
    deallocate(delta_v)
  end if

  call elsi_finalize(elsi_h)

  call MPI_Comm_Free(elsi_global_comm, ierr)

end subroutine elsi_finalize_scfloop

! This is a routine that is meant to be called by dhscf after
! computing the total potential V_scf.
! 
! It will find its minimum and maximum change between two invocations
! of dhscf (that is, between two buildings of the Hamiltonian from a DM: DM->H)
!
! When mixing the DM, the change in H between two solver steps (H->DM
! operation) (Delta_H below) is indeed related directly to the
! Delta_Vscf (see Lin Lin's paper).
! iscf = 1 DM(0) -> H(0) -> DM_out(1) -~> DM_mix(1)  
! iscf = 2 DM_mix(1) -> H(1) -> DM_out(2) -~> DM_mix(2)  
! Delta_H = H(1) - H(0) = (Delta_Vscf)

! When mixing H there is an extra mixing step that destroys this correspondence.
! iscf = 1 DM(0) -> H(0) -> DM_out(1) -> H_out(1) -~> H_mix(1)  
! iscf = 2 H_mix(1) -> DM_out(2) -> H_out(2) -~> H_mix(2)  
! Delta_H = H_mix(1) - H(0) /= (Delta_Vscf)
! However, for the simple case of linear mixing:
! H_mix(1) = alpha*H_out(1) + (1-alpha)*H(0), it can be seen that
! Delta_H = H_mix(1) - H(0) = alpha*(H_out(1)-H(0))
! so Delta_H in this case is a "damped" (Delta_Vscf)
! 
! It is then heuristically useful to employ the (Delta_Vscf)
! bracketing shifts, as described in the paper, for both cases. In the
! case of mixing H, the bracketing shift will be more conservative.
! (Could there be "Pulay" mixing sequences for which this does not
! hold?)

! Regarding the mechanics: For now this is confined to a single
! scf loop. At iscf=1, the v_old and delta_v arrays are allocated,
! and V_old is simply filled with the current V. In the corresponding
! solver step, we still do not need to re-bracket mu, and we use
! an initial bracket.
! For subsequent scf steps, Delta_V and dv_min=min(Delta_V) and
! dv_max=max(Delta_V) are computed with dhscf data. The solver calls
! check the shifts and update the bracket.
! At the end of the scf loop, the cleaning operations include the
! deallocation of V_old and delta_V.

! The two operations (computation of Delta_V and re-bracketing) are
! now decoupled, and this routine could be passed directly to dhscf as
! "things to do after computing V_scf", a possible general handler
! that could serve many other cases.  (the check for PEXSI_SOLVER
! below could be elided if the routine is enabled as handler only in
! that case)

subroutine elsi_save_potential(n_pts, n_spin, v_scf, comm)

  use m_mpi_utils, only: globalize_min, globalize_max

  integer,  intent(in) :: n_pts
  integer,  intent(in) :: n_spin
  real(dp), intent(in) :: v_scf(n_pts,n_spin)
  integer , intent(in) :: comm  ! The Siesta communicator used for grid operations
                                ! since this routine is called from dhscf...

  real(dp) :: tmp

  if (which_solver == PEXSI_SOLVER) then

    if (.not. allocated(v_old)) then
      allocate(v_old(n_pts,n_spin))
      allocate(delta_v(n_pts,n_spin))

      v_old = v_scf

    else

      delta_v = v_scf - v_old
      v_old = v_scf

      ! Get minimum and maximum of change of total potential
      tmp = minval(delta_v)
      call globalize_min(tmp, dv_min, comm=comm)

      tmp = maxval(delta_v)
      call globalize_max(tmp, dv_max, comm=comm)

      if (ionode) then
         print *, " * (dhscf) min and max Delta-V: ", dv_min/eV, dv_max/eV
      endif

    end if

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

#if SIESTA__ELSI || SIESTA__PEXSI

subroutine elsi_complex_solver(iscf, n_basis, n_basis_l, n_spin, nnz_l, numh, row_ptr, &
     col_idx, qtot, temp, ham, ovlp, DM, ef, ets, &
     nkpnt, kpt_n, kpt, weight, kpt_comm, Get_EDM_Only)

  use fdf,         only: fdf_get
  use m_mpi_utils, only: globalize_sum
  use parallel,    only: BlockSize
#ifdef MPI
  use mpi_siesta
#endif
  use class_Distribution
  use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
  use alloc, only: de_alloc  ! to deallocate some pointers in transfer matrices
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
  real(dp), intent(out)   :: ef        ! Fermi energy
  real(dp), intent(out)   :: ets       ! Entropy/k, dimensionless
  integer,  intent(in)    :: nkpnt     ! number of k-points
  integer,  intent(in)    :: kpt_n
  real(dp), intent(in)    :: kpt(3:)
  real(dp), intent(in)    :: weight
  integer,  intent(in)    :: kpt_comm

  integer :: ierr
  integer :: n_state
  integer :: nnz_g

  real(dp) :: energy

  integer, allocatable, dimension(:) :: row_ptr2

  integer :: elsi_Spatial_comm, elsi_Spin_comm

  type(distribution) :: dist_global
  type(distribution) :: dist_spin(2)

  type(aux_matrix) :: pkg_global, pkg_spin   ! Packages for transfer

  integer :: my_spin

  integer :: my_no_l
  integer :: my_nnz_l
  integer :: my_nnz
  integer, allocatable  :: my_row_ptr2(:)
  integer  :: i, ih, ispin, spin_rank, global_rank

  integer, pointer  :: my_col_idx(:)
  complex(dp), pointer :: my_S(:)
  complex(dp), pointer :: my_H(:)
  complex(dp), allocatable, target :: my_DM(:)

  integer :: date_stamp

  logical :: Get_EDM_Only

  external :: timer

#ifndef MPI
  call die("This ELSI solver interface needs MPI")
#endif

  call timer("elsi-complex-solver", 1)

  ! Initialization
  if (iscf == 1) then

    ! Get ELSI options
    call elsi_get_opts()

    ! Number of states to solve when calling an eigensolver
    n_state = min(n_basis, n_basis/2+5)

    ! Now we have all ingredients to initialize ELSI
    call elsi_init(elsi_h, which_solver, MULTI_PROC, SIESTA_CSC, n_basis, &
      qtot, n_state)

    ! Output
    if (ionode) then
      call elsi_set_output(elsi_h, out_level)
      call elsi_set_output_log(elsi_h, out_json)
      call elsi_set_write_unit(elsi_h, 6)
    endif

    ! Possible ill-conditioning of S
    call elsi_set_illcond_check(elsi_h, illcond_check)
    call elsi_set_illcond_tol(elsi_h, illcond_tol)

    ! Broadening
    call elsi_set_mu_broaden_scheme(elsi_h, which_broad)
    call elsi_set_mu_broaden_width(elsi_h, temp)
    call elsi_set_mu_mp_order(elsi_h, mp_order)

    ! Solver settings
    call elsi_set_elpa_solver(elsi_h, elpa_flavor)
    call elsi_set_elpa_n_single(elsi_h, elpa_n_single)
    call elsi_set_elpa_autotune(elsi_h, elpa_autotune)
    call elsi_set_elpa_gpu(elsi_h, elpa_gpu)

    call elsi_set_omm_flavor(elsi_h, omm_flavor)
    call elsi_set_omm_n_elpa(elsi_h, omm_n_elpa)
    call elsi_set_omm_tol(elsi_h, omm_tol)

    if (pexsi_tasks_per_pole /= ELSI_NOT_SET) then
      call elsi_set_pexsi_np_per_pole(elsi_h, pexsi_tasks_per_pole)
    end if

    call elsi_set_pexsi_n_mu(elsi_h, pexsi_n_mu)
    call elsi_set_pexsi_n_pole(elsi_h, pexsi_n_pole)
    call elsi_set_pexsi_inertia_tol(elsi_h, pexsi_inertia_tol)

    call elsi_set_pexsi_np_symbo(elsi_h, pexsi_tasks_symbolic)
    call elsi_set_pexsi_temp(elsi_h, temp)

    if (sips_n_slice /= ELSI_NOT_SET) then
      call elsi_set_sips_n_slice(elsi_h, sips_n_slice)
    end if

    call elsi_set_sips_n_elpa(elsi_h, sips_n_elpa)
    call elsi_set_sips_ev_min(elsi_h, -10.0_dp)
    call elsi_set_sips_ev_max(elsi_h,  10.0_dp)

    call elsi_set_ntpoly_method(elsi_h, ntpoly_method)
    call elsi_set_ntpoly_filter(elsi_h, ntpoly_filter)
    call elsi_set_ntpoly_tol(elsi_h, ntpoly_tol)

 endif   ! iscf == 1

 if ( (which_solver == PEXSI_SOLVER) .and. &
      .not. Get_EDM_only) then
    ! Set the proper bounds for the chemical potential
    if (iscf == 1) then
       call elsi_set_pexsi_mu_min(elsi_h, pexsi_initial_mu_min)
       call elsi_set_pexsi_mu_max(elsi_h, pexsi_initial_mu_max)
    else
       call elsi_get_pexsi_mu_min(elsi_h, mu_min)
       call elsi_get_pexsi_mu_max(elsi_h, mu_max)
       if (ionode) then
          print *, "*-- current mu_min, mu_max:", mu_min/eV, mu_max/eV
       endif
       mu_min = mu_min+dv_min
       mu_max = mu_max+dv_max
       ! Adjust chemical potential range for PEXSI
       call elsi_set_pexsi_mu_min(elsi_h, mu_min)
       call elsi_set_pexsi_mu_max(elsi_h, mu_max)
       if (ionode) then
          print *, "*-- updated mu_min, mu_max:", mu_min/eV, mu_max/eV
       endif

    endif   ! iscf == 1
 end if

      !print *, global_rank, "| ", " Entering elsi_complex_solver"

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

       call mpi_comm_rank( elsi_global_Comm, global_rank, ierr )

       ! MPI logic for spin polarization

       ! Split the communicator in spins and get distribution objects
       ! for the data redistribution needed
       ! Note that dist_spin is an array
       call get_spin_comms_and_dists(kpt_comm,kpt_comm, &  !! **** kpt_comm as global?
            blocksize, n_spin, &
            dist_global,dist_spin, elsi_spatial_comm, elsi_spin_comm)

       ! Find out which spin team we are in, and tag the spin we work on
       call mpi_comm_rank( elsi_Spin_Comm, spin_rank, ierr )
       my_spin = spin_rank+1  ! {1,2}

      !print *, global_rank, "| ", "spin ", my_spin, " After spin splitting"

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
                                             pkg_spin,dist_spin(ispin),kpt_Comm)

          call timer("redist_orbs_fwd", 2)

          if (my_spin == ispin) then  ! Each team gets their own data

             !nrows = pkg_spin%norbs          ! or simply 'norbs'
             my_no_l = pkg_spin%no_l
             my_nnz_l    = pkg_spin%nnzl
             call MPI_AllReduce(my_nnz_l,my_nnz,1,MPI_integer,MPI_sum,elsi_Spatial_Comm,ierr)
             ! generate off-by-one row pointer
             allocate(my_row_ptr2(my_no_l+1))
             my_row_ptr2(1) = 1
             do ih = 1,my_no_l
                my_row_ptr2(ih+1) = my_row_ptr2(ih) + pkg_spin%numcols(ih)
             enddo

             my_col_idx => pkg_spin%cols
             my_S => pkg_spin%complex_vals(1)%data
             my_H => pkg_spin%complex_vals(2)%data

             allocate(my_DM(my_nnz_l))
          endif

          ! Clean pkg_global
          nullify(pkg_global%complex_vals(1)%data)
          nullify(pkg_global%complex_vals(2)%data)
          deallocate(pkg_global%complex_vals)
          nullify(pkg_global%numcols)
          nullify(pkg_global%cols)

       enddo

       !print *, global_rank, "| ", "spin ", my_spin, "Done spin transfers"

       call elsi_set_csc(elsi_h, my_nnz, my_nnz_l, my_no_l, my_col_idx, my_row_ptr2)
       deallocate(my_row_ptr2)

       call elsi_set_csc_blk(elsi_h, BlockSize)
       call elsi_set_spin(elsi_h, n_spin, my_spin)
       call elsi_set_kpoint(elsi_h, nkpnt, kpt_n, weight)
       call elsi_set_mpi(elsi_h, elsi_Spatial_comm)
       call elsi_set_mpi_global(elsi_h, elsi_global_comm)

    endif  ! n_spin

  call timer("elsi-solver", 1)

  if (n_spin == 1) then
     if (.not.Get_EDM_Only) then
       call elsi_dm_complex_sparse(elsi_h, ham, ovlp, DM, energy)
       call elsi_get_entropy(elsi_h, ets)
     else
       call elsi_get_edm_complex_sparse(elsi_h, DM)
     endif
  else
     ! Solve DM, and get (at every step for now) EDM, Fermi energy, and entropy
     ! Energy and entropy already summed over spins
     if (.not.Get_EDM_Only) then
       call elsi_dm_complex_sparse(elsi_h, my_H, my_S, my_DM, energy)
       call elsi_get_entropy(elsi_h, ets)
     else
       call elsi_get_edm_complex_sparse(elsi_h, my_DM)
     endif
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

           deallocate(pkg_spin%complex_vals)
           allocate(pkg_spin%complex_vals(1))
           pkg_spin%complex_vals(1)%data => my_DM(1:my_nnz_l)

        endif

        ! pkg_global is clean now
        call timer("redist_orbs_bck", 1)
        call redistribute_spmatrix(n_basis,pkg_spin,dist_spin(ispin) &
                                          ,pkg_global,dist_global,kpt_Comm)
        call timer("redist_orbs_bck", 2)

        ! Clean pkg_spin
        if (my_spin == ispin) then
           ! Each team deallocates during "its" spin cycle
           deallocate(my_DM)

           nullify(pkg_spin%complex_vals(1)%data)    ! formerly pointing to DM
           deallocate(pkg_spin%complex_vals)
           ! allocated in the direct transfer
           call de_alloc(pkg_spin%numcols,"pkg_spin%numcols","elsi_solver")
           call de_alloc(pkg_spin%cols,   "pkg_spin%cols",   "elsi_solver")
        endif


        ! In future, pkg_global%vals(1) could be pointing to DM (or EDM),
        ! and the 'redistribute' routine check whether the vals arrays are
        ! associated, to use them instead of allocating them.
        DM(:,ispin)  = pkg_global%complex_vals(1)%data(:)
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
        deallocate(pkg_global%complex_vals)

     enddo

     call MPI_Comm_Free(elsi_Spatial_comm, ierr)
     call MPI_Comm_Free(elsi_Spin_comm, ierr)

  else    ! n_spin == 1

       ! Nothing else to do
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

# endif

end module m_elsi_interface
