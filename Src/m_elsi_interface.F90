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
!        -  MPI.Nprocs.SIESTA is not working
!        -  Spin (and k-point?)
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

! real*8 eo(maxo,nspin,nk)  : Eigenvalues 
! real*8 qo(maxo,nspin,nk)  : Occupations of eigenstates
! real*8 Dnew(maxnd,spin%DM)  : Output Density Matrix
! real*8 Enew(maxnd,spin%EDM)  : Output Energy-Density Matrix
! real*8 ef                      : Fermi energy
! real*8 Entropy                 : Electronic entropy


     real(dp), intent(inout) :: H(:,:), S(:)    ! Note!
     integer, intent(in) ::  iscf, maxnh, no_u, no_l, no_s, nkpnt
     integer, intent(in) ::  neigwanted, nspin
     integer, intent(in) ::  indxuo(no_s), listh(maxnh), numh(no_l), listhptr(no_l)

      real(dp), intent(out) ::  Dscf(maxnh,nspin), ef, Escf(maxnh,nspin), Entropy
      real(dp), intent(out) ::  eo(no_u,nspin,nkpnt), qo(no_u,nspin,nkpnt)
      real(dp), intent(in)  ::  kpoint(3,nkpnt), qtot, temp, kweight(nkpnt), occtol,xijo(3,maxnh)

      logical :: gamma, using_aux_cell

      external die

      gamma = ((nkpnt == 1) .and. (sum(abs(kpoint(:,1))) == 0.0_dp))
      using_aux_cell =  (no_s /= no_u)

      if (gamma) then
         if  (.not. using_aux_cell) then
            call elsi_solver(iscf, no_u, no_l, nspin, &
                     maxnh, listhptr, listh, qtot, temp, &
                     H, S, Dscf, Escf, ef, Entropy)
         else
            call die("Cannot do gamma-point ELSI with supercell yet")
            ! Fold arrays
            ! call routine
            ! UNFOLD !!!
         endif
      else
         ! Always using aux cell...
         call die("Cannot do k-points with ELSI yet")
         ! Fold arrays
         !call elsi_solver_complex( )
         ! UNFOLD !!!
      endif

end subroutine getdm_elsi

! This version uses separate distributions for Siesta (setup_H et al) and ELSI
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

  implicit none

  integer,  intent(in)    :: iscf      ! SCF step counter
  integer,  intent(in)    :: n_basis   ! Global basis
  integer,  intent(in)    :: n_basis_l ! Local basis
  integer,  intent(in)    :: n_spin
  integer,  intent(in)    :: nnz_l     ! Local nonzero
  integer,  intent(in)    :: row_ptr(n_basis_l)
  integer,  intent(in)    :: col_idx(nnz_l)
  real(dp), intent(in)    :: qtot
  real(dp), intent(in)    :: temp
  real(dp), intent(inout) :: ham(nnz_l,n_spin)
  real(dp), intent(inout) :: ovlp(nnz_l)
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

  external :: timer

#ifndef MPI
  call die("This ELSI solver interface needs MPI")
#endif

  ! Global communicator is a duplicate of passed communicator
  call MPI_Comm_Dup(MPI_Comm_DFT, elsi_global_comm, ierr)
  call MPI_Comm_Rank(elsi_global_comm, mpirank, ierr)

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
    case ("elpa")
      which_solver = ELPA_SOLVER
    case ("omm")
      which_solver = OMM_SOLVER
    case ("pexsi")
      which_solver = PEXSI_SOLVER
    case ("sips")
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

    ! Sparsity pattern
    call globalize_sum(nnz_l, nnz_g, comm=elsi_global_comm)

    allocate(row_ptr2(n_basis_l+1))

    row_ptr2(1:n_basis_l) = row_ptr(1:n_basis_l)+1
    row_ptr2(n_basis_l+1) = nnz_l+1

    call elsi_set_csc(elsi_h, nnz_g, nnz_l, n_basis_l, col_idx, row_ptr2)
    call elsi_set_csc_blk(elsi_h, BlockSize)

    deallocate(row_ptr2)

    ! MPI
!    call elsi_set_spin(elsi_h, n_spin, elsi_spin_id)
    call elsi_set_mpi(elsi_h, elsi_global_comm)
    call elsi_set_mpi_global(elsi_h, elsi_global_comm)

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
  end if

  call timer("elsi-solver", 1)

  ! Solve DM, and get (at every step for now) EDM, Fermi energy, and entropy
  call elsi_dm_real_sparse(elsi_h, ham, ovlp, dm, energy)
  call elsi_get_edm_real_sparse(elsi_h, edm)
  call elsi_get_mu(elsi_h, ef)
  call elsi_get_entropy(elsi_h, ets)

  ets = ets/temp

  call timer("elsi-solver", 2)

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
