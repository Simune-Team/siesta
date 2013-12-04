module m_pexsi_solver
    use precision, only  : dp

  public :: pexsi_solver

  real(dp), save :: prevDmax  ! For communication of max diff in DM in scf loop
                              ! which is used in the heuristics for N_el tolerance
  public :: prevDmax

CONTAINS

! This version uses separate distributions for Siesta (setup_H et al) and PEXSI.
!
  subroutine pexsi_solver(iscf, no_u, no_l, nspin,  &
       maxnh, numh, listhptr, listh, H, S, qtot, DM, EDM, &
       ef, Entropy, temp, delta_Ef)

    use fdf
    use parallel, only   : SIESTA_worker, BlockSize
    use parallel, only   : SIESTA_Group, SIESTA_Comm
    use m_mpi_utils, only: globalize_sum, globalize_max
    use m_mpi_utils, only: broadcast
    use units,       only: Kelvin, eV
    use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
    use class_Dist
    use alloc,             only: re_alloc, de_alloc
    use siesta_options,    only: dDtol
#ifdef MPI
    use mpi_siesta
#endif
    use m_pexsi_interface, only: f_ppexsi_solve_interface

    implicit          none

    integer, intent(in)  :: iscf  ! scf step number
    integer, intent(in)  :: maxnh, no_u, no_l, nspin
    integer, intent(in), target  :: listh(maxnh), numh(no_l), listhptr(no_l)
    real(dp), intent(in), target :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: qtot
    real(dp), intent(out), target:: DM(maxnh,nspin), EDM(maxnh,nspin)
    real(dp), intent(out)        :: ef  ! Fermi energy
    real(dp), intent(out)        :: Entropy ! Entropy/k, dimensionless
    real(dp), intent(in)         :: temp   ! Electronic temperature
    real(dp), intent(in)         :: delta_Ef  ! Estimated shift in E_fermi

#ifndef MPI
    call die("PEXSI needs MPI")
#else

    integer :: PEXSI_Comm, World_Comm
    integer :: PEXSI_Group, World_Group

    integer :: ispin, maxnhtot, ih, nnzold, i, pexsiFlag

    integer  :: ordering, isInertiaCount, numInertiaCounts, numMinICountShifts, numNodesTotal
    integer  :: muIter
    real(dp) :: muZeroT

    real(dp), save :: mu
    real(dp), save :: muMin0, muMax0
    real(dp), save :: muMinInertia, muMaxInertia
    real(dp), save :: muLowerEdge, muUpperEdge
    real(dp), save :: muSolverInput, muMinSolverInput, muMaxSolverInput
    real(dp), save :: previous_pexsi_temperature

    real(dp)       :: safe_width_ic, safe_width_solver
    real(dp)       :: safe_dDmax_NoInertia, safe_dDmax_Ef_inertia
    real(dp)       :: safe_dDmax_Ef_solver
    logical        :: do_inertia_count
    logical, save  :: first_call = .true.
    real(dp)       :: eBandStructure, eBandH, on_the_fly_tolerance

    integer        :: info, infomax

!Lin variables
integer :: nrows, nnz, nnzLocal, numColLocal
integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
!
real(dp), pointer, dimension(:) :: &
        HnzvalLocal=>null(), SnzvalLocal=>null(),  &
        DMnzvalLocal => null() , EDMnzvalLocal => null(), &
        FDMnzvalLocal => null()
!
real(dp), pointer, dimension(:) :: muList=>null(), &
                                   numElectronList=>null(), &
                                   numElectronDrvList=>null(), &
                                   shiftList=>null(), inertiaList=>null()
logical  :: PEXSI_worker
integer  :: numPole, nptsInertia
real(dp) :: temperature, numElectronExact, numElectron, gap, deltaE
real(dp) :: muInertia
real(dp) :: muMinPEXSI, muMaxPEXSI
integer  :: muMaxIter
integer  :: npPerPole
integer  :: npSymbFact
integer  :: mpirank, ierr
integer  :: isSIdentity
integer  :: inertiaMaxIter, inertiaIter, inertiaMaxNumRounds
logical  :: inertiaExpertDriver
logical  :: use_annealing
real(dp) :: annealing_preconditioner, temp_factor
real(dp) :: annealing_target_factor
real(dp) :: pexsi_temperature
real(dp) :: inertiaNumElectronTolerance, &
            inertiaMinNumElectronTolerance, &
            inertiaEnergyTolerance, &
            inertiaMuTolerance, &
            PEXSINumElectronToleranceMin, &
            PEXSINumElectronToleranceMax, &
            PEXSINumElectronTolerance
real(dp) :: lateral_expansion_solver, lateral_expansion_inertia
real(dp) :: free_bs_energy

!------------

real(dp) :: buffer1

external         :: timer
character(len=6) :: msg

type(aux_matrix) :: m1, m2
type(Dist)       :: dist1, dist2
integer          :: pbs, norbs, scf_step


! "SIESTA_Worker" means a processor which is in the Siesta subset.
! NOTE:  fdf calls will assign values to the whole processor set,
! but some other variables will have to be re-broadcast (see examples
! below)

World_Comm = true_MPI_Comm_World

if (SIESTA_worker) then
   ! deal with intent(in) variables
   norbs = no_u
   scf_step = iscf
endif
call broadcast(norbs,comm=World_Comm)
call broadcast(scf_step,comm=World_Comm)
call broadcast(prevDmax,comm=World_Comm)
call broadcast(dDtol,comm=World_Comm)

!  Find rank in global communicator
call mpi_comm_rank( World_Comm, mpirank, ierr )

call newDistribution(BlockSize,SIESTA_Group,dist1,TYPE_BLOCK_CYCLIC,"bc dist")

! Group and Communicator for first-pole team of PEXSI workers
!
npPerPole  = fdf_get("PEXSI.np-per-pole",4)
call MPI_Comm_Group(World_Comm, World_Group, Ierr)
call MPI_Group_incl(World_Group, npPerPole,   &
                    (/ (i,i=0,npPerPole-1) /),&
                    PEXSI_Group, Ierr)
call MPI_Comm_create(World_Comm, PEXSI_Group,&
                     PEXSI_Comm, Ierr)

PEXSI_worker = (mpirank < npPerPole)

! Number of processors for symbolic factorization
! Only relevant for PARMETIS/PT_SCOTCH
npSymbFact = fdf_get("PEXSI.np-symbfact",npPerPole)


pbs = norbs/npPerPole
call newDistribution(pbs,PEXSI_Group,dist2,TYPE_PEXSI,"px dist")


if (SIESTA_worker) then
   call timer("pexsi", 1)

   ispin = 1
   if (nspin /=1) then
      call die("Spin polarization not yet supported in PEXSI")
   endif

   numElectronExact = qtot 

   ! Note that the energy units for the PEXSI interface are arbitrary, but
   ! H, the interval limits, and the temperature have to be in the
   ! same units. Siesta uses Ry units.

   temperature      = temp

   if (mpirank==0) write(6,"(a,f10.2)") &
               "Electronic temperature (K): ", temperature/Kelvin

   m1%norbs = norbs
   m1%no_l  = no_l
   m1%nnzl  = sum(numH(1:no_l))
   m1%numcols => numH
   m1%cols    => listH
   allocate(m1%vals(2))
   m1%vals(1)%data => S(:)
   m1%vals(2)%data => H(:,ispin)

endif  ! SIESTA_worker

call timer("redist_orbs_fwd", 1)
call redistribute_spmatrix(norbs,m1,dist1,m2,dist2,World_Comm)
call timer("redist_orbs_fwd", 2)

if (PEXSI_worker) then

   nrows = m2%norbs          ! or simply 'norbs'
   numColLocal = m2%no_l
   nnzLocal    = m2%nnzl
   call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,PEXSI_Comm,ierr)

  call re_alloc(colptrLocal,1,numColLocal+1,"colptrLocal","pexsi_solver")
  colptrLocal(1) = 1
  do ih = 1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + m2%numcols(ih)
  enddo

  rowindLocal => m2%cols
  SnzvalLocal => m2%vals(1)%data
  HnzvalLocal => m2%vals(2)%data

  call re_alloc(DMnzvalLocal,1,nnzLocal,"DMnzvalLocal","pexsi_solver")
  call re_alloc(EDMnzvalLocal,1,nnzLocal,"EDMnzvalLocal","pexsi_solver")
  call re_alloc(FDMnzvalLocal,1,nnzLocal,"FDMnzvalLocal","pexsi_solver")

  call memory_all("after setting up H+S for PEXSI (PEXSI_workers)",PEXSI_comm)

endif ! PEXSI worker

call memory_all("after setting up H+S for PEXSI",World_comm)


isSIdentity = 0

numPole          = fdf_get("PEXSI.num-poles",20)
gap              = fdf_get("PEXSI.gap",0.0_dp,"Ry")

! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.
deltaE           = fdf_get("PEXSI.delta-E",3.0_dp,"Ry")


! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
PEXSINumElectronToleranceMin = fdf_get("PEXSI.num-electron-tolerance-lower-bound",0.01_dp)
PEXSINumElectronToleranceMax = fdf_get("PEXSI.num-electron-tolerance-upper-bound",0.5_dp)
! muMaxIter should be 1 or 2 later when combined with SCF.
muMaxIter        = fdf_get("PEXSI.mu-max-iter",10)

! How to expand the intervals in case of need.
lateral_expansion_solver = fdf_get("PEXSI.lateral-expansion-solver",1.0_dp*eV,"Ry")
lateral_expansion_inertia = fdf_get("PEXSI.lateral-expansion-inertia",3.0_dp*eV,"Ry")


! Initial guess of chemical potential and containing interval
! When using inertia counts, this interval can be wide.
! Note that mu, muMin0 and muMax0 are saved variables
if (first_call) then
   mu = fdf_get("PEXSI.mu",-0.60_dp,"Ry")
   ! Lower/Upper bound for the chemical potential.
   muMin0           = fdf_get("PEXSI.mu-min",-1.0_dp,"Ry")
   muMax0           = fdf_get("PEXSI.mu-max", 0.0_dp,"Ry")

   ! start with largest tolerance
   ! (except if overriden by user)
   PEXSINumElectronTolerance = fdf_get("PEXSI.num-electron-tolerance",&
                                       PEXSINumElectronToleranceMax)
   first_call = .false.
else
!
!  Here we could also check whether we are in the first scf iteration
!  of a multi-geometry run...
!
   ! Use a moving tolerance, based on how far DM_out was to DM_in
   ! in the previous iteration (except if overriden by user)

   on_the_fly_tolerance = Max(PEXSINumElectronToleranceMin, &
                              Min(prevDmax*1.0, PEXSINumElectronToleranceMax))
   PEXSINumElectronTolerance =  fdf_get("PEXSI.num-electron-tolerance",&
                                        on_the_fly_tolerance)
endif

! Arrays for reporting back information about the PEXSI iterations
call re_alloc(muList,1,muMaxIter,"muList","pexsi_solver")
call re_alloc(numElectronList,1,muMaxIter,"numElectronList","pexsi_solver")
call re_alloc(numElectronDrvList,1,muMaxIter,"numElectronDrvList","pexsi_solver")

!-----------------------------------------------------------------------------
! Use inertia counts?
isInertiaCount = fdf_get("PEXSI.inertia-count",1)
! For how many scf steps?
numInertiaCounts = fdf_get("PEXSI.inertia-counts",3)


! Maximum number of iterations for computing the inertia
! in a given scf step (until a proper bracket is obtained)
inertiaMaxIter   = fdf_get("PEXSI.inertia-max-iter",5)
inertiaMaxNumRounds  = inertiaMaxIter

! Call the inertia-count routine one step at a time, and perform
! convergence checks in the caller. Note that "inertiaMaxIter" will
! still be honored, in the sense that the total number of rounds will
! be capped by it.

inertiaExpertDriver = fdf_get("PEXSI.inertia-expert-driver",.false.)
if (inertiaExpertDriver) then
   inertiaMaxIter = 1
endif

! Stop inertia count if Ne(muMax) - Ne(muMin) < inertiaNumElectronTolerance
! This is the only stopping criterion available in non-expert mode
! Note that this number should grow with the size of the system

inertiaNumElectronTolerance = fdf_get("PEXSI.inertia-num-electron-tolerance",20)

! If the electron tolerance is too low the bracketting will not work.
! Set a minimum tolerance 
inertiaMinNumElectronTolerance = fdf_get("PEXSI.inertia-min-num-electron-tolerance",10)

inertiaNumElectronTolerance = max(inertiaNumElectronTolerance,inertiaMinNumElectronTolerance)


! Stop inertia count if (muMax - muMin) < inertiaEnergyTolerance
! Useful for metals only. By default, use a very small tolerance to deactivate this
! criterion. Reasonable values otherwise might be 0.1-0.2 eV.
!
inertiaEnergyTolerance = fdf_get("PEXSI.inertia-energy-tolerance",1.0e-5_dp*eV,"Ry")

! Stop inertia count if mu has not changed much from iteration to iteration.
! By default, use a very small tolerance to deactivate this
! criterion. Reasonable values otherwise might be 0.1-0.2 eV.
!
inertiaMuTolerance = fdf_get("PEXSI.inertia-mu-tolerance",1.0e-8_dp*eV,"Ry")

! Since we use the processor teams corresponding to the different poles, 
! the number of points in the energy interval should be a multiple
! of numNodesTotal/npPerPole.
! We can avoid serializing the calculations if we use only the
! available teams in a single step, subject to a minimum number:

! Minimum number of sampling points for inertia counts
numMinICountShifts = fdf_get("PEXSI.inertia-min-num-shifts", 10)

call mpi_comm_size( World_Comm, numNodesTotal, ierr )
nptsInertia = numNodesTotal/npPerPole
do
   if (nptsInertia < numMinICountShifts) then
      nptsInertia = nptsInertia + numNodesTotal/npPerPole
   else
      exit
   endif
enddo

! Arrays for reporting back information about the integrated DOS
! computed by the inertia count method.
call re_alloc(shiftList,1,nptsInertia,"shiftList","pexsi_solver")
call re_alloc(inertiaList,1,nptsInertia,"inertiaList","pexsi_solver")

! Ordering flag:
!   1: Use METIS
!   0: Use PARMETIS/PTSCOTCH
ordering = fdf_get("PEXSI.ordering",1)
!
! Broadcast these to the whole processor set, just in case
! (They were set only by the Siesta workers)
!
call MPI_Bcast(nrows,1,MPI_integer,0,World_Comm,ierr)
call MPI_Bcast(nnz,1,MPI_integer,0,World_Comm,ierr)
call MPI_Bcast(numElectronExact,1,MPI_double_precision,0,World_Comm,ierr)
call MPI_Bcast(temperature,1,MPI_double_precision,0,World_Comm,ierr)
call MPI_Bcast(delta_Ef,1,MPI_double_precision,0,World_Comm,ierr)
!
safe_width_ic = fdf_get("PEXSI.safe-width-ic-bracket",4.0_dp*eV,"Ry")
safe_width_solver = fdf_get("PEXSI.safe-width-solver-bracket",2.0_dp*eV,"Ry")
safe_dDmax_NoInertia = fdf_get("PEXSI.safe-dDmax-no-inertia",0.05)
safe_dDmax_Ef_Inertia = fdf_get("PEXSI.safe-dDmax-ef-inertia",0.1)
safe_dDmax_Ef_solver = fdf_get("PEXSI.safe-dDmax-ef-solver",0.05)

use_annealing = fdf_get("PEXSI.use-annealing",.false.)
if (use_annealing) then
   annealing_preconditioner = fdf_get("PEXSI.annealing-preconditioner",1.0_dp)
!   By default, the temperature goes to the target at a level 10 times dDtol
   annealing_target_factor = fdf_get("PEXSI.annealing-target-factor",10.0_dp)

   if (scf_step > 1 ) then

      ! Examples for target_factor = 10, dDtol=0.0001:
      ! prevDmax=0.1, preconditioner=1, factor=3
      ! prevDmax=0.1, preconditioner=2, factor=5
      ! prevDmax=0.1, preconditioner=3, factor=7
      ! prevDmax<=0.001, factor = 1
      ! prevDmax<0.001, factor = 1

      temp_factor = (log10(prevDmax/(annealing_target_factor*dDtol)))
      temp_factor = 1 + annealing_preconditioner * max(0.0_dp, temp_factor)

      pexsi_temperature = temp_factor * temperature
      if (pexsi_temperature > previous_pexsi_temperature) then
         if (mpirank==0) write(6,"(a,f10.2)") &
              "Will not raise PEXSI temperature to: ", &
              pexsi_temperature/Kelvin
         pexsi_temperature = previous_pexsi_temperature
      endif
      previous_pexsi_temperature = pexsi_temperature
   else
      ! No heuristics for now for first step
      previous_pexsi_temperature = huge(1.0_dp)
      pexsi_temperature = temperature
      !   Keep in mind for the future if modifying T at the 1st step
      !      previous_pexsi_temperature = pexsi_temperature
   endif
else
      pexsi_temperature = temperature
endif
if (mpirank==0) write(6,"(a,f10.2)") &
     "Current PEXSI temperature (K): ", pexsi_temperature/Kelvin

!
!  Possible inertia-count stage
!
do_inertia_count = .false.
if (isInertiaCount .ne. 0) then
  if (scf_step .le. numInertiaCounts) then
     if (mpirank == 0) write(6,"(a,i4)") "&o Inertia-count step scf_step<numIC ", scf_step
     do_inertia_count = .true.
  endif
  if (numInertiaCounts < 0) then
     if (scf_step <= -numInertiaCounts) then
        if (mpirank == 0) write(6,"(a,i4)") "&o Inertia-count step scf_step<-numIC ", scf_step
        do_inertia_count = .true.
     else if (prevDmax > safe_dDmax_NoInertia) then
        if (mpirank == 0) write(6,"(a,i4)") "&o Inertia-count step as prevDmax > safe_Dmax ", scf_step
        do_inertia_count = .true.
     endif
  endif
endif

if (do_inertia_count) then

 ! Proper bracketting
 if (scf_step > 1) then
   if (prevDmax < safe_dDmax_Ef_inertia) then
      ! Shift brackets using estimate of Ef change from previous iteration
      !
      if (mpirank == 0) write(6,"(a)") "&o Inertia-count bracket shifted by Delta_Ef"
      muMin0 = muMinInertia + delta_Ef
      muMax0 = muMaxInertia + delta_Ef
      mu     = mu + delta_Ef
   else
      if (mpirank == 0) write(6,"(a)") "&o Inertia-count safe bracket"
      muMin0 = min(muLowerEdge - 0.5*safe_width_ic, muMinInertia)
      muMax0 = max(muUpperEdge + 0.5*safe_width_ic, muMaxInertia)
   endif
 endif

  call do_inertia()
  call memory_all("after inertia-count ",World_comm)

endif

!
!  do actual solve
!
  ! Use the starting values, or the output of previous inertia-count iterations
  ! (see above)

  if (do_inertia_count) then
     if (mpirank == 0) write(6,"(a)") "&o Solver bracket from Inertia count"
     muSolverInput = muInertia
     muMinSolverInput = muMinInertia
     muMaxSolverInput = muMaxInertia
  else
     if (scf_step > 1) then
        if (prevDmax < safe_dDmax_Ef_solver) then
           ! Shift brackets using estimate of Ef change from previous iteration
           if (mpirank == 0) write(6,"(a)") "&o Solver bracket shifted by delta_Ef"
           muSolverInput = mu + delta_Ef
           muMinSolverInput = muMinSolverInput + delta_Ef
           muMaxSolverInput = muMaxSolverInput + delta_Ef
        else
           if (mpirank == 0) write(6,"(a)") "&o Safe Solver bracket"
           muSolverInput = mu 
           muMinSolverInput = min(mu - 0.5*safe_width_solver, muMinSolverInput)
           muMaxSolverInput = max(mu + 0.5*safe_width_solver, muMaxSolverInput)
        endif
     else
        if (mpirank == 0) write(6,"(a)") "&o Solver bracket from initial values"
        muSolverInput = mu
        muMinSolverInput = muMin0
        muMaxSolverInput = muMax0
     endif

  endif

solver_loop: do

   if(mpirank == 0) then
     write (6,"(a,2f9.4,a,f9.4,a,f9.5)") 'Calling PEXSI   (eV): [', &
                                              muMinSolverInput/eV, &
                                              muMaxSolverInput/eV, &
                                              "] estimated mu: ", muSolverInput/eV, &
                                              ' Tol: ', PEXSINumElectronTolerance
   endif

   call timer("pexsi-solver", 1)

   call f_ppexsi_solve_interface(&
        ! input parameters
        nrows,&
        nnz,&
        nnzLocal,&
        numColLocal,&
        colptrLocal,&
        rowindLocal,&
        HnzvalLocal,&
        isSIdentity,&
        SnzvalLocal,&
        pexsi_temperature,&
        numElectronExact,&
        muSolverInput,&
        muMinSolverInput,&
        muMaxSolverInput,&
        gap,&
        deltaE,&
        numPole,&
        muMaxIter,&
        PEXSINumElectronTolerance,&
        ordering,&
        npPerPole,&
        npSymbFact,&
        World_Comm,&
! output parameters
        DMnzvalLocal,&
        EDMnzvalLocal,&
        FDMnzvalLocal,&
        mu,&
        numElectron,&
        muMinPEXSI,&
        muMaxPEXSI,&
        muIter,&
        muList,&
        numElectronList,&
        numElectronDrvList,&
        info)

   call memory_all("after solver ",World_comm)
   call timer("pexsi-solver", 2)

  ! Make sure that all processors report info=1 when any of them
  ! raises it...
  ! This should be done inside the routine
   call globalize_max(info,infomax,comm=World_Comm)
   info = infomax

   if (info /= 0) then
      write(msg,"(i6)") info 
      call die("Error in pexsi solver routine. Info: " // msg)
   endif


   if (abs(numElectron-numElectronExact) > PEXSInumElectronTolerance) then

      ! The solver did not converge.
      ! We need to do an inertia-count step...
      ! check which side of the interval we need to expand
      !
      if ((numElectron-numElectronExact) > 0) then
         muMin0 = muMin0 - lateral_expansion_solver 
      else
         muMax0 = muMax0 + lateral_expansion_solver 
      endif
      if (mpirank == 0) then
         write (6,"(a,f10.4,a,f8.4,a,f9.4,a)") 'The PEXSI solver did not converge. Nel - Nel_exact: ', &
                                          (numElectron - numElectronExact), &
                                          ". Tolerance: ", PEXSInumElectronTolerance, &
                                          '. Expanding interval by ', lateral_expansion_solver/eV, ' eV.'

         write(6,"(a,i3)") " #&s Number of solver iterations: ", muIter
         write(6,"(a3,2a12,a20)") "it", "mu", "N_e", "dN_e/dmu"
         do i = 1, muIter
            write(6,"(i3,2f12.4,g20.5,2x,a2)") i, muList(i)/eV, &
                 numElectronList(i), numElectronDrvList(i)*eV, "&s"
         end do
      endif

      call do_inertia()

      muSolverInput = muInertia
      muMinSolverInput = muMinInertia
      muMaxSolverInput = muMaxInertia

      !  will take another iteration of the solver
   else
      exit solver_loop
   endif

enddo solver_loop

!------------ End of solver step

   if (mpirank == 0) then
            write(6,"(a,i3)") " #&s Number of solver iterations: ", muIter
      write(6,"(a3,2a12,a20)") "it", "mu", "N_e", "dN_e/dmu"
      do i = 1, muIter
         write(6,"(i3,2f12.4,g20.5,2x,a2)") i, muList(i)/eV, &
              numElectronList(i), numElectronDrvList(i)*eV, "&s"
      end do
   endif

if (PEXSI_worker) then

   free_bs_energy = 0.0_dp
   eBandStructure = 0.0_dp
   eBandH = 0.0_dp
   do i = 1,nnzLocal
      free_bs_energy = free_bs_energy + SnzvalLocal(i) * &
           ( FDMnzvalLocal(i) )
      eBandStructure = eBandStructure + SnzvalLocal(i) * &
           ( EDMnzvalLocal(i) )
      eBandH = eBandH + HnzvalLocal(i) * &
           ( DMnzvalLocal(i) )
   enddo

   call de_alloc(FDMnzvalLocal,"FDMnzvalLocal","pexsi_solver")
   call de_alloc(colPtrLocal,"colPtrLocal","pexsi_solver")

   ! These operations in PEXSI group now
   call globalize_sum( free_bs_energy, buffer1, comm=PEXSI_comm )
   ! Note that FDM has an extra term: -mu*N
   free_bs_energy = buffer1 + mu*numElectron
   call globalize_sum( eBandStructure, buffer1, comm=PEXSI_comm )
   eBandStructure = buffer1
   call globalize_sum( eBandH, buffer1, comm=PEXSI_comm )
   eBandH = buffer1

   if( mpirank == 0 ) then
      write(*, *) "mu (eV)       = ", mu/eV
      write(*, *) "muMinSolverInput  = ", muMinSolverInput/eV
      write(*, *) "muMaxSolverInput  = ", muMaxSolverInput/eV
      write(*, *) "muMinPEXSI    = ", muMinPEXSI/eV
      write(*, *) "muMaxPEXSI    = ", muMaxPEXSI/eV
      write(*, *) "muZeroT (eV)  = ", muZeroT/eV
      write(*, *) "numElectron   = ", numElectron
      write(*, *) "eBandS (eV) = ", eBandStructure/eV
      write(*, *) "eBandH (eV) = ", eBandH/eV
      write(*, *) "freeBandEnergy (eV) = ", (free_bs_energy)/eV
      write(*, *) "eBandS (Ry) = ", eBandStructure
      write(*, *) "eBandH (Ry) = ", eBandH
      write(*, *) "freeBandEnergy (Ry) = ", (free_bs_energy)

   endif

   ef = mu
   Entropy = - (free_bs_energy - ebandStructure) / temp

   call de_alloc(m2%vals(1)%data,"m2%vals(1)%data","pexsi_solver")
   call de_alloc(m2%vals(2)%data,"m2%vals(2)%data","pexsi_solver")

   m2%vals(1)%data => DMnzvalLocal(1:nnzLocal)
   m2%vals(2)%data => EDMnzvalLocal(1:nnzLocal)

endif ! PEXSI_worker

call de_alloc(muList,            "muList",            "pexsi_solver")
call de_alloc(numElectronList,   "numElectronList",   "pexsi_solver")
call de_alloc(numElectronDrvList,"numElectronDrvList","pexsi_solver")

! Prepare m1 to receive the results
if (SIESTA_worker) then
   nullify(m1%vals(1)%data)    ! formerly pointing to S
   nullify(m1%vals(2)%data)    ! formerly pointing to H
   deallocate(m1%vals)
   nullify(m1%numcols)         ! formerly pointing to numH
   nullify(m1%cols)            ! formerly pointing to listH
endif

call timer("redist_orbs_bck", 1)
call redistribute_spmatrix(norbs,m2,dist2,m1,dist1,World_Comm)
call timer("redist_orbs_bck", 2)

if (PEXSI_worker) then
   call de_alloc(DMnzvalLocal, "DMnzvalLocal", "pexsi_solver")
   call de_alloc(EDMnzvalLocal,"EDMnzvalLocal","pexsi_solver")

   nullify(m2%vals(1)%data)    ! formerly pointing to DM
   nullify(m2%vals(2)%data)    ! formerly pointing to EDM
   deallocate(m2%vals)
   call de_alloc(m2%numcols,"m2%numcols","pexsi_solver") ! allocated in the direct transfer
   call de_alloc(m2%cols,   "m2%cols",   "pexsi_solver")
endif

! We assume that the root node is common to both communicators
if (SIESTA_worker) then
   call broadcast(ef,comm=SIESTA_Comm)
   call broadcast(Entropy,comm=SIESTA_Comm)
   ! In future, m1%vals(1,2) could be pointing to DM and EDM,
   ! and the 'redistribute' routine check whether the vals arrays are
   ! associated, to use them instead of allocating them.
   DM(:,ispin)  = m1%vals(1)%data(:)    
   EDM(:,ispin) = m1%vals(2)%data(:)    
   ! Check no_l
   if (no_l /= m1%no_l) then
      call die("Mismatch in no_l")
   endif
   ! Check listH
   if (any(listH(:) /= m1%cols(:))) then
      call die("Mismatch in listH")
   endif

   call de_alloc(m1%vals(1)%data,"m1%vals(1)%data","pexsi_solver")
   call de_alloc(m1%vals(2)%data,"m1%vals(2)%data","pexsi_solver")
   deallocate(m1%vals)
   call de_alloc(m1%numcols,"m1%numcols","pexsi_solver") ! allocated in the direct transfer
   call de_alloc(m1%cols,   "m1%cols",   "pexsi_solver")

   call timer("pexsi", 2)

endif


call de_alloc(shiftList,"shiftList","pexsi_solver")
call de_alloc(inertiaList,"shiftList","pexsi_solver")

call delete(dist1)
call delete(dist2)
if (PEXSI_worker) then
   call MPI_Comm_Free(PEXSI_Comm, ierr)
   call MPI_Group_Free(PEXSI_Group, ierr)
endif

#endif 

CONTAINS 

!----------
! Improve bracketting for mu with inertia counts,
! starting with the interval (muMin0, muMax0)
! Results:  muMinInertia, muMaxInertia
!           muLowerEdge, muUpperEdge
!           muInertia
!
subroutine do_inertia()

    use m_pexsi_interface, only: f_ppexsi_inertiacount_interface
    use m_convergence, only: converger_t
    use m_convergence, only: reset, set_tolerance, is_converged, add_value
    use m_interpolate, only: interpolate

    logical        :: interval_problem, one_more_round
    logical        :: bad_lower_bound
    logical        :: bad_upper_bound
    real(dp)       :: inertia_electron_width, inertia_energy_width
    real(dp)       :: inertia_original_electron_width
    integer        :: nInertiaRounds
    real(dp)       :: numElectronMax, numElectronMin
    type(converger_t)  ::  conv_mu

   nInertiaRounds = 0

refine_interval: do

   if ( (mpirank == 0) .and. (nInertiaRounds == 0) ) then
     write (6,"(a,2f9.4,a,a,i4)") 'Calling inertiaCount: [', &
                                      muMin0/eV, muMax0/eV, "] (eV)", &
                                      " Nshifts: ", nptsInertia
   endif 


   call timer("pexsi-inertia-ct", 1)

   call f_ppexsi_inertiacount_interface(&
   ! input parameters
        nrows,&
        nnz,&
        nnzLocal,&
        numColLocal,&
        colptrLocal,&
        rowindLocal,&
        HnzvalLocal,&
        isSIdentity,&
        SnzvalLocal,&
        pexsi_temperature,&
        numElectronExact,&
        muMin0,&
        muMax0,&
        nptsInertia,&
        inertiaMaxIter,&
        inertiaNumElectronTolerance,&
        ordering,&
        npPerPole,&
        npSymbFact,&
        World_Comm,&
! output parameters
        muMinInertia,&
        muMaxInertia,&
        muLowerEdge,&
        muUpperEdge,&
        inertiaIter,&
        shiftList,&
        inertiaList,&
        info)

   muInertia    = (muLowerEdge + muUpperEdge) / 2d0

   call timer("pexsi-inertia-ct", 2)

   call globalize_max(info,infomax,comm=World_Comm)
   info = infomax

   interval_problem = .false.
   if (info /= 0) then

      if(mpirank == 0) then
         bad_lower_bound = (inertiaList(1) > (numElectronExact - 0.1)) 
         bad_upper_bound = (inertiaList(nptsInertia) < (numElectronExact + 0.1)) 
      endif

      call broadcast(bad_lower_bound,comm=World_Comm)
      call broadcast(bad_upper_bound,comm=World_Comm)

      if (bad_lower_bound) then
         interval_problem =  .true.
         muMin0 = muMin0 - lateral_expansion_inertia ! 0.5
         if (mpirank==0) then
            write (6,"(a,2f12.4,a,2f10.4)") 'Wrong inertia-count interval (lower end). Counts: ', &
                                 inertiaList(1), inertiaList(nptsInertia), &
                                 ' New interval: ', muMin0/eV, muMax0/eV
         endif
      endif
      if (bad_upper_bound) then
         interval_problem =  .true.
         muMax0 = muMax0 + lateral_expansion_inertia ! 0.5
         if (mpirank==0) then
            write (6,"(a,2f12.4,a,2f10.4)") 'Wrong inertia-count interval (upper end). Counts: ', &
                                 inertiaList(1), inertiaList(nptsInertia), &
                                 ' New interval: ', muMin0/eV, muMax0/eV
         endif
      endif

      if (interval_problem) then
          ! do nothing more, stay in loop
	  nInertiaRounds = 0
      else
         write(msg,"(i6)") info 
         call die("Non-interval error in inertia count routine. Info: " // msg)
      endif

   else

      ! "Expert mode" of operation. At every inertia-count iteration,
      ! we check the energy and electron number widths and the convergence of mu,
      ! and decide whether to go one more round
      !
      ! By driving the routine one-step at a time, the only major difference
      ! would be extra calls to the symbolic factorization...

      if (.not. inertiaExpertDriver) exit refine_interval

       nInertiaRounds = nInertiaRounds + 1

      if (mpirank==0) then
         inertia_energy_width = (muMaxInertia - muMinInertia)
         ! Note that this is the width of the starting interval...
         inertia_original_electron_width = (inertiaList(nptsInertia) - inertiaList(1))
         ! Compute the actual electron width
         numElectronMax = interpolate(shiftList,inertiaList,muMaxInertia)
         numElectronMin = interpolate(shiftList,inertiaList,muMinInertia)
         inertia_electron_width = (numElectronMax - numElectronMin)

         write (6,"(a,2f9.4,a,f9.4,3(a,f10.3))") ' -- new bracket (eV): [', &
              muMinInertia/eV, muMaxInertia/eV,  &
              "] estimated mu: ", muInertia/eV, &
              " Nel width: ", inertia_electron_width, &
              " (Base: ", inertia_original_electron_width, &
	      " ) E width: ", inertia_energy_width/eV

         if (nInertiaRounds == 1) then
           call reset(conv_mu)
           call set_tolerance(conv_mu,inertiaMuTolerance)
         endif
         call add_value(conv_mu, muInertia)

         !
         one_more_round = .true.
         if (inertia_original_electron_width < inertiaNumElectronTolerance) then
            write (6,"(a)") 'Leaving inertia loop: electron tolerance'
            one_more_round = .false.
         endif
         if (inertia_electron_width < inertiaMinNumElectronTolerance) then
            write (6,"(a)") 'Leaving inertia loop: minimum workable electron tolerance'
            one_more_round = .false.
         endif
         if (inertia_energy_width < inertiaEnergyTolerance) then
            write (6,"(a,f12.6)") 'Leaving inertia loop: energy tolerance: ', inertiaEnergyTolerance/eV
            one_more_round = .false.
         endif
         if (is_converged(conv_mu)) then
            write (6,"(a,f12.6)") 'Leaving inertia loop: mu tolerance: ', inertiaMuTolerance/eV
            one_more_round = .false.
         endif
	 if (nInertiaRounds == inertiaMaxNumRounds) then
            write (6,"(a)") 'Leaving inertia loop: too many rounds'
            one_more_round = .false.
         endif
      endif
      call broadcast(one_more_round,comm=World_Comm)
      
      if (one_more_round) then
         ! stay in loop
	 muMin0 = muMinInertia
	 muMax0 = muMaxInertia
      else
         exit refine_interval
      endif
   endif

enddo refine_interval


   if(mpirank == 0) then
     write (*,"(/,a)") 'PEXSI inertia count executed'
     write (*,"(a,f10.4)") ' mu (eV)=', muInertia/eV
     write (*,"(a,f10.4)") ' mu lowerEdge (eV):', muLowerEdge/eV
     write (*,"(a,f10.4)") ' mu upperEdge (eV):', muUpperEdge/eV
     write (*,"(a,f10.4)") ' muMin (eV):', muMinInertia/eV
     write (*,"(a,f10.4)") ' muMax (eV):', muMaxInertia/eV

     write(*,"(/,a)") "Cumulative DOS by inertia count:"
     do i=1, nptsInertia
        write(*,"(f10.4,f10.4)") shiftList(i)/eV, inertiaList(i)
     enddo
  end if
end subroutine do_inertia

end subroutine pexsi_solver
end module m_pexsi_solver
