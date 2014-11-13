module m_pexsi_solver
    use precision, only  : dp

  implicit none

  public :: pexsi_solver

  real(dp), save :: prevDmax  ! For communication of max diff in DM in scf loop
                              ! which is used in the heuristics for N_el tolerance
  public :: prevDmax


CONTAINS

! This version uses separate distributions for Siesta (setup_H et al) and PEXSI.
! It uses the simple KSDFT driver
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
use f_ppexsi_interface
use iso_c_binding
use m_pexsi, only: plan, pexsi_initialize_scfloop

#ifdef TRACING_SOLVEONLY
      use extrae_module
#endif

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

type(f_ppexsi_options) :: options

integer :: numTotalPEXSIIter
integer :: numTotalInertiaIter
real(dp) :: totalEnergyH
real(dp) :: totalEnergyS
real(dp) :: totalFreeEnergy


    integer :: PEXSI_Comm, World_Comm
    integer :: PEXSI_Group, World_Group

    integer :: ispin, maxnhtot, ih, nnzold, i, pexsiFlag

    integer  :: isInertiaCount, numInertiaCounts

    real(dp), save :: muMin0, muMax0, mu
    real(dp), save :: muMinInertia, muMaxInertia
    real(dp), save :: previous_pexsi_temperature
    logical, save  :: first_call = .true.

    real(dp)       :: bs_energy, eBandH, on_the_fly_tolerance

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
real(dp) :: temperature, numElectronExact, numElectron
integer  :: npPerPole
integer  :: mpirank, ierr
integer  :: isSIdentity

real(dp) :: pexsi_temperature, two_kT

real(dp) :: PEXSINumElectronToleranceMin, &
            PEXSINumElectronToleranceMax, &
            PEXSINumElectronTolerance
real(dp) :: free_bs_energy

!------------

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


! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
PEXSINumElectronToleranceMin =  &
         fdf_get("PEXSI.num-electron-tolerance-lower-bound",0.01_dp)
PEXSINumElectronToleranceMax =  &
         fdf_get("PEXSI.num-electron-tolerance-upper-bound",0.5_dp)


! Initial guess of chemical potential and containing interval
! When using inertia counts, this interval can be wide.
! Note that mu, muMin0 and muMax0 are saved variables
if (first_call) then
   ! Lower/Upper bound for the chemical potential.
   muMin0           = fdf_get("PEXSI.mu-min",-1.0_dp,"Ry")
   muMax0           = fdf_get("PEXSI.mu-max", 0.0_dp,"Ry")
   mu = fdf_get("PEXSI.mu",-0.60_dp,"Ry")

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

   call get_on_the_fly_tolerance(prevDmax,on_the_fly_tolerance)

   ! Override if tolerance is explicitly specified in the fdf file
   PEXSINumElectronTolerance =  fdf_get("PEXSI.num-electron-tolerance",&
                                        on_the_fly_tolerance)
endif

!-----------------------------------------------------------------------------

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
call get_current_temperature(pexsi_temperature)
!
!  Set guard smearing for later use
!
two_kT = 2.0_dp * pexsi_temperature

!                                                                             
!  New interface.
if (iscf == 1) then
   call pexsi_initialize_scfloop(World_Comm,npPerPole,mpirank)
endif
!
call f_ppexsi_set_default_options( options )

options%muPEXSISafeGuard = fdf_get("PEXSI.mu-pexsi-safeguard",0.05_dp,"Ry")
options%maxPEXSIIter = fdf_get("PEXSI.mu-max-iter",10)

isSIdentity = 0

options%numPole  = fdf_get("PEXSI.num-poles",40)
options%gap      = fdf_get("PEXSI.gap",0.0_dp,"Ry")

! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.
options%deltaE     = fdf_get("PEXSI.delta-E",3.0_dp,"Ry") ! Lin: 10 Ry...

! Ordering flag:
!   1: Use METIS
!   0: Use PARMETIS/PTSCOTCH
options%ordering = fdf_get("PEXSI.ordering",1)

! Number of processors for symbolic factorization
! Only relevant for PARMETIS/PT_SCOTCH
options%npSymbFact = fdf_get("PEXSI.np-symbfact",npPerPole)

options%verbosity = fdf_get("PEXSI.verbosity",1)

options%temperature = pexsi_temperature
options%numElectronPEXSITolerance = PEXSINumElectronTolerance
!
call f_ppexsi_load_real_symmetric_hs_matrix(&
      plan,&
      options,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,&
      rowindLocal,&
      HnzvalLocal,&
      isSIdentity,&
      SnzvalLocal,&
      info) 

if (mpirank == 0) then
   print *, "Info in load_real_sym_hs_matrix: ", info
endif

if (iscf == 1) then
   ! This is only needed for inertia-counting
   call f_ppexsi_symbolic_factorize_real_symmetric_matrix(&
        plan, &
        options,&
        info)
   if (mpirank == 0) then
      print *, "Info in real symb_fact in iscf==1: ", info
   endif

   call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
        plan, &
        options,&
        info)
   if (mpirank == 0) then
      print *, "Info in complex symb_fact in iscf==1: ", info
   endif
endif
options%isSymbolicFactorize = 0 ! We do not need it anymore
!
!  do actual solve
!
call timer("pexsi-solver", 1)

! Use inertia counts?
isInertiaCount = fdf_get("PEXSI.inertia-count",1)
! For how many scf steps?
numInertiaCounts = fdf_get("PEXSI.inertia-counts",3)

if (need_inertia_counting()) then
   options%isInertiaCount = 1
   ! Stop inertia count if mu has not changed much from iteration to iteration.
   options%muInertiaTolerance =  &
             fdf_get("PEXSI.inertia-mu-tolerance",0.05_dp,"Ry")
   ! One-sided expansion of interval if correct mu falls outside it
   options%muInertiaExpansion =  &
             fdf_get("PEXSI.lateral-expansion-inertia",3.0_dp*eV,"Ry") 

   call get_bracket_for_inertia_count( )  
   options%muMin0 = muMin0
   options%muMax0 = muMax0
   if (mpirank == 0) then
     write (6,"(a,2f9.4,a)") 'Calling inertiaCount: [', &
                                      muMin0/eV, muMax0/eV, "] (eV)"
   endif

else
   ! Will call the PEXSI solver directly
   options%isInertiaCount = 0
   call get_bracket_for_solver()
   options%muMin0 = muMin0
   options%muMax0 = muMax0
   options%mu0       = mu
   if(mpirank == 0) then
     write (6,"(a,2f9.4,a,f9.4,a,f9.5)") 'Calling solver directly (eV): [', &
                                      muMin0/eV, &
                                      muMax0/eV, &
                                     "] prev. mu: ", mu/eV, &
                                     ' Tol: ', PEXSINumElectronTolerance
   endif

endif
      
call f_ppexsi_dft_driver(&
  plan,&
  options,&    ! includes muMin0, muMax0
  numElectronExact,& ! in
  mu,&             ! out
  numElectron,&    ! out
  muMinInertia,&   ! out
  muMaxInertia,&   ! out
  numTotalInertiaIter,&  ! out
  numTotalPEXSIIter,&    ! out
  info)

!! AG
!  if (numElectron not near numElectronExact) ...

if( info .ne. 0 ) then
	call mpi_finalize( ierr )
	call exit(info)
endif

if( PEXSI_worker ) then
  call f_ppexsi_retrieve_real_symmetric_dft_matrix(&
    plan,&
    DMnzvalLocal,&
    EDMnzvalLocal,&
    FDMnzvalLocal,&
    totalEnergyH,&
    totalEnergyS,&
    totalFreeEnergy,&
    info)

!!$  ! These are the "band-structure" (free)-energies  
!!$  if( mpirank == 0 ) then
!!$    write(*,*) "Output from the main program."
!!$    write(*,*) "Total energy (H*DM)         = ", totalEnergyH
!!$    write(*,*) "Total energy (S*EDM)        = ", totalEnergyS
!!$    ! This is computed as the trace of (S*FDM)     
!!$    write(*,*) "Total free energy           = ", totalFreeEnergy
!!$  endif

endif

!------------ End of solver step

   if (mpirank == 0) then
      write(6,"(a,i3)") " #&s Number of solver iterations: ", numTotalPEXSIIter
      write(6,"(a,i3)") " #&s Number of inertia iterations: ", numTotalInertiaIter
      write(6,"(a,f12.4)") " #&s muMinInertia: ", muMinInertia
      write(6,"(a,f12.4)") " #&s muMaxInertia: ", muMaxInertia
      write(6,"(a,f12.5,f12.4,2x,a2)") "mu, N_e:", mu/eV, &
              numElectron, "&s"
   endif

if (PEXSI_worker) then

   free_bs_energy = totalFreeEnergy
   bs_energy = totalEnergyS
   eBandH = totalEnergyH

   call de_alloc(FDMnzvalLocal,"FDMnzvalLocal","pexsi_solver")
   call de_alloc(colPtrLocal,"colPtrLocal","pexsi_solver")

   if( mpirank == 0 ) then
      write(*, *) "mu (eV)       = ", mu/eV
      write(*, *) "muMinInertia    = ", muMinInertia/eV
      write(*, *) "muMaxInertia    = ", muMaxInertia/eV
      write(*, *) "numElectron   = ", numElectron
      write(*, *) "eBandS (eV) = ", bs_energy/eV
      write(*, *) "eBandH (eV) = ", eBandH/eV
      write(*, *) "freeBandEnergy (eV) = ", (free_bs_energy)/eV
      write(*, *) "eBandS (Ry) = ", bs_energy
      write(*, *) "eBandH (Ry) = ", eBandH
      write(*, *) "freeBandEnergy (Ry) = ", (free_bs_energy)

   endif

   ef = mu
   ! Note that we use the S*EDM version of the band-structure energy
   ! to estimate the entropy, by comparing it to S*FDM This looks
   ! consistent, but note that the EDM is not used in Siesta to
   ! estimate the total energy, only the DM (via the density) (that
   ! is, the XC and Hartree correction terms to Ebs going into Etot
   ! are estimated using the DM)

   Entropy = - (free_bs_energy - bs_energy) / temp

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

   call timer("pexsi-solver", 2)
   call timer("pexsi", 2)

endif


call de_alloc(shiftList,"shiftList","pexsi_solver")
call de_alloc(inertiaList,"shiftList","pexsi_solver")

call delete(dist1)
call delete(dist2)

! Step 3. Clean up */

! We cannot finalize now if we are going to reuse
! the plan in subsequent iterations...
! We need an extra module to take care of this

if (PEXSI_worker) then
   call MPI_Comm_Free(PEXSI_Comm, ierr)
   call MPI_Group_Free(PEXSI_Group, ierr)
endif


#endif 

CONTAINS
!
! This routine encodes the heuristics to compute the
! tolerance dynamically.
!
subroutine get_on_the_fly_tolerance(dDmax,tolerance)
real(dp), intent(in)  :: dDmax
real(dp), intent(out) :: tolerance

real(dp) :: tolerance_preconditioner
real(dp) :: tolerance_target_factor, tolerance_exp
real(dp), save :: previous_tolerance
logical :: new_algorithm

new_algorithm = fdf_get("PEXSI.dynamical-tolerance",.false.)
!
!
if (new_algorithm) then

!   By default, the tolerance goes to the (minimum) target 
!   at a level 5 times dDtol

   tolerance_target_factor = fdf_get("PEXSI.tolerance-target-factor",5.0_dp)

!
!  This can range in a (0.5,2.0) interval, approximately

   tolerance_preconditioner = fdf_get("PEXSI.tolerance-preconditioner",1.0_dp)

   if (scf_step > 1 ) then

      tolerance_exp = log10(dDmax/(tolerance_target_factor*dDtol))
      ! 
  !   range = log10(PEXSINumElectronToleranceMax/PEXSINumElectronToleranceMin)
      tolerance_exp = max(tolerance_exp,0.0_dp)*tolerance_preconditioner
      tolerance = PEXSINumElectronToleranceMin * 10.0_dp**tolerance_exp
      tolerance = min(tolerance,PEXSINumElectronToleranceMax)

      if (tolerance > previous_tolerance) then
         if (mpirank==0) write(6,"(a,f10.2)") &
              "Will not raise PEXSI solver tolerance to: ", &
              tolerance
         tolerance = previous_tolerance
      endif
      previous_tolerance = tolerance
   else
      ! No heuristics for now for first step
      ! Note that this should really change in MD or geometry optimization
      previous_tolerance = huge(1.0_dp)
      tolerance = PEXSINumElectronToleranceMax

   endif
else
   tolerance = Max(PEXSINumElectronToleranceMin, &
                              Min(dDmax*1.0, PEXSINumElectronToleranceMax))
endif

if (mpirank==0) write(6,"(a,f10.2)") &
     "Current PEXSI solver tolerance: ", tolerance

end subroutine get_on_the_fly_tolerance

!------------------------------------------------------------------
! This function will determine whether an initial inertia-counting
! stage is needed, based on user input and the level of convergence
!
! Variables used through host association for now:
!
!      isInertiaCount
!      numInertiaCounts
!      scf_step
!      prevDmax, safe_dDmax_NoInertia
!
! Some logging output is done, so this function is not pure.

function need_inertia_counting() result(do_inertia_count)
logical :: do_inertia_count

real(dp) :: safe_dDmax_NoInertia

safe_dDmax_NoInertia = fdf_get("PEXSI.safe-dDmax-no-inertia",0.05)

do_inertia_count = .false.

if (isInertiaCount .ne. 0) then
  if (scf_step .le. numInertiaCounts) then
     if (mpirank == 0) write(6,"(a,i4)") "&o Inertia-count step scf_step<numIC \
", scf_step
     do_inertia_count = .true.
  endif
  if (numInertiaCounts < 0) then
     if (scf_step <= -numInertiaCounts) then
        if (mpirank == 0) write(6,"(a,i4)") "&o Inertia-count step scf_step<-nu\
mIC ", scf_step
        do_inertia_count = .true.
     else if (prevDmax > safe_dDmax_NoInertia) then
        if (mpirank == 0) write(6,"(a,i4)") "&o Inertia-count step as prevDmax \
> safe_Dmax ", scf_step
        do_inertia_count = .true.
     endif
  endif
endif
end function need_inertia_counting

!---------------------------------------------------------------
!  Chooses the proper interval for the call to the driver
!  in case we need a stage of inertia counting  
!
subroutine get_bracket_for_inertia_count()

 real(dp)       :: safe_width_ic
 real(dp)       :: safe_dDmax_Ef_inertia

 safe_width_ic = fdf_get("PEXSI.safe-width-ic-bracket",4.0_dp*eV,"Ry")
 safe_dDmax_Ef_Inertia = fdf_get("PEXSI.safe-dDmax-ef-inertia",0.1)

 ! Proper bracketing                                                           
 if (scf_step > 1) then
   if (prevDmax < safe_dDmax_Ef_inertia) then
      ! Shift brackets using estimate of Ef change from previous iteration      
      !                                                                         
      if (mpirank == 0) write(6,"(a)") "&o Inertia-count bracket shifted by Delta_Ef"
      ! This might be risky, if the final interval of the previous iteration    
      ! is too narrow. We should broaden it by o(kT)                            
      ! The usefulness of delta_Ef is thus debatable...                         

      muMin0 = muMinInertia + delta_Ef - two_kT
      muMax0 = muMaxInertia + delta_Ef + two_kT
   else
      ! Use a large enough interval around the previous estimation of           
      ! mu (the gap edges are not available...)  
      if (mpirank == 0) write(6,"(a)") "&o Inertia-count safe bracket"
!      muMin0 = min(muLowerEdge - 0.5*safe_width_ic, muMinInertia)
      muMin0 = min(mu - 0.5*safe_width_ic, muMinInertia)
!      muMax0 = max(muUpperEdge + 0.5*safe_width_ic, muMaxInertia)
      muMax0 = max(mu + 0.5*safe_width_ic, muMaxInertia)
   endif
 else
    if (mpirank == 0) write(6,"(a)") "&o Inertia-count called with iscf=1 parameters"
 endif
end subroutine get_bracket_for_inertia_count

subroutine get_bracket_for_solver()

    real(dp)       :: safe_width_solver
    real(dp)       :: safe_dDmax_Ef_solver

safe_width_solver = fdf_get("PEXSI.safe-width-solver-bracket",2.0_dp*eV,"Ry")
safe_dDmax_Ef_solver = fdf_get("PEXSI.safe-dDmax-ef-solver",0.05)

! Do nothing for now
! No setting of  muMin0 and muMax0 yet, pending clarification of flow

  if (scf_step > 1) then
     if (prevDmax < safe_dDmax_Ef_solver) then
        if (mpirank == 0) write(6,"(a)") "&o Solver mu shifted by delta_Ef"
        mu = mu + delta_Ef
     endif
     ! Always provide a safe bracket around mu, in case we need to fallback
     ! to executing a cycle of inertia-counting
     if (mpirank == 0) write(6,"(a)") "&o Safe solver bracket around mu"
     muMin0 = mu - safe_width_solver
     muMax0 = mu + safe_width_solver
  else
     if (mpirank == 0) write(6,"(a)") "&o Solver called with iscf=1 parameters"
     ! do nothing. Keep mu, muMin0 and muMax0 as they are inherited
  endif
end subroutine get_bracket_for_solver
!------------------------------------------------------
! If using the "annealing" feature, this routine computes
! the current temperature to use in the PEXSI solver
!
subroutine get_current_temperature(pexsi_temperature)
  real(dp), intent(out) :: pexsi_temperature

 logical  :: use_annealing
 real(dp) :: annealing_preconditioner, temp_factor
 real(dp) :: annealing_target_factor

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
end subroutine get_current_temperature

end subroutine pexsi_solver

end module m_pexsi_solver
