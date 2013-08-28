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
       ef, Entrop, temp)

    use fdf
    use parallel, only   : SIESTA_worker, BlockSize
    use parallel, only   : SIESTA_Group, SIESTA_Comm
    use m_mpi_utils, only: globalize_sum, globalize_max
    use m_mpi_utils, only: broadcast
    use units,       only: Kelvin, eV
    use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
    use class_Dist
    use alloc,             only: re_alloc, de_alloc
#ifdef MPI
    use mpi_siesta
#endif
    implicit          none

    integer, intent(in)  :: iscf  ! scf step number
    integer, intent(in)  :: maxnh, no_u, no_l, nspin
    integer, intent(in), target  :: listh(maxnh), numh(no_l), listhptr(no_l)
    real(dp), intent(in), target :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: qtot
    real(dp), intent(out), target:: DM(maxnh,nspin), EDM(maxnh,nspin)
    real(dp), intent(out)        :: ef  ! Fermi energy
    real(dp), intent(out)        :: Entrop ! Entropy/k, dimensionless
    real(dp), intent(in)         :: temp   ! Electronic temperature

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
    real(dp)       :: muSolverInput, muMinSolverInput, muMaxSolverInput
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
real(dp) :: muInertia, muMinInertia, muMaxInertia, muLowerEdge, muUpperEdge
real(dp) :: muMinPEXSI, muMaxPEXSI
integer  :: muMaxIter
integer  :: npPerPole
integer  :: npSymbFact
integer  :: mpirank, ierr
integer  :: isSIdentity
integer  :: inertiaMaxIter, inertiaIter
real(dp) :: inertiaNumElectronTolerance, &
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

interface
 ! subroutine f_ppexsi_solve_interface
   include "pexsi_solve.h"
 end subroutine f_ppexsi_solve_interface
end interface

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

   if (mpirank==0) write(6,"(a,g12.5,a,f10.2)") &
                          "Electronic temperature: ", temperature, &
                          ". In Kelvin:", temperature/Kelvin

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

endif ! PEXSI worker

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

! How to expand the intervals in case of need. Default: ~ 3 eV
lateral_expansion_solver = fdf_get("PEXSI.lateral-expansion-solver",0.2_dp,"Ry")
lateral_expansion_inertia = fdf_get("PEXSI.lateral-expansion-inertia",0.2_dp,"Ry")


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
   !
   ! Here we have to decide what to do with the previous 
   ! step's muMin0 and muMax0:
   !
   ! - Use a rigid shift based on Tr(DM*DeltaH)
   ! - Do nothing and just inherit the interval
   ! - Do something else.
   ! 
   ! For now, do nothing, but as a safeguard 
   ! use numInertiaCounts > 1 or 2... below
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

! Stop inertia count if Ne(muMax) - Ne(muMin) < inertiaNumElectronTolerance
inertiaNumElectronTolerance = fdf_get("PEXSI.inertia-num-electron-tolerance",10)

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

!
!  Possible inertia-count stage
!
if ((isInertiaCount .ne. 0) .and. (scf_step .le. numInertiaCounts)) then

  call do_inertia()

  muSolverInput = muInertia
  muMinSolverInput = muMinInertia
  muMaxSolverInput = muMaxInertia

  ! re-define the starting interval for successive steps
  ! this is the most reasonable starting point, pending
  ! further refining (see above)

  muMin0 = muMinInertia
  muMax0 = muMaxInertia

else !no inertia count

  ! Use the starting values, or the output of previous inertia-count iterations
  ! (see above)

  muSolverInput = mu
  muMinSolverInput = muMin0
  muMaxSolverInput = muMax0

end if

!
!  do actual solve
!

solver_loop: do

   if(mpirank == 0) then
     write (6,"(a,f12.5,a,f12.5,a,f12.5,a,a,f8.5)") 'Calling PEXSI solver. mu: ', &
                                              muSolverInput/eV, &
                                              ' in [',muMinSolverInput/eV, &
                                              ',', muMaxSolverInput/eV, "] (eV)", &
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
        temperature,&
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
      if (mpirank == 0) then
         write(6,"(a,2f14.6,a,f14.6)") &
            "The PEXSI solver did not converge well. Nel, Nel_exact: ", &
            numElectron, numElectronExact, &
            " Requested tol: ", PEXSInumElectronTolerance
         write(6,"(a)") "You need a better interval. Computing inertia-counts..."
      endif
      ! we need to do an inertia-count step...
      ! check which side of the interval we need to expand
      ! Expand by 0.2 Ry ~ 3 eV
      !
      if ((numElectron-numElectronExact) > 0) then
         muMin0 = muMin0 - lateral_expansion_solver !0.2_dp 
      else
         muMax0 = muMax0 + lateral_expansion_solver ! 0.2_dp
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

if (fdf_get("PEXSI.InheritSolverInterval",.false.)) then
   ! re-define the starting interval		       
   muMin0 = muMinPEXSI
   muMax0 = muMaxPEXSI
endif

!------------ End of solver step

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
      write(*,*) "Number of mu iterations: ", muIter
      write(*,"(a3,2a12,a20)") "it", "mu", "N_e", "dN_e/dmu"
      do i = 1, muIter
         write(*,"(i3,2f12.4,g20.5)") i, muList(i)/eV, &
                 numElectronList(i), numElectronDrvList(i)*eV
      end do
   endif

   ef = mu
   Entrop = - (free_bs_energy - ebandStructure) / temp

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
   call broadcast(Entrop,comm=SIESTA_Comm)
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

    logical        :: interval_problem
    logical        :: bad_lower_bound
    logical        :: bad_upper_bound

interface
 ! subroutine f_ppexsi_inertiacount_interface
   include "pexsi_inertia.h"
 end subroutine f_ppexsi_inertiacount_interface
end interface


search_interval: do

   if(mpirank == 0) then
     write (6,"(a,f12.5,a,f12.5,a,a,i4)") 'Calling inertia_count with interval: [', &
                                      muMin0/eV, ",", muMax0/eV, "] (eV)", &
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
        temperature,&
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

   call timer("pexsi-inertia-ct", 2)

   call globalize_max(info,infomax,comm=World_Comm)
   info = infomax

   interval_problem = .false.
   if (info /= 0) then
!       call mpi_gather(inertiaList(1),1,MPI_Double_Precision, lower_inertia(:),1, &&
!                       MPI_Double_Precision,0,World_Comm, ierr)
!       call mpi_gather(inertiaList(nptsInertia),1,MPI_Double_Precision, upper_inertia(:),1, &&
!                       MPI_Double_Precision,0,World_Comm, ierr)
!      write(6,"(a,i4,2f12.4,4x,2f12.4)") "Node, lb, ub, lin, uin:", &
!               mpirank, shiftlist(1), shiftlist(nptsInertia), inertialist(1), inertialist(nptsInertia)
      if(mpirank == 0) then
         write(msg,"(i6)") info 
         write (6,"(/,a)") 'PEXSI inertia count ended in error. Info: ' // msg
         write (6,"(a,2f12.4)") 'Lower, Upper counts: ', inertiaList(1), &
                                                    inertiaList(nptsInertia)
         bad_lower_bound = (inertiaList(1) > (numElectronExact - 0.1)) 
         bad_upper_bound = (inertiaList(nptsInertia) < (numElectronExact + 0.1)) 
      endif

      call broadcast(bad_lower_bound,comm=World_Comm)
      call broadcast(bad_upper_bound,comm=World_Comm)
      call broadcast(inertiaList(1),comm=World_Comm)
      call broadcast(inertiaList(nptsInertia),comm=World_Comm)

      if (bad_lower_bound) then
         interval_problem =  .true.
         muMin0 = muMin0 - lateral_expansion_inertia ! 0.5
         if (mpirank==0) then
            write(6,"(a,f10.4)") "At the lower end, inertiaList: ", &
                                  inertiaList(1)
            write (6,"(a,2f10.4)") 'Updating the interval: ', &
                                  muMin0/eV, muMax0/eV
         endif
      endif
      if (bad_upper_bound) then
         interval_problem =  .true.
         muMax0 = muMax0 + lateral_expansion_inertia ! 0.5
         if (mpirank==0) then
            write(6,"(a,f10.4)") "At the upper end, inertiaList: ", &
                                  inertiaList(nptsInertia)
            write (6,"(a,2f10.4)") 'Updating the interval: ', &
                                  muMin0/eV, muMax0/eV
         endif
      endif

      if (interval_problem) then
          ! do nothing more
      else
         write(msg,"(i6)") info 
         call die("Non-interval error in inertia count routine. Info: " // msg)
      endif

   else
      exit search_interval
   endif

enddo search_interval

   muInertia    = (muLowerEdge + muUpperEdge) / 2d0

   if(mpirank == 0) then
     write (6,"(/,a)") 'PEXSI inertia count executed'
     write (6,"(a,f10.4)") ' mu (eV)=', muInertia/eV
     write (6,"(a,f10.4)") ' mu lowerEdge (eV):', muLowerEdge/eV
     write (6,"(a,f10.4)") ' mu upperEdge (eV):', muUpperEdge/eV
     write (6,"(a,f10.4)") ' muMin (eV):', muMinInertia/eV
     write (6,"(a,f10.4)") ' muMax (eV):', muMaxInertia/eV

     write(6,"(/,a)") "Cumulative DOS by inertia count:"
     do i=1, nptsInertia
        write(6,"(f10.4,f10.4)") shiftList(i)/eV, inertiaList(i)
     enddo
  end if
end subroutine do_inertia

end subroutine pexsi_solver
end module m_pexsi_solver
