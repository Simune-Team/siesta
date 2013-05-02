module m_pexsi_solver

  public :: pexsi_solver

CONTAINS

! This version uses the standard PEXSI distribution (by tricking Siesta into
! using it)
!
  subroutine pexsi_solver(iscf, no_u, no_l, nspin,  &
       maxnh, numh, listhptr, listh, H, S, qtot, DM, EDM, &
       ef, freeEnergyCorrection, temp)

    use precision, only  : dp
    use fdf
    use parallel, only   : worker, ionode
    use m_mpi_utils, only: globalize_sum
    use units,       only: Kelvin, eV
#ifdef MPI
    use mpi_siesta
#endif
    implicit          none

    integer, intent(in)  :: iscf  ! scf step number
    integer, intent(in)  :: maxnh, no_u, no_l, nspin
    integer, intent(in), target  :: listh(maxnh), numh(*), listhptr(*)
    real(dp), intent(in), target :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: qtot
    real(dp), intent(out), target:: DM(maxnh,nspin), EDM(maxnh,nspin)
    real(dp), intent(out)        :: ef  ! Fermi energy
    real(dp), intent(out)        :: freeEnergyCorrection
    real(dp), intent(in)         :: temp   ! Electronic temperature

#ifndef MPI
    call die("PEXSI needs MPI")
#else

    integer :: siesta_comm 
    integer :: MPIerror, stat(MPI_STATUS_SIZE), count
    integer :: bs, norbs_slack, nnz_slack
    integer :: ispin, maxnhtot, ih, nnzold, i

    integer  :: ordering, isInertiaCount, numInertiaCounts
    integer  :: muIter
    real(dp) :: muZeroT

    real(dp), save :: mu
    real(dp), save :: muMin0, muMax0
    logical, save  :: first_call = .true.
    real(dp)       :: eBandStructure, eBandH

    integer        :: info

!Lin variables
integer :: nrows, nnz, nnzLocal, numColLocal
integer, pointer, dimension(:) ::  tmpi => null()
integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
!
real(dp), pointer, dimension(:) :: &
        HnzvalLocal=>null(), SnzvalLocal=>null(),  &
        DMnzvalLocal => null() , EDMnzvalLocal => null(), &
        FDMnzvalLocal => null()
!
real(dp), allocatable, dimension(:) :: muList, &
                                       numElectronList, &
                                       numElectronDrvList, &
                                       shiftList, inertiaList
integer  :: numPole
real(dp) :: temperature, numElectronExact, numElectron, gap, deltaE
real(dp) :: muInertia, muMinInertia, muMaxInertia, muLowerEdge, muUpperEdge
real(dp) :: muMinPEXSI, muMaxPEXSI
integer  :: muMaxIter
integer  :: npPerPole, nprow, npcol
integer  :: mpirank, mpisize, ierr
integer  :: isSIdentity
integer  :: inertiaMaxIter, inertiaIter
real(dp) :: inertiaNumElectronTolerance, PEXSINumElectronTolerance
!------------

real(dp) :: buffer1

external      :: timer

interface
 ! subroutine f_ppexsi_inertiacount_interface
   include "pexsi_inertia.h"
 end subroutine f_ppexsi_inertiacount_interface
end interface

interface
 ! subroutine f_ppexsi_solve_interface
   include "pexsi_solve.h"
 end subroutine f_ppexsi_solve_interface
end interface

! "Worker" means a processor which is in the Siesta subset (which
! comprises the first npPerPole processors in this implementation)
! NOTE:  fdf calls will assign values to the whole processor set,
! but some other variables will have to be re-broadcast (see examples
! below)

! NOTE: Some comments in the code below are placeholders for a new
! version with independent Siesta/PEXSI distributions (work in progress)

if (worker) then
   siesta_comm = mpi_comm_world
   call timer("pexsi", 1)

   ispin = 1
   if (nspin /=1) then
      call die("Spin polarization not yet supported in PEXSI")
   endif

   call mpi_comm_rank( SIESTA_COMM, mpirank, ierr )
   call mpi_comm_size( SIESTA_COMM, mpisize, ierr )

   npPerPole = mpisize
   numElectronExact = qtot 

   ! Note that the energy units for the PEXSI interface are arbitrary, but
   ! H, the interval limits, and the temperature have to be in the
   ! same units. Siesta uses Ry units.

   temperature      = temp

   if (IOnode) write(6,"(a,g12.5,a,f10.2)") &
                          "Electronic temperature: ", temperature, &
                          ". In Kelvin:", temperature/Kelvin

   call MPI_Barrier(Siesta_comm,ierr)

   nrows = no_u
!  numColLocal = num_local_elements(dist2,no_u,myrank2)
   numColLocal = no_l

   ! Figure out the communication needs
   ! call analyze_comms(comms)
   ! call do_transfers(comms,numh,numh_pexsi, &
   !                   g1,g2,mpi_comm)

 !   nnzLocal = sum(numh_pexsi(1:numColLocal))
 ! call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,PEXSI_comm,MPIerror)
   nnzLocal = sum(numh(1:numColLocal))
   call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,Siesta_comm,MPIerror)


!  print *, "Node ", mpirank, ": no_l, numColLocal, maxnh, nnzLocal: ", &
!           no_l, numColLocal, maxnh, nnzLocal

  allocate(colptrLocal(1:numColLocal+1))
  colptrLocal(1) = 1
  do ih = 1,numColLocal
!    colptrLocal(ih+1) = colptrLocal(ih) + numh_pexsi(ih)
     colptrLocal(ih+1) = colptrLocal(ih) + numh(ih)
  enddo

   allocate(rowindLocal(1:nnzLocal))
   allocate(HnzvalLocal(1:nnzLocal))
   allocate(SnzvalLocal(1:nnzLocal))
   allocate(DMnzvalLocal(1:nnzLocal))
   allocate(EDMnzvalLocal(1:nnzLocal))
   allocate(FDMnzvalLocal(1:nnzLocal))

!  call generate_commsnnz(comms,commsnnz)
!  call do_transfers(commsnnz,listh,rowindLocal, &
!                        g1, g2, mpi_comm)
!  call do_transfers(commsnnz,H,HnzvalLocal, &
!                        g1, g2, mpi_comm)
!  call do_transfers(commsnnz,S,SnzvalLocal, &
!                        g1, g2, mpi_comm)

   rowindLocal(1:nnzLocal) = listh(1:nnzLocal)
   HnzvalLocal(1:nnzLocal) = H(1:nnzLocal,ispin)
   SnzvalLocal(1:nnzLocal) = S(1:nnzLocal)

endif ! worker

isSIdentity = 0

numPole          = fdf_get("PEXSI.num-poles",20)
gap              = fdf_get("PEXSI.gap",0.0d0)

! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.
deltaE           = fdf_get("PEXSI.delta-E",3.0d0)

! Initial guess of chemical potential and containing interval
! When using inertia counts, this interval can be wide.
! Note that muMin0 and muMax0 are saved variables
if (first_call) then
   mu = fdf_get("PEXSI.mu",-0.60_dp)
   ! Lower/Upper bound for the chemical potential.
   muMin0           = fdf_get("PEXSI.mu-min",-1.0d0)
   muMax0           = fdf_get("PEXSI.mu-max", 0.0d0)
   first_call = .false.
endif

! Maximum number of iterations for computing the inertia
inertiaMaxIter   = fdf_get("PEXSI.inertia-max-iter",3)

! muMaxIter should be 1 or 2 later when combined with SCF.
muMaxIter        = fdf_get("PEXSI.mu-max-iter",10)

! Arrays for reporting back information about the PEXSI iterations
allocate( muList( muMaxIter ) )
allocate( numElectronList( muMaxIter ) )
allocate( numElectronDrvList( muMaxIter ) )

! Arrays for reporting back information about the integrated DOS
! computed by the inertia count method. Since we use the processor
! teams corresponding to the different poles, the number of points
! in the energy interval is "numPole"

allocate( shiftList( numPole ) )
allocate( inertiaList( numPole ) )

! Stop inertia count if Ne(muMax) - Ne(muMin) < inertiaNumElectronTolerance
inertiaNumElectronTolerance = fdf_get("PEXSI.inertia-num-electron-tolerance",1d-1)

! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
PEXSINumElectronTolerance = fdf_get("PEXSI.num-electron-tolerance",1d-1)

! Number of processors used for each pole. At the moment use mpisize.
! Later can be changed to 
!npPerPole        = fdf_get("PEXSI.np-per-pole",mpisize)

! Ordering flag
ordering = fdf_get("PEXSI.ordering",1)
isInertiaCount = fdf_get("PEXSI.inertia-count",1)
numInertiaCounts = fdf_get("PEXSI.inertia-counts",1)
!
! Broadcast these to the whole processor set, just in case
!
call MPI_Bcast(npPerPole,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(nrows,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(nnz,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(numElectronExact,1,MPI_double_precision,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(temperature,1,MPI_double_precision,0,true_MPI_COMM_world,ierr)

call MPI_Bcast(ordering,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(isInertiaCount,1,MPI_integer,0,true_MPI_COMM_world,ierr)

!
!  do inertia count
!

if ((isInertiaCount .ne. 0) .and. (iscf .le. numInertiaCounts)) then

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
        numPole,&
        inertiaMaxIter,&
        inertiaNumElectronTolerance,&
        ordering,&
        npPerPole,&
        true_MPI_COMM_WORLD,&
! output parameters
        muMinInertia,&
        muMaxInertia,&
        muLowerEdge,&
        muUpperEdge,&
        inertiaIter,&
        shiftList,&
        inertiaList,&
        info)

   if (info /= 0) then
      call die("Error in inertia count routine")
   endif

   muInertia    = (muLowerEdge + muUpperEdge) / 2d0

   if( mpirank == 0 ) then
     write (*,*) 'PEXSI inertia count executed, mu=', muInertia, &
            'lowerEdge:', muLowerEdge, 'upperEdge:', muUpperEdge, &
            'muMin:', muMinInertia, 'muMax:', muMaxInertia
     write(6,"(a)") "Cumulative DOS by inertia count:"
     do i=1, numPole
        write(6,"(f10.4,f10.4)") shiftList(i)/eV, inertiaList(i)
     enddo
  end if

else !no inertia count

  muInertia = mu
  muMinInertia = muMin0
  muMaxInertia = muMax0

end if

!
!  do actual solve
!

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
        muInertia,&
        muMinInertia,&
        muMaxInertia,&
        gap,&
        deltaE,&
        numPole,&
        muMaxIter,&
        PEXSINumElectronTolerance,&
        ordering,&
        npPerPole,&
        true_MPI_COMM_WORLD,&
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

   if (info /= 0) then
      call die("Error in pexsi solver routine")
   endif

! save the mu-range for the next run
! If we do not do this, we will use always the latest inertia-count
! guess, or the initial (muMin0, muMax0) bracket.

!    muMin0 = muMinPEXSI
!    muMax0 = muMaxPEXSI

!  call reverse_comms(commsnnz,commsnnz_reverse)
!  call do_transfers(commsnnz_reverse,DMnzvalLocal,DM, &
!                        g2, g1, mpi_comm)
!  call do_transfers(commsnnz_reverse,EDMnzvalLocal,EDM, &
!                        g2, g1, mpi_comm)
!  allocate(FDM(1:maxnh))
!  call do_transfers(commsnnz_reverse,FDMnzvalLocal,FDM, &
!                        g2, g1, mpi_comm)

if (worker) then

   ! Watch out here. We have to assume that rank 0 is common to both
   ! subgroups
   ! if (mpirank == 0) then
   !  ef = mu
   !  call MPI_Bcast(ef, .... siesta subgroup)
   ! endif

   ef = mu

   ! this will disappear
   DM(1:nnzLocal,ispin) = DMnzvalLocal(1:nnzLocal)
   EDM(1:nnzLocal,ispin) = EDMnzvalLocal(1:nnzLocal)

   freeEnergyCorrection = 0.0_dp
   eBandStructure = 0.0_dp
   eBandH = 0.0_dp
   do i = 1,nnzLocal
      ! Watch out for FDM nzvalLocal (see above)
      freeEnergyCorrection = freeEnergyCorrection + SnzvalLocal(i) * &
           ( FDMnzvalLocal(i) - EDMnzvalLocal(i) )
      eBandStructure = eBandStructure + SnzvalLocal(i) * &
           ( EDMnzvalLocal(i) )
      eBandH = eBandH + HnzvalLocal(i) * &
           ( DMnzvalLocal(i) )
   enddo

   ! These operations in siesta group (implicit communicator)
   call globalize_sum( freeEnergyCorrection, buffer1 )
   freeEnergyCorrection = buffer1 + mu*numElectron
   call globalize_sum( eBandStructure, buffer1 )
   eBandStructure = buffer1
   call globalize_sum( eBandH, buffer1 )
   eBandH = buffer1

   if( mpirank == 0 ) then
      write(*, *) "mu            = ", mu
      write(*, *) "muMinInertia  = ", muMinInertia
      write(*, *) "muMaxInertia  = ", muMaxInertia
      write(*, *) "muMinPEXSI    = ", muMinPEXSI
      write(*, *) "muMaxPEXSI    = ", muMaxPEXSI
      write(*, *) "mu (eV)       = ", mu/eV
      write(*, *) "muZeroT (eV)  = ", muZeroT/eV
      write(*, *) "numElectron   = ", numElectron
      write(*, *) "eBandStructure (Ry) = ", eBandStructure
      write(*, *) "eBandStructure (eV) = ", eBandStructure/eV
      write(*, *) "eBandH (eV) = ", eBandH/eV
      write(*, *) "freeEnergy (eV) = ", (eBandStructure + freeEnergyCorrection)/eV
      write(*,*) "Number of mu iterations: ", muIter
      write(*,"(a3,2a12,a20)") "it", "mu", "N_e", "dN_e/dmu"
      do i = 1, muIter
         write(*,"(i3,2f12.4,g20.5)") i, muList(i)/eV, &
                 numElectronList(i), numElectronDrvList(i)*eV
      end do
   endif

  deallocate(rowindLocal)
  deallocate(HnzvalLocal)
  deallocate(SnzvalLocal)
  deallocate(DMnzvalLocal)
  deallocate(EDMnzvalLocal)
  deallocate(FDMnzvalLocal)
  deallocate(colPtrLocal)
! deallocate(FDM)

    call timer("pexsi", 2)
 endif ! worker

 deallocate( shiftList )
 deallocate( inertiaList )

#endif 

end subroutine pexsi_solver
end module m_pexsi_solver
