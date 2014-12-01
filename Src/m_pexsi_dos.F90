module m_pexsi_DOS
    use precision, only  : dp

  public :: pexsi_DOS

CONTAINS

! This version uses separate distributions 
! for Siesta (setup_H et al) and PEXSI.
!
  subroutine pexsi_DOS(no_u, no_l, nspin,  &
       maxnh, numh, listhptr, listh, H, S, qtot, ef)

    use fdf
    use parallel, only   : SIESTA_worker, BlockSize
    use parallel, only   : SIESTA_Group, SIESTA_Comm
    use m_mpi_utils, only: globalize_max
    use m_mpi_utils, only: broadcast
    use units,       only: eV
    use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
    use class_Dist
    use alloc,             only: re_alloc, de_alloc
#ifdef MPI
    use mpi_siesta
#endif

    use iso_c_binding
    use f_ppexsi_interface, only: f_ppexsi_options
    use f_ppexsi_interface, only: f_ppexsi_plan_finalize
    use f_ppexsi_interface, only: f_ppexsi_plan_initialize
    use f_ppexsi_interface, only: f_ppexsi_inertia_count_real_symmetric_matrix
    use f_ppexsi_interface, only: f_ppexsi_load_real_symmetric_hs_matrix
    use f_ppexsi_interface, only: f_ppexsi_set_default_options
    use f_ppexsi_interface, &
          only: f_ppexsi_symbolic_factorize_real_symmetric_matrix

    implicit          none

    integer, intent(in)  :: maxnh, no_u, no_l, nspin
    integer, intent(in), target  :: listh(maxnh), numh(no_l), listhptr(no_l)
    real(dp), intent(in), target :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in)        :: qtot ! Total number of electrons
    real(dp), intent(in)        :: ef  ! Fermi energy

#ifndef MPI
    call die("PEXSI needs MPI")
#else

    integer :: PEXSI_Comm, World_Comm
    integer :: PEXSI_Group, World_Group

    integer :: ispin, maxnhtot, ih, i

    integer  :: ordering, numNodesTotal


    real(dp)   :: e1, e2, delta
    integer    :: j, ncalls, npoints, lun
    real(dp), allocatable :: edos(:)
    real(dp),  allocatable :: intdos(:)
    real(dp)   :: emin, emax, deltaE

    integer        :: info, infomax

!Lin variables
integer :: nrows, nnz, nnzLocal, numColLocal

integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
real(dp), pointer, dimension(:) :: HnzvalLocal=>null(), SnzvalLocal=>null()

logical  :: PEXSI_worker
integer  :: nShifts

integer  :: npPerPole
integer  :: npSymbFact
integer  :: mpirank, ierr
integer  :: isSIdentity
logical  :: ef_reference

type(f_ppexsi_options) :: options
integer(c_intptr_t)    :: plan
  integer :: numProcRow, numProcCol
  integer :: outputFileIndex


!------------
external         :: timer
character(len=6) :: msg

type(aux_matrix) :: m1, m2
type(Dist)       :: dist1, dist2
integer          :: pbs, norbs


! "SIESTA_Worker" means a processor which is in the Siesta subset.
! NOTE:  fdf calls will assign values to the whole processor set,
! but some other variables will have to be re-broadcast (see examples
! below)

World_Comm = true_MPI_Comm_World

if (SIESTA_worker) then
   ! deal with intent(in) variables
   norbs = no_u
endif
call broadcast(norbs,comm=World_Comm)

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

! -- Prepare plan
numProcRow = sqrt(dble(npPerPole))
numProcCol = numProcRow

if ((numProcRow * numProcCol) /= npPerPole) then
  call die("not perfect square")
endif

outputFileIndex = mpirank

   plan = f_ppexsi_plan_initialize(&
      World_Comm,&
      numProcRow,&
      numProcCol,&
      outputFileIndex,&
      info) 

if (mpirank == 0) then
   print *, "Info in plan_initialize: ", info
endif
! -- 

pbs = norbs/npPerPole
call newDistribution(pbs,PEXSI_Group,dist2,TYPE_PEXSI,"px dist")


if (SIESTA_worker) then
   call timer("pexsi_dos", 1)

   ispin = 1
   if (nspin /=1) then
      call die("Spin polarization not yet supported in PEXSI")
   endif

   ! Note that the energy units for the PEXSI interface are arbitrary, but
   ! H, the interval limits, and the temperature have to be in the
   ! same units. Siesta uses Ry units.


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
! Will allocate m2
call redistribute_spmatrix(norbs,m1,dist1,m2,dist2,World_Comm)
call timer("redist_orbs_fwd", 2)

if (PEXSI_worker) then

   nrows = m2%norbs          ! or simply 'norbs'
   numColLocal = m2%no_l
   nnzLocal    = m2%nnzl
   call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,PEXSI_Comm,ierr)

  call re_alloc(colptrLocal,1,numColLocal+1,"colptrLocal","pexsi_DOS")
  colptrLocal(1) = 1
  do ih = 1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + m2%numcols(ih)
  enddo

  rowindLocal => m2%cols
  SnzvalLocal => m2%vals(1)%data
  HnzvalLocal => m2%vals(2)%data

endif ! PEXSI worker

isSIdentity = 0

! Since we use the processor teams corresponding to the different poles, 
! the number of points in the energy interval should be a multiple
! of numNodesTotal/npPerPole.

call mpi_comm_size( World_Comm, numNodesTotal, ierr )
nShifts = numNodesTotal/npPerPole

!
! Broadcast these to the whole processor set, just in case
! (They were set only by the Siesta workers)
!
call MPI_Bcast(nrows,1,MPI_integer,0,World_Comm,ierr)
call MPI_Bcast(nnz,1,MPI_integer,0,World_Comm,ierr)
call MPI_Bcast(ef,1,MPI_double_precision,0,World_Comm,ierr)

  ! Absolute range in principle
  emin = fdf_get("PEXSI.DOS.emin",-1.0_dp,"Ry")
  emax = fdf_get("PEXSI.DOS.emax",+1.0_dp,"Ry")

  ! The range is around fermi level by default
  ef_reference = fdf_get("PEXSI.DOS.Ef.Reference",.true.)
  if (ef_reference) then
     emin = emin + ef
     emax = emax + ef
  endif

  npoints = fdf_get("PEXSI.DOS.npoints",200)

  ! Make it a multiple of the number of samples per call
  npoints = ceiling(dble(npoints)/nShifts) * nShifts

  allocate (edos(npoints),intdos(npoints))

    if(mpirank == 0) then
       write (6,"(a,f12.5,a,f12.5,a,a,i4)")  &
                    'Calling inertia_count for DOS: [', &
                    emin/eV, ",", emax/eV, "] (eV)", &
                    " Nshifts: ", npoints
    endif

call f_ppexsi_set_default_options( options )
! Ordering flag:
!   1: Use METIS
!   0: Use PARMETIS/PTSCOTCH
options%ordering = fdf_get("PEXSI.ordering",1)

! Number of processors for symbolic factorization
! Only relevant for PARMETIS/PT_SCOTCH
options%npSymbFact = fdf_get("PEXSI.np-symbfact",npPerPole)

options%verbosity = fdf_get("PEXSI.verbosity",1)

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

   call f_ppexsi_symbolic_factorize_real_symmetric_matrix(&
        plan, &
        options,&
        info)
   if (mpirank == 0) then
      print *, "Info in real symb_fact in iscf==1: ", info
   endif

    
    delta = (emax-emin)/(npoints-1)
    do j = 1, npoints
       edos(j) = emin + (j-1)*delta
    enddo

    call timer("pexsi-raw-inertia-ct", 1)
    call f_ppexsi_inertia_count_real_symmetric_matrix(&
      plan,&
      options,&
      npoints,&
      edos,&
      intdos,&
      info) 

    call timer("pexsi-raw-inertia-ct", 2)

    if(mpirank == 0) then
       if (info /= 0) then
          write(6,*) "DOS call Info : ", info
          call die("Error exit from raw-inertia-count routine")
       endif
      call pxfflush(6)
    endif	
    
!    if(mpirank == 0) then

 if (mpirank == 0) then
    call io_assign(lun)
    open(unit=lun,file="PEXSI_INTDOS",form="formatted",status="unknown", &
         position="rewind",action="write")
    write(lun,"(2f15.6,i6,a)") ef/eV, qtot, npoints, &
                "# (Ef, qtot, npoints) / npoints lines: E(eV), IntDos(E)"
    do j=1,npoints
       write(lun,"(f15.6,f15.2)") edos(j)/eV, intdos(j)
    enddo
 endif

deallocate(edos,intdos)

   if (SIESTA_worker) then
      call timer("pexsi_dos", 2)
     ! deallocate(m1%vals)
   endif

call delete(dist1)
call delete(dist2)

call f_ppexsi_plan_finalize( plan, info )

if (PEXSI_worker) then

    call de_alloc(colptrLocal,"colptrLocal","pexsi_DOS")

    call de_alloc(m2%numcols,"m2%numcols","m_pexsi_dos")
    call de_alloc(m2%cols,"m2%cols","m_pexsi_dos")
    do j=1,size(m2%vals)
       call de_alloc(m2%vals(j)%data,"m2%vals(j)%data","m_pexsi_dos")
    enddo
    deallocate(m2%vals)

endif

if (PEXSI_worker) then
   call MPI_Comm_Free(PEXSI_Comm, ierr)
   call MPI_Group_Free(PEXSI_Group, ierr)
endif

#endif 


end subroutine pexsi_DOS
end module m_pexsi_DOS
