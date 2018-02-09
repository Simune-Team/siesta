!
! Alberto Garcia, Feb. 2018
!
! ELSI DM-based interface to Siesta. It uses the sparse matrices from
! Siesta (re-distributed into the PEXSI CSC format by custom routines
! taken from the legacy PEXSI interface), and obtains the DM and EDM
! matrices in sparse form.
!
! The elsi_solver routine is in principle able to perform (spin-polarized)
! calculations for real matrices (i.e., at the Gamma point). 
!
! It drives ELPA now, but the same ideas can be used (with minimal changes)
! for the PEXSI solver.
!
! The structure of the solver routine is such that it will detect when
! it is being run for the first scf step, so it can perform any needed
! initialization. The same idea can be used for the diagonalization
! mode, with (more extensive) appropriate changes.
!
! The module also exports the 'elsi_finalize_scfloop' routine, to be called
! from the appropriate place. Some variables are kept at the module level
! (communicators, distribution objects...) for this.

! Usage: Compile Siesta with -DSIESTA__ELSI
!        Define
!             SolutionMethod ELSI
!        in the fdf file
! Optionally, use
!            MPI.Nprocs.SIESTA NpSiesta
! to request fewer nodes for Siesta non-solver operations
! (which might have worse scalability than the solver)
!
module m_elsi_interface

#ifdef SIESTA__ELSI

 use precision, only  : dp
 use class_Distribution
 use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
 use elsi
 
 implicit none

 integer, parameter   :: ELPA = 1           ! solver
 integer, parameter   :: MULTI_PROC = 1     ! parallel_mode
 integer, parameter   :: PEXSI_CSC = 1      ! distribution
 integer, parameter   :: FERMI = 1          ! broadening (no Methfessel-Paxton yet)
 
 type(elsi_handle)    :: elsi_h
 
 type(Distribution)   :: dist1
 type(Distribution), allocatable, target   :: dist2_spin(:)

 integer :: ELSI_Spatial_Comm, ELSI_Spin_Comm, World_Comm
 integer :: ELSI_Spatial_Group, World_Group
 integer :: nspin, spin_rank
 logical :: ELSI_worker
 type(aux_matrix), allocatable, target :: m1_spin(:)

 public :: elsi_solver, elsi_finalize_scfloop
 private
 
CONTAINS

! This version uses separate distributions for Siesta 
! (setup_H et al) and ELSI operations. ELSI uses the CSC sparse distribution.
!
subroutine elsi_solver(iscf, no_u, no_l, nspin_in,  &
     maxnh, numh, listhptr, listh, H, S, qtot, DM, EDM, &
     ef, Entropy, temp)

    use fdf
    use parallel, only   : SIESTA_worker, BlockSize
    use parallel, only   : SIESTA_Group, SIESTA_Comm
    use m_mpi_utils, only: globalize_sum, globalize_max
    use m_mpi_utils, only: broadcast
    use units,       only: Kelvin, eV
    use alloc,             only: re_alloc, de_alloc
#ifdef MPI
    use mpi_siesta
#endif
    use elsi

  implicit          none

  integer, intent(in)  :: iscf  ! scf step number
  integer, intent(in)  :: maxnh, no_u, no_l, nspin_in
  integer, intent(in), target  :: listh(maxnh), numh(no_l), listhptr(no_l)
  real(dp), intent(in), target :: H(maxnh,nspin_in), S(maxnh)
  real(dp), intent(in) :: qtot
  real(dp), intent(out), target:: DM(maxnh,nspin_in), EDM(maxnh,nspin_in)
  real(dp), intent(out)        :: ef  ! Fermi energy
  real(dp), intent(out)        :: Entropy ! Entropy/k, dimensionless
  real(dp), intent(in)         :: temp   ! Electronic temperature

integer        :: ih, i
integer        :: info
logical        :: write_ok
!------------
external         :: timer
integer          :: mpirank, ierr
!
real(dp)  :: temperature, numElectronExact
integer   :: norbs, scf_step
!
type(Distribution), pointer :: dist2
integer, allocatable :: elsi_ranks_in_world(:)
integer  :: nworkers_SIESTA
integer, allocatable :: siesta_ranks_in_world(:)
integer, allocatable :: ELSI_ranks_in_World_Spin(:,:)
integer :: numNodesTotal
integer :: npSpatial
!
integer  :: pbs, color, spatial_rank
type(aux_matrix) :: m2
type(aux_matrix), pointer :: m1
integer :: nrows, nnz, nnzLocal, numColLocal
integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
!
real(dp), pointer, dimension(:) :: &
        HnzvalLocal=>null(), SnzvalLocal=>null(),  &
        DMnzvalLocal => null() , EDMnzvalLocal => null()
!
integer :: ispin, elsi_spin

!
real(dp)       :: bs_energy, eBandH
real(dp)       :: energy, mu
real(dp)       :: buffer1
!  --------  for serial compilation
#ifndef MPI
  call die("This ELSI solver interface needs MPI")
#endif
!
! Our global communicator is a duplicate of the passed communicator
!
call MPI_Comm_Dup(MPI_Comm_DFT, World_Comm, ierr)
call mpi_comm_rank( World_Comm, mpirank, ierr )

! NOTE:  fdf calls will assign values to the whole processor set,
! but some other variables will have to be re-broadcast (see examples
! below)

call timer("elsi", 1)  

if (SIESTA_worker) then

   ! rename some intent(in) variables, which are only
   ! defined for the Siesta subset

   norbs = no_u
   nspin = nspin_in
   scf_step = iscf
   numElectronExact = qtot 

   ! Note that the energy units for the PEXSI interface are arbitrary, but
   ! H, the interval limits, and the temperature have to be in the
   ! same units. Siesta uses Ry units.

   temperature      = temp

   if (mpirank==0) write(6,"(a,f10.2)") &
               "Electronic temperature (K): ", temperature/Kelvin
endif
!
call broadcast(norbs,comm=World_Comm)
call broadcast(scf_step,comm=World_Comm)
call broadcast(numElectronExact,World_Comm)
call broadcast(temperature,World_Comm)
call broadcast(nspin,World_Comm)

! Initialization
if (scf_step == 1) then
   
   call MPI_Comm_Group(World_Comm,World_Group, ierr)
   call MPI_Group_Size(SIESTA_Group, nworkers_SIESTA, ierr)
   allocate(siesta_ranks_in_world(nworkers_SIESTA))
   call MPI_Group_translate_ranks( SIESTA_Group, nworkers_SIESTA, &
        (/ (i,i=0,nworkers_SIESTA-1) /), &
        World_Group, siesta_ranks_in_world, ierr )
   call newDistribution(dist1,World_Comm,siesta_ranks_in_world, &
                     TYPE_BLOCK_CYCLIC,BlockSize,"bc dist")
   deallocate(siesta_ranks_in_world)
   call MPI_Barrier(World_Comm,ierr)

   call mpi_comm_size( World_Comm, numNodesTotal, ierr )
   npSpatial = numNodesTotal / nspin

! "Row" communicator for independent PEXSI operations on each spin
! The name refers to "spatial" degrees of freedom.
   color = mod(mpirank,nspin)    ! {0,1} for nspin = 2, or {0} for nspin = 1
   call MPI_Comm_Split(World_Comm, color, mpirank, ELSI_Spatial_Comm, ierr)

! "Column" communicator for spin reductions
   color = mpirank/nspin       
   call MPI_Comm_Split(World_Comm, color, mpirank, ELSI_Spin_Comm, ierr)

   call mpi_comm_rank( ELSI_Spatial_Comm, spatial_rank, ierr )
   call mpi_comm_rank( ELSI_Spin_Comm, spin_rank, ierr )
   ELSI_worker = (spatial_rank < npSpatial) ! Could be spin up or spin down

!  CSC (PEXSI) distribution blocksize 
   pbs = norbs/npSpatial

   allocate(elsi_ranks_in_world(npSpatial))
   call MPI_Comm_Group(World_Comm, World_Group, Ierr)
   call MPI_Comm_Group(ELSI_Spatial_Comm, ELSI_Spatial_Group, Ierr)
   call MPI_Group_translate_ranks( ELSI_Spatial_Group, npSpatial, &
        (/ (i,i=0,npSpatial-1) /), &
        World_Group, elsi_ranks_in_world, ierr )

! Include the actual world ranks in the distribution object

   allocate (ELSI_ranks_in_World_Spin(npSpatial,nspin))
   call MPI_AllGather(elsi_ranks_in_world,npSpatial,MPI_integer,&
        Elsi_Ranks_in_World_Spin(1,1),npSpatial, &
        MPI_integer,ELSI_Spin_Comm,ierr)

! Create distributions known to all nodes
   allocate(dist2_spin(nspin))
   do ispin = 1, nspin
      call newDistribution(dist2_spin(ispin), World_Comm, &
           Elsi_Ranks_in_World_Spin(:,ispin),  &
           TYPE_PEXSI, pbs, "px dist")
   enddo
   deallocate(elsi_ranks_in_world,Elsi_Ranks_in_World_Spin)
   call MPI_Barrier(World_Comm,ierr)

   allocate(m1_spin(nspin))
   
endif  ! First scf step -- communicator setup

elsi_spin = spin_rank+1  ! {1,2}
! This is done serially on the Siesta side, each time
! filling in the structures in one ELSI set

do ispin = 1, nspin

   m1 => m1_spin(ispin)

   if (SIESTA_worker) then
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

   ! Note that we cannot simply wrap this in a elsi_spin test, as
   ! there are Siesta nodes in both spin sets.
   ! We must discriminate the PEXSI workers by the distribution info
   dist2 => dist2_spin(ispin)
   call redistribute_spmatrix(norbs,m1,dist1,m2,dist2,World_Comm)
   
   call timer("redist_orbs_fwd", 2)

   if (ELSI_worker .and. (elsi_spin == ispin) ) then

      nrows = m2%norbs          ! or simply 'norbs'
      numColLocal = m2%no_l
      nnzLocal    = m2%nnzl
      call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,&
           ELSI_Spatial_Comm,ierr)

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

   endif ! ELSI worker
enddo


! Make these available to all
! (Note that the values are those on process 0, which is in the spin=1 set
! In fact, they are only needed for calls to the interface, so the broadcast
! could be over ELSI_Spatial_Comm only.

call MPI_Bcast(nrows,1,MPI_integer,0,World_Comm,ierr)
call MPI_Bcast(nnz,1,MPI_integer,0,World_Comm,ierr)


if (scf_step == 1) then
   ! Now we have all the ingredients to initialize the ELSI interface

   ! Use ELPA for now. PEXSI would be easy now
   call elsi_init(elsi_h, ELPA, MULTI_PROC, PEXSI_CSC, &
                  norbs, NumElectronExact,norbs/2 + 5)
   call elsi_set_csc(elsi_h, nnz, nnzLocal, numColLocal, &
        rowindLocal, colptrLocal)
   
   call elsi_set_spin(elsi_h, nspin, elsi_spin)
   call elsi_set_mpi(elsi_h, ELSI_Spatial_Comm)
   call elsi_set_mpi_global(elsi_h, World_Comm)

   call elsi_set_output(elsi_h, 3)
   call elsi_set_write_unit(elsi_h, 6)
   !
   ! Is there more documentation about the meaning of the
   ! width? I am just using k_B*T (temp is in units of energy,
   ! but Siesta uses rydberg...)
   !
   call elsi_set_mu_broaden_scheme(elsi_h, FERMI)
   call elsi_set_mu_broaden_width(elsi_h, temperature)
endif

! Solve for the DM, and get (at every step for now) the EDM and mu

call timer("elsi-solver", 1)

call elsi_dm_real_sparse(elsi_h, HnzvalLocal, SnzvalLocal, &
                                 DMnzvalLocal, energy)
call elsi_get_edm_real_sparse(elsi_h, EDMnzvalLocal)
call elsi_get_mu(elsi_h, mu)

call timer("elsi-solver", 2)

!!$   if (nspin == 2) then
!!$      ! The matrices have to be divided by two...   ???
!!$      DMnzvalLocal(:) = 0.5_dp * DMnzvalLocal(:)
!!$      EDMnzvalLocal(:) = 0.5_dp * EDMnzvalLocal(:)
!!$   endif

! Double check the energy values

if (ELSI_worker) then

   bs_energy = 0.0_dp
   eBandH = 0.0_dp
   do i = 1,nnzLocal
      bs_energy = bs_energy + SnzvalLocal(i) * &
           ( EDMnzvalLocal(i) )
      eBandH = eBandH + HnzvalLocal(i) * &
           ( DMnzvalLocal(i) )
   enddo

   ! First, reduce over the spatial communicator

   call globalize_sum( bs_energy, buffer1, comm=ELSI_Spatial_Comm )
   bs_energy = buffer1
   call globalize_sum( eBandH, buffer1, comm=ELSI_Spatial_Comm )
   eBandH = buffer1

   ! Now, reduce over both spins

   call globalize_sum( bs_energy, buffer1, comm=ELSI_Spin_Comm )
   bs_energy = buffer1
   call globalize_sum( eBandH, buffer1, comm=ELSI_Spin_Comm )
   eBandH = buffer1

   ! This output block will be executed only if World's root node is
   ! in one of the leading pole groups. This might not be so

   if ((mpirank == 0)) then
      write(6, "(a,f12.4)") "#&s Tr(S*EDM) (eV) = ", bs_energy/eV
      write(6,"(a,f12.4)") "#&s Tr(H*DM) (eV) = ", eBandH/eV
      write(6,"(a,f12.4)") "#&s ELSI energy (eV) = ", (energy)/eV
   endif

   ef = mu
   ! Note that we use the S*EDM version of the band-structure energy
   ! to estimate the entropy, by comparing it to S*FDM This looks
   ! consistent, but note that the EDM is not used in Siesta to
   ! estimate the total energy, only the DM (via the density) (that
   ! is, the XC and Hartree correction terms to Ebs going into Etot
   ! are estimated using the DM)

   ! How do we get the Entropy from ELSI?
   ! It should be there for the eigenvector solvers
   ! With the new PEXSI, it seems that FDM is gone, since its pole
   ! structure is different with the new method... What do we do?
   
   Entropy = 0.0_dp   ! ???? No info about Free energy??

   ! ef and Entropy are now known to the leading-pole processes
endif ! ELSI_worker


do ispin = 1, nspin

   m1 => m1_spin(ispin)

   if (ELSI_worker .and. (elsi_spin == ispin)) then
      ! Prepare m2 to transfer

      call de_alloc(colPtrLocal,"colPtrLocal","elsi_solver")

      call de_alloc(m2%vals(1)%data,"m2%vals(1)%data","elsi_solver")
      call de_alloc(m2%vals(2)%data,"m2%vals(2)%data","elsi_solver")

      m2%vals(1)%data => DMnzvalLocal(1:nnzLocal)
      m2%vals(2)%data => EDMnzvalLocal(1:nnzLocal)

   endif

   ! Prepare m1 to receive the results
   if (SIESTA_worker) then
      nullify(m1%vals(1)%data)    ! formerly pointing to S
      nullify(m1%vals(2)%data)    ! formerly pointing to H
      deallocate(m1%vals)
      nullify(m1%numcols)         ! formerly pointing to numH
      nullify(m1%cols)            ! formerly pointing to listH
   endif

   call timer("redist_orbs_bck", 1)
   dist2 => dist2_spin(ispin)
   call redistribute_spmatrix(norbs,m2,dist2,m1,dist1,World_Comm)
   call timer("redist_orbs_bck", 2)

   if (ELSI_worker .and. (elsi_spin == ispin)) then
      call de_alloc(DMnzvalLocal, "DMnzvalLocal", "elsi_solver")
      call de_alloc(EDMnzvalLocal,"EDMnzvalLocal","elsi_solver")

      nullify(m2%vals(1)%data)    ! formerly pointing to DM
      nullify(m2%vals(2)%data)    ! formerly pointing to EDM
      deallocate(m2%vals)
      ! allocated in the direct transfer
      call de_alloc(m2%numcols,"m2%numcols","elsi_solver")
      call de_alloc(m2%cols,   "m2%cols",   "elsi_solver")
   endif

   ! We assume that Siesta's root node also belongs to one of the
   ! ELSI communicators.
   
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

      call de_alloc(m1%vals(1)%data,"m1%vals(1)%data","elsi_solver")
      call de_alloc(m1%vals(2)%data,"m1%vals(2)%data","elsi_solver")
      deallocate(m1%vals)
      ! allocated in the direct transfer
      call de_alloc(m1%numcols,"m1%numcols","elsi_solver") 
      call de_alloc(m1%cols,   "m1%cols",   "elsi_solver")

   endif
enddo
call timer("elsi", 2)
end subroutine elsi_solver

subroutine elsi_finalize_scfloop()
    ! Free communicators, etc

    integer :: ispin, ierr
    
    call elsi_finalize(elsi_h)

    call delete(dist1)
    do ispin = 1, nspin  
       call delete(dist2_spin(ispin))
    enddo
    deallocate(dist2_spin)
    deallocate(m1_spin)

    call MPI_Comm_Free(ELSI_Spatial_Comm, ierr)
    call MPI_Comm_Free(ELSI_Spin_Comm, ierr)
    call MPI_Comm_Free(World_Comm, ierr)

    call MPI_Group_Free(ELSI_Spatial_Group, ierr)
    call MPI_Group_Free(World_Group, ierr)
  end subroutine elsi_finalize_scfloop
  

#endif  /* SIESTA__ELSI */

end module m_elsi_interface
