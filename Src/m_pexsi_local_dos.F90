      MODULE m_pexsi_local_dos
      private
      public :: pexsi_local_dos

      CONTAINS

      subroutine pexsi_local_dos( )
      use m_energies

      use sparse_matrices
      USE siesta_options
      use siesta_geom
      use atomlist,       only: indxuo, indxua           
      use atomlist,       only: qtot, no_u, no_l
      use atomlist,       only: iphorb                   
      use atomlist,       only: datm, no_s, iaorb        
      use fdf
      use files,          only : slabel     
      use parallel,       only:  SIESTA_worker
      use files,          only : label_length            
      use m_ntm
      use m_forces,       only: fa
      use m_eo
      use m_spin,         only: nspin, qs, efs
      use m_gamma
      use m_dhscf,        only: dhscf
      implicit none

      character(len=label_length+5), external :: paste

      integer :: dummy_iscf = 1
      character(len=label_length+5) :: fildos

      real(dp)  :: dummy_str(3,3), dummy_strl(3,3)  ! for dhscf call
      real(dp)  :: dummy_dipol(3)

      real(dp)  :: factor, g2max
      real(dp)  :: energy, broadening

      ! Find local density of states with Selected Inversion
      ! Only the first pole group participates 
      ! in the computation of the selected
      ! inversion for a single shift.
      ! In this version, this set coincides with the set of siesta SIESTA_workers.

         energy = fdf_get('PEXSI.LocalDOS.Energy',0.0_dp,"Ry")
         broadening = fdf_get('PEXSI.LocalDOS.Broadening',0.01_dp,"Ry")

         ! Note that we re-use Dscf, so it will be obliterated
         call get_LDOS_SI( no_u, no_l, nspin,  &
              maxnh, numh, listh, H, S,  &
              Dscf, energy, broadening)

      if (SIESTA_worker) then
         !Find the LDOS in the real space mesh
         fildos = paste( slabel, '.LDSI' )
         g2max = g2cut
         call dhscf( nspin, no_s, iaorb, iphorb, no_l, &
                   no_u, na_u, na_s, isa, xa_last, indxua,  &
                   ntm, 0, 0, 0, fildos, ' ', ' ', ' ', ' ', ' ', &
                   maxnh, numh, listhptr, listh, Dscf, Datm, maxnh, H, &
                   Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc, &
                   dummy_dipol, dummy_str, fa, dummy_strl )
                    ! next to last argument is dummy here,
                    ! as no forces are calculated
                    ! todo: make all these optional

      endif

    END subroutine pexsi_local_dos

    subroutine get_LDOS_SI( no_u, no_l, nspin,  &
         maxnh, numh, listh, H, S,  &
         LocalDOSDM, energy, broadening)

      use precision, only  : dp
      use fdf
      use units,       only: eV
      use sys,         only: die
      use m_mpi_utils, only: globalize_max, broadcast
      use parallel, only   : SIESTA_worker, BlockSize
      use parallel, only   : SIESTA_Group, SIESTA_Comm
      use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
      use class_Dist
      use alloc,             only: re_alloc, de_alloc
#ifdef MPI
    use mpi_siesta
#endif
    use m_pexsi_interface, only: f_ppexsi_localdos_interface

    implicit          none

    integer, intent(in)          :: maxnh, no_u, no_l, nspin
    integer, intent(in), target  :: listh(maxnh), numh(no_l)
    real(dp), intent(in), target :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in)         :: energy, broadening
    real(dp), intent(out)        :: LocalDOSDM(maxnh,nspin)


    integer :: PEXSI_Comm, World_Comm
    integer :: PEXSI_Group, World_Group

    integer :: ispin, maxnhtot, ih, nnzold, i, pbs, norbs, npPerPole
    integer :: npSymbFact
    logical :: PEXSI_worker

    integer  :: ordering
    integer  :: info, infomax

    type(aux_matrix) :: m1, m2
    type(Dist)       :: dist1, dist2

!Lin variables
integer :: nrows, nnz, nnzLocal, numColLocal
integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
!
real(dp), pointer, dimension(:) :: &
        HnzvalLocal=>null(), SnzvalLocal=>null(),  &
        DMnzvalLocal => null()
!
integer  :: mpirank, ierr
integer  :: isSIdentity
!------------

external         :: timer
character(len=6) :: msg



#ifndef MPI
    call die("PEXSI-LDOS needs MPI")
#else

! "SIESTA_Worker" means a processor which is in the Siesta subset.
! NOTE:  fdf calls will assign values to the whole processor set,
! but some other variables will have to be re-broadcast (see examples
! below)

World_Comm = true_MPI_Comm_World

if (SIESTA_worker) then
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

! Number of processors for symbolic factorization
! Only relevant for PARMETIS/PT_SCOTCH
npSymbFact = fdf_get("PEXSI.np-symbfact",npPerPole)

pbs = norbs/npPerPole
call newDistribution(pbs,PEXSI_Group,dist2,TYPE_PEXSI,"px dist")

if (SIESTA_worker) then
   call timer("pexsi-ldos", 1)

   ispin = 1
   if (nspin /=1) then
      call die("Spin polarization not yet supported in PEXSI")
   endif

   ! Note that the energy units for the PEXSI interface are arbitrary, but
   ! H and the reference energy for the LDOS and the broadening
   ! have to be in the same units. Siesta uses Ry units.

   call MPI_Barrier(Siesta_comm,ierr)

   m1%norbs = norbs
   m1%no_l  = no_l
   m1%nnzl  = sum(numH(1:no_l))
   m1%numcols => numH
   m1%cols    => listH
   allocate(m1%vals(2))
   m1%vals(1)%data => S(:)
   m1%vals(2)%data => H(:,ispin)

endif  ! SIESTA_worker

call redistribute_spmatrix(norbs,m1,dist1,m2,dist2,World_Comm)

if (PEXSI_worker) then

   nrows = m2%norbs          ! or simply 'norbs'
   numColLocal = m2%no_l
   nnzLocal    = m2%nnzl
   call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,PEXSI_Comm,ierr)

  call re_alloc(colptrLocal,1,numColLocal+1,"colptrLocal","pexsi_ldos")
  colptrLocal(1) = 1
  do ih = 1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + m2%numcols(ih)
  enddo

  rowindLocal => m2%cols
  SnzvalLocal => m2%vals(1)%data
  HnzvalLocal => m2%vals(2)%data

  call re_alloc(DMnzvalLocal,1,nnzLocal,"DMnzvalLocal","pexsi_ldos")

  isSIdentity = 0

  ! Ordering flag
  ordering = fdf_get("PEXSI.ordering",1)
!
! No need to broadcast these to the whole processor set.
! (They were set by the PEXSI workers, which are the only
!  ones making the interface call)
!
!call MPI_Bcast(nrows,1,MPI_integer,0,World_Comm,ierr)
!call MPI_Bcast(nnz,1,MPI_integer,0,World_Comm,ierr)

   if(mpirank == 0) then
     write (*,*) 'Calling PEXSI LDOS routine...'
     write(*, *) "energy (eV)       = ", energy/eV
     write(*, *) "broadening (eV)   = ", broadening/eV
   endif 

! (Note that only the first-pole team does this)

	call f_ppexsi_localdos_interface(&
		nrows,&
		nnz,&
		nnzLocal,&
		numColLocal,&
		colptrLocal,&
		rowindLocal,&
		HnzvalLocal,&
		isSIdentity,&
		SnzvalLocal,&
		Energy,&
		Broadening,&
		ordering,&
                npSymbFact,&
		PEXSI_Comm,&
		DMnzvalLocal,&   ! LDOS/ partial DM
		info)

  ! Make sure that all processors report info=1 when any of them
  ! raises it...
  ! This should be done inside the routine
   call globalize_max(info,infomax,comm=PEXSI_comm)
   info = infomax

   if (info /= 0) then
      write(msg,"(i6)") info 
      call die("Error in pexsi LDOS routine. Info: " // msg)
   endif

   call de_alloc(colPtrLocal,"colPtrLocal","pexsi_ldos")

   deallocate(m2%vals(1)%data)
   deallocate(m2%vals(2)%data)
   deallocate(m2%vals)
   allocate(m2%vals(1))
   m2%vals(1)%data => DMnzvalLocal(1:nnzLocal)

endif ! PEXSI_worker

! Now we should not re-allocate m1...
! or better:
if (SIESTA_worker) then
   nullify(m1%vals(1)%data)    ! formerly pointing to S
   nullify(m1%vals(2)%data)    ! formerly pointing to H
   deallocate(m1%vals)
   nullify(m1%numcols)         ! formerly pointing to numH
   nullify(m1%cols)            ! formerly pointing to listH
endif

call redistribute_spmatrix(norbs,m2,dist2,m1,dist1,World_Comm)

if (PEXSI_worker) then
   call de_alloc(DMnzvalLocal,"DMnzvalLocal","pexsi_ldos")

   nullify(m2%vals(1)%data)    ! formerly pointing to DM
   deallocate(m2%vals)
   deallocate(m2%numcols)      ! allocated in the direct transfer
   deallocate(m2%cols)         !  "
endif


if (SIESTA_worker) then

   LocalDOSDM(:,ispin)  = m1%vals(1)%data(:)  
   ! Check no_l
   if (no_l /= m1%no_l) then
      call die("Mismatch in no_l")
   endif
   ! Check listH
   if (any(listH(:) /= m1%cols(:))) then
      call die("Mismatch in listH")
   endif

   deallocate(m1%vals(1)%data)
   deallocate(m1%vals)
   deallocate(m1%cols)
   deallocate(m1%numcols)

   call timer("pexsi-ldos", 2)

endif

call delete(dist1)
call delete(dist2)
if (PEXSI_worker) then
   call MPI_Comm_Free(PEXSI_Comm, ierr)
   call MPI_Group_Free(PEXSI_Group, ierr)
endif

#endif ! mpi

end subroutine get_LDOS_SI

END module m_pexsi_local_dos
