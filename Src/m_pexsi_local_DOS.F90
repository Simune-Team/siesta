      MODULE m_pexsi_local_DOS
      private
      public :: pexsi_local_DOS

      CONTAINS

      subroutine pexsi_local_DOS( )
      use m_energies

      use sparse_matrices
      USE siesta_options
      use siesta_geom
      use atomlist,       only: indxuo, indxua           
      use atomlist,       only: qtot, no_u, no_l
      use atomlist,       only: iphorb                   
      use atomlist,       only: datm, no_s, iaorb        
      use fdf
      use sys,            only: die                      
      use files,          only : slabel     
      use parallel,       only:  worker
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
      ! In this version, this set coincides with the set of siesta workers.
      if (worker) then

         energy = fdf_get('PEXSI.LocalDOS.Energy',0.0_dp,"eV")
         broadening = fdf_get('PEXSI.LocalDOS.Broadening',0.5_dp,"eV")

         ! Note that we re-use Dscf, so it will be obliterated
         call get_LDOS_SI( no_u, no_l, nspin,  &
              maxnh, numh, listh, H, S,  &
              Dscf, energy, broadening)
      endif

      if (worker) then
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

    END subroutine pexsi_local_DOS

    subroutine get_LDOS_SI( no_u, no_l, nspin,  &
         maxnh, numh, listh, H, S,  &
         LocalDOSDM, energy, broadening)

      use precision, only  : dp
      use fdf
      use units,       only: eV
      use m_mpi_utils, only: globalize_max
#ifdef MPI
    use mpi_siesta
#endif
    implicit          none

    integer, intent(in)          :: maxnh, no_u, no_l, nspin
    integer, intent(in)          :: listh(maxnh), numh(*)
    real(dp), intent(in)         :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in)         :: energy, broadening
    real(dp), intent(out)        :: LocalDOSDM(maxnh,nspin)

      interface
         ! subroutine f_ppexsi_localdos_interface
         include "pexsi_localdos.h"
      end interface

    integer :: siesta_comm, single_pole_comm
    integer :: ispin, maxnhtot, ih, nnzold, i

    integer  :: ordering
    integer  :: info, infomax

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


!!if (worker) then
   siesta_comm = mpi_comm_world
   Single_Pole_Comm = Siesta_Comm

!  Find rank
   call mpi_comm_rank( Single_Pole_Comm, mpirank, ierr )

   call timer("pexsi-ldos", 1)

   ispin = 1
   if (nspin /=1) then
      call die("Spin polarization not yet supported in PEXSI")
   endif

   ! Note that the energy units for the PEXSI interface are arbitrary, but
   ! H and the reference energy for the LDOS and the broadening
   ! have to be in the same units. Siesta uses Ry units.

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
   call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,Siesta_comm,ierr)

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

!!!!!endif ! worker

isSIdentity = 0

! Ordering flag
ordering = fdf_get("PEXSI.ordering",1)
!
! Broadcast these to the whole processor set, just in case
! (They were set only by the Siesta workers)
!
!call MPI_Bcast(nrows,1,MPI_integer,0,Single_Pole_Comm,ierr)
!call MPI_Bcast(nnz,1,MPI_integer,0,Single_Pole_Comm,ierr)

   if(mpirank == 0) then
     write (*,*) 'Calling PEXSI LDOS routine...'
     write(*, *) "energy (eV)       = ", energy/eV
     write(*, *) "broadening (eV)   = ", broadening/eV
   endif 

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
		Single_Pole_Comm,&
		DMnzvalLocal,&   ! LDOS/ partial DM
		info)

  ! Make sure that all processors report info=1 when any of them
  ! raises it...
  ! This should be done inside the routine
   call globalize_max(info,infomax,comm=single_pole_comm)
   info = infomax

   if (info /= 0) then
      write(msg,"(i6)") info 
      call die("Error in pexsi LDOS routine. Info: " // msg)
   endif

!!if (worker) then

   ! this will disappear
   LocalDOSDM(1:nnzLocal,ispin) = DMnzvalLocal(1:nnzLocal)

  deallocate(rowindLocal)
  deallocate(HnzvalLocal)
  deallocate(SnzvalLocal)
  deallocate(DMnzvalLocal)
  deallocate(colPtrLocal)

  call timer("pexsi-ldos", 2)
!! endif ! worker

#endif ! mpi

end subroutine get_LDOS_SI

END module m_pexsi_local_DOS
