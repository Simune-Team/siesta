! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! ##################################################################
! ##              PDOS curves resolved into atoms                 ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##  Adapted and streamlined by Nick Papior Andersen             ##
! ##################################################################

subroutine AtomPDOS(uTOTDOS,uORBDOS,LowdinFlag,spin_F, &
     noBufL_P_L,noD, &
     na_u,IsoAt1,IsoAt2,lasto, &
     Energy,wGF, &
     GF,GFRGF,S)

  use precision,       only : dp
  use sys,             only : die
  use parallel,        only : Node, Nodes, IOnode
  use m_lowdin,        only : Lowdin_Trans
  use m_tbt_out,       only : out_AtomPDOS_Tot, out_AtomPDOS_Orb
  use units,           only : Pi, eV
#ifdef MPI
  use mpi_siesta, only : MPI_Comm_World
  use mpi_siesta, only : MPI_ISend,MPI_IRecv
  use mpi_siesta, only : MPI_Wait,MPI_Status_Size
  use mpi_siesta, only : MPI_Double_Precision
#endif

  implicit none
! ************************
! * INPUT variables      *
! ************************
  integer, intent(in)     :: uTOTDOS,uORBDOS  ! units for writing out info
  logical, intent(in)     :: LowdinFlag       ! if .true. then Loewdin Charges
  real(dp),intent(in)     :: spin_F           ! The spin factor of the calculation
  integer, intent(in)     :: noBufL_P_L       ! number of orbitals in the left  device region
  integer, intent(in)     :: noD              ! number of orbitals in the device region
  integer, intent(in)     :: na_u             ! Atoms in the unit cell
  integer, intent(in)     :: IsoAt1,IsoAt2    ! the projection region for the atoms
  integer, intent(in)     :: lasto(0:na_u)    ! The lasto orbital for the atoms
  real(dp), intent(in)    :: Energy           ! energy for the analysis
  real(dp), intent(in)    :: wGF              ! the weight for the energy point
! 'rho = S^1.2 x G x S^1.2' else
! Mulliken charges 'rho = G x S'
  complex(dp), intent(in) :: GF(noD,noD)        !
  complex(dp), intent(in) :: GFRGF(noD,noD)     !
! The variable can not be changed
  complex(dp), intent(inout) :: S(noD,noD)         ! overlap matrix in case of Loewdin S=S^1.2
 
! ************************
! * LOCAL variables      *
! ************************
  real(dp), parameter :: r1dPi = 1.0_dp/Pi ! Local Pi factor
  integer             :: i, no
  integer             :: ia, io ! loops
  real(dp), dimension(IsoAt1:IsoAt2) :: PDOSTot, PDOSLeft, PDOSRight ! The dos regions

  real(dp),    dimension(:,:),  allocatable :: PDOST
  complex(dp), dimension (:,:), allocatable :: GFS
  complex(dp), dimension (:,:), allocatable :: GFRS

#ifdef MPI
  real(dp), dimension(:,:),   allocatable :: bufE
  real(dp), dimension(:,:,:), allocatable :: buf
  integer :: req, status(MPI_Status_Size)
  integer :: iNode, MPIerror
#endif

  call timer('AtomPDOS',1)

  ! Retrieve the maximum number of orbitals on one atom
  ! This is found from the atoms that are to be calculated
  no = maxval(lasto(IsoAt1:IsoAt2)-lasto(IsoAt1-1:IsoAt2-1))
! Allocate used arrays
  allocate(PDOST(no,IsoAt1:IsoAt2))
  call memory('A','D',no*(IsoAt2-IsoAt1+1),'AtomPDOS')
  allocate(GFS(noD,noD))
  call memory('A','Z',noD*noD,'AtomPDOS')
  allocate(GFRS(noD,noD))
  call memory('A','Z',noD*noD,'AtomPDOS')

  GFS  = dcmplx(0.0_dp,0.0_dp)
  GFRS = dcmplx(0.0_dp,0.0_dp)

  if ( LowdinFlag ) then
     GFS  = GF
     call Lowdin_Trans(.false.,noD,GFS,S)
     GFRS = GFRGF
     call Lowdin_Trans(.false.,noD,GFRS,S)
  else
     call zgemm('N','N',noD,noD,noD, &
          dcmplx(1.0_dp,0.0_dp),GF,noD,S,noD, &
          dcmplx(0.0_dp,0.0_dp),GFS,noD) 

     call zgemm('N','N',noD,noD,noD, &
          dcmplx(1.0_dp,0.0_dp),GFRGF,noD,S,noD, &
          dcmplx(0.0_dp,0.0_dp),GFRS,noD) 
  end if

! Start calculating the AtomPDOS
  do ia = IsoAt1, IsoAt2
! Initialize for the ia'th atom
     PDosLeft(ia)  = 0.0_dp
     PDosRight(ia) = 0.0_dp
     PDosTot(ia)   = 0.0_dp
     PDosT(:,ia)   = 0.0_dp
     
     no = lasto(ia-1) - noBufL_P_L

     do i = 1 , lasto(ia) - lasto(ia-1)
        io = i + no

        PDOST  (i,ia) =               - spin_F * r1dPi * dimag( GFS(io,io))
        PDOSTot  (ia) = PDOSTot  (ia) - spin_F * r1dPi * dimag( GFS(io,io))
        PDOSRight(ia) = PDOSRight(ia) - spin_F * r1dPi * dimag(GFRS(io,io))
     end do

     PDOSLeft(ia) = PDOSTot(ia) - PDOSRight(ia)
  end do

! Clean up to prepare for MPI routines
  call memory('D','Z',noD*noD,'AtomPDOS')
  deallocate(GFS)
  call memory('D','Z',noD*noD,'AtomPDOS')
  deallocate(GFRS)
  
#ifdef MPI
  no = size(PDOST,1)

  if ( Node == 0 ) then
     allocate(buf(3+no,IsoAt2-IsoAt1+1,0:Nodes-1))
     allocate(bufE(2,0:Nodes-1))
  else
     allocate(buf(3+no,IsoAt2-IsoAt1+1,1))
     allocate(bufE(2,1))
  end if

  if ( Node == 0 ) then
     bufE(1,Node)         = Energy
     bufE(2,Node)         = wGF
     buf(1,:,Node)        = PDOSTot(:)
     buf(2,:,Node)        = PDOSLeft(:)
     buf(3,:,Node)        = PDOSRight(:)
     buf(4:4+no-1,:,Node) = PDOST(:,:)
  else
     bufE(1,1)            = Energy
     bufE(2,1)            = wGF
     buf(1,:,1)           = PDOSTot(:)
     buf(2,:,1)           = PDOSLeft(:)
     buf(3,:,1)           = PDOSRight(:)
     buf(4:4+no-1,:,1)    = PDOST(:,:)
  end if

  if ( Node == 0 .and. Nodes > 1 ) then
     do iNode = 1 , Nodes - 1
        call MPI_IRecv(buf(1,1,iNode),(no+3)*(IsoAt2-IsoAt1+1), &
             MPI_Double_Precision, iNode, iNode, MPI_Comm_World, req, MPIerror)
        call MPI_Wait(req,status,MPIerror)
        call MPI_IRecv(bufE(1,iNode),2, &
             MPI_Double_Precision, iNode, iNode, MPI_Comm_World, req, MPIerror)
        call MPI_Wait(req,status,MPIerror)
     end do
  else if ( Node /= 0 ) then
     call MPI_ISend(buf(1,1,1),(no+3)*(IsoAt2-IsoAt1+1), &
          MPI_Double_Precision, 0, Node, MPI_Comm_World, req, MPIerror)
     call MPI_Wait(req,status,MPIerror)
     call MPI_ISend(bufE(1,1),2, &
          MPI_Double_Precision, 0, Node, MPI_Comm_World, req, MPIerror)
     call MPI_Wait(req,status,MPIerror)
  end if

  do iNode = 0 , Nodes - 1
     do ia = 1 , IsoAt2 - IsoAt1 + 1
        no = lasto(IsoAt1+ia-1) - lasto(IsoAt1+ia-2)
        call out_AtomPDOS_Tot(uTOTDOS,ia+IsoAt1-1,bufE(1,iNode),bufE(2,iNode), &
             buf(1,ia,iNode),buf(2,ia,iNode),buf(3,ia,iNode))
        call out_AtomPDOS_Orb(uORBDOS,ia+IsoAt1-1,bufE(1,iNode),bufE(2,iNode), &
             buf(4,ia,iNode),no)
     end do
  end do

  deallocate(buf,bufE)

#else
  do ia = IsoAt1 , IsoAt2
     no = lasto(ia) - lasto(ia-1)
     call out_AtomPDOS_Tot(uTOTDOS,ia,Energy,wGF, &
          PDOSTot(ia),PDOSLeft(ia),PDOSRight(ia))
     call out_AtomPDOS_Orb(uORBDOS,ia,Energy,wGF, &
          PDOST(1,ia),no)
  end do
#endif

  call memory('D','D',no*(IsoAt2-IsoAt1+1),'AtomPDOS')
  deallocate(PDOST)

  call timer('AtomPDOS',2)

end subroutine AtomPDOS


