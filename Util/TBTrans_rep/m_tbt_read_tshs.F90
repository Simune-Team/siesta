! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_tbt_read_tshs

contains
  subroutine tbt_read_tshs(HSfile,no_s,no_u,nspin, &
       ucell, na_u, xa, lasto, &
       maxnh , numh , listhptr , listh , xij , indxuo, &
       H, S, &
       Gamma, Ef)

    use precision,    only : dp
    use parallel,     only : Node, Nodes, IONode
    use m_ts_io,      only : ts_iohs

#ifdef MPI
    use mpi_siesta
#endif

    implicit none

! **************************
! * INPUT variables        *
! **************************
    character(len=*) :: HSfile    ! file name of TSHS

! **************************
! * OUTPUT variables       *
! **************************
  integer, intent(out)  :: no_s            ! # orbitals in supercell cell local to this processor
  integer, intent(out)  :: no_u            ! # orbitals in unit cell global
  integer, intent(out)  :: nspin           ! spins in the system
  real(dp), intent(out) :: ucell(3,3)      ! The unit cell that the transformation is happening in
  integer, intent(out)  :: na_u            ! The number of atoms in the unit cell
  real(dp), pointer, intent(out) :: xa(:,:)      ! The atomic coordinates in the unit cell
  integer, pointer, intent(out)  :: lasto(:)     ! The last orbital of each atom referenced in the unit cell
  integer, intent(out)  :: maxnh           ! Maximum number of orbitals interacting
  integer, pointer, intent(out)  :: numh(:)      ! # of nonzero elements of each row of hamiltonian matrix
  integer, pointer, intent(out)  :: listhptr(:)  ! Pointer to each row (-1) of the hamiltonian matrix
  integer, pointer, intent(out)  :: listh(:)     ! Nonzero hamiltonian-matrix element column
!                                                    indexes for each matrix row
  real(dp), pointer, intent(out) :: xij(:,:)     ! Vectors between orbital centers (sparse)
!                                                    (not only gamma point)
  integer, pointer, intent(out)  :: indxuo(:)    ! Index of equivalent orbital in unit cell
!                                                    Unit cell orbitals must be the first in
!                                                    orbital lists, i.e. indxuo.le.no_l, with
!                                                    no_l the number of orbitals in unit cell
  real(dp), pointer, intent(out) :: H(:,:)        ! Hamiltonian in sparse form
  real(dp), pointer, intent(out) :: S(:)        ! Overlap in sparse form
  logical, intent(in)  :: Gamma                ! Gamma point calculation (i.e. one k-point)
  real(dp), intent(out) :: Ef                   ! Fermi energy

! ***********************
! * LOCAL variables     *
! ***********************
! >>>> Related to the Electrode TSHS
    integer,  dimension(:), pointer :: isa ! The atomic species
    logical :: ts_gamma ! Read gamma from file
    real(dp) :: kdispl(3)
    integer  :: kscell(3,3)
    real(dp) :: qtot,Temp ! total charges, electronic temperature
    integer :: dummy1,dummy2 ! dummy variables
    character(8) :: iotask

#ifdef MPI
    integer :: MPIerror
#endif

    if ( IONode ) then
       iotask = 'READ'
       call ts_iohs(iotask,Gamma, .false., no_u, no_s, nspin, &
            indxuo, maxnh, numh, listhptr, listh, H, S, qtot, Temp, &
            xij, len(HSfile), HSfile, na_u, lasto, isa, Ef, &
            ucell, kscell, kdispl, ts_gamma, xa, dummy1,dummy2,check_kcell=.false.)
       ! Deallocate isa, not needed anymore
       call memory('D','I',na_u,'tbt_readHS')
       deallocate(isa)
    end if

! Do communication of variables
#ifdef MPI
    call MPI_Bcast(no_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(no_s,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(na_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nspin,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(maxnh,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(Ef,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(ucell(1,1),3*3,MPI_Double_Precision,0, &
         MPI_Comm_World,MPIerror)
    if ( .not. IONode ) then
       if ( .not. Gamma ) then 
          allocate(indxuo(no_s))
          call memory('A','I',no_s,'read_tshs')
          allocate(xij(3,maxnh))
          call memory('A','D',maxnh*3,'read_tshs')
       end if
       allocate(xa(3,na_u))
       call memory('A','D',3*na_u,'read_tshs')
       allocate(lasto(0:na_u))
       call memory('A','I',1+na_u,'read_tshs')
       allocate(numh(no_u))
       call memory('A','I',no_u,'read_tshs')
       allocate(listhptr(no_u))
       call memory('A','I',no_u,'read_tshs')
       allocate(listh(maxnh))
       call memory('A','I',maxnh,'read_tshs')
       allocate(H(maxnh,nspin))
       call memory('A','D',maxnh*nspin,'read_tshs')
       allocate(S(maxnh))
       call memory('A','D',maxnh,'read_tshs')
    end if
    call MPI_Bcast(xa(1,1),3*na_u,MPI_Double_Precision,0, &
         MPI_Comm_World,MPIerror)
    if ( .not. Gamma ) then
       call MPI_Bcast(indxuo,no_s,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(xij(1,1),3*maxnh,MPI_Double_Precision,0, &
            MPI_Comm_World,MPIerror)
    end if
    call MPI_Bcast(lasto(0),1+na_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(numh,no_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(listhptr,no_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(listh,maxnh,MPI_Integer,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(H(1,1),maxnh*nspin,MPI_Double_Precision,0, &
         MPI_Comm_World,MPIerror)
    call MPI_Bcast(S,maxnh,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
#endif

  end subroutine tbt_read_tshs
end module m_tbt_read_tshs
