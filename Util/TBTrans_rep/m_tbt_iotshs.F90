! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_tbt_iotshs
!
! Routines that are used for Input and Output of files 
! This has been adapted to be used in TBtrans forms
!
!=============================================================================
! CONTAINS:
!          1) tbt_read_tshs

  implicit none

  public :: tbt_read_tshs

  private

contains

  subroutine tbt_read_tshs(HSfile,no_s,no_u,nspin, &
       ucell, na_u, xa, lasto, &
       maxnh , numh , listhptr , listh , xij , indxuo, &
       H, S, &
       Gamma, Ef)

! *********************************************************************
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible
! Modified by M.Paulsson 2009 to:
! 1: To include information of which FC step for phonon calculations
! 2: To only save the overlap matrix if onlyS flag is set
!    (Used for e-ph coupling calculations)
! 3: File format changed to unify Copenhagen/Barcelona Transiesta vers.
! 4: Smaller files by writing arrays directly instead of element wise
! *************************** INPUT **********************************
! character(*) HSfile         : file name to read
! *************************** OUTPUT *********************************
! integer no_s                : Number of basis orbitals per supercell
! integer no_u                : Number of basis orbitals per unit cell
! integer nspin               : Spin polarization (1 or 2)
! real*8  ucell(3,3)          : unit cell of TSHS file
! integer na_u                : Number of atoms per unit cell
! real*8  xa(3,na_u)          : Atomic coordinates in the unitcell
! integer lasto(0:na_u)       : The last orbital of the atom in the unit cell
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(no_u)          : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(no_u)      : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! integer indxuo(no_s)        : Index of orbitals in supercell
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! logical Gamma               : Where it a Gamma calculation
! real*8  Ef                  : The Fermi level of the system

!
!  Modules
!
    use precision,     only : dp
    use parallel,      only : Node, Nodes, IONode
    use sys,           only : die
    use m_tbt_kpoints, only : kscell, kdispl

#ifdef MPI
    use mpi_siesta
#endif

! **************************
! * INPUT variables        *
! **************************
    character(len=*),intent(in) :: HSfile    ! file name of TSHS

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
    logical, intent(out)  :: Gamma           ! Gamma point calculation (i.e. one k-point)
    real(dp), intent(out) :: Ef              ! Fermi energy

! **************************
! * LOCAL variables        *
! **************************
    integer, allocatable :: isa(:)
    logical :: onlySfile
    real(dp) :: qtot, Temp
    integer :: iu ! File unit ID
    integer, pointer  :: iaorb(:)     ! The equivalent atomic index for a given orbital index


    ! Loop counters
    integer :: i, j, ih, is

! Variables concerning the k-point sampling in the TranSIESTA run
    integer, dimension(3,3) :: ts_kscell_file
    real(dp), dimension(3)  :: ts_kdispl_file 
    logical :: ts_gamma_scf_file

    integer :: iTmp ! Used for temporary reading
    logical :: found, same

#ifdef MPI
    integer :: MPIerror
#endif


! Check if input file exists
    inquire(file=HSfile, exist=found)
    if ( .not. found ) call die("Could not find file: "//trim(HSfile))

    if ( IONode ) then
! Open file
       call io_assign( iu )
       open( iu, file=HSfile, form='unformatted', status='old' )      

! Read dimensions
       read(iu) na_u, no_u, no_s, nspin, maxnh

! Allocate arrays that are going to be read now
       allocate(xa(3,na_u))
       allocate(isa(na_u)) 
       call memory('A','D',3*na_u,'read_tshs')
       call memory('A','I',na_u,'read_tshs')
! Read Geometry information
       read(iu) xa
       read(iu) isa   
       read(iu) ucell  
       deallocate(isa) 
       call memory('D','I',na_u,'read_tshs')

! Read k-point sampling information
       read(iu) Gamma
       read(iu) onlySfile
       read(iu) ts_gamma_scf_file
       read(iu) ts_kscell_file
       read(iu) ts_kdispl_file
! actually iStep and ia1
       read(iu) iTmp, iTmp

! Check whether file is compatible from a gamma point of view
       if (onlySfile) then
          write(*,*) 'TSHS file does not contain Hamiltonian'
          call die('read_tshs: onlyS flag, no H in file')
       endif

! Check the k-point sampling
! In this case as we are in TBTrans we do not worry that they are not the same.
! In fact it can often be increased in accuracy while calculating the transmission.
! However, we write out that we have encountered a non-matching k-point sampling...
       same = .true.
       do i = 1,3
          same = same .and. ( ts_kdispl_file(i) == kdispl(i) )
       end do
       if ( .not. same ) then
          write(*,*) "WARNING: TSHS and TBTrans does not have same k-displacements." 
          write(*,'(a,2F8.4)') 'TSHS k displacements    :', (ts_kdispl_file(j),j=1,2)
          write(*,'(a,2F8.4)') 'TBTrans k displacements :', (kdispl(j),j=1,2)
       end if
       same = .true.
       do j = 1,3
          do i = 1,3
             same = same .and. ( ts_kscell_file(i,j) == kscell(i,j) )
          end do
       end do
       if ( .not. same ) then
          write(*,*) "WARNING: TSHS and TBTrans does not have same k-cell." 
          do i = 1,3
             write(*,'(a,/,3(2x,i3))') 'TSHS k cell    :', (ts_kscell_file(i,j),j=1,3)
          end do
          do i = 1,3
             write(*,'(a,/,3(2x,i3))') 'TBTrans k cell :', (kscell(i,j),j=1,3)
          end do
       end if

! Read sparse listings
       allocate(lasto(0:na_u))
       call memory('A','I',1+na_u,'read_tshs')
       read(iu) lasto


       if (.not.Gamma) then
! Allocate arrays that are going to be read now
          allocate(indxuo(no_s))
          call memory('A','I',no_s,'read_tshs')
          read(iu) (indxuo(ih),ih=1,no_s)
       endif

! Allocate local array for global numh
       allocate(numh(no_u))
       call memory('A','I',no_u,'read_tshs')

! Read numh and send to appropriate Node
       read(iu) numh


! Read Electronic Structure Information
       read(iu) qtot,temp
       read(iu) Ef

! Create listhptr
       call memory('A','I',no_u,'read_tshs')
       allocate(listhptr(no_u))
       listhptr(1) = 0

       do ih = 2,no_u
          listhptr(ih) = listhptr(ih-1) + numh(ih-1)
       enddo

! Read listh
! Allocate lish
       allocate(listh(maxnh))
       call memory('A','I',maxnh,'read_tshs')

       do ih = 1,no_u
          read(iu) listh(listhptr(ih)+1:listhptr(ih)+numh(ih))
       enddo

! Read Overlap matrix
! Allocate S
       allocate(S(maxnh))
       call memory('A','D',maxnh,'read_tshs')
       do ih = 1,no_u
          read(iu) S(listhptr(ih)+1:listhptr(ih)+numh(ih))
       enddo

! Read Hamiltonian
! Allocate H
       allocate(H(maxnh,nspin))
       call memory('A','D',maxnh*nspin,'read_tshs')
       do is = 1,nspin
          do ih = 1,no_u
             read(iu) H(listhptr(ih)+1:listhptr(ih)+numh(ih),is)
          enddo
       enddo

       if (.not.Gamma) then
! Read interorbital vectors for K point phasing
! Allocate xij
          allocate(xij(3,maxnh))
          call memory('A','D',3*maxnh,'read_tshs')
          do ih = 1,no_u
             read(iu) (xij(i,listhptr(ih)+1:listhptr(ih)+numh(ih)),i=1,3)
          enddo
       endif

! Close file
       call io_close( iu )
    end if

! Do communication of variables
#ifdef MPI
    call MPI_Bcast(Gamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(na_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(no_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(no_s,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
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

end module m_tbt_iotshs
