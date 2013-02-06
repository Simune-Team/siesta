module m_ts_io
!
! Routines that are used for Input and Output of files 
!
!=============================================================================
! CONTAINS:
!          1) ts_read_TSHS_na
!          2) ts_iohs


  implicit none

  public :: ts_read_TSHS_na
  public :: ts_read_TSHS_lasto

  public :: ts_read_TSHS
#ifndef TBTRANS
  public :: ts_write_TSHS
#endif
  public :: filenameTSHS

  private

contains

  ! Reads in the number of atoms in the electrode. This is for easy operation
  ! In the options reading phase. We need the size of the electrodes to determine 
  ! the number of atoms in the electrode.
  subroutine ts_read_TSHS_na(TSHS,na_u)
    use parallel, only : IONode
#ifdef MPI
    use mpi_siesta
#endif
! ***********************
! * INPUT variables     *
! ***********************
    character(len=200), intent(in) :: TSHS

! ***********************
! * OUTPUT variables    *
! ***********************    
    integer, intent(out)           :: na_u

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: uTSHS,tmp(4)
#ifdef MPI
    integer :: MPIerror
#endif
    
    if ( IONode ) then
       call io_assign(uTSHS)
       open(file=TSHS,unit=uTSHS,form='unformatted')
       read(uTSHS) na_u, tmp(1:4) !na_u, no_u, no_s, Enspin, maxnh
       call io_close(uTSHS)
    end if

#ifdef MPI
    call MPI_Bcast(na_u,1,MPI_INTEGER,0,MPI_Comm_World,MPIerror)
#endif

  end subroutine ts_read_TSHS_na


  ! Reads in the orbitals in the electrode. This is for easy operation
  subroutine ts_read_TSHS_lasto(TSHS,na_u,lasto)
    use parallel, only : IONode
#ifdef MPI
    use mpi_siesta
#endif
! ***********************
! * INPUT variables     *
! ***********************
    character(len=200), intent(in) :: TSHS
    integer,            intent(in) :: na_u

! ***********************
! * OUTPUT variables    *
! ***********************    
    integer,           intent(out) :: lasto(0:na_u)

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: uTSHS
#ifdef MPI
    integer :: MPIerror
#endif
    
    if ( IONode ) then
       call io_assign(uTSHS)
       open(file=TSHS,unit=uTSHS,form='unformatted')
       read(uTSHS) ! na_u, no_u, no_s, Enspin, maxnh
       read(uTSHS) ! xa
       read(uTSHS) ! isa   
       read(uTSHS) ! ucell  
       read(uTSHS) ! gammaonfile   
       read(uTSHS) ! onlySfile
       read(uTSHS) ! ts_gamma_scf_file       
       read(uTSHS) ! ts_kscell_file
       read(uTSHS) ! ts_kdispl_file  
       read(uTSHS) ! istep, ia1
       read(uTSHS) lasto

       call io_close(uTSHS)
    end if

#ifdef MPI
    call MPI_Bcast(lasto,na_u+1,MPI_INTEGER,0,MPI_Comm_World,MPIerror)
#endif

  end subroutine ts_read_TSHS_lasto

! TBTrans will not need to write anything... (maybe in the future...)
#ifndef TBTRANS
! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! First we give the parameters that MUST be the same in the files.
! Also we supply the system label which leverages the need for the 
! files module
  subroutine ts_write_tshs(filename, &
       onlyS, Gamma, TSGamma, &
       ucell, na_u, no_l, no_u, no_s, maxnh, nspin,  &
       kscell, kdispl, &
       xa, iza, lasto, &
       numh, listhptr, listh, xij, indxuo, &
       H, S, Ef, &
       Qtot, Temp, &
       istep, ia1)

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
! logical       Gamma         : Is only gamma point used?
! logical       TSGamma       : Is only TS gamma point used?
! logical       onlyS         : Should only overlap matrix be saved?
! ******************** INPUT or OUTPUT (depending on task) ***********
! integer no_u                : Number of basis orbitals per unit cell
! integer no_s                : Number of basis orbitals per supercell
! integer Enspin              : Spin polarization (1 or 2)
! integer indxuo(no_s)        : Index of orbitals in supercell
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(maxnh,Enspin)     : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! TSS Begin
! ********************* ADDED ARGUMENTS FOR TRANSIESTA ****************
! integer fnlength            : file name length
! character(fnlength)  fname  : file nema for input or output
! integer na_u                : Number of atoms per unit cell
! integer istep, ia1          : Force constant step and atom number
! logical check_kcell        : Do a check of the kscell and kdispl
! TSS End
! *************************** UNITS ***********************************
! Units should be consistent between task='read' and 'write'
! *********************************************************************

!
!  Modules
!
    use precision,    only : dp
    use parallel,     only : IONode

#ifdef MPI
    use m_glob_sparse
#endif

    implicit none

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(in) :: onlyS
    logical, intent(in) :: Gamma, TSGamma
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, no_l, no_u, no_s, maxnh, nspin
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: iza(na_u)
    integer, intent(in) :: numh(no_l), listhptr(no_l)
    integer, intent(in) :: listh(maxnh)
    real(dp), intent(in) :: xij(3,maxnh)
    integer, intent(in) :: indxuo(no_s)
    integer, intent(in) :: lasto(0:na_u)
    real(dp), intent(in) :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: Ef
    integer, intent(in) :: kscell(3,3)
    real(dp), intent(in) :: kdispl(3)
    real(dp), intent(in) :: Qtot,Temp
    integer, intent(in) :: istep, ia1
    
! ************************
! * LOCAL variables      *
! ************************
    integer :: iu
    integer :: ispin, i,j
    integer :: maxnhg

#ifdef MPI
    integer,  allocatable :: numhg(:), listhptrg(:), listhg(:)
    real(dp), allocatable :: xijg(:,:), Mg(:)
#endif

    external :: io_assign, io_close

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE ts_io_write' )
#endif

#ifdef MPI
    call glob_sparse_arrays(no_l,no_u,no_s,maxnh, &
         numh ,listhptr ,listh ,xij , Gamma,&
         numhg,listhptrg,maxnhg,listhg,xijg)
    allocate(Mg(maxnhg))
    call memory('A','D',maxnhg,'globArrays')

    call glob_sparse_matrix(no_l,no_u,no_s, &
         maxnh,  numh , listhptr , S , &
         maxnhg, numhg, listhptrg, Mg)
#else
    maxnhg = maxnh
#endif

    if (IONode) then
! Open file
       call io_assign( iu )
       open( iu, file=filename, form='unformatted', status='unknown' )

! Write Dimensions Information
       write(iu) na_u, no_u, no_s, nspin, maxnhg
       
! Write Geometry information
       write(iu) xa
       write(iu) iza
       write(iu) ucell

! Write k-point samplung information
       write(iu) Gamma
! Is this only an S containing object?
       write(iu) onlyS
       write(iu) TSGamma
       write(iu) kscell
       write(iu) kdispl
       write(iu) istep, ia1

       write(iu) lasto

       if (.not.Gamma) then
          write(iu) indxuo
       endif

#ifdef MPI
       write(iu) numhg
#else
       write(iu) numh
#endif

! Write Electronic Structure Information
       write(iu) Qtot,Temp
       write(iu) Ef

! Write listh
       do i = 1 , no_u
#ifdef MPI
          write(iu) listhg(listhptrg(i)+1:listhptrg(i)+numhg(i))
#else
          write(iu) listh(listhptr(i)+1:listhptr(i)+numh(i))
#endif
       end do

! Write Overlap matrix
       do i = 1 , no_u
#ifdef MPI
          write(iu) Mg(listhptrg(i)+1:listhptrg(i)+numhg(i))
#else
          write(iu) S(listhptr(i)+1:listhptr(i)+numh(i))
#endif
       end do

    end if
    
    if (.not. onlyS) then
! Write Hamiltonian	 
       do ispin = 1 , nspin	 
#ifdef MPI
          call glob_sparse_matrix(no_l,no_u,no_s, &
               maxnh,  numh , listhptr , H(:,ispin), &
               maxnhg, numhg, listhptrg, Mg)
#endif
          if ( IONode ) then
             do i = 1 , no_u
#ifdef MPI
                write(iu) Mg(listhptrg(i)+1:listhptrg(i)+numhg(i))
#else
                write(iu) H(listhptr(i)+1:listhptr(i)+numh(i),ispin)
#endif
             end do
          end if
       end do
    end if  ! onlyS

    if (IONode) then
       if (.not.Gamma) then
          do i = 1 , no_u
#ifdef MPI
             write(iu) (xijg(j,listhptrg(i)+1:listhptrg(i)+numhg(i)),j=1,3)
#else
             write(iu) (xij(j,listhptr(i)+1:listhptr(i)+numh(i)),j=1,3)
#endif
          end do
       end if

! Close file
       call io_close( iu )
    end if

#ifdef MPI
    call glob_sparse_arrays_dealloc(no_u, Gamma, &
         maxnhg, &
         numhg,listhptrg,listhg,xijg)
    call glob_sparse_matrix_dealloc(maxnhg, Mg)
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS ts_io_write' )
#endif

  end subroutine ts_write_tshs

#endif

  subroutine ts_read_tshs(filename, &
       onlyS, Gamma, TSGamma, ucell, na_u, no_l, no_u, no_s, maxnh, nspin,  &
       kscell, kdispl, &
       xa, iza, lasto, &
       numh, listhptr, listh, xij, indxuo, &
       H, S, Ef, &
       Qtot, Temp, &
       istep, ia1, &
       Bcast)

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
! logical       Gamma         : Is only gamma point used?
! logical       TSGamma       : Is only TS gamma point used?
! logical       onlyS         : Should only overlap matrix be saved?
! ******************** INPUT or OUTPUT (depending on task) ***********
! integer no_u                : Number of basis orbitals per unit cell
! integer no_s                : Number of basis orbitals per supercell
! integer Enspin              : Spin polarization (1 or 2)
! integer indxuo(no_s)        : Index of orbitals in supercell
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(maxnh,Enspin)     : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! TSS Begin
! ********************* ADDED ARGUMENTS FOR TRANSIESTA ****************
! integer fnlength            : file name length
! character(fnlength)  fname  : file nema for input or output
! integer na_u                : Number of atoms per unit cell
! integer istep, ia1          : Force constant step and atom number
! logical check_kcell        : Do a check of the kscell and kdispl
! TSS End
! *************************** UNITS ***********************************
! Units should be consistent between task='read' and 'write'
! *********************************************************************

!
!  Modules
!
    use precision,    only : dp
    use parallel,     only : IONode
    use sys,          only : die
#ifdef MPI
    use mpi_siesta
#endif

    implicit none

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(out) :: onlyS
    logical, intent(out) :: Gamma, TSGamma
    real(dp), intent(out) :: ucell(3,3)
    integer, intent(out) :: na_u, no_l, no_u, no_s, maxnh, nspin
    real(dp), pointer, intent(out) :: xa(:,:) ! (3,na_u)
    integer, pointer, intent(out) :: iza(:) !(na_u)
    integer, pointer, intent(out) :: numh(:), listhptr(:) ! (no_u)
    integer, pointer, intent(out) :: listh(:) !(maxnh)
    real(dp), pointer, intent(out) :: xij(:,:) ! (3,maxnh)
    integer, pointer, intent(out) :: indxuo(:) ! (no_s)
    integer, pointer, intent(out) :: lasto(:) ! (0:na_u) 
    real(dp), pointer, intent(out) :: H(:,:), S(:) !(maxnh,nspin),(maxnh)
    real(dp), intent(out) :: Ef
    integer, intent(out) :: kscell(3,3)
    real(dp), intent(out) :: kdispl(3)
    real(dp), intent(out) :: Qtot,Temp
    ! These have to be set before entrance (makes it possible to read
    ! in FCrun TSHS files...
    integer, intent(out) :: istep, ia1
    ! If true it will broadcast every information within the code...
    logical, optional, intent(in) :: Bcast
    
! ************************
! * LOCAL variables      *
! ************************
    integer :: iu
    integer :: ispin,i,j,all_I(8)
    integer :: maxnhg
    logical :: lBcast, exist
#ifdef MPI
    integer :: MPIerror
#endif

    external :: io_assign, io_close

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE ts_io_read' )
#endif

    ! Determine whether to broadcast afterwards
    lBcast = .false.
    if ( present(Bcast) ) then
       lBcast = Bcast
    end if

    if (IONode) then
       inquire(file=filename,exist=exist)
       if ( .not. exist ) then
          call die('ERROR: Could not read '//trim(filename)//'.')
       end if
! Open file
       call io_assign( iu )
       open( iu, file=filename, form='unformatted', status='old' )

! Read Dimensions Information
       read(iu) na_u, no_u, no_s, nspin, maxnh
       no_l = no_u
       
! Read Geometry information
       nullify(xa,iza)
       allocate(xa(3,na_u),iza(na_u)) 
       call memory('A','D',3*na_u,'iohs')
       call memory('A','I',na_u,'iohs')
       read(iu) xa
       read(iu) iza
       read(iu) ucell

! Read k-point samplung information
       read(iu) Gamma
! Is this only an S containing object?
       read(iu) onlyS
       read(iu) TSGamma
       read(iu) kscell
       read(iu) kdispl
       read(iu) istep, ia1

       nullify(lasto)
       allocate(lasto(0:na_u))
       call memory('A','I',1+na_u,'iohs')
       read(iu) lasto

       if (.not.Gamma) then
          nullify(indxuo)
          allocate(indxuo(1:no_s))
          call memory('A','I',no_s,'iohs')
          read(iu) indxuo
       endif

       nullify(numh,listhptr)
       allocate(numh(no_u))
       call memory('A','I',no_u,'iohs')
       read(iu) numh

       call memory('A','I',no_u,'iohs')
       allocate(listhptr(no_u))
       listhptr(1) = 0
       do i = 2 , no_u
          listhptr(i) = listhptr(i-1) + numh(i-1)
       end do

! Read Electronic Structure Information
       read(iu) Qtot,Temp
       read(iu) Ef

! Read listh
       nullify(listh)
       allocate(listh(maxnh))
       call memory('A','I',maxnh,'iohs')
       do i = 1 , no_u
          read(iu) listh(listhptr(i)+1:listhptr(i)+numh(i))
       end do

! Read Overlap matrix
       nullify(S)
       allocate(S(maxnh))
       call memory('A','D',maxnh,'iohs')
       do i = 1 , no_u
          read(iu) S(listhptr(i)+1:listhptr(i)+numh(i))
       end do

       if (.not. onlyS) then
          nullify(H)
          allocate(H(maxnh,nspin))
          call memory('A','D',maxnh*nspin,'iohs')
! Read Hamiltonian	 
          do ispin = 1 , nspin	 
             do i = 1 , no_u
                read(iu) H(listhptr(i)+1:listhptr(i)+numh(i),ispin)
             end do
          end do
       end if  ! onlyS

       if (.not.Gamma) then
          nullify(xij)
          allocate(xij(3,maxnh))
          call memory('A','D',3*maxnh,'iohs')
          do i = 1 , no_u
             read(iu) (xij(j,listhptr(i)+1:listhptr(i)+numh(i)),j=1,3)
          end do
       end if

! Close file
       call io_close( iu )
    end if

    if ( lBcast ) then
#ifdef MPI
       all_I(1:8) = (/na_u, no_l, no_u, no_s, maxnh, nspin,istep,ia1 /)
       call MPI_Bcast(all_I(1),8,MPI_Integer,0,MPI_Comm_World,MPIerror)
       na_u = all_I(1)
       no_l = all_I(2)
       no_u = all_I(3)
       no_s = all_I(4)
       maxnh = all_I(5)
       nspin = all_I(6)
       istep = all_I(7)
       ia1 = all_I(8)
       call MPI_Bcast(Qtot,1,DAT_double,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(Temp,1,DAT_double,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(Ef,1,DAT_double,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(ucell(1,1),9,DAT_double,0, &
            MPI_Comm_World,MPIerror)
       call MPI_Bcast(kscell(1,1),9,MPI_Integer,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(kdispl(1),3,DAT_Double,0,MPI_Comm_World,MPIerror)

       if ( .not. IONode ) then
          if ( .not. Gamma ) then 
! they behave as dummy arrays in case of Gamma == .true.
             allocate(indxuo(no_s))
             call memory('A','I',no_s,'iohs')
             allocate(xij(3,maxnh))
             call memory('A','D',maxnh*3,'iohs')
          end if
          allocate(xa(3,na_u))
          call memory('A','D',3*na_u,'iohs')
          allocate(iza(na_u))
          call memory('A','I',na_u,'iohs')
          allocate(lasto(0:na_u))
          call memory('A','I',1+na_u,'iohs')
          allocate(numh(no_u))
          call memory('A','I',no_u,'iohs')
          allocate(listhptr(no_u))
          call memory('A','I',no_u,'iohs')
          allocate(listh(maxnh))
          call memory('A','I',maxnh,'iohs')
          allocate(H(maxnh,nspin))
          call memory('A','D',maxnh*nspin,'iohs')
          allocate(S(maxnh))
          call memory('A','D',maxnh,'iohs')
       end if
       call MPI_Bcast(onlyS,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(Gamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(TSGamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(xa(1,1),3*na_u,DAT_Double,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(iza,na_u,MPI_Integer,0,MPI_Comm_World,MPIerror)
       if ( .not. Gamma ) then
          call MPI_Bcast(indxuo,no_s,MPI_Integer,0, MPI_Comm_World,MPIerror)
          call MPI_Bcast(xij(1,1),3*maxnh,DAT_Double,0, &
               MPI_Comm_World,MPIerror)
       end if
       call MPI_Bcast(lasto(0),1+na_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(numh,no_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(listhptr,no_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(listh,maxnh,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(H(1,1),maxnh*nspin,DAT_Double,0, &
            MPI_Comm_World,MPIerror)
       call MPI_Bcast(S,maxnh,DAT_Double,0, MPI_Comm_World,MPIerror)
#endif
    end if

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS ts_io_read' )
#endif

  end subroutine ts_read_tshs


  function filenameTSHS(slabel,istep,ia1,onlyS)
    character(len=*), intent(in) :: slabel
    integer, intent(in) :: istep, ia1
    logical, intent(in) :: onlyS
    character(len=255) :: filenameTSHS
    integer :: fL
    
    interface
       function paste(s1,s2)
         character(len=*), intent(in) :: s1,s2
         character(len=255) :: paste
       end function paste
    end interface

    ! Initialize...
    filenameTSHS = ' '

    ! This is the criteria for not doing an FCrun
    ! There will never be an atom denoted by index 0
    if ( ia1 /= 0 ) then
       write(filenameTSHS,'(2i4)') ia1, istep
       filenameTSHS = paste( slabel, filenameTSHS(1:8) )
    else
       filenameTSHS = slabel
    end if
    
    fL = len_trim(filenameTSHS)
    if ( onlyS ) then
       if ( filenameTSHS(fL-5:fL) .ne. '.onlyS' ) then
          filenameTSHS = paste( filenameTSHS , '.onlyS' )
       end if
    else
       if ( filenameTSHS(fL-4:fL) .ne. '.TSHS' ) then
          filenameTSHS = paste( filenameTSHS , '.TSHS' )
       end if
    end if

  end function filenameTSHS

end module m_ts_io
