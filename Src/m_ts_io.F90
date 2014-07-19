module m_ts_io

  use precision, only : dp
  use parallel, only : IONode

  implicit none

  public :: ts_read_TSHS_opt
  public :: ts_read_TSHS
  public :: ts_write_TSHS
  public :: fname_TSHS
  public :: TSHS_version

  private

contains

  subroutine ts_read_TSHS_opt(TSHS,DUMMY,na_u,no_u,no_s,nspin,maxnh, &
       xa, ucell, nsc, Qtot, Temp, Ef, &
       Gamma,Gamma_TS,kscell,kdispl,OnlyS,lasto, &
       Bcast)
#ifdef MPI
    use mpi_siesta
#endif
! ***********************
! * INPUT variables     *
! ***********************
    character(len=*), intent(in) :: TSHS

! ***********************
! * OUTPUT variables    *
! *********************** 
    integer, optional :: DUMMY ! MUST NEVER BE PASSED
    integer, intent(out), optional :: na_u, no_u, no_s, nspin, maxnh, lasto(:), kscell(3,3), nsc(3)
    real(dp), intent(out), optional :: xa(:,:), ucell(3,3), Qtot, Temp, Ef, kdispl(3)
    logical, intent(out), optional :: Gamma, Gamma_TS, OnlyS
    logical, intent(in), optional :: Bcast

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: uTSHS, tmp(5), version, tkscell(3,3), tnsc(3)
    real(dp), allocatable :: txa(:,:)
    real(dp) :: rtmp(3), tucell(3,3)
    logical :: fGamma, bool(3)
#ifdef MPI
    integer :: buffer_size, ipos
    character(len=1), allocatable :: buffer(:)
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    if ( present(DUMMY) ) call die('ts_read_TSHS_opt: Arguments has to be &
         &named. Please correct sources.')

    version = tshs_version(tshs)

    if ( IONode ) then

       select case ( version )
       case ( 0 , 1 )
          ! do nothing
       case default
          call die('Unsupported TSHS version file')
       end select

       call io_assign(uTSHS)
       open(file=trim(TSHS),unit=uTSHS,form='unformatted')
       if ( version /= 0 ) then
          read(uTSHS) ! version number
       end if
       read(uTSHS) tmp(1:5) !na_u, no_u, no_s, Enspin, maxnh
       allocate(txa(3,tmp(1)))
       if ( present(na_u) ) na_u = tmp(1)
       if ( present(no_u) ) no_u = tmp(2)
       if ( present(no_s) ) no_s = tmp(3)
       if ( present(nspin) ) nspin = tmp(4)
       if ( present(maxnh) ) maxnh = tmp(5)
       if ( version == 0 ) then
          read(uTSHS) txa
          read(uTSHS) ! iza
          read(uTSHS) tucell
          tnsc = 0
       else if ( version == 1 ) then
          read(uTSHS) tnsc ! corrected nsc
          read(uTSHS) tucell, txa
       end if
       if ( present(nsc) ) nsc = tnsc
       if ( present(ucell) ) ucell = tucell
       if ( present(xa) )    xa    = txa
       if ( version == 0 ) then
          read(uTSHS) fGamma ! SIESTA_Gamma
          if ( present(Gamma) ) Gamma = fGamma
          if ( present(OnlyS) ) then
             read(uTSHS) OnlyS
          else
             read(uTSHS) ! OnlyS
          end if
          
          if ( present(Gamma_TS) ) then
             read(uTSHS) Gamma_TS
          else
             read(uTSHS) ! Gamma_TS
          end if
          if ( present(kscell) ) then
             read(uTSHS) kscell
          else
             read(uTSHS) ! ts_kscell_file
          end if
          if ( present(kdispl) ) then
             read(uTSHS) kdispl
          else
             read(uTSHS) ! ts_kdispl_file  
          end if
       else if ( version == 1 ) then
          read(uTSHS) fGamma, bool(1:2)
          if ( present(Gamma) ) Gamma = fGamma
          if ( present(Gamma_TS) ) Gamma_TS = bool(1)
          if ( present(OnlyS) ) OnlyS = bool(2)
          read(uTSHS) tkscell, rtmp
          if ( present(kscell) ) kscell = tkscell
          if ( present(kdispl) ) kdispl = rtmp
          read(uTSHS) rtmp(1:3)
          if ( present(Ef) )   Ef   = rtmp(1)
          if ( present(Qtot) ) Qtot = rtmp(2)
          if ( present(Temp) ) Temp = rtmp(3)
       end if
       read(uTSHS) ! istep, ia1

       if ( present(lasto) ) then
          if ( size(lasto) /= tmp(1)+1 ) call die('ts_read_TSHS: Wrong size of lasto')
          read(uTSHS) lasto
       else
          read(uTSHS) ! lasto
       end if

       if ( version == 0 .and. .not. fGamma ) then
          read(uTSHS) ! indxuo
       end if

       read(uTSHS) ! numh

       if ( version == 0 ) then
          read(uTSHS) rtmp(1:2)
          if ( present(Qtot) ) Qtot = rtmp(1)
          if ( present(Temp) ) Temp = rtmp(2)

          if ( present(Ef) ) then
             read(uTSHS) Ef
          else
             read(uTSHS) ! Ef
          end if
       end if
       
       deallocate(txa)
       call io_close(uTSHS)

    end if

#ifdef MPI
    if ( present(Bcast) ) then
       ! if we do not request broadcasting, then return...
       if ( .not. Bcast ) return
    end if

    ! Broadcast na_u (for easy reference)
    call MPI_Bcast(tmp(1),1,MPI_Integer,0,MPI_Comm_World,MPIerror)

#ifdef MPI_OLD
    if ( present(na_u) ) &
         call MPI_Bcast(na_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(no_u) ) &
         call MPI_Bcast(no_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(no_s) ) &
         call MPI_Bcast(no_s,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(nspin) ) &
         call MPI_Bcast(nspin,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(maxnh) ) &
         call MPI_Bcast(maxnh,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(nsc) ) &
         call MPI_Bcast(nsc,3,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(xa) ) &
         call MPI_Bcast(xa(1,1),3*tmp(1),MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(ucell) ) &
         call MPI_Bcast(ucell(1,1),9,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(Gamma) ) &
         call MPI_Bcast(Gamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
    if ( present(OnlyS) ) &
         call MPI_Bcast(OnlyS,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
    if ( present(Gamma_TS) ) &
         call MPI_Bcast(Gamma_TS,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
    if ( present(kscell) ) &
         call MPI_Bcast(kscell,9,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(kdispl) ) &
         call MPI_Bcast(kdispl,3,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(lasto) ) &
         call MPI_Bcast(lasto(1),tmp(1)+1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(Qtot) ) &
         call MPI_Bcast(Qtot,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(Temp) ) &
         call MPI_Bcast(Temp,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(Ef) ) &
         call MPI_Bcast(Ef,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
#else

    ! this should be more than enough...
    buffer_size = 8 * (tmp(1) * 6 + 103)
    allocate(buffer(buffer_size))
    ! position of data in buffer...
    ipos = 0


    if ( IONode ) then
       if ( present(na_u) ) & !  4
            call MPI_Pack(na_u,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(no_u) ) & !  4
            call MPI_Pack(no_u,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(no_s) ) & !  4
            call MPI_Pack(no_s,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(nspin) ) & !  4
            call MPI_Pack(nspin,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(maxnh) ) & !  4
            call MPI_Pack(maxnh,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(xa) ) &    ! 8 * 3 * na_u
            call MPI_Pack(xa(1,1),3*tmp(1),MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(ucell) ) & ! 8 * 3 * 3
            call MPI_Pack(ucell(1,1),9,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(nsc) ) &   ! 4
            call MPI_Pack(nsc(1),3,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Gamma) ) & ! 4
            call MPI_Pack(Gamma,1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(kscell) ) &! 4 * 3 * 3
            call MPI_Pack(kscell(1,1),9,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(kdispl) ) &! 8 * 3
            call MPI_Pack(kdispl(1),3,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(OnlyS) ) & ! 4
            call MPI_Pack(OnlyS,1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Gamma_TS) ) &! 4
            call MPI_Pack(Gamma_TS,1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(lasto) ) & ! 4 * (na_u+1)
            call MPI_Pack(lasto(1),tmp(1)+1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Qtot) ) &  ! 8
            call MPI_Pack(Qtot,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Temp) ) &  ! 8
            call MPI_Pack(Temp,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Ef) ) &    ! 8
            call MPI_Pack(Ef,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)

       if ( ipos >= buffer_size .or. ipos < 0 .or. MPIerror /= MPI_Success ) then
          call die('Error in estimating the buffer-size for the &
               &TSHS reading. Please contact the developers')
       end if

    end if

    call MPI_Bcast(buffer,buffer_size,MPI_Packed, &
         0, MPI_Comm_World, MPIerror)

    if ( .not. IONode ) then
       if ( present(na_u) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            na_u,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(no_u) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            no_u,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(no_s) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            no_s,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(nspin) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            nspin,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(maxnh) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            maxnh,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(xa) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            xa(1,1),3*tmp(1),MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(ucell) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            ucell(1,1),9,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(nsc) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            nsc(1),3,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(Gamma) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Gamma,1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(kscell) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            kscell(1,1),9,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(kdispl) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            kdispl(1),3,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(OnlyS) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            OnlyS,1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(Gamma_TS) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Gamma_TS,1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(lasto) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            lasto(1),tmp(1)+1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(Qtot) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Qtot,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(Temp) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Temp,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(Ef) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Ef,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
    end if       

    deallocate(buffer)

#endif

#endif

  end subroutine ts_read_TSHS_opt

  subroutine ts_read_TSHS(filename, onlyS, Gamma, TSGamma, &
       ucell, nsc, na_u, no_l, no_u, no_s, maxnh, nspin,  &
       kscell, kdispl, &
       xa, lasto, &
       numh, listhptr, listh, xij, &
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

    use sys,          only : die
#ifdef MPI
    use mpi_siesta
#endif
    use geom_helper,  only : iaorb, ucorb

    implicit none

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(out) :: onlyS
    logical, intent(out) :: Gamma, TSGamma
    real(dp), intent(out) :: ucell(3,3)
    integer, intent(out) :: na_u, no_l, no_u, no_s, maxnh, nspin, nsc(3)
    real(dp), pointer, intent(out) :: xa(:,:) ! (3,na_u)
    integer, pointer, intent(out) :: numh(:), listhptr(:) ! (no_u)
    integer, pointer, intent(out) :: listh(:) !(maxnh)
    real(dp), pointer, intent(out) :: xij(:,:) ! (3,maxnh)
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
    integer :: iu, version, n_s
    integer :: ispin,i,j, all_I(8), ind, io, ia, ja, tm(3), is
    integer, allocatable :: offsets(:,:)
    integer, allocatable :: indxuo(:)
    logical :: lBcast, exist
#ifdef MPI
    integer :: MPIerror
#endif

    external :: io_assign, io_close

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE ts_io_read' )
#endif

    nullify(xa,lasto,numh,listhptr,listh,S,H,xij)

    ! Determine whether to broadcast afterwards
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    ! Get file version
    version = tshs_version(filename)

    nsc(:)  = 0

    if ( IONode ) then
       select case ( version )
       case ( 0 , 1 )
          ! do nothing
       case default
          call die('Unsupported TSHS version file')
       end select
       inquire(file=filename,exist=exist)
       if ( .not. exist ) then
          call die('ERROR: Could not read '//trim(filename)//'.')
       end if

       ! Open file
       call io_assign( iu )
       open( iu, file=filename, form='unformatted', status='old' )
       
       ! Read Dimensions Information
       if ( version /= 0 ) then
          read(iu) version
       end if
       read(iu) na_u, no_u, no_s, nspin, maxnh
       
       no_l = no_u
       
       ! Read Geometry information
       allocate(xa(3,na_u)) 
       call memory('A','D',3*na_u,'iohs')
       call memory('A','I',na_u,'iohs')
       if ( version == 0 ) then
          read(iu) xa
          read(iu) ! iza
          read(iu) ucell
       else if ( version == 1 ) then
          read(iu) nsc ! original nsc, and the corrected nsc
          read(iu) ucell, xa
       end if

       ! Read k-point sampling information
       if ( version == 0 ) then
          read(iu) Gamma
          read(iu) onlyS
          read(iu) TSGamma
          read(iu) kscell
          read(iu) kdispl
       else if ( version == 1 ) then
          read(iu) Gamma, TSGamma, onlyS
          read(iu) kscell, kdispl
          read(iu) Ef, Qtot, Temp
       end if
       read(iu) istep, ia1

       allocate(lasto(0:na_u))
       call memory('A','I',1+na_u,'iohs')
       read(iu) lasto

       if ( version == 0 ) then
          allocate(indxuo(no_s))
          read(iu) indxuo
          do i = 1 , no_s
             if ( indxuo(i) /= ucorb(i,no_u) ) &
                  call die('Error in indxuo, not recognized')
          end do
          deallocate(indxuo)
       end if

       allocate(numh(no_u))
       call memory('A','I',no_u,'iohs')
       read(iu) numh

       call memory('A','I',no_u,'iohs')
       allocate(listhptr(no_u))
       listhptr(1) = 0
       do i = 2 , no_u
          listhptr(i) = listhptr(i-1) + numh(i-1)
       end do

       if ( version == 0 ) then
          ! Read Electronic Structure Information
          read(iu) Qtot,Temp
          read(iu) Ef
       end if

       ! Read listh
       allocate(listh(maxnh))
       call memory('A','I',maxnh,'iohs')
       do i = 1 , no_u
          read(iu) listh(listhptr(i)+1:listhptr(i)+numh(i))
       end do

       ! Read Overlap matrix
       allocate(S(maxnh))
       call memory('A','D',maxnh,'iohs')
       do i = 1 , no_u
          read(iu) S(listhptr(i)+1:listhptr(i)+numh(i))
       end do

       if ( .not. onlyS ) then
          allocate(H(maxnh,nspin))
          call memory('A','D',maxnh*nspin,'iohs')
          do ispin = 1 , nspin	 
             do i = 1 , no_u
                read(iu) H(listhptr(i)+1:listhptr(i)+numh(i),ispin)
             end do
          end do
       end if  ! onlyS

       if ( .not. Gamma ) then

          allocate(xij(3,maxnh))
          call memory('A','D',3*maxnh,'iohs')

          if ( version == 0 ) then
             do i = 1 , no_u
                read(iu) (xij(j,listhptr(i)+1:listhptr(i)+numh(i)),j=1,3)
             end do

          else if ( version == 1 ) then
             
             ! Number of supercells
             n_s = no_s / no_u
             allocate(offsets(3,0:n_s-1))
             read(iu) offsets

             do io = 1 , no_u
                ia = iaorb(io,lasto)
                do j = 1 , numh(io)
                   ind = listhptr(io) + j
                   ja = iaorb(ucorb(listh(ind),no_u),lasto)

                   is = ( listh(ind) - 1 ) / no_u
                   tm(:) = offsets(:,is)

                   xij(:,ind) = ucell(:,1) * tm(1) &
                        + ucell(:,2) * tm(2) &
                        + ucell(:,3) * tm(3) &
                        + xa(:,ja) - xa(:,ia)

                end do
             end do
   
             ! clean-up
             deallocate(offsets)

          end if
       end if

       ! Close file
       call io_close( iu )
    end if

    if ( lBcast ) then
#ifdef MPI
       all_I(1:8) = (/na_u, no_l, no_u, no_s, maxnh, nspin,istep,ia1 /)
       call MPI_Bcast(all_I,8,MPI_Integer,0,MPI_Comm_World,MPIerror)
       na_u = all_I(1)
       no_l = all_I(2)
       no_u = all_I(3)
       no_s = all_I(4)
       maxnh = all_I(5)
       nspin = all_I(6)
       istep = all_I(7)
       ia1 = all_I(8)
       call MPI_Bcast(Gamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(TSGamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(Qtot,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(nsc,3,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(Temp,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(Ef,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(ucell(1,1),9,MPI_Double_Precision,0, &
            MPI_Comm_World,MPIerror)
       call MPI_Bcast(kscell(1,1),9,MPI_Integer,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(kdispl(1),3,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(onlyS,1,MPI_Logical,0,MPI_Comm_World,MPIerror)

       if ( .not. IONode ) then
          if ( .not. Gamma ) then 
             ! they behave as dummy arrays in case of Gamma == .true.
             allocate(xij(3,maxnh))
             call memory('A','D',maxnh*3,'iohs')
          end if
          allocate(xa(3,na_u))
          call memory('A','D',3*na_u,'iohs')
          allocate(lasto(0:na_u))
          call memory('A','I',1+na_u,'iohs')
          allocate(numh(no_u))
          call memory('A','I',no_u,'iohs')
          allocate(listhptr(no_u))
          call memory('A','I',no_u,'iohs')
          allocate(listh(maxnh))
          call memory('A','I',maxnh,'iohs')
          if ( .not. onlyS ) then
             allocate(H(maxnh,nspin))
             call memory('A','D',maxnh*nspin,'iohs')
          end if
          allocate(S(maxnh))
          call memory('A','D',maxnh,'iohs')
       end if
       call MPI_Bcast(xa(1,1),3*na_u,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
       if ( .not. Gamma ) then
          call MPI_Bcast(xij(1,1),3*maxnh,MPI_Double_Precision,0, &
               MPI_Comm_World,MPIerror)
       end if
       call MPI_Bcast(lasto(0),1+na_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(numh,no_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(listh,maxnh,MPI_Integer,0, MPI_Comm_World,MPIerror)
       if ( .not. onlyS ) then
          call MPI_Bcast(H(1,1),maxnh*nspin,MPI_Double_Precision,0, &
               MPI_Comm_World,MPIerror)
       end if
       call MPI_Bcast(S,maxnh,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
#endif
    end if

    if ( .not. IONode ) then
       listhptr(1) = 0
       do i = 2 , no_u
          listhptr(i) = listhptr(i-1) + numh(i-1)
       end do
    end if

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS ts_io_read' )
#endif

  end subroutine ts_read_TSHS

  subroutine ts_write_TSHS(filename, &
       onlyS, Gamma, TSGamma, &
       ucell, nsc, na_u, no_l, no_u, no_s, maxnh, nspin,  &
       kscell, kdispl, &
       xa, lasto, &
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
! *********************************************************************

    use m_hs_matrix, only : nsc_to_offsets, offset2idx, list_col_correct
    use m_hs_matrix, only : set_HS_available_transfers
    use geom_helper, only : ucorb

#ifdef MPI
    use m_glob_sparse
#endif

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(in) :: onlyS
    logical, intent(in) :: Gamma, TSGamma
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: nsc(3), na_u, no_l, no_u, no_s, maxnh, nspin
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: numh(no_l), listhptr(no_l)
    integer, intent(in) :: listh(maxnh)
    real(dp), intent(in) :: xij(3,maxnh)
    integer, intent(in) :: indxuo(no_s)
    integer, intent(in) :: lasto(0:na_u)
    real(dp), intent(in) :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: Ef
    integer, intent(in) :: kscell(3,3)
    real(dp), intent(in) :: kdispl(3)
    real(dp), intent(in) :: Qtot, Temp
    integer, intent(in) :: istep, ia1
    
! ************************
! * LOCAL variables      *
! ************************
    integer :: iu
    integer :: ispin, i, n_s
    integer :: maxnhg
    integer :: lnsc(3), lno_s, tms(2,3)
    integer, pointer :: offsets(:,:)
    integer, allocatable :: listhg(:)
#ifdef MPI
    integer,  allocatable :: numhg(:), listhptrg(:)
    real(dp), allocatable :: xijg(:,:), Mg(:)
#endif

    external :: io_assign, io_close

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE ts_io_write' )
#endif

    ! This setting is (I am afraid) not constant
    ! with system. I had suspected that n_s ALWAYS would 
    ! equal no_s / no_u, but this seems not to be the case.
    n_s = no_s / no_u 
    if ( mod(no_s,no_u) /= 0 ) call die('Error in supercell orbitals, no_s')
    if ( n_s /= product(nsc) ) call die('Error in supercell orbitals, nsc')

    ! Check that indxuo is constructed in a sane format
    do i = 1 , no_s
       if ( indxuo(i) /= ucorb(i,no_u) ) &
            call die('Error in indxuo, could not understand the format, &
            &please consult the developers.')
    end do

#ifdef MPI
    call glob_sparse_arrays(no_l,no_u,no_s,maxnh, &
         numh ,listhptr ,listh ,xij , Gamma,&
         numhg,listhptrg,maxnhg,listhg,xijg)

    call set_HS_available_transfers(ucell,na_u,xa,lasto,no_u,maxnhg, &
         xijg,numhg,listhptrg,listhg,tms)

    allocate(Mg(maxnhg))
    call memory('A','D',maxnhg,'globArrays')

    call glob_sparse_matrix(no_l,no_u,no_s, &
         maxnh,  numh , listhptr , S , &
         maxnhg, numhg, listhptrg, Mg)
#else

    call set_HS_available_transfers(ucell,na_u,xa,lasto,no_u,maxnh, &
         xij,numh,listhptr,listh,tms)

    maxnhg = maxnh
    allocate(listhg(maxnhg))
    listhg = listh
#endif

    ! Calculate new supercells
    do i = 1 , 3
       lnsc(i) = 2 * maxval(abs(tms(:,i))) + 1
    end do
    lno_s = product(lnsc) * no_u

    ! Always calculate offsets
    ! If lnsc is then 1, then fine! :)
    call nsc_to_offsets(lnsc,offsets)

    ! Correct the listhg array
#ifdef MPI
    call list_col_correct(ucell,na_u, no_u,maxnhg, &
         lasto, xa, numhg, listhptrg, listhg, xijg, lnsc)
#else
    call list_col_correct(ucell,na_u, no_u,maxnhg, &
         lasto, xa, numh, listhptr, listhg, xij, lnsc)
#endif

    if ( IONode ) then

       ! Open file
       call io_assign( iu )
       open( iu, file=filename, form='unformatted', status='unknown' )

       ! Write file version
       write(iu) 1 ! This is version one of the file format

       ! Write Dimensions Information
       write(iu) na_u, no_u, lno_s, nspin, maxnhg
       write(iu) lnsc ! The "corrected" supercells
       
       ! Write Geometry information
       write(iu) ucell, xa

       ! Write k-point sampling information
       write(iu) Gamma, TSGamma, onlyS
       write(iu) kscell, kdispl

       ! Write Electronic Structure Information
       write(iu) Ef, Qtot, Temp

       ! The phonon-calculation
       write(iu) istep, ia1

       write(iu) lasto

#ifdef MPI
       write(iu) numhg
#else
       write(iu) numh
#endif

       ! Write listh
       do i = 1 , no_u
#ifdef MPI
          write(iu) listhg(listhptrg(i)+1:listhptrg(i)+numhg(i))
#else
          write(iu) listhg(listhptr(i)+1:listhptr(i)+numh(i))
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
    
    if ( .not. onlyS ) then
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

    if ( IONode ) then

       if ( .not. Gamma ) then
          write(iu) offsets
       end if

       ! Close file
       call io_close( iu )

    end if

#ifdef MPI
    call glob_sparse_arrays_dealloc(no_u, Gamma, &
         maxnhg, numhg, listhptrg, listhg, xijg)
    call glob_sparse_matrix_dealloc(maxnhg, Mg)
#else
    deallocate(listhg)
#endif

    deallocate(offsets)
       
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS ts_io_write' )
#endif

  end subroutine ts_write_TSHS

  function fname_TSHS(slabel,istep,ia1,onlyS) result(fname)
    character(len=*), intent(in) :: slabel
    integer, intent(in) :: istep, ia1
    logical, intent(in) :: onlyS
    character(len=255) :: fname
    integer :: fL

    ! Initialize...
    fname = ' '

    ! This is the criteria for not doing an FCrun
    ! There will never be an atom denoted by index 0
    if ( ia1 /= 0 ) then
       write(fname,'(2i4)') ia1, istep
       fname = trim(slabel)//fname(1:8)
    else
       fname = slabel
    end if
    
    fL = len_trim(fname)
    if ( onlyS ) then
       if ( fname(fL-5:fL) .ne. '.onlyS' ) then
          fname = trim(fname)//'.onlyS'
       end if
    else
       if ( fname(fL-4:fL) .ne. '.TSHS' ) then
          fname = trim(fname)//'.TSHS'
       end if
    end if

  end function fname_TSHS

  function TSHS_version(fname) result(version)
    character(len=*), intent(in) :: fname
    integer :: version
    integer :: iu
    integer :: na_u, no_u, no_s, nspin, maxnh, err
    
    external :: io_assign, io_close

    if ( .not. IONode ) return

    ! Open file
    call io_assign( iu )
    open( iu, file=fname, form='unformatted', status='unknown' )

    read(iu,iostat=err) na_u, no_u, no_s, nspin, maxnh
    if ( err == 0 ) then
       ! we can successfully read 5 integers
       version = 0
    else
       backspace(iu)
       read(iu,iostat=err) version
    end if

    call io_close(iu)

  end function TSHS_version

end module m_ts_io
