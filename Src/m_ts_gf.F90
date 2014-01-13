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
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
module m_ts_GF
!
! Routines that are used for reading the Green's function files.
! It is very closely related to m_ts_electrode (consider moving it to there)
! 
! This module constitutes the routines that are needed to ensure the
! correct format of the GF files. 
! In order to limit the size of the m_ts_electrode file this has been created.
!
! A call to check_green is used to ensure the correct format of the GF file.
! It checks as much information as is available and dies if non-conforming
! entities exists. This routine should only be called by IONode!
!
! A call to read_Green will read in the header of the GF file and do "basic"
! checks against array sizes. Thus a call to check_green is adviced before 
! read_Green!
! This routine will also distribute the arrays in an MPI run.

! In both routines i have added a prefix "c_" which clarifies what is to
! be checked and what is the read in parameter of the file.
!
!=============================================================================
! CONTAINS:
!          1) do_Green
!          2) read_Green
!          3) check_Green


  implicit none

  public :: do_Green, read_Green, check_Green

  private

contains

  subroutine do_Green(El, optReUseGF, &
       ucell,nkpnt,kpoint,kweight, &
       RemUCellDistance,xa_EPS, &
       CalcDOS)
    
    use precision,  only : dp
    use parallel  , only : IONode
    use sys ,       only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Bcast, MPI_Integer, MPI_Logical
#endif
    use m_ts_cctype
    use m_ts_electype
    use m_ts_electrode, only : create_Green

    use m_ts_contour_eq
    use m_ts_contour_neq

    implicit none
    
! ***********************
! * INPUT variables     *
! ***********************
    type(Elec), intent(inout)     :: El
    logical, intent(in)           :: optReUseGF ! Should we re-use the GF files if they exists?    
    integer, intent(in)           :: nkpnt ! Number of k-points
    real(dp), intent(in)          :: kpoint(3,nkpnt) ! k-points
    real(dp), intent(in)          :: kweight(nkpnt) ! weights of kpoints
    logical, intent(in)           :: RemUCellDistance ! Whether to remove the unit cell distance in the Hamiltonian.
    real(dp), intent(in)          :: xa_Eps ! coordinate precision check
    real(dp), dimension(3,3)      :: ucell ! The unit cell of the CONTACT
    logical, intent(in)           :: CalcDOS

! ***********************
! * LOCAL variables     *
! ***********************
    complex(dp), allocatable :: ZBulkDOS(:,:) ! DOS at energy points
    integer :: uGF, i, iE, NEn
    logical :: errorGF, exist, ReUseGF
    character(len=NAME_LEN) :: rGFtitle
    complex(dp), allocatable :: ce(:)
    type(ts_c_idx) :: c
#ifdef MPI
    integer :: MPIerror
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE do_Green' )
#endif
    
! check the file for existance
    inquire(file=trim(GFfile(El)),exist=exist)
    
    ReUseGF = optReUseGF
! If it does not find the file, calculate the GF
    if ( exist ) then
       if (IONode ) then
          write(*,*) "Electrode Green's function file: '"//&
               trim(GFfile(El))//"' already exist."
          if ( .not. ReUseGF ) then
             write(*,*)"Green's function file '"//&
                  trim(GFfile(El))//"' is requested overwritten."
          end if
       end if
    else
       ReUseGF = .false.
    end if

#ifdef MPI
    call MPI_Bcast(ReUseGF,1,MPI_Logical,0, &
         MPI_Comm_World,MPIerror)
#endif
   
    errorGF = .false.

    ! we need to create all the contours
    NEn = N_Eq_E() + N_nEq_E()
    allocate(ce(NEn))
    iE = 0
    do i = 1 , N_Eq_E()
       c = Eq_E(i)
       ce(i) = c%e
    end do
    iE = N_Eq_E()
    do i = 1 , N_nEq_E()
       c = nEq_E(i)
       ce(iE+i) = c%e
    end do
       
    ! We return if we should not calculate it
    if ( .not. ReUseGF ) then

       call create_Green(El, &
            ucell,nkpnt,kpoint,kweight, &
            RemUCellDistance, &
            NEn,ce, &
            CalcDOS,ZBulkDOS)

    else

       ! Check that the Green's functions are correct!
       ! This is needed as create_Green returns if user requests not to
       ! overwrite an already existing file.
       ! This check will read in the number of orbitals and atoms in the
       ! electrode surface Green's function.
       ! Check the GF file
       if ( IONode ) then
          call io_assign(uGF)
          open(file=GFfile(El),unit=uGF,form='UNFORMATTED')
          
          ! Read in the title and rewind
          read(uGF) rGfTitle
          rewind(uGF)
          
          call check_Green(uGF,El, &
               ucell,nkpnt,kpoint,kweight, &
               NEn, ce, &
               RemUCellDistance, xa_Eps, errorGF)
          
          write(*,'(/,4a,/)') "Using GF-file '",trim(GFfile(El)), &
               "' with title: '",trim(rGfTitle)//"'"
          
          call io_close(uGF)
       endif

    end if

!     Check the error in the GF file
#ifdef MPI
    call MPI_Bcast(errorGF,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
#endif
    if ( errorGF ) &
         call die("Error in GFfile: "//trim(GFFile(El))//". Please move or delete")

    deallocate(ce)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS do_Green' )
#endif

  end subroutine do_Green


! This routine requires a call to check_Green before.
! It will return the k-points of the electrode Green's function file

! ##################################################################
! ##            Read-in header of Greens function file            ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                                                              ##
! ## Changed to F90 by Jose-Luis Mozos, jlm@icmab.es              ##
! ## Changed to only read header by Nick P. Andersen              ##
! ##    will only check against integer information and Ef shift. ##
! ##################################################################
  subroutine read_Green(funit,El,c_nkpar,c_NEn, c_RemUCell)
    
    use precision, only : dp
    use parallel,  only : IONode
    use sys ,      only : die
#ifdef MPI
    use mpi_siesta, only: MPI_Double_Precision => MPI_double_precision
    use mpi_siesta, only: MPI_logical, MPI_comm_world, MPI_Bcast
    use mpi_siesta, only: MPI_integer
#endif
    use m_ts_electype
    real(dp) , parameter :: EPS = 1d-7
    
! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in)  :: funit ! unit of gf-file
    type(Elec), intent(in) :: El
    integer, intent(in)  :: c_nkpar,c_NEn
    logical, intent(in)  :: c_RemUCell

! ***********************
! * LOCAL variables     *
! ***********************
    character(len=FILE_LEN) :: curGFfile  ! Name of the GF file
    character(len=NAME_LEN) :: curGFtitle ! title of the GF file

    integer :: nspin,nkpar,na,no,NA1,NA2,NA3,NEn
    real(dp) :: mu ! The Fermi energy shift due to a voltage
    real(dp), allocatable :: xa(:,:)
    integer, allocatable :: lasto(:)
    real(dp) :: ucell(3,3)
    logical :: errorGf , RemUCell
#ifdef MPI
    integer :: MPIerror
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE read_Green' )
#endif

    errorGF = .false.

    
    io_read: if ( IONode ) then

       ! Retrieve name of file currently reading
       inquire(unit=funit,name=curGFfile)

       read(funit) curGFtitle
       ! read electrode information
       read(funit) nspin, ucell
       read(funit) na, no ! used atoms and used orbitals
       read(funit) ! xa, lasto
       read(funit) NA1,NA2,NA3
       read(funit) mu
       ! read contour information
       read(funit) RemUCell
       read(funit) nkpar
       read(funit) ! kpoints, kweight
       read(funit) NEn
       read(funit)! ce

       ! Check unit cell distances..
       if ( RemUCell .neqv. c_RemUCell ) then
          write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
          if ( RemUCell ) then
             write(*,*)"The GF file has no inner unit cell distances. You have requested &
               &that they are preserved!"
          else
             write(*,*)"The GF file has inner unit cell distances. You have requested &
               &that they are not preserved!"
          end if
          write(*,'(2(a,l2))')"Found: ",RemUCell,", expected: ",c_RemUCell
          errorGF = .true.
       end if

       ! Check Fermi shift
       if ( dabs(El%mu%mu-mu) > EPS ) then
          write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
          write(*,*)"The chemical shift in the electrode does not match the &
               &required shift!"
          write(*,'(2(a,f12.6))')"Found: ",mu,", expected: ",El%mu%mu
          errorGF = .true.
       end if

       ! Check # of energy points
       if (NEn .ne. c_NEn) then
          write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: NEn=',NEn,' expected:', c_NEn
          errorGF = .true.
       end if

       ! Check # of atoms
       if (na .ne. UsedAtoms(El)) then
          write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: na_u=',na,' expected:', UsedAtoms(El)
          errorGF = .true.
       end if

       ! Check # of q-points
       if (Rep(El) .ne. NA1*NA2*NA3) then
          write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: unexpected no. q-points'
          errorGF = .true.
       end if

       ! Check # of k-points
       if ( nkpar .ne. c_nkpar ) then
          write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
          write(*,*) 'read_Green: Unexpected number of kxy-points'
          write(*,*) 'read_Green: ERROR: nkpt=',nkpar,' expected:', c_nkpar
          errorGF = .true.
       end if

       ! Check # of spin
       if (nspin .ne. Spin(El) ) then
          write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: nspin=',nspin,' expected:', Spin(El)
          errorGF = .true.
       end if

       ! Check # of orbitals
       if (no .ne. UsedOrbs(El)) then
          write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: no=',no,' expected:', UsedOrbs(El)
          errorGF = .true.
       end if

    end if io_read

    if ( errorGF ) then
       call die("Error in reading GFfile: "//trim(curGFfile))
    end if

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS read_Green' )
#endif

  end subroutine read_Green



! ##################################################################
! ##    Check header of Greens function file against settings     ##
! ##                            By                                ##
! ##      Nick Papior Andersen, nickpapior@gmail.com              ##
! ##                                                              ##
! ## Checks information an returns number of atoms and orbitals   ##
! ##################################################################
  subroutine check_Green(funit,El, &
       c_ucell,c_nkpar,c_kpar,c_wkpar, &
       c_NEn,c_ce, &
       c_RemUCell, xa_Eps, errorGF)

    use precision, only: dp
    use units,     only: Ang
    use m_ts_cctype
    use m_ts_electype

    real(dp) , parameter :: EPS = 1d-7

! ***********************
! * INPUT variables     *
! ***********************
! file for reading, Green's function file
    integer, intent(in)        :: funit
    type(Elec), intent(in)     :: El
    real(dp), intent(in)       :: c_ucell(3,3) ! Unit cell of the CONTACT
    ! k-point information
    integer, intent(in)        :: c_nkpar
    real(dp), intent(in)       :: c_kpar(3,c_nkpar) , c_wkpar(c_nkpar)
! Energy point on the contour used 
    integer, intent(in)        :: c_NEn
    complex(dp), intent(in)    :: c_ce(c_NEn)
    logical, intent(in)        :: c_RemUCell ! Should the Green's function file have the inner cell distances or not?
    real(dp), intent(in)       :: xa_Eps
! ***********************
! * OUTPUT variables    *
! ***********************
! Return whether it is a correct Green's function file
    logical, intent(out)       :: errorGF

! ***********************
! * LOCAL variables     *
! ***********************
    character(NAME_LEN) :: curGFtitle ! Title, currently not used
    real(dp) :: mu ! The energy shift in the Fermi energy
    integer :: nspin, na, no, nkpar ! spin, # of atoms, # of orbs, # k-points
    integer :: NA1,NA2,NA3 ! # repetitions in x, # repetitions in y
    real(dp), allocatable :: xa(:,:) ! electrode atomic coordinates
    integer, allocatable :: lasto(:) ! the electrode orbitals of the atoms
    real(dp), allocatable :: kpar(:,:) ! k-points
    real(dp), allocatable :: wkpar(:) ! k-point weights

    integer :: NEn ! # energy points on the contour
    complex(dp), allocatable :: ce(:)

! Helpers..
    character(200) :: curGFfile
    real(dp) :: ucell(3,3)
    integer :: iEn
    integer :: i,j,k,iq, iaa, ia
    real(dp) :: c_xa_o(3), xa_o(3)
    real(dp) :: wqbtmp,qbtmp(3), ktmp(3), kpt(3)
    logical :: localErrorGf, eXa, RemUCell

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE check_Green' )
#endif

    ! Initialize it to be an error unless the entire routine is runned through.
    ! Upon normal exit it will be changed to .FALSE.
    localErrorGf = errorGF
    errorGF = .true.
    
    ! Retrieve name of file currently reading
    inquire(unit=funit,name=curGFfile)

    ! Read in header of the file
    read(funit) curGFtitle

    ! Read in electrode information
    read(funit) nspin, ucell
    read(funit) na,no
    allocate(xa(3,na),lasto(na))
    read(funit) xa,lasto
    read(funit) NA1,NA2,NA3
    read(funit) mu

    if ( Spin(El) /= nspin ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Number of spin is wrong!"
       write(*,'(2(a,i2))') "Found: ",nspin,", expected: ",Spin(El)
       localErrorGf = .true.
    end if
    if ( any(abs(unitcell(El)-ucell) > EPS) ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Unit-cell is not consistent!"
       write(*,*) "Found (Ang):"
       write(*,'(3(3(tr1,f10.5),/))') ucell/Ang
       write(*,*) "Expected (Ang):"
       write(*,'(3(3(tr1,f10.5),/))') unitcell(El)/Ang
       localErrorGf = .true.
    end if
    if ( UsedAtoms(El) /= na ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Number of atoms is wrong!"
       write(*,'(2(a,i2))') "Found: ",na,", expected: ",UsedAtoms(El)
       localErrorGf = .true.
    end if
    if ( UsedOrbs(El) /= no ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Number of orbitals is wrong!"
       write(*,'(2(a,i2))') "Found: ",no,", expected: ",UsedOrbs(El)
       localErrorGf = .true.
    end if
    
    ! Initialize error parameter
    eXa = .false.
    ! We check elsewhere that the electrode is consistent with
    ! the FDF input
    do ia = 1 , min(na,UsedAtoms(El)) ! in case it is completely wrong
       do i = 1 , 3
          eXa= eXa .or. &
               abs(xa(i,ia)-El%xa_used(i,ia)) > xa_Eps
       end do
    end do
    if ( eXa ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Atomic coordinates are wrong:"
       write(*,'(1x,a,t35,a)') &
            "Structure of GF electrode","| Electrode:"
       write(*,'(t3,3a10,''  |'',3a10)') &
            "X (Ang)","Y (Ang)","Z (Ang)", &
            "X (Ang)","Y (Ang)","Z (Ang)"
       do ia = 1, na
          write(*,'(t3,3f10.5,''  |'',3f10.5)') &
               xa(:,ia)/Ang, El%xa_used(:,ia)/Ang
       end do
       localErrorGf = .true.
    end if
    deallocate(xa)

    if ( any(lasto - El%lasto_used /= 0) ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Number of orbitals on used atoms is wrong!"
       write(*,'(a,1000i3)') "Found    lasto: ",lasto
       write(*,'(a,1000i3)') "Expected lasto: ",El%lasto_used
       localErrorGf = .true.
    end if
    deallocate(lasto)

    if ( RepA1(El)/=NA1 .or. RepA2(El)/=NA2 .or. RepA3(El)/=NA3 ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Number of repetitions is wrong!"
       write(*,'(2(a,i3))') "Found NA1: ",NA1,", expected NA1: ",RepA1(El)
       write(*,'(2(a,i3))') "Found NA2: ",NA2,", expected NA2: ",RepA2(El)
       write(*,'(2(a,i3))') "Found NA3: ",NA3,", expected NA3: ",RepA3(El)
       localErrorGf = .true.
    end if
    if ( abs(El%mu%mu-mu) > EPS ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"The chemical shift in the electrode does not match the &
            &required shift!"
       write(*,'(2(a,f12.6))')"Found: ",mu,", expected: ",El%mu%mu
       localErrorGf = .true.
    end if


    ! Read in general information about the context
    read(funit) RemUCell
    read(funit) nkpar
    allocate(kpar(3,nkpar),wkpar(nkpar))
    read(funit) kpar,wkpar

    if ( RemUCell .neqv. c_RemUCell ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       if ( RemUCell ) then
          write(*,*)"The GF file has no inner unit cell distances. You have requested &
               &that they are preserved!"
          write(*,'(2(a,l2))')"Found: ",RemUCell,", expected: ",c_RemUCell
       end if
       localErrorGf = .true.
    end if
    if ( c_nkpar /= nkpar ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Number of k-points is wrong!"
       write(*,'(2(a,i4))') "Found: ",nkpar,", expected: ",c_nkpar
       localErrorGf = .true.
    end if
    
    ! Check k-points
    do i = 1 , min(c_nkpar,nkpar)
       ! As the k-points are in the Electrode unit cell
       ! we need to compare that with those of the CONTACT cell!
       ! The advantage of this is that the GF files can be re-used for
       ! the same system with different lengths between the electrode layers.
       call kpoint_convert(ucell,kpar(:,i),ktmp,1)
       ktmp(1) = ktmp(1) * real(NA1,dp)
       ktmp(2) = ktmp(2) * real(NA2,dp)
       ktmp(3) = ktmp(3) * real(NA3,dp)
       call kpoint_convert(c_ucell,ktmp,kpt,-1)
       if ( dabs(c_kpar(1,i)-kpt(1)) > EPS .or. &
            dabs(c_kpar(2,i)-kpt(2)) > EPS .or. &
            dabs(c_kpar(3,i)-kpt(3)) > EPS ) then
          write(*,*)"k-points are not the same:"
          do j = 1 , min(c_nkpar,nkpar)
             call kpoint_convert(ucell,kpar(:,i),ktmp,1)
             ktmp(1) = ktmp(1) * real(NA1,dp)
             ktmp(2) = ktmp(2) * real(NA2,dp)
             ktmp(3) = ktmp(3) * real(NA3,dp)
             call kpoint_convert(c_ucell,ktmp,kpt,-1)
             write(*,'(3f12.5,a,3f12.5)') c_kpar(:,j),'  :  ',kpt(:)
          end do
          localErrorGf = .true.
          exit
       end if
       if ( dabs(c_wkpar(i)-wkpar(i)) > EPS ) then
          write(*,*)"k-point weights are not the same:"
          do j = 1 , c_nkpar
             write(*,'(f12.5,a,f12.5)') c_wkpar(j),'  :  ',wkpar(j)
          end do
          localErrorGf = .true.
          exit
       end if
    end do
    deallocate(kpar,wkpar)

    ! Read in information about the contour
    read(funit) NEn
    allocate(ce(NEn))
    read(funit) ce

    ! Check energy points
    if ( c_NEn /= NEn ) then
       write(*,*)"ERROR: Green's function file: "//trim(curGFfile)
       write(*,*)"Number of energy points is not as expected!"
       write(*,'(2(a,i4))') "Found: ",NEn,", expected: ",c_NEn
       localErrorGf = .true.
    end if
    do iEn = 1 , NEn
       if ( cdabs(ce(iEn)-c_ce(iEn)) > EPS ) then
          write(*,*) ' Warning: contours differ by >', EPS
       end if
       if ( cdabs(ce(iEn)-c_ce(iEn)) > 10.d0*EPS ) then
          write(*,*) ' ERROR  : contours differ by >', 10.d0*EPS
          localErrorGf = .true.
       end if
    end do
    deallocate(ce)

    errorGF = localErrorGf
    
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS check_Green' )
#endif
    
  end subroutine check_Green

end module m_ts_GF
