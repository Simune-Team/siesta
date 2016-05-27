!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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

  subroutine do_Green(tElec, HSFile, GFFile, GFTitle, &
       ElecValenceBandBot, optReUseGF, &
       nkpnt,kpoint,kweight, &
       NBufAt,NUsedAtoms,NA1,NA2, RemUCellDistance, &
       ucell,xa,nua,NEn,contour,chem_shift,ZBulkDOS,nspin)
    
    use precision,  only : dp
    use parallel  , only : IONode
    use sys ,       only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Bcast, MPI_Integer, MPI_Logical
#endif
    use m_ts_cctype
    use m_ts_electrode, only : create_Green

    implicit none
    
! ***********************
! * INPUT variables     *
! ***********************
    character(len=1), intent(in) :: tElec   ! 'L' for Left electrode, 'R' for right
    character(len=*), intent(in) :: HSFile  ! The electrode TSHS file
    character(len=*), intent(in) :: GFFile  ! The electrode GF file to be saved to
    character(len=*), intent(in) :: GFTitle ! The title to be written in the GF file
    logical, intent(in)          :: ElecValenceBandBot ! Whether or not to calculate electrodes valence bandbottom
    logical, intent(in)          :: optReUseGF ! Should we re-use the GF files if they exists?    
    integer, intent(in)          :: nkpnt ! Number of k-points
    real(dp),dimension(3,nkpnt),intent(in) :: kpoint ! k-points
    real(dp),dimension(nkpnt),intent(in) :: kweight ! weights of kpoints
    integer, intent(in)            :: NBufAt,NA1,NA2 ! Buffer/Rep a1/Rep a2
    logical, intent(in)            :: RemUCellDistance ! Whether to remove the unit cell distance in the Hamiltonian.
    integer, intent(in)            :: NUsedAtoms ! Needs update here
    integer, intent(in)            :: nua ! Full system count of atoms in unit cell
    real(dp), dimension(3,3)       :: ucell ! The unit cell of the CONTACT
    real(dp), intent(in)           :: xa(3,nua) ! Coordinates in the system for the TranSIESTA routine
    integer, intent(in)            :: nspin ! spin in system
    integer, intent(in)            :: NEn ! Number of energy points
    type(ts_ccontour), intent(in)  :: contour(NEn) ! contour for GF
    real(dp), intent(in)           :: chem_shift ! the Fermi-energy we REQUIRE the electrode

! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), intent(out)       :: ZBulkDOS(NEn,nspin) ! DOS at energy points

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: uGF
    logical :: errorGF, exist, ReUseGF
#ifdef MPI
    integer :: MPIerror
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE do_Green' )
#endif
    
! check the file for existance
    inquire(file=GFfile,exist=exist)
    
    ReUseGF = optReUseGF
! If it does not find the file, calculate the GF
    if ( exist ) then
       if (IONode ) then
          write(*,*) "Electrode Green's function file: '"//&
               trim(GFFile)//"' already exist."
          if ( .not. ReUseGF ) then
             write(*,*)"Green's function file '"//&
                  trim(GFFile)//"' is requested overwritten."
          end if
       end if
    else
       ReUseGF = .false.
    end if

    ! We return if we should not calculate it
    if ( .not. ReUseGF ) then
       ! Create the GF file
       call create_Green(tElec,HSFile, GFFile, GFTitle, &
            ElecValenceBandBot, &
            nkpnt,kpoint,kweight, &
            NBufAt,NUsedAtoms,NA1,NA2, RemUCellDistance, &
            ucell,xa,nua,NEn,contour,chem_shift,ZBulkDOS,nspin)
    end if

!
!     Check that the Green's functions are correct!
!     This is needed as create_Green returns if user requests not to
!     overwrite an already existing file.
!     This check will read in the number of orbitals and atoms in the
!     electrode surface Green's function.
!
    errorGF = .false.
! Check the GF file
    if(IONode) then
       call io_assign(uGF)
       open(file=GFFile,unit=uGF,form='UNFORMATTED')

       if ( tElec == 'L' ) then
          call check_Green(uGF,chem_shift,ucell, &
               NUsedAtoms*NA1*NA2,xa(1,NBufAt+1), &
               nspin,nkpnt,kpoint, &
               kweight,NEn,contour,NA1,NA2,RemUCellDistance,errorGF)
       else if ( tElec == 'R' ) then
          call check_Green(uGF,chem_shift,ucell, &
               NUsedAtoms*NA1*NA2,xa(1,nua-NBufAt-NUsedAtoms*NA1*NA2+1), &
               nspin,nkpnt,kpoint, &
               kweight,NEn,contour,NA1,NA2,RemUCellDistance,errorGF)
       end if
       call io_close(uGF)
    endif

!     Check the error in the GF file
#ifdef MPI
    call MPI_Bcast(errorGF,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
#endif
    if ( errorGF ) &
         call die("Error in GFfile: "//trim(GFFile)//". Please move or delete")

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
  subroutine read_Green(funit,print_title,c_EfShift,c_nkpar,c_NEn,c_nua,c_NA1,c_NA2, c_RemUCell, &
       c_no,c_nspin, &
       nkpar,kpar,wkpar,nq,wq,qb)
    
    use precision, only : dp
    use parallel,  only : IONode
    use sys ,      only : die
#ifdef MPI
    use mpi_siesta, only: DAT_double => MPI_double_precision
    use mpi_siesta, only: MPI_logical, MPI_comm_world, MPI_Bcast
    use mpi_siesta, only: MPI_integer
#endif
    real(dp) , parameter :: EPS = 1d-7
    
! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in)  :: funit ! unit of gf-file
    logical, intent(in)  :: print_title
    real(dp), intent(in) :: c_EfShift ! The fermi shift in the electrode
    integer, intent(in)  :: c_nkpar,c_NEn,c_nua,c_NA1,c_NA2,c_no,c_nspin
    logical, intent(in)  :: c_RemUCell

! ***********************
! * OUTPUT variables    *
! ***********************    
    integer, intent(out) :: nkpar,nq
    real(dp),allocatable, intent(inout) :: kpar(:,:), wkpar(:)
    real(dp),allocatable, intent(inout) :: qb(:,:), wq(:)

! ***********************
! * LOCAL variables     *
! ***********************
    character(200) :: GFfile ! Name of the GF file
    integer :: NEn,nua,NA1,NA2,no,nspin
    real(dp) :: EfShift            ! The Fermi energy shift due to a voltage
    real(dp), dimension(:,:), allocatable :: xa
    complex(dp), dimension(:), allocatable :: contour,wgf
    real(dp) :: ucell(3,3)
    character(len=200) :: GFtitle
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
       inquire(unit=funit,name=GFfile)

       read(funit) GFtitle
       if ( print_title ) then
          write(*,'(1x,a)')   "Reading GF file, with title:"
          write(*,'(1x,a)')   "  "//trim(GFfile)
          write(*,'(1x,a,/)') "Title: '"//trim(GFtitle)//"'"
       end if
       read(funit) EfShift,NEn
       read(funit) RemUCell
       read(funit) nua,NA1,NA2,nkpar,nq
       read(funit) nspin,ucell
       allocate(xa(3,nua))
       read(funit) xa
       deallocate(xa) ! We expect check_Green to catch this...

! Check unit cell distances..
       if ( RemUCell .neqv. c_RemUCell ) then
          write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
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
       if ( dabs(c_EfShift-EfShift) > EPS ) then
          write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
          write(*,*)"The chemical shift in the electrode does not match the &
               &required shift!"
          write(*,'(2(a,f12.6))')"Found: ",EfShift,", expected: ",c_EfShift
          errorGF = .true.
       end if

! Check # of energy points
       if (NEn .ne. c_NEn) then
          write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
          write(*,*) 'read_Green: ERROR: NEn=',NEn,' expected:', c_NEn
          errorGF = .true.
       end if

! Check # of atoms
       if (nua .ne. c_nua) then
          write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
          write(*,*) 'read_Green: ERROR: nua=',nua,' expected:', c_nua
          errorGF = .true.
       end if

! Check # of q-points
       if (c_NA1*c_NA2 .ne. NA1*NA2) then
          write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
          write(*,*) 'read_Green: ERROR: unexpected no. q-points'
          errorGF = .true.
       end if

! Check # of k-points
       if ( nkpar .ne. c_nkpar ) then
          write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
          write(*,*) 'read_Green: Unexpected number of kxy-points'
          write(*,*) 'read_Green: ERROR: nkpt=',nkpar,' expected:', c_nkpar
          errorGF = .true.
       end if

! Check # of spin
       if (nspin .ne. c_nspin) then
          write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
          write(*,*) 'read_Green: ERROR: nspin=',nspin,' expected:', c_nspin
          errorGF = .true.
       end if

       allocate(contour(NEn))
       allocate(wgf(NEn))
       call memory('A','Z',NEn*2,'read_green')
       allocate(kpar(3,nkpar))
       allocate(wkpar(nkpar))
       call memory('A','D',nkpar*4,'transiesta')
       allocate(qb(3,nq))
       allocate(wq(nq))
       call memory('A','D',nq*4,'transiesta')

       read(funit) contour,wgf,kpar,wkpar,qb,wq
       call memory('D','Z',size(contour)*2,'read_green')
       deallocate(contour)
       deallocate(wgf)

       read(funit) no

! Check # of orbitals
       if (no .ne. c_no) then
          write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
          write(*,*) 'read_Green: ERROR: no=',no,' expected:', c_no
          errorGF = .true.
       end if

    end if io_read

    ! Do broadcast etc

#ifdef MPI
    call MPI_Bcast(errorGF,1,MPI_LOGICAL,0,MPI_Comm_World,MPIerror)
#endif

    if ( errorGF ) then
       call die("Error in reading GFfile: "//trim(GFfile))
    end if

#ifdef MPI
    call MPI_Bcast(nkpar,1,MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nq,1,MPI_integer,0,MPI_Comm_World,MPIerror)
    
    if(.not. IONode) then
       allocate(kpar(3,nkpar))
       allocate(wkpar(nkpar))
       call memory('A','D',4*nkpar,'transiesta')
       allocate(qb(3,nq))
       allocate(wq(nq))
       call memory('A','D',4*nq,'transiesta')
    endif
    
    call MPI_Bcast(qb(1,1),3*nq,DAT_double,0,MPI_Comm_World ,MPIerror)
    call MPI_Bcast(wq,nq,DAT_double,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(kpar(1,1),3*nkpar,DAT_double,0,MPI_Comm_World ,MPIerror)
    call MPI_Bcast(wkpar,nkpar,DAT_double,0,MPI_Comm_World,MPIerror)
#endif

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
  subroutine check_Green(funit,c_EfShift,c_ucell,c_nua,c_xa,c_nspin,c_nkpar,c_kpar,c_wkpar, &
       c_NEn,c_contour, &
       c_NA1,c_NA2,c_RemUCell,errorGF)

    use precision, only: dp
    use units,     only: Ang
    use m_ts_cctype

    real(dp) , parameter :: EPS = 1d-7
    real(dp) , parameter :: EPS_xa = 1d-4

! ***********************
! * INPUT variables     *
! ***********************
! file for reading, Green's function file
    integer, intent(in)        :: funit
    real(dp), intent(in)       :: c_EfShift ! The Fermi energy shift of the electrode
    real(dp), intent(in)       :: c_ucell(3,3) ! Unit cell of the CONTACT
    integer, intent(in)        :: c_nua ! number of atoms in the electrode which match the expanded electrode
    real(dp), intent(in)       :: c_xa(3,c_nua) ! Atomic coordinates of the electrode in Transiesta
! spin of system
    integer, intent(in)        :: c_nspin
! k-point information
    integer, intent(in)        :: c_nkpar
    real(dp), intent(in)       :: c_kpar(3,c_nkpar) , c_wkpar(c_nkpar)
! Energy point on the contour used 
    integer, intent(in)        :: c_NEn
    type(ts_ccontour), intent(in) :: c_contour(c_NEn)
! We cannot check for number of atoms in the unit cell.
! TODO Add this so that it is possible.
! Repetition information
    integer, intent(in)        :: c_NA1,c_NA2
    logical, intent(in)        :: c_RemUCell ! Should the Green's function file have the inner cell distances or not?
! ***********************
! * OUTPUT variables    *
! ***********************
! Return whether it is a correct Green's function file
    logical, intent(out)       :: errorGF

! ***********************
! * LOCAL variables     *
! ***********************
    character(200) :: GFtitle ! Title, currently not used
    real(dp) :: EfShift ! The energy shift in the Fermi energy

    integer :: NA1,NA2 ! # repetitions in x, # repetitions in y
    integer :: nua,no,nkpar,nspin ! # of atoms, # of orbs, # k-points, # spin
    real(dp), dimension (:,:), allocatable :: xa ! electrode atomic coordinates
    real(dp), dimension (:,:), allocatable :: kpar ! k-points
    real(dp), dimension (:), allocatable :: wkpar ! k-point weights
    integer :: nqb
    real(dp), dimension (:,:), allocatable :: qb ! q-points
    real(dp), dimension (:), allocatable :: wqb ! q-point weights

    integer :: NEn ! # energy points on the contour
    complex(dp), dimension(:), allocatable :: contour,wGF ! Energies and weights

! Helpers..
    character(200) :: GFfile
    real(dp) :: ucell(3,3)
    integer :: iEn
    integer :: i,j,iq, iaa, ia
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
    inquire(unit=funit,name=GFfile)

! Read in header of the file
    read(funit) GFtitle
    read(funit) EfShift,NEn
! Check Fermi shift
    if ( dabs(c_EfShift-EfShift) > EPS ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"The chemical shift in the electrode does not match the &
            &required shift!"
       write(*,'(2(a,f12.6))')"Found: ",EfShift,", expected: ",c_EfShift
       localErrorGf = .true.
    end if
! Check energy points
    if ( c_NEn /= NEn ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of energy points is not as expected!"
       write(*,'(2(a,i4))') "Found: ",NEn,", expected: ",c_NEn
       localErrorGf = .true.
    end if
    read(funit) RemUCell
! Check unit cell distances..
    if ( RemUCell .neqv. c_RemUCell ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       if ( RemUCell ) then
          write(*,*)"The GF file has no inner unit cell distances. You have requested &
               &that they are preserved!"
          write(*,'(2(a,l2))')"Found: ",RemUCell,", expected: ",c_RemUCell
          localErrorGf = .true.
       end if
    end if

    ! Read in integers (also the returned number of atoms)
    read(funit) nua,NA1,NA2,nkpar,nqb
    if ( c_NA1 /= NA1 .or. c_NA2 /= NA2 .or. c_NA1*c_NA2 /= nqb ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of repetitions is wrong!"
       write(*,'(2(a,i3))') "Found NA1: ",NA1,", expected NA1: ",c_NA1
       write(*,'(2(a,i3))') "Found NA2: ",NA2,", expected NA2: ",c_NA2
       localErrorGf = .true.
    end if
    if ( c_nua /= nua*NA1*NA2 ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of atoms is wrong!"
       write(*,'(2(a,i2))') "Found: ",nua,", expected: ",c_nua/NA1/NA2
       localErrorGf = .true.
    end if
    if ( c_nkpar /= nkpar ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of k-points is wrong!"
       write(*,'(2(a,i4))') "Found: ",nkpar,", expected: ",c_nkpar
       localErrorGf = .true.
    end if

    read(funit) nspin,ucell
    if ( c_nspin /= nspin ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of spin is wrong!"
       write(*,'(2(a,i2))') "Found: ",nspin,", expected: ",c_nspin
       localErrorGf = .true.
    end if

    ! Read in electrode coordinates
    ! We know that NA[12] == c_NA[12] 
    allocate(xa(3,nua))
    read(funit) xa 
    ! Check electrode coordinates
    ! Save origo of System electrode
    c_xa_o(:) = c_xa(:,1)
    ! Save origo of electrode 
    xa_o(:) = xa(:,1)

    ! Initialize error parameter
    eXa = .false.
    iaa = 1
    do ia = 1 , nua
       do j=0,NA2-1
          do i=0,NA1-1
             eXa=eXa.or.abs(xa(1,ia)-xa_o(1) + &
                  ucell(1,1)*i+ucell(1,2)*j - &
                  c_xa(1,iaa)+c_xa_o(1)) > EPS_xa
             eXa=eXa.or.abs(xa(2,ia)-xa_o(2) + &
                  ucell(2,1)*i+ucell(2,2)*j - &
                  c_xa(2,iaa)+c_xa_o(2)) > EPS_xa
             eXa=eXa.or.abs(xa(3,ia)-xa_o(3) - &
                  c_xa(3,iaa)+c_xa_o(3)) > EPS_xa
             iaa = iaa + 1
          end do
       end do
    end do
    if ( eXa ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Atomic coordinates are wrong:"
       write(*,'(1x,a,t35,a)') &
            "Structure of GF electrode","| System electrode:"
       write(*,'(t3,3a10,''  |'',3a10)') &
            "X (Ang)","Y (Ang)","Z (Ang)", &
            "X (Ang)","Y (Ang)","Z (Ang)"
       iaa = 1
       do ia = 1, nua
          do j=0,NA2-1
             do i=0,NA1-1
                write(*,'(t3,3f10.5,''  |'',3f10.5)') &
                     (xa(1,ia)-xa_o(1)+ucell(1,1)*i+ucell(1,2)*j)/Ang, &
                     (xa(2,ia)-xa_o(2)+ucell(2,1)*i+ucell(2,2)*j)/Ang, &
                     (xa(3,ia)-xa_o(3))/Ang, &
                     (c_xa(1,iaa)-c_xa_o(1))/Ang, &
                     (c_xa(2,iaa)-c_xa_o(2))/Ang, &
                     (c_xa(3,iaa)-c_xa_o(3))/Ang
                iaa = iaa + 1
             end do
          end do
       end do
       localErrorGf = .true.
    end if
    deallocate(xa)

! Read in and check contour, k-point, q-point information
! First allocate all arrays

! Energy contour
    allocate(contour(NEn))
    allocate(wGF(NEn))
    call memory('A','Z',NEn*2,'check_GF')

! k-points
    allocate(kpar(3,nkpar))
    allocate(wkpar(nkpar))
    call memory('A','D',nkpar*4,'check_GF')

! q-points used for expanding (b for units of reciprocal vectors)
    allocate(qb(3,nqb))
    allocate(wqb(nqb))
    call memory('A','D',nqb*4,'check_GF')

    read(funit) contour,wgf,kpar,wkpar,qb,wqb

! Check contours
    do iEn=1,NEn
       if ( cdabs(contour(iEn)-c_contour(iEn)%c) > EPS ) then
          write(*,*) ' Warning: contours differ by >', EPS
       end if
       if ( cdabs(contour(iEn)-c_contour(iEn)%c) > 10.d0*EPS ) then
          write(*,*) ' ERROR  : contours differ by >', 10.d0*EPS
          localErrorGf = .true.
       end if
       if (cdabs(wgf(iEn)-c_contour(iEn)%w) > EPS ) then 
          write(*,*) ' ERROR: contour weights differ by >',EPS
          localErrorGf = .true.
       end if
    end do
    call memory('D','Z',NEn*2,'check_GF')
    deallocate(contour)
    deallocate(wgf)

! Check k-points
    do i = 1 , c_nkpar
       ! As the k-points are in the Electrode unit cell
       ! we need to compare that with those of the CONTACT cell!
       ! The advantage of this is that the GF files can be re-used for
       ! the same system with different lengths between the electrode layers.
       call kpoint_convert(ucell,kpar(:,i),ktmp,1)
       ktmp(1) = ktmp(1) * real(NA1,dp)
       ktmp(2) = ktmp(2) * real(NA2,dp)
       call kpoint_convert(c_ucell,ktmp,kpt,-1)
       if ( dabs(c_kpar(1,i)-kpt(1)) > EPS .or. &
            dabs(c_kpar(2,i)-kpt(2)) > EPS .or. &
            dabs(c_kpar(3,i)-kpt(3)) > EPS ) then
          write(*,*)"k-points are not the same:"
          do j = 1 , c_nkpar
             call kpoint_convert(ucell,kpar(:,i),ktmp,1)
             ktmp(1) = ktmp(1) * real(NA1,dp)
             ktmp(2) = ktmp(2) * real(NA2,dp)
             call kpoint_convert(c_ucell,ktmp,kpt,-1)
             write(*,'(3f12.5,a,3f12.5)') c_kpar(:,j),'  :  ',kpt(:)
          end do
          localErrorGf = .true.
       end if
       if ( dabs(c_wkpar(i)-wkpar(i)) > EPS ) then
          write(*,*)"k-point weights are not the same:"
          do j = 1 , c_nkpar
             write(*,'(f12.5,a,f12.5)') c_wkpar(j),'  :  ',wkpar(j)
          end do
          localErrorGf = .true.
       end if
    end do
    call memory('D','D',nkpar*4,'check_GF')
    deallocate(kpar,wkpar)

! Check q-points, have already checked size:
    iq = 0
    qbtmp(3) = 0._dp
    do i=1,NA1
       do j=1,NA2
          iq = iq+1
          qbtmp(1)=(1.0_dp*(i-1))/real(NA1,dp)
          qbtmp(2)=(1.0_dp*(j-1))/real(NA2,dp)
          wqbtmp  = 1.0_dp       /real(NA1*NA2,dp)
          if ( dabs(qbtmp(1)-qb(1,iq)) > EPS .or. &
               dabs(qbtmp(2)-qb(2,iq)) > EPS ) then
             write(*,*)"Expansion q-points are not the same:"
             write(*,'(3f12.5,a,3f12.5)') qbtmp(:),'  :  ',qb(:,iq)
             localErrorGf = .true.
          end if
          if ( dabs(qbtmp(1)-qb(1,iq)) > EPS .or. &
               dabs(qbtmp(2)-qb(2,iq)) > EPS ) then
             write(*,*)"Expansion q-point weights are not the same:"
             write(*,'(f12.5,a,f12.5)') wqbtmp,'  :  ',wqb(iq)
             localErrorGf = .true.
          end if
       end do
    end do
    call memory('D','D',nqb*4,'check_GF')
    deallocate(qb,wqb)

! TODO, need to check this in some fashion
    read(funit) no

    errorGF = localErrorGf

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS check_Green' )
#endif

  end subroutine check_Green


end module m_ts_GF
