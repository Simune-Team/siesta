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
!          1) read_Green
!          2) check_Green


  implicit none

  public :: read_Green, check_Green

  private

contains

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
  subroutine read_Green(funit,c_EfShift,c_nkpar,c_NEn,c_nua,c_NA1,c_NA2, &
       c_no,c_nspin, &
       nkpar,kpar,wkpar,nq,wq,qb)
    
    use precision, only : dp
    use parallel,  only : IONode
    use sys ,      only : die
#ifdef MPI
    use mpi_siesta
#endif
    real(dp) , parameter :: EPS = 1d-7
    
! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in)  :: funit ! unit of gf-file
    real(dp), intent(in) :: c_EfShift ! The fermi shift in the electrode
    integer, intent(in)  :: c_nkpar,c_NEn,c_nua,c_NA1,c_NA2,c_no,c_nspin

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
    complex(dp), dimension(:), allocatable :: contour,wgf
    real(dp) :: ucell(3,3)
    character(len=200) :: GFtitle
    logical :: errorGf 
#ifdef MPI
    integer :: MPIerror
#endif

#ifdef DEBUG
    call write_debug( 'PRE read_Green' )
#endif

    errorGF = .false.

    
    io_read: if ( IONode ) then

       ! Retrieve name of file currently reading
       inquire(unit=funit,name=GFfile)

       read(funit) GFtitle
       read(funit) EfShift,NEn
       read(funit) nua,NA1,NA2,nkpar,nq
       read(funit) nspin,ucell

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

#ifdef DEBUG
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
  subroutine check_Green(funit,c_EfShift,c_ucell,c_nspin,c_nkpar,c_kpar,c_wkpar, &
       c_NEn,c_contour,c_wgf, &
       c_NA1,c_NA2,errorGF, &
       nua,no)

    use precision, only: dp

    real(dp) , parameter :: EPS = 1d-7

! ***********************
! * INPUT variables     *
! ***********************
! file for reading, Green's function file
    integer, intent(in)        :: funit
    real(dp), intent(in)       :: c_EfShift ! The Fermi energy shift of the electrode
    real(dp), intent(in)       :: c_ucell(3,3) ! Unit cell of the CONTACT
! spin of system
    integer, intent(in)        :: c_nspin
! k-point information
    integer, intent(in)        :: c_nkpar
    real(dp), intent(in)       :: c_kpar(3,c_nkpar) , c_wkpar(c_nkpar)
! Energy point on the contour used 
    integer, intent(in)        :: c_NEn
    complex(dp), intent(in)    :: c_contour(c_NEn),c_wgf(c_NEn)
! We cannot check for number of atoms in the unit cell.
! TODO Add this so that it is possible.
! Repetition information
    integer, intent(in)        :: c_NA1,c_NA2
! ***********************
! * OUTPUT variables    *
! ***********************
! Return whether it is a correct Green's function file
    logical, intent(out)       :: errorGF
! Return number of atoms in the Green's function file
    integer, intent(out)       :: nua
! Return the number of orbitals used in the unit cell of this Green's function
    integer, intent(out)       :: no

! ***********************
! * LOCAL variables     *
! ***********************
    character(200) :: GFtitle ! Title, currently not used
    real(dp) :: EfShift ! The energy shift in the Fermi energy

    integer :: NA1,NA2 ! # repetitions in x, # repetitions in y
    integer :: nkpar,nspin ! # k-points, # spin
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
    integer :: i,j,iq
    real(dp) :: wqbtmp,qbtmp(3), ktmp(3), kpt(3)
    logical :: localErrorGf

#ifdef DEBUG
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
       return
    end if
! Check energy points
    if ( c_NEn /= NEn ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of energy points is not as expected!"
       write(*,'(2(a,i4))') "Found: ",NEn,", expected: ",c_NEn
       return
    end if

    ! Read in integers (also the returned number of atoms)
    read(funit) nua,NA1,NA2,nkpar,nqb
    if ( c_NA1 /= NA1 .or. c_NA2 /= NA2 .or. c_NA1*c_NA2 /= nqb ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of repetitions is wrong!"
       write(*,'(2(a,i3))') "Found NA1: ",NA1,", expected NA1: ",c_NA1
       write(*,'(2(a,i3))') "Found NA2: ",NA2,", expected NA2: ",c_NA2
       return
    end if
    if ( c_nkpar /= nkpar ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of k-points is wrong!"
       write(*,'(2(a,i4))') "Found: ",nkpar,", expected: ",c_nkpar
       return
    end if

    read(funit) nspin,ucell
    if ( c_nspin /= nspin ) then
       write(*,*)"ERROR: Green's function file: "//TRIM(GFfile)
       write(*,*)"Number of spin is wrong!"
       write(*,'(2(a,i2))') "Found: ",nspin,", expected: ",c_nspin
       return
    end if

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
       if ( cdabs(contour(iEn)-c_contour(iEn)) > EPS ) then
          write(*,*) ' Warning: contours differ by >', EPS
       end if
       if ( cdabs(contour(iEn)-c_contour(iEn)) > 10.d0*EPS ) then
          write(*,*) ' ERROR  : contours differ by >', 10.d0*EPS
          return
       end if
       if (cdabs(wgf(iEn)-c_wgf(iEn)) > EPS ) then 
          write(*,*) ' ERROR: contour weights differ by >',EPS
          return
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
          return
       end if
       if ( dabs(c_wkpar(i)-wkpar(i)) > EPS ) then
          write(*,*)"k-point weights are not the same:"
          do j = 1 , c_nkpar
             write(*,'(f12.5,a,f12.5)') c_wkpar(j),'  :  ',wkpar(j)
          end do
          return
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
             return
          end if
          if ( dabs(qbtmp(1)-qb(1,iq)) > EPS .or. &
               dabs(qbtmp(2)-qb(2,iq)) > EPS ) then
             write(*,*)"Expansion q-point weights are not the same:"
             write(*,'(f12.5,a,f12.5)') wqbtmp,'  :  ',wqb(iq)
             return
          end if
       end do
    end do
    call memory('D','D',nqb*4,'check_GF')
    deallocate(qb,wqb)

! TODO, need to check this in some fashion
    read(funit) no

    errorGF = localErrorGf

#ifdef DEBUG
    call write_debug( 'POS check_Green' )
#endif

  end subroutine check_Green


end module m_ts_GF
