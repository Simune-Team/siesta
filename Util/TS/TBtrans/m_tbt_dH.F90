! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---


! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication

! Enable the feature of manipulating with the
! Hamiltonian.

! This is easily done by letting the user
! create a NetCDF file which contain dH elements
! which will enter the NEGF equations through:
!   H = H + dH
! We note that the user can input several different
!  dH 
! corresponding to different situations.
! If dH is complex, one can use dH as an onsite
! self-energy change, thus introducing spurious
! effects.

! There are 4 levels of usage:
!   1. constant dH
!      This will enter in a common equation
!   2. k-dependent dH
!      This will only enter for the corresponding
!      k-point.
!         H_k = H_k + dH_k' with dH_k' = 0 for k/=k'
!   3. Energy
!      This will as 2. only enter for specific
!      energy points.
!         H(E) = H(E) + dH(E') with dH(E') = 0 for E/=E'
!   4. k and E dependent dH
!         H_k(E) = H_k(E) + dH_k'(E') with dH_k'(E') = 0 for E/=E' and k/=k'

! Note that the highest level has precedence above the others,
! so specifying both 1. and 4. will only be the equivalent of
! using level 4.

! This is particularly handy for TB calculations
! but can prove just as useful for regular DFT
! calculations as one can use an already existing
! SIESTA.nc file, and then do several "case" studies
! on such a file, thus actually having FULL control
! over EVERYTHING.

module m_tbt_dH

  use precision, only : dp

  use class_zSpData1D
  use m_tbt_save, only : tNodeE

  implicit none
 
  private

  type :: tTBTdH
     ! Designator of the current level of the 
     ! quantities in the dH designator
     integer :: lvl = -1

     !> The spin component that should be read in
     integer :: ispin = 1
     
     ! As this is a reciprocal cell, k-point,
     ! we initialize it immediately
     ! This should, for level 2 calculations 
     ! speed things up as we do not need to read the 
     ! sparse matrices always.
     real(dp) :: bkpt(3) = 2.12345_dp
     ! To ease the handling of different elements
     ! we only keep the complex quantity
     ! YES, we could limit this to only the real
     ! part for purely real designators.
     ! However, I do not suspect that users
     ! will use this to change extreme amounts
     ! of elements. (200.000 elements take up 3 MB)
     type(zSpData1D) :: dH
  end type tTBTdH

  ! The dH global variable
  type(tTBTdH), save :: dH

#ifdef NCDF_4

  ! File name that contains all information.
  character(len=250), save :: fname_dH
  
  ! Logical to determine whether the
  ! file should be read independently of the nodes
  ! or whether, there should be an IO-node.
  logical, save :: cdf_r_parallel = .false.

  ! To ease the looking up of the level 2-4 elements.
  !  Level 1:
  logical :: has_lvl1 = .false.
  logical :: is_real1 = .true.

  !  Level 2:
  logical :: is_real2 = .true.
  integer, save :: n_k2 = 0
  real(dp), allocatable, save :: bkpt2(:,:)

  !  Level 3:
  logical :: is_real3 = .true.
  integer, save :: n_E3 = 0
  real(dp), allocatable, save :: E3(:)

  !  Level 4:
  logical :: is_real4 = .true.
  integer, save :: n_k4 = 0
  real(dp), allocatable, save :: bkpt4(:,:)
  integer, save :: n_E4 = 0
  real(dp), allocatable, save :: E4(:)

#endif

  ! The algorithm for insertion
  ! If it is 0 it will be assumed that the user
  ! has few dH elements compared to the full structure
  ! If it is 1 it will be assumed that the user
  ! has many dH elements, on the same order of the full
  ! structure.
  integer, save :: insert_algo = 0

  ! Whether this is used or not
  logical :: use_dH = .false.

  public :: tTBTdH, dH, use_dH

  public :: init_dH_options, print_dH_options

#ifdef NCDF_4
  public :: read_next_dH, clean_dH
  public :: read_Sp_dH
#endif
  public :: add_zdH_TriMat, add_zdH_Mat

contains

  subroutine init_dH_options( )

    use parallel, only : Node, Nodes
    use fdf
    use m_os, only : file_exist

#ifdef NCDF_4
    use netcdf_ncdf, ncdf_parallel => parallel
#endif
    
#ifdef MPI
    use mpi_siesta, only : MPI_Bcast, MPI_Comm_World
    use mpi_siesta, only : MPI_Integer, MPI_Logical
    use mpi_siesta, only : MPI_Double_Precision
#endif

#ifdef NCDF_4
    type(hNCDF) :: ndH, grp
#endif
    character(len=20) :: char
    logical :: exists
#ifdef MPI
    integer :: MPIerror
#endif

    use_dH = .false.

#ifdef NCDF_4

    ! Initialize dH to not have any values
    dH%lvl = 0

    ! just set a file-name that should never be created by any user :).
    fname_dH = fdf_get('TBT.dH','NONE1234567890')

    ! If the file exists, use it
    if ( .not. file_exist(fname_dH, Bcast = .true.) ) then
       fname_dH = ' '
       return
    end if

    ! Tell tbtrans to use it
    use_dH = .true.

    ! Ok, the file exists, lets see if all can see the file...
    cdf_r_parallel = .false.
    if ( file_exist(fname_dH, all = .true.) ) then
       
       ! We can check whether we should read in parallel
       ! We default it to read parallelly
       cdf_r_parallel = fdf_get('TBT.dH.Parallel',.true.)

    end if

    ! Insertion algorithm (we default to the user having
    ! very few elements to change)
    insert_algo = 0
    char = fdf_get('TBT.dH.Algorithm','sparse')
    if ( leqi(char,'sparse') ) then
       ! Will loop on the sparsity pattern
       insert_algo = 0
    else if ( leqi(char,'block') .or. leqi(char,'region') ) then
       insert_algo = 1
    end if

    ! The user cannot decide if only one core
    if ( Nodes == 1 ) cdf_r_parallel = .true.

    ! Read in options
    if ( cdf_r_parallel ) then
       call ncdf_open(ndH,fname_dH, mode = NF90_SHARE , &
            parallel = .true. )
    else
       call ncdf_open(ndH,fname_dH, mode = NF90_NOWRITE )
    end if

    ! Read in quantities
    ! Check if level 4 exists
    exists = .false. ! required for cdf_r_parallel == .false.
    call ncdf_inq_grp(ndH,'LEVEL-4',exist = exists)
    if ( exists ) then

       call ncdf_open_grp(ndH,'LEVEL-4',grp)

       ! check the real
       call ncdf_inq_var(grp,'dH',exist=is_real4)

       ! Read in quantities for level 4
       call ncdf_inq_dim(grp, 'ne' , len = n_E4)
       call ncdf_inq_dim(grp,'nkpt', len = n_k4)
       
       if ( n_E4 == 0 .or. n_k4 == 0 ) then

          ! do nothing, we just skip dH for level 4
          if ( Node == 0 ) then
             write(*,'(a)')'tbtrans: Level 4 exists, but has no E or k points. &
                  &It has been discarded.'
          end if

          n_E4 = 0
          n_k4 = 0

       else

          ! Write out that we are using level 4
          if ( Node == 0 ) then
             write(*,'(a,2(i0,a))')'tbtrans: Level 4 will be used, and has ', &
                  n_k4,' k-points and ',n_E4,' energy points.'
          end if

          allocate(E4(n_E4))
          call ncdf_get_var(grp,'E',E4)
          allocate(bkpt4(3,n_k4))
          call ncdf_get_var(grp,'kpt',bkpt4)

       end if
       
    end if

#ifdef MPI
    if ( .not. cdf_r_parallel ) then

       call MPI_Bcast(exists,1,MPI_Logical, 0, &
            MPI_Comm_World, MPIerror)

    if ( exists ) then

       call MPI_Bcast(is_real4,1,MPI_Logical, 0, &
            MPI_Comm_World, MPIerror)

       call MPI_Bcast(n_E4,1,MPI_Integer, 0, &
            MPI_Comm_World, MPIerror)

       if ( n_E4 > 0 ) then
          if ( Node /= 0 ) allocate(E4(n_E4))
          call MPI_Bcast(E4,n_E4,MPI_Double_Precision, 0, &
               MPI_Comm_World, MPIerror)
          call MPI_Bcast(n_k4,1,MPI_Integer, 0, &
               MPI_Comm_World, MPIerror)
          if ( Node /= 0 ) allocate(bkpt4(3,n_k4))
          call MPI_Bcast(bkpt4(1,1),n_k4*3,MPI_Double_Precision, 0, &
               MPI_Comm_World, MPIerror)
       end if

    end if

    end if
#endif

    ! Check if level 3 exists
    exists = .false.
    call ncdf_inq_grp(ndH,'LEVEL-3',exist = exists)
    if ( exists ) then

       call ncdf_open_grp(ndH,'LEVEL-3',grp)

       ! check the real
       call ncdf_inq_var(grp,'dH',exist=is_real3)

       ! Read in quantities for level 3
       call ncdf_inq_dim(grp, 'ne' , len = n_E3)
       
       if ( n_E3 == 0 ) then

          ! do nothing, we just skip dH for level 3
          if ( Node == 0 ) then
             write(*,'(a)')'tbtrans: Level 3 exists, but has no E points. &
                  &It has been discarded.'
          end if

          n_E3 = 0

       else

          ! Write out that we are using level 3
          if ( Node == 0 ) then
             write(*,'(a,i0,a)')'tbtrans: Level 3 will be used, and has ', &
                  n_E3,' energy points.'
          end if

          allocate(E3(n_E3))
          call ncdf_get_var(grp,'E',E3)

       end if
       
    end if

#ifdef MPI
    if ( .not. cdf_r_parallel ) then

       call MPI_Bcast(exists,1,MPI_Logical, 0, &
            MPI_Comm_World, MPIerror)

    if ( exists ) then

       call MPI_Bcast(is_real3,1,MPI_Logical, 0, &
            MPI_Comm_World, MPIerror)

       call MPI_Bcast(n_E3,1,MPI_Integer, 0, &
            MPI_Comm_World, MPIerror)

       if ( n_E3 > 0 ) then
          if ( Node /= 0 ) allocate(E3(n_E3))
          call MPI_Bcast(E3,n_E3,MPI_Double_Precision, 0, &
               MPI_Comm_World, MPIerror)
       end if

    end if

    end if
#endif

    ! Check if level 2 exists
    exists = .false.
    call ncdf_inq_grp(ndH,'LEVEL-2',exist = exists)
    if ( exists ) then

       call ncdf_open_grp(ndH,'LEVEL-2',grp)

       ! check the real
       call ncdf_inq_var(grp,'dH',exist=is_real2)

       ! Read in quantities for level 2
       call ncdf_inq_dim(grp,'nkpt', len = n_k2)
       
       if ( n_k2 == 0 ) then

          ! do nothing, we just skip dH for level 2
          if ( Node == 0 ) then
             write(*,'(a)')'tbtrans: Level 2 exists, but has no points. &
                  &It has been discarded.'
          end if

          n_k2 = 0

       else

          ! Write out that we are using level 2
          if ( Node == 0 ) then
             write(*,'(a,i0,a)')'tbtrans: Level 2 will be used, and has ', &
                  n_k2,' k-points.'
          end if

          allocate(bkpt2(3,n_k2))
          call ncdf_get_var(grp,'kpt',bkpt2)

       end if
       
    end if

#ifdef MPI
    if ( .not. cdf_r_parallel ) then

       call MPI_Bcast(exists,1,MPI_Logical, 0, &
            MPI_Comm_World, MPIerror)

    if ( exists ) then

       call MPI_Bcast(is_real2,1,MPI_Logical, 0, &
            MPI_Comm_World, MPIerror)

       call MPI_Bcast(n_k2,1,MPI_Integer, 0, &
            MPI_Comm_World, MPIerror)

       if ( n_k2 > 0 ) then
          if ( Node /= 0 ) allocate(bkpt2(3,n_k2))
          call MPI_Bcast(bkpt2(1,1),n_k2*3,MPI_Double_Precision, 0, &
               MPI_Comm_World, MPIerror)
       end if

    end if
       
    end if
#endif

    ! Check if level 1 exists
    exists = .false.
    call ncdf_inq_grp(ndH,'LEVEL-1',exist = exists)
    if ( exists ) then
       
       has_lvl1 = .true.
       call ncdf_open_grp(ndH,'LEVEL-1',grp)

       ! check the real
       call ncdf_inq_var(grp,'dH',exist=is_real1)

       ! Write out that we are using level 1
       if ( Node == 0 ) then
          write(*,'(a)')'tbtrans: Level 1 will be used.'
       end if

    end if

#ifdef MPI
    if ( .not. cdf_r_parallel ) then

       call MPI_Bcast(exists,1,MPI_Logical, 0, &
            MPI_Comm_World, MPIerror)

    if ( exists ) then

       has_lvl1 = .true.
       call MPI_Bcast(is_real1,1,MPI_Logical, 0, &
            MPI_Comm_World, MPIerror)

    end if

    end if
#endif

    call ncdf_close(ndH)

#endif

  end subroutine init_dH_options

  subroutine print_dH_options( )

    use parallel, only : IONode

    character(len=*), parameter :: f10='(''tbt: '',a,t53,''='',tr4,a)'
    character(len=*), parameter :: f11='(''tbt: '',a)'
    character(len=*), parameter :: f1 ='(''tbt: '',a,t53,''='',tr4,l1)'

    if ( .not. IONode ) return

#ifdef NCDF_4
    if ( len_trim(fname_dH) == 0 ) then
       write(*,f11)'No delta-Hamiltonian'
       return
    end if
    
    write(*,f10)'User selected dH file', trim(fname_dH)
    write(*,f1) 'Reading in parallel', cdf_r_parallel
    if ( insert_algo == 0 ) then
       write(*,f10)'Inner loop algorithm', 'sparse'
    else if ( insert_algo == 1 ) then
       write(*,f10)'Inner loop algorithm', 'region'
    end if

    if ( .not. cdf_r_parallel .and. n_E3+n_k4+n_E4 > 0 ) then
       write(*,'(a)')'tbtrans: WARNING --- Using level 3 or 4 you must, at least, &
            &have all energy-points in the dH file.'
       write(*,'(a)')'tbtrans: WARNING --- This restriction can be circumventet &
            &if you can use a parallel read (this is highly advised).'
    end if
#else
    write(*,f11)'delta-Hamiltonian not enabled (NetCDF4)'
#endif
    
  end subroutine print_dH_options

#ifdef NCDF_4

  subroutine read_Sp_dH(no_u,sp)

    use class_Sparsity
    use class_OrbitalDistribution

    use m_sparsity_handling, only : Sp_union
    use netcdf_ncdf, ncdf_parallel => parallel

#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif

    use m_ncdf_io, only : cdf_r_Sp

    ! Read in the sparsity pattern
    integer, intent(in) :: no_u
    type(Sparsity), intent(inout) :: sp

    ! Temporary sparsity pattern
    type(Sparsity) :: sp_tmp
    type(OrbitalDistribution) :: fdist

    type(hNCDF) :: ndH, grp

    character(len=7) :: igrp
    logical :: has_level(4)
    integer :: i 

    call delete(sp)

    ! Return if the file is not used
    if ( .not. use_dH ) return

#ifdef MPI
    call newDistribution(no_u,MPI_Comm_Self,fdist,name='TBT-fake dist')
#else
    call newDistribution(no_u,-1           ,fdist,name='TBT-fake dist')
#endif

    if ( cdf_r_parallel ) then
       call ncdf_open(ndH,fname_dH, mode = IOR(NF90_SHARE,NF90_NOWRITE) , &
            parallel = .true. )
    else
       call ncdf_open(ndH,fname_dH, mode = NF90_NOWRITE )
    end if

    ! Assign whether they are existing or not
    has_level(1) = has_lvl1
    has_level(2) = n_k2 > 0
    has_level(3) = n_E3 > 0
    has_level(2) = n_k4 > 0

    do i = 1 , 4
       if ( .not. has_level(i) ) cycle
       
       write(igrp,'(a,i0)') 'LEVEL-',i

       call ncdf_open_grp(ndH,igrp,grp)
       call cdf_r_Sp(grp,no_u,sp_tmp, 'SpdH', Bcast = .not. cdf_r_parallel )
       call Sp_union(fdist,sp_tmp,sp,sp)

    end do

    call delete(fdist)
    call delete(sp_tmp)

    call ncdf_close(ndH)
    
  end subroutine read_Sp_dH

  ! Updates the dH type to the current energy-point
  subroutine read_next_dH(no_u,bkpt,nE)

    use parallel, only : Node, Nodes

#ifdef MPI
    use mpi_siesta, only : MPI_Gather
    use mpi_siesta, only : MPI_Comm_World, MPI_Integer
#endif

    ! Input variables
    integer, intent(in) :: no_u
    real(dp), intent(in) :: bkpt(3)
    type(tNodeE), intent(in) :: nE

    integer :: ik, iE, iN
    integer, allocatable :: nlvl(:)

#ifdef MPI
    integer :: MPIerror
#endif

    if ( .not. use_dH ) return

#ifdef TBTRANS_TIMING
    call timer('read-dH',1)
#endif

    ! Everybody figures out which level they correspond to
    ik = 0
    iE = 0
    allocate(nlvl(0:Nodes-1))
    nlvl(Node) = 0
    if ( n_k4 > 0 ) then
       ik = idx_k(bkpt,bkpt4)
       if ( ik > 0 ) then
          iE = idx_E(nE%E(Node),E4)
          if ( iE == 0 ) ik = 0
       end if
       if ( ik + iE /= 0 ) nlvl(Node) = 4
    end if
    if ( nlvl(Node) == 0 .and. n_E3 > 0 ) then
       iE = idx_E(nE%E(Node),E3)
       if ( iE /= 0 ) nlvl(Node) = 3
    end if
    if ( nlvl(Node) == 0 .and. n_k2 > 0 ) then
       ik = idx_k(bkpt,bkpt2)
       if ( ik /= 0 ) nlvl(Node) = 2
    end if
    if ( nlvl(Node) == 0 .and. has_lvl1 ) then
       nlvl(Node) = 1
    end if
    
    !print *,Node,nlvl(node),ik,bkpt,iE,nE%E(Node) * 13.60580_dp

#ifdef MPI
    if ( .not. cdf_r_parallel ) then
       ! Gather levels on IO
       call MPI_Gather(nlvl(Node),1,MPI_Integer, &
            nlvl(0),1,MPI_Integer,0,MPI_Comm_World, MPIerror)
    end if
#endif

    ! In case there is only one IO node
    if ( .not. cdf_r_parallel ) then

       ! We can only all read the same level
       ! This is because a change in sparsity pattern
       ! might prohibit the Bcast mechanism for the sparsity
       ! patterns.
       if ( sum(nlvl - nlvl(0)) /= 0 ) then
          write(*,'(a)')'Error in using dH functionality'
          write(*,'(a)')'When using non-parallel reading of dH file you must &
               &ensure that at each iteration each core will use the same level.'
          write(*,'(a)')'For (easy) full functionality please see if you can place the dH-file &
               &so that all MPI-cores can see it.'
          call die('Differing level designation and non-MPI IO, please see output...')
       end if

    end if

    if ( nlvl(Node) /= dH%lvl ) then

       call clean_dH( )

    end if

    if ( nlvl(Node) == 4 ) then

       call sub_read_dH(dH,4,is_real4,ik,iE)

    else if ( nlvl(Node) == 3 ) then

       ! We should not have two same energy-points
       ! consecutively
       call sub_read_dH(dH,3,is_real3,0,iE)

    else if ( nlvl(Node) == 2 ) then
   
       if ( dH%lvl == 2 .and. &
            sum(abs(dH%bkpt(:) - bkpt(:))) < 0.00001_dp ) then

          ! We do not need to read anything...
          ! It already has the correct dH

       else

          dH%bkpt(:) = bkpt(:)
          call sub_read_dH(dH,2,is_real2,ik,0)

       end if

    else if ( nlvl(Node) == 1 ) then
       
       if ( dH%lvl == 1 ) then

          ! We do not need to read anything...
          ! It already has the correct dH

       else

          call sub_read_dH(dH,1,is_real1,0,0)

       end if

    end if

    ! Change dH to the current level
    dH%lvl = nlvl(Node)

    deallocate(nlvl)

#ifdef TBTRANS_TIMING
    call timer('read-dH',2)
#endif

  contains

    subroutine sub_read_dH(dH, lvl, is_real, ik, iE)

      use parallel, only : Node

      use class_Sparsity
      use class_OrbitalDistribution
      use netcdf_ncdf, ncdf_parallel => parallel
      use m_ncdf_io, only : cdf_r_Sp

#ifdef MPI
      use mpi_siesta, only : MPI_Send, MPI_Recv
      use mpi_siesta, only : MPI_Comm_World, MPI_Comm_Self, MPI_Status_Size
      use mpi_siesta, only : MPI_Integer
      use mpi_siesta, only : MPI_Double_Precision, MPI_Double_Complex
#endif

      type(tTBTdH), intent(inout) :: dH
      integer, intent(in) :: lvl
      logical, intent(in) :: is_real
      integer, intent(in) :: ik, iE

      real(dp), allocatable :: rH(:)
      complex(dp), pointer :: zH(:)

      type(hNCDF) :: grp
      character(len=7) :: igrp

      type(Sparsity) :: sp
      type(OrbitalDistribution) :: fdist
      
      integer :: nnz, oE
      integer :: start(4)

#ifdef MPI
      ! We must figure out which level each node
      ! lives on
      integer :: MPIerror, status(Mpi_status_size)
#endif

      ! At this point we know which levels we should read
      ! on each node
      write(igrp,'(a,i0)') 'LEVEL-',lvl
      
      !print *,Node,lvl,ik,iE
      
      if ( cdf_r_parallel ) then
         call ncdf_open(grp,fname_dH, mode = IOR(NF90_SHARE,NF90_NOWRITE) , &
              group = igrp , parallel = .true. )
      else
         call ncdf_open(grp,fname_dH, mode = NF90_NOWRITE , group = igrp )
      end if
      
      ! First read in sparsity pattern (the user can have
      ! different sparsity patterns for each level)
      if ( dH%lvl /= lvl ) then
         
         call cdf_r_Sp(grp,no_u,sp, 'SpdH', Bcast = .not. cdf_r_parallel )
         
#ifdef MPI
         call newDistribution(no_u,MPI_Comm_Self,fdist,name='TBT-fake dist')
#else
         call newDistribution(no_u,-1           ,fdist,name='TBT-fake dist')
#endif
          
         ! Create the data container
         call newzSpData1D(sp,fdist,dH%dH,name='TBT dH')
         
         call delete(sp)
         call delete(fdist)
         
      end if
      
      zH  => val(dH%dH)
      nnz =  size(zH)

      ! If the dH file is not read by NF90_SHARE
      ! We must have the IO-node to read and distribute
      if ( is_real ) then
         allocate(rH(nnz))
      end if
      
      start(:)    =  1
      start(2)    =  dH%ispin
      if ( lvl == 1 ) then
      else if ( lvl == 2 ) then
         start(3) = ik
      else if ( lvl == 3 ) then
         start(3) = iE
      else if ( lvl == 4 ) then
         start(3) = iE
         start(4) = ik
      end if
         
#ifdef MPI
      if ( .not. cdf_r_parallel ) then
      if ( lvl == 1 .or. lvl == 2 ) then
         if ( is_real ) then
            call ncdf_get_var(grp,'dH',rH, start = start )
            call MPI_Bcast(rH,nnz,MPI_Double_Precision, 0, &
                 MPI_Comm_World, MPIerror)
!$OMP parallel workshare default(shared)
            zH(:) = rH(:)
!$OMP end parallel workshare
            deallocate(rH)
         else
            call ncdf_get_var(grp,'dH',zH, start = start )
            call MPI_Bcast(zH,nnz,MPI_Double_Complex, 0, &
                 MPI_Comm_World, MPIerror)
         end if

         return

      else 
      if ( Node == 0 ) then
         do iN = 1 , Node - 1
            ! Retrieve the energy point (the k-point IS the same)
            call MPI_Recv(oE,1,MPI_Integer, iN, iN, &
                 MPI_Comm_World, status, MPIerror)
            start(3) = oE
            if ( is_real ) then
               call ncdf_get_var(grp,'dH',rH, start = start )
               call MPI_Send(rH,nnz,MPI_Double_Precision, iN, iN, &
                    MPI_Comm_World, MPIerror)
            else
               call ncdf_get_var(grp,'dH',zH, start = start )
               call MPI_Send(zH,nnz,MPI_Double_Complex, iN, iN, &
                    MPI_Comm_World, MPIerror)
            end if
         end do
         ! Get back the original starting place
         start(3) = iE
      else
         call MPI_Send(iE,1,MPI_Integer, 0, Node, &
              MPI_Comm_World, MPIerror)
         if ( is_real ) then
            call MPI_Recv(rH,nnz,MPI_Double_Precision, 0, Node, &
                 MPI_Comm_World, status, MPIerror)
         else
            call MPI_Recv(zH,nnz,MPI_Double_Complex, 0, Node, &
                 MPI_Comm_World, status, MPIerror)
         end if
      end if
      end if
      end if
#endif
      
      if ( is_real ) then
         call ncdf_get_var(grp,'dH',rH, start = start )
!$OMP parallel workshare default(shared)
         zH(:) = rH(:)
!$OMP end parallel workshare
         deallocate(rH)
      else
         call ncdf_get_var(grp,'dH',zH, start = start )
      end if

      call ncdf_close(grp)

    end subroutine sub_read_dH
    
    function idx_k(bk,fbk) result(i)
      real(dp), intent(in) :: bk(3), fbk(:,:)
      integer :: i
      do i = 1 , size(fbk,dim=2)
         if ( abs( fbk(1,i) - bk(1) ) + &
              abs( fbk(2,i) - bk(2) ) + &
              abs( fbk(3,i) - bk(3) ) < 0.0001_dp ) then
            return
         end if
      end do
      i = 0
    end function idx_k

    function idx_E(E,fE) result(i)
      real(dp), intent(in) :: E, fE(:)
      integer :: i
      do i = 1 , size(fE)
         if ( abs(fE(i) - E) < 7.349806700083788e-06_dp ) then
            return
         end if
      end do
      i = 0
    end function idx_E

    subroutine correct_idx(ik,iE)
      integer, intent(inout), optional :: ik, iE
      if ( present(iE) .and. present(ik) ) then
         if ( iE == 0 ) then
            ik = 0
         else if ( ik == 0 ) then
            iE = 0
         end if
      end if
    end subroutine correct_idx

  end subroutine read_next_dH

  subroutine clean_dH( )

    dH%lvl  = 0
    dH%bkpt = 2.12425_dp
    call delete(dH%dH)

  end subroutine clean_dH

#endif

  ! Add the dH to the tri-diagonal matrix
  subroutine add_zdH_TriMat( zdH, GFinv_tri , r )

    use class_zTriMat
    use class_Sparsity
    use m_region

    use intrinsic_missing, only : SFIND, MODP

    type(zSpData1D), intent(inout) :: zdH
    type(zTriMat), intent(inout) :: GFinv_tri
    type(tRgn), intent(in) :: r

    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: dH(:), GFinv(:)

    integer :: idx, iu, ju, ind, jo, no

    sp => spar(zdH)
    dH => val(zdH)

    call attach(sp, nrows_g=no, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    Gfinv => val(Gfinv_tri)

    if ( insert_algo == 0 ) then

!$OMP parallel do default(shared), private(iu,jo,ind,ju,idx)
       do ju = 1, r%n
          jo = r%r(ju) ! get the orbital in the big sparsity pattern
          if ( l_ncol(jo) /= 0 ) then

          ! Loop on entries here...
          do ind = l_ptr(jo) + 1 , l_ptr(jo) + l_ncol(jo)
             iu = rgn_pivot(r,modp(l_col(ind),no))
             ! Check whether this element should be added
             if ( iu == 0 ) cycle

             idx = index(Gfinv_tri,ju,iu)
          
             GFinv(idx) = GFinv(idx) - dH(ind)
          end do

          end if
       end do
!$OMP end parallel do

    else if ( insert_algo == 1 ) then

!$OMP parallel do default(shared), private(iu,jo,ind,ju,idx)
       do ju = 1, r%n
          jo = r%r(ju) ! get the orbital in the big sparsity pattern
          if ( l_ncol(jo) /= 0 ) then

          ! Loop on entries here...
          do iu = 1 , r%n

             ! Check if the orbital exists in the region
             ! We are dealing with a UC sparsity pattern.
             ind = SFIND(l_col(l_ptr(jo)+1:l_ptr(jo) + l_ncol(jo)),r%r(iu))
             if ( ind == 0 ) cycle
             ind = l_ptr(jo) + ind
             
             idx = index(Gfinv_tri,ju,iu)
             
             GFinv(idx) = GFinv(idx) - dH(ind)
          end do
          
          end if
       end do
!$OMP end parallel do
       
    end if

  end subroutine add_zdH_TriMat

  ! Add the dH to the tri-diagonal matrix
  subroutine add_zdH_Mat( zdH , r, off1,n1,off2,n2, M)

    use class_Sparsity
    use m_region
    use intrinsic_missing, only : SFIND

    type(zSpData1D), intent(inout) :: zdH
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    ! The sizes and offsets of the matrix
    integer, intent(in) :: off1, n1, off2, n2
    complex(dp), intent(inout) :: M(n1,n2)

    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: dH(:)

    integer :: iu, ju, ind, jo

    sp => spar(zdH)
    dH => val(zdH)

    call attach(sp, n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    if ( insert_algo == 0 ) then

!$OMP parallel do default(shared), private(iu,jo,ind,ju)
       do ju = 1 , n1
          jo = r%r(off1+ju) ! get the orbital in the sparsity pattern
          
          if ( l_ncol(jo) /= 0 ) then

          do ind = l_ptr(jo) + 1 , l_ptr(jo) + l_ncol(jo)
             iu = rgn_pivot(r,l_col(ind)) - off2
             if ( iu < 1 .or. n2 < iu ) cycle
             
             M(ju,iu) = M(ju,iu) - dH(ind)
          end do

          end if

       end do
!$OMP end parallel do

    else if ( insert_algo == 1 ) then


!$OMP parallel do default(shared), private(iu,jo,ind,ju)
       do ju = 1 , n1
          jo = r%r(off1+ju) ! get the orbital in the sparsity pattern
          
          if ( l_ncol(jo) /= 0 ) then

          ! Loop on entries here...
          do iu = 1 , n2

             ! Check if the orbital exists in the region
             ! We are dealing with a UC sparsity pattern.
             ind = SFIND(l_col(l_ptr(jo)+1:l_ptr(jo) + l_ncol(jo)),r%r(off2+iu))
             if ( ind == 0 ) cycle
             ind = l_ptr(jo) + ind
             
             M(ju,iu) = M(ju,iu) - dH(ind)

          end do

          end if

       end do
!$OMP end parallel do

    end if

  end subroutine add_zdH_Mat
  
end module m_tbt_dH
