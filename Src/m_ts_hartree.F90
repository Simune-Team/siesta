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
! This code segment has been improved or fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
!
module m_ts_hartree

! Module for fixing the Hartree potential so that the potential fluctuations
! does not go wild.
! This is necessary to get a stable SCF solution
!
! Created and copyrighted by: Nick Papior Andersen, 2014
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the author.
  
  use precision, only : dp
  use m_ts_electype
  use m_ts_tdir, only : ts_tdir

  implicit none
  
  private
  save

  ! The idea is to have sub routines in this module to do
  ! various Hartree potential fixes
  public :: ts_init_hartree_fix
  public :: ts_hartree_fix

  ! The electrode that provides the basin of the constant potential
  type(Elec), pointer :: El => null()

  ! Whether we should employ the electrode
  ! basal plane or not
  logical, public :: elec_basal_plane = .false.

contains

  subroutine ts_init_hartree_fix(ucell,na_u,xa,meshG,nsm, N_Elec, Elecs)

    use units, only : Ang
    use intrinsic_missing, only: VNORM
    use m_ts_mesh, only : meshl, offset_i, offset_r, dMesh, dL
    use parallel, only : IONode
#ifdef MPI
    use mpi_siesta, only : MPI_AllReduce, MPI_Sum
    use mpi_siesta, only : MPI_Comm_World, MPI_integer
#endif

! ***********************
! * INPUT variables     *
! ***********************
    real(dp),   intent(in) :: ucell(3,3)
    integer,    intent(in) :: na_u
    real(dp),   intent(in) :: xa(3,na_u)
    integer,    intent(in) :: meshG(3), nsm
    integer,    intent(in) :: N_Elec
    type(Elec), intent(in), target :: Elecs(N_Elec)

    integer :: iE
    real(dp) :: area, tmp
    integer :: i1, i2, i3, nlp
    real(dp) :: ll(3), llZ(3), llYZ(3)
    real(dp), external :: volcel
#ifdef MPI
    integer :: MPIerror
#endif

    ! We now were to put the Hartree correction 
    if ( .not. elec_basal_plane ) return

    ! Easy determination of largest basal plane of electrodes
    area = -1._dp
    do iE = 1 , N_Elec
       tmp = VOLCEL(Elecs(iE)%ucell)
       tmp = tmp / VNORM(Elecs(iE)%ucell(:,Elecs(iE)%t_dir))
       if ( tmp > area ) then
          area =  tmp
          El   => Elecs(iE)
       end if
    end do

    ! We check that we actually process something...
    nlp = 0
    do i3 = 0 , meshl(3) - 1
       llZ(:) = offset_r(:) + i3*dL(:,3)
       do i2 = 0 , meshl(2) - 1
          llYZ(:) = i2*dL(:,2) + llZ(:)
          do i1 = 0 , meshl(1) - 1
             ll(:) = i1*dL(:,1) + llYZ(:)
             if ( in_basal_Elec(El%p,ll,dMesh) ) then
                nlp = nlp + 1
             end if
          end do
       end do
    end do

#ifdef MPI
    call MPI_AllReduce(nlp,i1,1,MPI_integer,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    nlp = i1
#endif

    if ( IONode ) then
       write(*,*)
       write(*,'(3a)')   'transiesta: Using electrode: ',trim(El%Name),' for Hartree correction'
       write(*,'(a,i0)') 'transiesta: Number of points used: ',nlp
       if ( nlp == 0 ) then
          write(*,'(a)') 'transiesta: Basal plane of electrode '// &
               trim(El%name)//' might be outside of &
               &unit cell.'
          write(*,'(a)') 'transiesta: Please move structure so this point is &
               inside unit cell (Ang):'
          write(*,'(a,3(tr1,f13.5))') 'transiesta: Point (Ang):',&
               El%p%c/Ang
       end if
       write(*,*)
    end if

    if ( nlp == 0 ) then
       call die('The partitioning of the basal plane went wrong. &
            &No points are encapsulated.')
    end if

  end subroutine ts_init_hartree_fix

#ifdef INTEL_COMPILER_ERROR 
! This code segment will create an error with the
! intel compiler compiled at a high setting: TODO (not with Gfortran)
!  -m64 -O3 -xHost -fp-model source -fp-model except -ip -prec-div -prec-sqrt
! If I remove the pointer/target feature below it works beautifully!
! *** NOT GOOD ***
  integer, target :: i10, i20, i30
  integer, pointer :: iT
  if ( ts_tdir == 1 ) then
     i10 = 0
     iT => i10
  else if ( ts_tdir == 2 ) then
     i20 = 0
     iT => i20
  else if ( ts_tdir == 3 ) then
     i30 = 0
     iT => i30
  else
     call die('Hartree fix, not implemented')
  end if

  i10 = 0
  i20 = offset_i(2) - 1
  i30 = offset_i(3) - 1
  if ( iT <= 0 ) then
     imesh = 0
     i30 = offset_i(3) - 1
     do i3 = 0,meshl(3)-1
        i30 = i30 + 1
        i20 = offset_i(2) - 1
        do i2 = 0,meshl(2)-1
           i20 = i20 + 1
           do i10 = 0,meshl(1)-1
              imesh = imesh + 1
              if (iT.eq.0) then
                 nlp = nlp + 1
                 Vtot = Vtot + Vscf(imesh)
              end if
           end do
        end do
     end do
  end if
#endif

  ! Fix the potential
  subroutine ts_hartree_fix( ntpl , Vscf )
    use precision, only : grid_p
    use sys, only : die
#ifdef MPI
    use mpi_siesta, only : MPI_AllReduce, MPI_Sum
    use mpi_siesta, only : MPI_Comm_World, MPI_integer
    use mpi_siesta, only : MPI_double_precision
#endif
    use m_ts_mesh, only : meshl, offset_i, offset_r, dMesh, dL
    
    integer     , intent(in)    :: ntpl
    real(grid_p), intent(inout) :: Vscf(ntpl)

! Internal variables
    integer :: i1 , i2 , i3
    integer :: imesh, nlp
    integer :: i10, i20, i30
#ifdef MPI
    integer :: MPIerror
#endif
    real(dp) :: Vav, Vtot, temp
    real(dp) :: ll(3), llZ(3), llYZ(3)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE TS_VH_fix' )
#endif

    ! Initialize summation
    Vtot = 0._dp

    ! Initialize counters
    nlp   = 0
    imesh = 0

    if ( elec_basal_plane ) then
       ! This is an electrode averaging...
       do i3 = 0 , meshl(3) - 1
          llZ(:) = offset_r(:) + i3*dL(:,3)
          do i2 = 0 , meshl(2) - 1
             llYZ(:) = i2*dL(:,2) + llZ(:)
             do i1 = 0 , meshl(1) - 1
                ll(:) = i1*dL(:,1) + llYZ(:)
                imesh = imesh + 1
                if ( in_basal_Elec(El%p,ll,dMesh) ) then
                   nlp  = nlp + 1
                   Vtot = Vtot + Vscf(imesh)
                end if
             end do
          end do
       end do
    else if ( ts_tdir == 1 ) then
       do i3 = 1 , meshl(3)
          do i2 = 1 , meshl(2)
             i10 = offset_i(1) - 1
             do i1 = 1 , meshl(1)
                i10 = i10 + 1
                imesh = imesh + 1
                if ( i10 == 0 ) then
                   nlp  = nlp + 1
                   Vtot = Vtot + Vscf(imesh)
                end if
             end do
          end do
       end do
    else if ( ts_tdir == 2 ) then
       do i3 = 1 , meshl(3)
          i20 = offset_i(2) - 1
          do i2 = 1 , meshl(2)
             i20 = i20 + 1
             do i1 = 1 , meshl(1)
                imesh = imesh + 1
                if ( i20 == 0 ) then
                   nlp  = nlp + 1
                   Vtot = Vtot + Vscf(imesh)
                end if
             end do
          end do
       end do
    else if ( ts_tdir == 3 ) then
       i30 = offset_i(3) - 1
       do i3 = 1 , meshl(3)
          i30 = i30 + 1
          do i2 = 1 , meshl(2)
             do i1 = 1 , meshl(1)
                imesh = imesh + 1
                if ( i30 == 0 ) then
                   nlp  = nlp + 1
                   Vtot = Vtot + Vscf(imesh)
                end if
             end do
          end do
       end do
    else
       call die('Something went extremely wrong...Hartree Fix')
    end if

#ifdef MPI
    call MPI_AllReduce(Vtot,temp,1,MPI_double_precision,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    Vtot = temp
    call MPI_AllReduce(nlp,i1,1,MPI_integer,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    nlp = i1
#endif

    if ( nlp == 0 ) then
       call die('The partitioning of the basal plane went wrong. &
            &No points are encapsulated.')
    end if
    
    Vav = Vtot / real(nlp,dp)
    
    ! Align potential
    Vscf(1:ntpl) = Vscf(1:ntpl) - Vav

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS TS_VH_fix' )
#endif

  end subroutine ts_hartree_fix

end module m_ts_hartree
