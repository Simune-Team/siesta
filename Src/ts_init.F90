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
module m_ts_init

  implicit none
  private
  public :: ts_init

contains

  subroutine ts_init(wmix, kT, nspin, ucell, na_u, xa, lasto, no_u)
  ! Routine for initializing everything related to the Transiesta package.
  ! This is to comply with the SIESTA way of initializing at the beginning
  ! and make the code more managable.

  ! This routine does the following in that order:
  !  1. Read in ts_options
  !  2. Create Transiesta k-points
  !  3. Setup the contour path
  !  4. Create the GF related to the integration scheme
  
! Used modules
    use parallel, only : IONode
    use files, only : slabel

    use m_ts_gf,        only : do_Green
    use m_ts_electrode, only : init_Electrode_HS
    
    use m_ts_contour, only : contour_Eq, contour_EqL, contour_EqR, contour_nEq
    use m_ts_contour, only : sort_contour
    use m_ts_contour, only : NEn, contour
    use m_ts_kpoints, only : setup_ts_kpoint_grid
    use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight
    use m_ts_kpoints, only : ts_kscell, ts_kdispl
    use m_ts_cctype
    use m_ts_electype
    use m_ts_options ! Just everything (easier)
    use m_ts_method

    implicit none
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: wmix
    real(dp), intent(in) :: kT
    integer, intent(in)  :: nspin
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in)  :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: lasto(0:na_u)
    integer, intent(in)  :: no_u

! *********************
! * LOCAL variables   *
! *********************
    complex(dp), dimension(:,:), allocatable :: dos
    type(ts_ccontour), pointer :: c(:) => null()
    integer :: i, sNE, eNE
    integer :: nC, nTS
    logical :: RemUCellDistance

    RemUCellDistance = .false.

    ! Read in options for transiesta
    call read_ts_options( wmix, kT, ucell , na_u , xa, lasto )

    ! Setup the k-points, must be done after options reading (determine the transport
    ! direction)
    call setup_ts_kpoint_grid( ucell )

    ! If we actually have a transiesta run we need to process accordingly!
    if ( TSmode ) then

       ! Check the electrodes
       do  i = 1 , size(Elecs)
          call check_Elec(Elecs(i),nspin,na_u,xa,lasto,Elecs_xa_EPS, &
               kcell=ts_kscell,kdispl=ts_kdispl)
       end do

       ! initialize regions of the electrodes and device
       ! the number of LCAO orbitals on each atom will not change
       call ts_init_region_types(na_BufL,Elecs,na_BufR,na_u,lasto)

       ! here we will check for all the size requirements of Transiesta
       ! I.e. whether some of the optimizations can be erroneous or not

       ! Do a crude check of the sizes
       ! if the transiesta region is equal of size to or smaller 
       ! than the size of the combined electrodes, then the system
       ! is VERY WRONG...
       ! First calculate L/C/R sizes (we remember to subtract the buffer
       ! orbitals)
       nTS = no_u - no_BufL - no_BufR
       nC  = nTS  - sum(TotUsedOrbs(Elecs))
       if ( nC < 1 ) &
            call die("The contact region size is &
               &smaller than the electrode size. &
               &What have you done? Please correct this insanity...")

       if ( minval(TotUsedOrbs(Elecs)) < 2 ) &
            call die('We cannot perform sparse pattern on the electrode &
            &system.')

       ! Show every region of the Transiesta run
       call ts_show_regions(ucell,na_u,xa, &
            na_BufL,size(Elecs),Elecs,na_BufR)
       
       ! Create the contour lines
!       call setup_contour(IsVolt)

       ! Print out the contour path
!       call ts_print_contour(contour)
     
       ! Save the contour path to <slabel>.CONTOUR
!       call ts_io_contour(contour,slabel)

       ! We sort the contour to obtain the highest 
       ! numerical accuracy (simply sort by weight in ascending order)
       eNE = 0
!       if ( IsVolt ) then
!          c => contour_EqL()
!          eNE = eNE + size(c)
          ! 1) sort the left equilibrium points
!          call sort_contour(size(c),c)
!          c => contour_EqR()
!          eNE = eNE + size(c)
          ! 2) sort the right equilibrium points
!          call sort_contour(size(c),c)
!          c => contour_nEq()
!          eNE = eNE + size(c)
          ! 3) sort the non-equilibrium points
!          call sort_contour(size(c),c)
!       else
!          c => contour_Eq()
!          eNE = eNE + size(c)
          ! 1) sort the equilibrium points
!          call sort_contour(size(c),c)
!       end if
!       if ( eNE /= NEn ) call die('Wrong sorting of the contour')

       if ( .not. TS_Analyze ) then

          ! GF generation:
          do i = 1 , size(Elecs)

             ! initialize the electrode for Green's function calculation
             call init_Electrode_HS(Elecs(i),RemUCellDistance)

             call do_Green(Elecs(i), ReUseGF, &
                  ucell,ts_nkpnt,ts_kpoint,ts_kweight, &
                  RemUCellDistance, Elecs_xa_Eps, &
                  .false.)

             ! clean-up
             call delete(Elecs(i))

          end do

       end if

       ! Print out information in Green's function files
       ! Show the number of used atoms and orbitals
       if ( IONode ) then
          do i = 1 , size(Elecs)
             write(*,'(a30,i6,'' / '',i6)') &
                  trim(Name(Elecs(i)))//' : GF atoms    / Expanded atoms    : ', &
                  UsedAtoms(Elecs(i)), TotUsedAtoms(Elecs(i))
             write(*,'(a30,i6,'' / '',i6)') &
                  trim(Name(Elecs(i)))//' : GF orbitals / Expanded orbitals : ', &
                  UsedOrbs(Elecs(i)), TotUsedOrbs(Elecs(i))
          end do
          write(*,*) ! New line
       end if

    end if
    
  end subroutine ts_init

end module m_ts_init

