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

  subroutine ts_init(wmix, kT, nspin, ucell, na_u, xa, lasto, no_u, inicoor, fincoor )
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
    
    use m_ts_kpoints, only : setup_ts_kpoint_grid
    use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight
    use m_ts_kpoints, only : ts_kscell, ts_kdispl
    use m_ts_cctype
    use m_ts_electype
    use m_ts_options ! Just everything (easier)
    use m_ts_method

    use m_ts_global_vars, only : TSmode, TSinit
    use siesta_options, only : isolve, SOLVE_TRANSI

    use m_fixed, only : is_fixed

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
    integer, intent(in)  :: inicoor, fincoor
    
! *********************
! * LOCAL variables   *
! *********************
    complex(dp), dimension(:,:), allocatable :: dos
    integer :: i, sNE, eNE, ia
    integer :: nC, nTS

    if ( isolve .eq. SOLVE_TRANSI ) then
       TSmode = .true.
       ! If in TSmode default to initalization
       ! In case of 'DM.UseSaveDM TRUE' TSinit will be set accordingly
       TSinit = .true.
    end if

    ! Read in options for transiesta
    call read_ts_options( wmix, kT, ucell , na_u , xa, lasto )

    ! Setup the k-points, must be done after options reading 
    ! (determine the transport direction)
    call setup_ts_kpoint_grid( ucell )

    ! If we actually have a transiesta run we need to process accordingly!
    if ( .not. TSmode ) return

    ! Check the electrodes
    do i = 1 , N_Elec
       call check_Elec(Elecs(i),nspin,na_u,xa,lasto,Elecs_xa_EPS, &
            kcell=ts_kscell,kdispl=ts_kdispl)
    end do

    ! initialize regions of the electrodes and device
    ! the number of LCAO orbitals on each atom will not change
    call ts_init_region_types(na_BufL,Elecs,na_BufR,na_u,lasto)

    ! Check that an eventual CGrun will fix all electrodes and 
    ! buffer atoms
    if ( fincoor - inicoor > 0 ) then
       ! check fix
       do ia = 1 , na_BufL
          if ( .not. is_fixed(ia) ) then
             call die('All buffer atoms and electrode atoms *MUST* be &
                  &fixed while doing transiesta geometry optimizations. &
                  &Please correct left buffer.')
          end if
       end do
       do i = 1 , N_Elec
          do ia = Elecs(i)%idx_na , Elecs(i)%idx_na + TotUsedAtoms(Elecs(i)) - 1
             if ( .not. is_fixed(ia) ) then
                call die('All buffer atoms and electrode atoms *MUST* be &
                     &fixed while doing transiesta geometry optimizations. &
                     &Please correct electrodes.')
             end if
          end do
       end do
       do ia = na_u - na_BufR + 1 , na_u
          if ( .not. is_fixed(ia) ) then
             call die('All buffer atoms and electrode atoms *MUST* be &
                  &fixed while doing transiesta geometry optimizations. &
                  &Please correct right buffer.')
          end if
       end do
    end if

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
         na_BufL,N_Elec,Elecs,na_BufR)

    if ( .not. TS_Analyze ) then

       ! GF generation:
       do i = 1 , N_Elec

          ! initialize the electrode for Green's function calculation
          call init_Electrode_HS(Elecs(i),RemUCellDistance)

          call do_Green(Elecs(i), &
               ucell,ts_nkpnt,ts_kpoint,ts_kweight, &
               RemUCellDistance, Elecs_xa_Eps, &
               .false.)
          
          ! clean-up
          call delete(Elecs(i))
          
       end do
!call die('')
    else

       do i = 1 , N_Elec
          call delete(Elecs(i)) ! ensure clean electrode
          call read_Elec(Elecs(i),Bcast=.true.)
          
          ! print out the precision of the electrode (whether it extends
          ! beyond first principal layer)
          call check_Connectivity(Elecs(i))

          call delete(Elecs(i))
       end do
       
    end if

  end subroutine ts_init

end module m_ts_init

