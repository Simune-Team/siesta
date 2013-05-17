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

  subroutine ts_init(nspin, ucell, na_u, xa, lasto, no_u)
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

    use m_ts_gf,      only : do_Green
    use m_ts_contour, only : setup_contour, print_contour, io_contour
    use m_ts_contour, only : sort_contour
    use m_ts_contour, only : NEn, contour
    use m_ts_kpoints, only : setup_ts_kpoint_grid
    use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight
    use m_ts_electype
    use m_ts_options ! Just everything (easier)
    use m_ts_method

    implicit none
! *********************
! * INPUT variables   *
! *********************
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
    integer :: i, sNE, eNE
    integer :: nBL, nBR
    integer :: nL, nC, nR, nTS

    ! Read in options for transiesta
    call read_ts_options( ucell )

    ! Setup the k-points
    call setup_ts_kpoint_grid( ucell )

    ! If we actually have a transiesta run we need to process accordingly!
    if ( TSmode ) then

       ! Figure out the number of orbitals on the buffer atoms
       nBL = 0
       do i = 1 , NBufAtL
          nBL = nBL + lasto(i) - lasto(i-1)
       end do
       nBR = 0
       do i = na_u - NBufAtR+1 , na_u
          nBR = nBR + lasto(i) - lasto(i-1)
       end do

       ! here we will check for all the size requirements of Transiesta
       ! I.e. whether some of the optimizations can be erroneous or not

       ! First calculate L/C/R sizes (we remember to subtract the buffer
       ! orbitals)
       nTS = no_u - nBL - nBR
       nL  = TotUsedOrbs(ElLeft)
       nR  = TotUsedOrbs(ElRight)
       nC  = nTS  - nL  - nR
       if ( nC < 1 ) &
            call die('Transiesta central region does not &
            &exist. What have you done?')

       if ( ts_method /= TS_ORIGINAL ) then

          if ( nL < 2 .or. nR < 2 ) &
               call die('We cannot perform sparse pattern on the electrode &
               &system.')

          i = 0

          ! We must ensure that the Sigma[LR] have enough space to hold
          ! one line of the full matrix (we use them as work arrays
          ! when calculating the Gamma[LR]
          if ( nL ** 2 < nTS .or. nR ** 2 < nTS ) then
             call die('The current implementation requires that the &
                  &square of the orbitals in the electrodes are larger &
                  &than the dimension of the problem.')
          end if

          if ( ts_method == TS_SPARSITY_TRI .and. IsVolt ) then

             if ( nL > nC .or. nR > nC ) then
                write(*,'(a,2(i0,tr1,''/'',tr1),i0)') &
                     'Sizes are L/C/R: ',nL,nC,nR
                call die('Your system &
                     &has inappropriate sizes for memory limited &
                     &tri-diagonalization')
             end if

          end if

       end if

       ! Show every region of the Transiesta run
       call ts_show_regions(ucell,na_u,xa, &
            NBufAtL,ElLeft, ElRight,NBufAtR)
       
       ! Create the contour lines
       call setup_contour(IsVolt,Cmethod,VoltL,0.0d0,VoltR, &
            NCircle,Nline,Npol,NVolt, &
            0.d0, 0.d0, Ntransport, & ! Transport emin and emax
            CCEmin, GFEta, kt)

       ! Print out the contour path
       call print_contour()
     
       ! Save the contour path to <slabel>.CONTOUR
       call io_contour(slabel)

       ! We sort the contour to obtain the highest 
       ! numerical accuracy (simply sort by weight in ascending order)
       if ( IsVolt ) then
          sNE =       1
          eNE = sNE+NCircle+Npol+Nline-1
          ! 1) sort the left equilibrium points
          call sort_contour(eNE-sNE+1,contour(sNE:eNE))
          sNE = eNE + 1
          eNE = sNE+NCircle+Npol+Nline-1
          ! 2) sort the right equilibrium points
          call sort_contour(eNE-sNE+1,contour(sNE:eNE))
          sNE = eNE + 1
          eNE = sNE + NVolt - 1
          ! 3) sort the non-equilibrium points
          call sort_contour(eNE-sNE+1,contour(sNE:eNE))
       else
          sNE =       1
          eNE = sNE+NCircle+Npol+Nline-1
          ! 1) sort the equilibrium points
          call sort_contour(eNE-sNE+1,contour(sNE:eNE))
       end if
       if ( eNE /= NEn ) call die('Wrong sorting in the contour points.')

       ! GF generation:
       allocate(dos(NEn,nspin))
       call memory('A','Z',NEn*nspin,'transiesta')
     
       ! Create the Left GF file
       call do_Green('L',ElLeft, GFFileL, GFTitle, &
            ElecValenceBandBot, ReUseGF, &
            ts_nkpnt,ts_kpoint,ts_kweight, &
            NBufAtL, .false., & !For now TranSIESTA will only perform with inner-cell distances
            ucell,xa,na_u,NEn,contour,VoltL,.false.,dos,nspin)
       
       ! Create the Right GF file
       call do_Green('R',ElRight,GFFileR, GFTitle, &
            ElecValenceBandBot, ReUseGF, &
            ts_nkpnt,ts_kpoint,ts_kweight, &
            NBufAtR, .false., &
            ucell,xa,na_u,NEn,contour,VoltR,.false.,dos,nspin)
       
       call memory('D','Z',NEn*nspin,'transiesta')
       deallocate(dos)

       ! Print out information in Green's function files
       ! Show the number of used atoms and orbitals
       if ( IONode ) then
          write(*,'(/,a,i6,'' / '',i6)') &
               'Left : GF atoms    / Expanded atoms    : ',UsedAtoms(ElLeft), &
               TotUsedAtoms(ElLeft)
          write(*,'(a,i6,'' / '',i6)') &
               'Left : GF orbitals / Expanded orbitals : ',UsedOrbs(ElLeft), &
               TotUsedOrbs(ElLeft)
          write(*,'(a,i6,'' / '',i6)') &
               'Right: GF atoms    / Expanded atoms    : ',UsedAtoms(ElRight), &
               TotUsedAtoms(ElRight)
          write(*,'(a,i6,'' / '',i6)') &
               'Right: GF orbitals / Expanded orbitals : ',UsedOrbs(ElRight), &
               TotUsedOrbs(ElRight)

          write(*,*) ! New line
       end if
       
    end if
    
  end subroutine ts_init

end module m_ts_init

