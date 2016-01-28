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
module m_ts_init
  implicit none
  private
  public :: ts_init
contains
  subroutine ts_init(nspin, ucell, na_u, xa)
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
    use m_ts_contour, only : NEn, contour
    use m_ts_kpoints, only : setup_ts_kpoint_grid
    use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight
    use m_ts_options ! Just everything (easier)

    implicit none
! *********************
! * INPUT variables   *
! *********************
    integer, intent(in)  :: nspin
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in)  :: na_u
    real(dp), intent(in) :: xa(3,na_u)

! *********************
! * LOCAL variables   *
! *********************
    complex(dp), dimension(:,:), allocatable :: dos

    ! Read in options for transiesta
    call read_ts_options( ucell )

    ! Setup the k-points
    call setup_ts_kpoint_grid( ucell )

    ! If we actually have a transiesta run we need to process accordingly!
    if ( TSmode ) then

       ! Show every region of the Transiesta run
       call ts_show_regions(ucell,na_u,xa, &
            NBufAtL,NUsedAtomsL,NUsedAtomsR,NBufAtR, &
            NRepA1L,NRepA2L,NRepA1R,NRepA2R)
       
       ! Create the contour lines
       call setup_contour(IsVolt,Cmethod,VoltL,0.0d0,VoltR, &
            NCircle,Nline,Npol,NVolt, &
            0.d0, 0.d0, Ntransport, & ! Transport emin and emax
            CCEmin, GFEta, kt)

       ! Print out the contour path
       call print_contour()
     
       ! Save the contour path to <slabel>.CONTOUR
       call io_contour(slabel)


       ! GF generation:
       allocate(dos(NEn,nspin))
       call memory('A','Z',NEn*nspin,'transiesta')
     
       ! Create the Left GF file
       call do_Green('L',HSFileL, GFFileL, GFTitle, &
            ElecValenceBandBot, ReUseGF, &
            ts_nkpnt,ts_kpoint,ts_kweight, &
            NBufAtL,NUsedAtomsL,NRepA1L,NRepA2L, .false., & !For now TranSIESTA will only perform with inner-cell distances
            ucell,xa,na_u,NEn,contour,VoltL,dos,nspin)
       
       ! Create the Right GF file
       call do_Green('R',HSFileR,GFFileR, GFTitle, &
            ElecValenceBandBot, ReUseGF, &
            ts_nkpnt,ts_kpoint,ts_kweight, &
            NBufAtR,NUsedAtomsR,NRepA1R,NRepA2R, .false., &
            ucell,xa,na_u,NEn,contour,VoltR,dos,nspin)
       
       call memory('D','Z',NEn*nspin,'transiesta')
       deallocate(dos)

       ! Print out information in Green's function files
       ! Show the number of used atoms and orbitals
       if ( IONode ) then
          write(*,'(/,a,i6,'' / '',i6)') &
               'Left : GF atoms    / Expanded atoms    : ',NUsedAtomsL, &
               NUsedAtomsL*NRepA1L*NRepA2L
          write(*,'(a,i6,'' / '',i6)') &
               'Left : GF orbitals / Expanded orbitals : ',NUsedOrbsL, &
               NUsedOrbsL*NRepA1L*NRepA2L
          write(*,'(a,i6,'' / '',i6)') &
               'Right: GF atoms    / Expanded atoms    : ',NUsedAtomsR, &
               NUsedAtomsR*NRepA1R*NRepA2R
          write(*,'(a,i6,'' / '',i6)') &
               'Right: GF orbitals / Expanded orbitals : ',NUsedOrbsR, &
               NUsedOrbsR*NRepA1R*NRepA2R

          write(*,*) ! New line
       end if
       
    end if
    
  end subroutine ts_init
end module m_ts_init

