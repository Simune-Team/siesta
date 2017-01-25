! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!
! ##################################################################
! ##                           The                                ##
! ##     Localized Basis-set Nonequilibrium Greens function       ## 
! ##                      Transport Program                       ##
! ##                                                              ## 
! ##                            By                                ##
! ##         Nick Papior Andersen, nickpapior@gmail.com           ##
! ##       Nanotech, Technical University of Denmark (DTU)        ##
! ##                                                              ##
! ##################################################################
!
!
! Tight-binding transport program.
! Copyright by Nick Papior Andersen, 2012.
! The use of this program is allowed for non-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the author.
!
! I (the author) have gained permission by Mads Brandbyge to rewrite 
! the TBTrans utility to be streamlined.
!
! It is a complete rewrite with emphasis on re-using the already
! optimized routines in TranSIESTA to use the full potential of the 
! parallel environment it resides in.
!
! Therefore if there has been a lot of changes in the TranSIESTA 
! codes please update this as well.

! To see the current used files from SIESTA/TranSIESTA please
! see the Makefile as it will be the easiest to maintain.

! All codes are supplied in Fortran 90+ standards to move away from the 
! old F77 syntax.

program tbtrans

! ************************
! * SIESTA modules       *
! ************************
  use parallel,       only : Node, Nodes, IONode, operator(.PARCOUNT.)
  use units,          only : eV,Pi
  use sys,            only : die
  use precision,      only : dp
  use files,          only : slabel
  use m_timestamp,    only : timestamp
  use m_wallclock,    only : wallclock
  use alloc,          only : alloc_report

#ifdef MPI
  use mpi_siesta
#endif
  use m_hs_matrix,   only : set_hs_matrix, matrix_rem_left_right
  use m_hs_matrix,   only : matrix_symmetrize

! ************************
! * TranSIESTA modules   *
! ************************
  use m_ts_contour,   only : NEn, PNEn, contour
  use m_ts_gf,        only : do_Green, read_Green
  use m_ts_scattering,only : getSFE
  use m_ts_io,        only : ts_iohs
  use m_ts_electrode, only : create_Green

! ************************
! * TBtrans modules      *
! ************************
  use m_tbt_options, only : UseBulk, kT, GFEta
  use m_tbt_options, only : VoltFDF, VoltL, VoltR, IsVolt
  use m_tbt_options, only : NBufAtL, NRepA1L, NRepA2L, NUsedAtomsL
  use m_tbt_options, only : nuaL_GF => NUsedAtomsL
  use m_tbt_options, only : noL_GF  => NUsedOrbsL
  use m_tbt_options, only : GFFileL, HSFileL
  use m_tbt_options, only : NBufAtR, NRepA1R, NRepA2R, NUsedAtomsR
  use m_tbt_options, only : nuaR_GF => NUsedAtomsR
  use m_tbt_options, only : noR_GF  => NUsedOrbsR
  use m_tbt_options, only : GFFileR, HSFileR
  use m_tbt_options, only : HSFile, ElecValenceBandBot, ReUseGF
  use m_tbt_options, only : GFTitle
  use m_tbt_options, only : Emin, Emax, NeigCh
  use m_tbt_options, only : IsoAt1, IsoAt2
  use m_tbt_options, only : CalcIeig
  use m_tbt_options, only : CalcCOOP, AlignScat, CalcAtomPDOS
  use m_tbt_options, only : RemUCellDistances

  use m_tbt_kpoints, only : siesta_Gamma
  use m_tbt_kpoints, only : Gamma, spiral
  use m_tbt_kpoints, only : nkpnt, kpoint, kweight
  use m_tbt_kpoints, only : setup_tbt_kpoint_grid

  use m_tbt_out,     only : create_file, out_NEWLINE
  use m_tbt_out,     only : out_Trans, out_kpt_header
  use m_tbt_out,     only : out_EIG, out_TEIG
  use m_tbt_out,     only : out_DOS, out_Trans
  use m_tbt_out,     only : out_REGION, out_DEVICE

  use m_tbt_read_tshs,only: tbt_read_tshs

  implicit none

! *****************************************
! * Parameters used extensively           *
! *****************************************
  real(dp), parameter :: r1dPi = 1.0_dp/Pi

! *****************************************
! * Variables concerning the contour      *
! *****************************************
  complex(dp), dimension(:,:), allocatable :: ZBulkDos ! We need to have it complex to pass to create_Green


! *****************************************
! * Electrode variables                   *
! *****************************************
! Number of atoms in the electrode (after expansion)
  integer :: nuaL, nuaR
! Number of orbitals in the electrodes (after expansion)
  integer :: noL, noR
! The lasto shortened to only the electrodes
  integer, dimension(:), allocatable :: lastoL,lastoR

! Read-in bulk H, S, for a given q-points
  complex(dp), dimension(:,:,:), allocatable :: HAAL, SAAL
  complex(dp), dimension(:,:,:), allocatable :: HAAR, SAAR
! Self energies and scattering matrices
! Gamma=i*0.5*(Sigma - Sigma^dagger) 
  complex(dp), dimension(:,:),   allocatable :: SFEL
  complex(dp), dimension(:,:),   allocatable :: SFER

! Electrode k-points kpar-points and their weights:
! Notice that these are at the moment completely identical
! to the TranSIESTA k-points at time of creation.
! Therefore they should be the same as ts_kpoints until 
! implementation is completed to handle different electrode k-point
! sampling
  integer                                :: nkparL,nkparR
  real(dp), dimension (:,:), allocatable :: kparL ,kparR
  real(dp), dimension (:)  , allocatable :: wkparL,wkparR

! q-points and their weights:
! The notation of b tells the programmer that
! q is in units of b_i**-1 and NOT in Bohr**-1
  integer                                :: nqL,nqR
  real(dp), dimension (:,:), allocatable :: qLb,qRb
  real(dp), dimension (:),   allocatable :: wqL,wqR

! *****************************************
! * CONTACT variables                     *
! *****************************************
! Variables concerning the TSHS file, and the sparse format
  logical                           :: Gamma_Scat ! Is it a Gamma Calculation?
  integer                           :: na_u,no_u ! Unit cell atoms / orbitals
  integer                           :: maxnh,no_s,nspin ! Hamiltonian size / total orbitals / spins
  real(dp), dimension(:,:), pointer :: H ! Hamiltonian in SPARSE format
  real(dp), dimension(:), pointer   :: S ! Overlap in SPARSE format
  real(dp), dimension(:,:), pointer :: xa ! atomic coordinates
  real(dp), dimension(:,:), pointer :: xij ! differences with unitcell, differences with unitcell
  integer, dimension(:), pointer    :: lasto,numh,listhptr
  integer, dimension(:), pointer    :: listh,indxuo
  real(dp)                          :: Ef ! Efermi
! Specifics used for the calculation
  complex(dp), dimension(:,:), pointer :: Hk, Sk ! Hamilton and overlap in k-point
  complex(dp), dimension(:,:), pointer :: Hk_D, Sk_D ! Hamilton and overlap in k-point for the device region
  complex(dp), dimension(:,:), pointer :: Hk_iso, Sk_iso ! Hamilton and overlap in k-point for isolated region
#ifdef MPI
  complex(dp), dimension(:,:), pointer :: Hk_recv, Sk_recv ! Hamilton and overlap in k-point for isolated region
#endif
! As we employ a "double" parallel execution we use twice as much memory
! and then speed up the calculation by holding the next 'Nodes' number of k-point Hamiltonians
! in memory
#ifdef MPI
  complex(dp), dimension(:,:), pointer :: Hk_Node, Sk_Node ! Hamilton and overlap in k-point kept in node
#endif
  complex(dp), dimension(:,:), pointer :: GF, GFRGF
  ! Size of the contact without the buffer regions
  integer :: nou
  ! Size of the scattering region without the electrode and buffer orbitals 
  integer :: noD

! *****************************************
! * Buffer variables                      *
! *****************************************
  integer :: noBufL, noBufR

! *****************************************
! * Isolated region variables             *
! *****************************************
  integer :: Isoo1, Isoo2   ! Equivalent to IsoAt[12] to orbitals (with buffer)
  integer :: Isoo1C, Isoo2C ! Equivalent to IsoAt[12] to orbitals in the CONTACT region (without buffer)
  integer :: Isoo1D, Isoo2D ! Equivalent to IsoAt[12] to orbitals in the DEVICE region (without buffer and electrode)
  integer :: Isoo           ! Number of orbitals in the isolated region

! *****************************************
! * Transmission variables                *
! *****************************************
  complex(dp) :: ZEnergy, ZwGF
  real(dp) :: TotDOS, PDOS
  real(dp) :: TotTrans, Current
  complex(dp), allocatable :: tt(:,:)   ! Transmission matrix
  real(dp), allocatable :: eig(:)
  real(dp), allocatable :: TAv(:), TDOSAv(:)
  real(dp), allocatable :: PDOSAv(:)
  real(dp), allocatable :: TEig(:)      ! total eigentransmissions
  real(dp), allocatable :: TEigAv(:,:)  ! total eigentransmissions

! *****************************************
! * File units                            *
! *****************************************
  ! To write out....
  integer :: uT, uTAv, uTeig, uTeigAv, uDOSL, uDOSR, uIeig
  ! To read in GF files
  integer :: uGFL, uGFR
  ! For the COOP curves
  integer :: uC, uCL, uCR
  ! For the AtomPDOS curves
  integer :: uTOTDOS, uORBDOS

! *****************************************
! * LOCAL variables                       *
! *****************************************
  real(dp) :: EfermiL, EfermiR, Efermi0 ! Fermi levels in the CONTACT
  real(dp) :: dE ! The distance of the energy points
  real(dp) :: f_L, f_R ! Fermi function evaluation of the fermi levels at energies 'contour'
  real(dp) :: kpt(3), wkpt, kpt_Node(3)
  real(dp) :: ucell(3,3) ! The unit cell of the scattering region
  real(dp) :: rTmp, spin_F
  integer :: ierr
  logical :: errorGS
  complex(dp), allocatable :: dummyGAMMA(:,:)
  complex(dp), allocatable :: aux(:,:)
 
! loop counters
  integer :: i,j, ia, ia2,ikpt,iE,ispin

#ifdef MPI
  real(dp), allocatable :: buf_recv(:,:), buf_send(:,:)
  integer :: iNode ! The broadcast node for the Hamiltonian and overlap
  integer :: status(MPI_STATUS_SIZE), req
  integer :: MPIerror
#endif
! *********************************************
! * Start TBTRANS                             *
! *********************************************


! Initialize all options and variables
! This will NOT read in the k-points as we still need the unit cell !
  call tbt_init()

  EfermiL = VoltL
  EfermiR = VoltR
  EFermi0 = 0.0_dp ! (EfermiL + EfermiR)*.5_dp !Ry

! We will not redistribute them as that will be done automatically in the loops
  if ( NEn > 1 ) then
     dE = (Emax-Emin)/real(NEn,dp)
  else
     dE = 0.0_dp
  end if

! Make new line before reading in TSHS
  call out_NEWLINE(6)


! Read in the scattering region H and S
! From these H and S we create all subsequent Hk and Sk
  call tbt_read_tshs(HSFile,no_s,no_u,nspin, &
       ucell, na_u, xa, lasto, &
       maxnh , numh , listhptr , listh , xij , indxuo, &
       H, S, &
       siesta_Gamma, Ef)
  
  ! Write out system information for tbtrans
  if ( IONode ) then
     write(*,'(a)') 'Atomic coordinates and regions (Ang):'
     if ( NBufAtL > 0 ) &
        call out_REGION(6,na_u,xa,1, NBufAtL, 'Left buffer','/')
     call out_REGION(6,na_u,xa,NBufAtL+1,NBufAtL+nuaL_GF*NRepA1L*NRepA2L, &
          'Left electrode','#')
     call out_DEVICE(6,na_u,xa,NBufAtL+nuaL_GF*NRepA1L*NRepA2L+1, &
          na_u-NBufAtR-nuaR_GF*NRepA1R*NRepA2R, &
          IsoAt1,IsoAt2,'*')
     call out_REGION(6,na_u,xa,na_u-NBufAtR-nuaR_GF*NRepA1R*NRepA2R+1, &
          na_u-NBufAtR, 'Right electrode','#')
     if ( NBufAtR > 0 ) &
          call out_REGION(6,na_u,xa,na_u-NBufAtR+1, na_u, 'Right buffer','/')
     write(*,*) ''
  end if

  
! Now we can initialize the k-points....
  call setup_tbt_kpoint_grid(ucell)

! Make new line before Electrode Green's function check
  call out_NEWLINE(6)

! Create the spin-factor in DOS calculations
  if ( nspin == 1 ) then
     spin_F = 2.0_dp
  else
     spin_F = 1.0_dp
  end if


! #########################################################
! # Create the Green's functions for the contour points   #
! #

! Allocate for the bulkdos
  allocate(ZBulkDOS(NEn,nspin))
  call memory('A','Z',NEn*nspin,'tbtrans')
  ZBulkDOS(:,:) = 0.0_dp
  
! Create the Left GF file
  call do_Green('L',HSFileL, GFFileL, GFTitle, &
       ElecValenceBandBot, ReUseGF, &
       nkpnt,kpoint,kweight, &
       NBufAtL,NUsedAtomsL,NRepA1L,NRepA2L, RemUCellDistances, &
       ucell,xa,na_u,NEn,contour,EFermiL,ZBulkDOS,nspin)
  
  ! If we have created the new GF file we can write out the ZBulkDOS
  ! This is the ZBulkDOS for both spins....
  if ( IONode .and. sum(abs(ZBulkDOS(:,:))) > 0.0_dp ) then
     do ispin = 1 , nspin
        call create_file(slabel,'LDOS',ispin,nspin,uDOSL)
        write(uDOSL,'(a)') "# k-point averaged density of states in the &
             &left electrode weighted by the energy point weight"
        write(uDOSL,'(a,a9,tr1,a16)')"#","E [eV]", "DOS"
        do iE = 1 , NEn
           ZBulkDOS(iE,ispin) = -spin_F*r1dPi*dimag(ZBulkDOS(iE,ispin))
        end do
        call out_DOS(uDOSL,NEn,contour(:)%c,ZBulkDOS(:,ispin))
        call io_close(uDOSL)
        ZBulkDOS(:,ispin) = 0.0_dp
     end do
  end if

! Create the Right GF file
  call do_Green('R',HSFileR,GFFileR, GFTitle, &
       ElecValenceBandBot, ReUseGF, &
       nkpnt,kpoint,kweight, &
       NBufAtR,NUsedAtomsR,NRepA1R,NRepA2R, RemUCellDistances, &
       ucell,xa,na_u,NEn,contour,EFermiR,ZBulkDOS,nspin)

  ! If we have created the new GF file we can write out the ZBulkDOS
  ! This is the ZBulkDOS for both spins....
  if ( IONode .and. sum(abs(ZBulkDOS(:,:))) > 0.0_dp ) then
     do ispin = 1 , nspin
        call create_file(slabel,'RDOS',ispin,nspin,uDOSR)
        write(uDOSR,'(a)') "# k-point averaged density of states in the &
             &right electrode weighted by the energy point weight"
        write(uDOSR,'(a,a9,tr1,a16)')"#","E [eV]", "DOS"
        do iE = 1 , NEn
           ZBulkDOS(iE,ispin) = -spin_F*r1dPi*dimag(ZBulkDOS(iE,ispin))
        end do
        call out_DOS(uDOSR,NEn,contour(:)%c,ZBulkDOS(:,ispin))
        call io_close(uDOSR)
     end do
  end if
  
  call memory('D','Z',NEn*nspin,'tbtrans')
  deallocate(ZBulkDOS)

! Make the expansion of the electrodes in atoms and orbitals
  nuaL = nuaL_GF*NRepA1L*NRepA2L
  nuaR = nuaR_GF*NRepA1R*NRepA2R
  noL  = noL_GF *NRepA1L*NRepA2L
  noR  = noR_GF *NRepA1R*NRepA2R

! Print out information in Green's function files
! Show the number of used atoms and orbitals
  if ( IONode ) then
     write(*,'(/,a,i6,'' / '',i6)') &
          'Left : GF atoms    / Expanded atoms    : ',nuaL_GF,nuaL
     write(*,'(a,i6,'' / '',i6)') &
          'Left : GF orbitals / Expanded orbitals : ',noL_GF,noL
     write(*,'(a,i6,'' / '',i6)') &
          'Right: GF atoms    / Expanded atoms    : ',nuaR_GF,nuaR
     write(*,'(a,i6,'' / '',i6,/)') &
          'Right: GF orbitals / Expanded orbitals : ',noR_GF,noR
  end if


! expected no. states on Electrode atoms within the Green's function file: 
! Left
  allocate(lastoL(0:nuaL_GF)) 
  call memory('A','I',nuaL_GF+1,'tbtrans')

  lastoL(0)=0
  ia2=0
  do ia=NBufAtL+1,NBufAtL + nuaL, NRepA1L*NRepA2L
     ia2=ia2+1
     lastoL(ia2)=lastoL(ia2-1) + (lasto(ia) - lasto(ia-1))
  end do ! ia

  if(lastoL(nuaL_GF) .ne. noL_GF) then
     if(IONode) &
          write(*,*) 'ERROR: lastoL,noL_GF',lastoL,noL_GF
     call die('ERROR: Unexpected no. orbs. in L elec.')
  end if

! Right
  allocate(lastoR(0:nuaR_GF))    
  call memory('A','I',nuaR_GF+1,'tbtrans')

  lastoR(0)=0
  ia2=0
  do ia=na_u-(nuaR + NBufAtR)+1,na_u-NBufAtR, NRepA1R*NRepA2R
     ia2=ia2+1
     lastoR(ia2)=lastoR(ia2-1) + (lasto(ia) - lasto(ia-1))
  end do                 !ia
  
  if (lastoR(nuaR_GF) .ne. noR_GF) then
     if(IONode) &
          write(*,*) 'ERROR: lastoR,noR_GF',lastoR,noR_GF
     call die('ERROR: Unexpected no. orbs. in R elec.')
  end if
  
! #
! # Done with creating/ensuring existance of the GF files
! ##########################################################

! We have now completed use of all files.
! We will start by doing initial setup of the system, isolated system, 
! and where which orbitals belong...

  
! ###################################
! # Do buffer atoms
! #

! the first NBufAtL atoms will be removed
  noBufL = 0
  do ia = 1 , NBufAtL
     noBufL = noBufL+(lasto(ia)-lasto(ia-1))
  end do ! ia
! the last NBufAtR atoms will be removed
  noBufR = 0
  do ia = na_u - NBufAtR+1 , na_u
     noBufR = noBufR+(lasto(ia)-lasto(ia-1))
  end do ! ia

! # 
! # Done with finding the buffer atom regions
! #####################################



! No. states minus buffers
  nou = no_u - (noBufL+noBufR)

! Find the number of orbitals in the device region
! This means that it is the total number of orbitals,
! minus the buffers and electrodes
  noD = nou - (noL+noR)

  if ( noD < 1 ) then
     call die('Device orbital count < 0. Please check your input!')
  end if



! ###########################################
! # Find the isolated region in terms of orbitals
! #

! Retrive orbital in the CONTACT
  Isoo1 = lasto(IsoAt1-1) + 1 ! Orbital in full CONTACT
  Isoo2 = lasto(IsoAt2)       ! Orbital in full CONTACT
  Isoo1C = Isoo1 - noBufL     ! Orbital in CONTACT
  Isoo2C = Isoo2 - noBufL     ! Orbital in CONTACT
  Isoo1D = Isoo1C - noL       ! Orbital in DEVICE
  Isoo2D = Isoo2C - noL       ! Orbital in DEVICE
  Isoo = Isoo2 - Isoo1 + 1
  
  if (IOnode) then
     write(*,*) repeat('=',62) 
     write(*,'(a30,i4,a1,i4,a1)') &
          'Projection Region: atoms : [',IsoAt1,';',IsoAt2,']'
     write(*,'(a30,i5,a1,i5,a6,i5)') &
          'Projection Region: states: [',Isoo1,';',Isoo2,'] Tot:' &
          ,Isoo
  end if
! Error check for wrong setup
  if ( Isoo1 > Isoo2 ) then
     call die('Requested isolated region is wrong. &
          &Have you reversed the options?')
  end if

! #
! # The isolated region is now found
! ###########################################



! ###########################################
! # Read-in header of Green's functions
! # Prepare for the calculation
! # We read in the k-points that the electrode was generated with.
! # Furthermore we read in the expansion q-points
! #

  if (IONode) then
     call io_assign(uGFL)
     open(file=GFFileL,unit=uGFL,form='unformatted')
     call io_assign(uGFR)
     open(file=GFFileR,unit=uGFR,form='unformatted')
  end if

! Left (notice that nuaL_GF should equal NUsedAtomsL)
  call read_Green(uGFL,.true.,EfermiL,nkpnt,NEn,nuaL_GF, &
       NRepA1L,NRepA2L,RemUCellDistances, noL_GF,nspin, &
       nkparL,kparL,wkparL, &
       nqL,wqL,qLb)

! Right (notice that nuaR_GF should equal NUsedAtomsR)
  call read_Green(uGFR,.true.,EfermiR,nkpnt,NEn,nuaR_GF, &
       NRepA1R,NRepA2R,RemUCellDistances, noR_GF,nspin, &
       nkparR,kparR,wkparR, &
       nqR,wqR,qRb)


! #
! # done with reading header of GF files
! ###########################################



! ###########################################
! # Allocate used arrays
! #

  allocate(SFEL(noL,noL))
  allocate(HAAL(noL_GF,noL_GF,nqL))
  allocate(SAAL(noL_GF,noL_GF,nqL))
  call memory('A','Z',noL**2+noL_GF**2*nqL*2,'tbtrans')

  allocate(SFER(noR,noR))
  allocate(HAAR(noR_GF,noR_GF,nqR))
  allocate(SAAR(noR_GF,noR_GF,nqR))
  call memory('A','Z',noR**2+noR_GF**2*nqR*2,'tbtrans')

  ! We need a dummy array for passing to the getSFE routine
  ! as a dummy, we do not care of the contents!
  allocate(dummyGAMMA(max(noL,noR),max(noL,noR)))
  call memory('A','Z',max(noL,noR)**2,'tbtrans')

  if ( CalcIeig ) then
     allocate(eig(Isoo))
     call memory('A','D',nou,'tbtrans')
     allocate(aux(Isoo,Isoo))
     call memory('A','Z',Isoo*Isoo,'tbtrans')
  end if


! We have several Hamiltonians
! ... current calculation k-point
! ... This will always be pointed to array...
  nullify(Hk,Sk)
#ifndef MPI
  ! if sequential it needs to be allocated instead of pointing...
  allocate(Hk(nou,nou))
  allocate(Sk(nou,nou))
  call memory('A','Z',nou*nou*2,'tbtrans')
#endif

#ifdef MPI
! ... an array to recv the Hamiltonian over MPI
  nullify(Hk_recv,Sk_recv)
  allocate(Hk_recv(nou,nou))
  allocate(Sk_recv(nou,nou))
  call memory('A','Z',nou*nou*2,'tbtrans')

! ... an array to hold the k-point for future use
  nullify(Hk_Node,Sk_Node)
  allocate(Hk_Node(nou,nou))
  allocate(Sk_Node(nou,nou))
  call memory('A','Z',nou*nou*2,'tbtrans')
#endif

  if ( CalcIeig ) then
! ... an array to hold the k-point for isolated region
     nullify(Hk_iso,Sk_iso)
     allocate(Hk_iso(Isoo,Isoo))
     allocate(Sk_iso(Isoo,Isoo))
     call memory('A','Z',Isoo*Isoo*2,'tbtrans')
  end if

! ... an array to hold the k-point for the DEVICE
  nullify(Hk_D,Sk_D)
!  allocate(Hk_D(noD,noD)) ! Not used yet
!  call memory('A','Z',noD*noD,'tbtrans')
  allocate(Sk_D(noD,noD))
  call memory('A','Z',noD*noD,'tbtrans')

! allocate the Green's functions
  nullify(GF,GFRGF)
  allocate(GF(noD,noD))
  allocate(GFRGF(noD,noD))
  call memory('A','Z',noD*noD*2,'tbtrans')

! allocate the transmission matrix
  allocate(tt(noD,noD),TEig(NeigCh))
  call memory('A','Z',noD*noD,'tbtrans')
  call memory('A','D',NeigCh,'tbtrans')
! Data collecting arrays
  allocate(TAv(PNEn),TEigAv(PNEn,NeigCh))
  call memory('A','D',PNEn+PNEn*NeigCh,'tbtrans')
  allocate(TDOSAv(PNEn),PDOSAv(PNEn))
  call memory('A','D',PNEn*2,'tbtrans')


#ifdef MPI
  ! Allocate the buffers for send etc.
  allocate(buf_send(PNEn,NeigCh+3),buf_recv(PNEn,NeigCh+3))
#endif

! #
! # done with allocation
! ###########################################


! ###########################################
! ######   Transmission calculation    ######
! ##########     starts here      ###########
! ###########################################
  l_spin: do ispin = 1 , nspin
     
     ! Initialize the arrays as we are going to do a
     ! reduction
     TAv(:)      = 0.0_dp
     TEigAv(:,:) = 0.0_dp
     TDOSAv(:)   = 0.0_dp
     PDOSAv(:)   = 0.0_dp

     Current = 0.0_dp
     
     ! Create the files that are needed for the calculation
     call create_file(slabel,'TRANS',ispin,nspin,uT)
     call create_file(slabel,'AVTRANS',ispin,nspin,uTAv)
     if ( IONode ) then ! Write headers...
        write(uT,'(a10,3(tr1,a16))') "#   E [eV]","Trans [G0]","TotDOS","PDOS"
        write(uTAv,'(a)') "# Averaged transmission, total DOS and projected DOS"
        write(uTAv,'(a10,3(tr1,a16))') "#   E [eV]","Trans [G0]","TotDOS","PDOS"
     end if

     if (NEigch > 0) then
        call create_file(slabel,'TEIG',ispin,nspin,uTeig)
        call create_file(slabel,'AVTEIG',ispin,nspin,uTeigAv)
        if ( IONode ) then ! Write headers...
           write(uTeig,'(a10,tr1,a)') "#   E [eV]", &
                "Transmission eigenvalues [G0], descending order"
           write(uTeigAv,'(a)') &
                "# Averaged transmission eigenvalues, descending order"
           write(uTeigAv,'(a10,tr1,a)') "#   E [eV]", &
                "Transmission eigenvalues [G0], descending order"
        end if
     end if

     if ( CalcIeig ) then
        call create_file(slabel,'IEIG',ispin,nspin,uIeig)
        if ( IONode ) then ! Write headers...
           write(uIeig,'(a5,tr1,a)') "#    ", &
                "Hamiltonian Eigenvalues of isolated region [eV]"
        end if
     end if

     
     if ( CalcCOOP ) then
        call create_file(slabel,'COOP',ispin,nspin,uC)
        call create_file(slabel,'COOPL',ispin,nspin,uCL)
        call create_file(slabel,'COOPR',ispin,nspin,uCR)

        if ( IONode ) then ! Write headers...
           write(uC,'(a)') "# COOP between atoms"
           write(uC,'(a,2(a4,tr1),a10,3(tr1,a16))')"#","At1","At2","E [eV]", &
                "Total", "Left", "Right"
           write(uCL,'(a)') "# COOP between atoms and regions"
           write(uCL,'(a,a3,tr1,a10,3(tr1,a16))')"#","At.","E [eV]", &
                "L2L", "L", "L2R"
           write(uCR,'(a)') "# COOP between atoms and regions"
           write(uCR,'(a,a3,tr1,a10,3(tr1,a16))')"#","At.","E [eV]", &
                "R2R", "R", "R2L"
        end if
     end if

     if ( CalcAtomPDOS ) then
        call create_file(slabel,'TOTDOS',ispin,nspin,uTOTDOS)
        call create_file(slabel,'ORBDOS',ispin,nspin,uORBDOS)

        if ( IONode ) then ! Write headers...
           write(uTOTDOS,'(a)')"# Total Mulliken population on atoms"
           write(uTOTDOS,'(a4,tr1,a10,3(tr1,a16))')"# At.","E [eV]", &
                "Total", "Left", "Right"
           write(uORBDOS,'(a)')"# Mulliken population on orbitals"
           write(uORBDOS,'(a)')"# Orbitals are in same order as SIESTA Mulliken"
           write(uORBDOS,'(a4,tr1,a10,tr1,a)')"# At.","E [eV]", &
                "Mulliken pop. on orbital"
        end if
     end if

     l_kpt: do ikpt = 1 , nkpnt
        
        kpt(:) = kpoint(:,ikpt)
        wkpt   = kweight(ikpt)


        ! Write header to transport calculation routines...
        call out_kpt_header(uT,ikpt,kpt,wkpt,ucell)
        ! Also to STDOUT
        call out_kpt_header(6,ikpt,kpt,wkpt,ucell)
        if (NEigch > 0) then
           call out_kpt_header(uTeig,ikpt,kpt,wkpt,ucell)
        end if

        if ( CalcIeig ) then
           call out_kpt_header(uIeig,ikpt,kpt,wkpt,ucell)
        end if

        if ( CalcCOOP ) then
           call out_kpt_header(uC,ikpt,kpt,wkpt,ucell)
           call out_kpt_header(uCL,ikpt,kpt,wkpt,ucell)
           call out_kpt_header(uCR,ikpt,kpt,wkpt,ucell)
        end if

        if ( CalcAtomPDOS ) then
           call out_kpt_header(uTOTDOS,ikpt,kpt,wkpt,ucell)
           call out_kpt_header(uORBDOS,ikpt,kpt,wkpt,ucell)
        end if

        
! Check whether we are in the next loop for the Node loop
! Notice this will always succeed for Nodes == 1
        if ( MOD(ikpt-1,Nodes) == 0 ) then
! Retrieve the node specific k-point
           kpt_Node(:) = kpoint(:,min(ikpt+Node,nkpnt))

#ifdef MPI
           call set_HS_matrix(siesta_Gamma,ucell,na_u,no_u,no_s,maxnh, &
                xij,numh,listhptr,listh,indxuo,H(:,ispin),S, &
                kpt_Node, &
                Hk_Node,Sk_Node, &
                xa=xa, &
                RemUCellDistances=RemUCellDistances,lasto=lasto, &
                RemNFirstOrbitals=noBufL,RemNLastOrbitals=noBufR)
#else
           call set_HS_matrix(siesta_Gamma,ucell,na_u,no_u,no_s,maxnh, &
                xij,numh,listhptr,listh,indxuo,H(:,ispin),S, &
                kpt_Node, &
                Hk,Sk, &
                xa=xa, &
                RemUCellDistances=RemUCellDistances,lasto=lasto, &
                RemNFirstOrbitals=noBufL,RemNLastOrbitals=noBufR)
#endif

           
#ifdef MPI
           call matrix_rem_left_right(nou,Hk_Node,Sk_Node,noL,noR)
           call matrix_symmetrize(nou,Hk_Node,Sk_Node,Ef)
#else
           call matrix_rem_left_right(nou,Hk,Sk,noL,noR)
           call matrix_symmetrize(nou,Hk,Sk,Ef)
#endif

        end if

! Retrieve which node has the current Hamiltonian....
! Furthermore we associate the calculating arrays with the 
! correct array.
! If it is the iNode we simply point to the current [HS]k_Node array.
! Else we point to the recieving array.
#ifdef MPI 
        iNode = MOD(ikpt-1,Nodes)
        if ( Node == iNode ) then
           call MPI_Bcast(Hk_Node(1,1),nou*nou, &
                MPI_double_complex,iNode,MPI_Comm_World,MPIerror)
           call MPI_Bcast(Sk_Node(1,1),nou*nou, &
                MPI_double_complex,iNode,MPI_Comm_World,MPIerror)
           Hk => Hk_Node
           Sk => Sk_Node
        else
           call MPI_Bcast(Hk_recv(1,1),nou*nou, &
                MPI_double_complex,iNode,MPI_Comm_World,MPIerror)
           call MPI_Bcast(Sk_recv(1,1),nou*nou, &
                MPI_double_complex,iNode,MPI_Comm_World,MPIerror)
           Hk => Hk_recv
           Sk => Sk_recv
        end if
#endif 

! Copy over device region overlap, we can use this to faster calculate the COOP curves
! Besides the extra memory should be no problem
        do j = 1 , noD
           do i = 1 , noD
              ! Hk_D(i,j) = Hk(i+noL,j+noL) ! Not used for anything, yet
              Sk_D(i,j) = Sk(i+noL,j+noL)
           end do
        end do

        if ( CalcIeig ) then
! Copy over PDOS region, the Hk and Sk matrices are in sizes with the electrodes
           do j = 1 , Isoo
              do i = 1 , Isoo
                 Hk_iso(i,j) = Hk(i+Isoo1C-1,j+Isoo1C-1)
                 Sk_iso(i,j) = Sk(i+Isoo1C-1,j+Isoo1C-1)
              end do
           end do
           
! Calculate eigenvalues of Isolated region
           call cdiag(Hk_iso,Sk_iso,Isoo,Isoo,Isoo,eig,aux,Isoo,10,errorGS)
           if(errorGS) call die('ERROR in isolated diagonalization.') 
           call out_EIG(uIeig,Isoo,eig)
        end if

! We need to ensure that the full loop is always run by ALL nodes        
        l_E: do iE = Node + 1 , PNEn , Nodes
! Retrieve the node specific energy-point
           ZEnergy = contour(min(iE,NEn))%c
           ZwGF    = contour(min(iE,NEn))%w
           
! LEFT:
           call getSFE(UseBulk,uGFL,HAAL,SAAL,ZEnergy,ikpt, &
                nqL,qLb,wqL, &
                noL_GF, &
                nuaL_GF,lastoL,NRepA1L,NRepA2L, &
                noL,SFEL,dummyGAMMA, &
                min(Nodes,NEn-(iE-Node)+1),errorGS)
           
           if(errorGS) call die('ERROR in getSFE Left') 

! RIGHT:
           call getSFE(UseBulk,uGFR,HAAR,SAAR,ZEnergy,ikpt, &
                nqR,qRb,wqR, &
                noR_GF, &
                nuaR_GF,lastoR,NRepA1R,NRepA2R, &
                noR,SFER,dummyGAMMA, &
                min(Nodes,NEn-(iE-Node)+1),errorGS)

           if(errorGS) call die('ERROR in getSFE Right') 


! Calculate the transmission for the energy point and the k-point
           call transmission(UseBulk,nou,Hk,Sk, &
                noD,noL,SFEL,noR,SFER,ZEnergy, &
                GF,GFRGF,TotTrans,tt,ierr)

! Total the density of states
           TotDOS = 0.0_dp
           do j = 1 , noD
              do i = 1 , noD
                 TotDOS = TotDOS - r1dPi*dimag(dconjg(Sk_D(i,j))*GF(i,j))
              end do
           end do
           TotDOS = spin_F * TotDOS

! Find the "excluded" DOS and subtract from TotDOS
           PDOS = 0.0_dp
           do j = 1 , Isoo1D - 1
              do i = 1 , noD
                 PDOS = PDOS - r1dPi*dimag(Sk_D(i,j)*GF(i,j))
              end do
           end do
           do j = Isoo2D + 1, noD
              do i = 1 , noD
                 PDOS = PDOS - r1dPi*dimag(Sk_D(i,j)*GF(i,j))
              end do
           end do
           PDOS = TotDOS - spin_F * PDOS

           if ( ZwGF == dcmplx(0.0_dp,0.0_dp) ) then
              PDOS     = 0.0_dp
              TotDOS   = 0.0_dp
              TotTrans = 0.0_dp
           end if

           ! Sum up the averages
           TAv(iE)    = TAv(iE)    + TotTrans * wkpt
           TDOSAv(iE) = TDOSAv(iE) + TotDOS   * wkpt
           PDOSAv(iE) = PDOSAv(iE) + PDOS     * wkpt


           ! Do globalization of the TotTrans, TotDOS, PDOS
#ifdef MPI
           buf_send(1,1) = real(ZEnergy,dp)
           buf_send(2,1) = TotTrans
           buf_send(3,1) = TotDOS
           buf_send(4,1) = PDOS
           do iNode = 0 , Nodes - 1
              if ( IONode ) then
                 if ( iNode /= Node ) then
                    call MPI_IRecv(buf_recv(1,1),4,MPI_Double_Precision, &
                         iNode,iNode,MPI_Comm_World,req,MPIerror)
                    call MPI_Wait(req,status,MPIerror)
                 else
                    buf_recv(1:4,1) = buf_send(1:4,1)
                 end if
                 if ( NEn < iNode + iE ) cycle
                 call out_Trans(uT,buf_recv(1,1),buf_recv(2,1), &
                      buf_recv(3,1),buf_recv(4,1))
                 ! Write to STDOUT
                 call out_Trans(6, buf_recv(1,1),buf_recv(2,1), &
                      buf_recv(3,1),buf_recv(4,1))
              else if ( iNode == Node ) then
                 call MPI_ISend(buf_send(1,1),4,MPI_Double_Precision, &
                      0,iNode,MPI_Comm_World,req,MPIerror) 
                 call MPI_Wait(req,status,MPIerror)
              end if
           end do
#else
           call out_Trans(uT,real(ZEnergy,dp),TotTrans,TotDOS,PDOS)
           ! Write to STDOUT
           call out_Trans(6, real(ZEnergy,dp),TotTrans,TotDOS,PDOS)
#endif

           ! Do the COOP curve
           if ( CalcCOOP ) then
              call COOP(uC,uCL,uCR, spin_F, &
                   IsoAt1,IsoAt2, &
                   noBufL,noL,noD, &
                   na_u,lasto,GF,GFRGF,Sk_D, &
                   iE,dreal(ZEnergy))
           end if

           if ( CalcAtomPDOS ) then
              call AtomPDOS(uTOTDOS,uORBDOS,.false.,spin_F, &
                   noBufL+noL,noD, &
                   na_u, IsoAt1, IsoAt2, lasto, &
                   real(ZEnergy,dp),real(ZwGF,dp), &
                   GF, GFRGF, Sk_D)
           end if

           ! Calculate the eigenchannels of the device
           if (NEigch > 0) then
              TEig(:) = 0.0_dp
              if ( ZwGF /= dcmplx(0.0_dp,0.0_dp) ) then
                 call tt_eig(noD,tt,NEigch,TEig)
              end if

              ! Calculate the mean transmission eigenvalues
              do i = 1 , NEigch
                 TEigAv(iE,i) = TEigAv(iE,i)+TEig(i)*wkpt
              end do

#ifdef MPI
              do iNode = 0 , Nodes - 1
                 if ( IONode ) then
                    if ( iNode /= Node ) then
                       call MPI_IRecv(TEig,NeigCh,MPI_Double_Precision, &
                            iNode,iNode,MPI_Comm_World,req,MPIerror)
                       call MPI_Wait(req,status,MPIerror)
                    end if
                    if ( NEn < iNode + iE ) cycle
                    call out_TEIG(uTeig,real(contour(iE+iNode)%c,dp), &
                         NeigCh,TEig)
                 else if ( iNode == Node ) then
                    call MPI_ISend(TEig,NeigCh,MPI_Double_Precision, &
                         0,iNode,MPI_Comm_World,req,MPIerror) 
                    call MPI_Wait(req,status,MPIerror)
                 end if
              end do
#else
              call out_TEIG(uTeig,real(ZEnergy,dp), &
                   NeigCh,TEig)
#endif
           end if

           ! Calculate the actual current
           f_L = fermi(dreal(ZEnergy), EFermiL, kT)
           f_R = fermi(dreal(ZEnergy), EFermiR, kT)

           Current = Current + (f_L-f_R)*TotTrans*dE*wkpt
        end do l_E

        call out_NEWLINE(uT)
        call out_NEWLINE(6)
        if (NEigch > 0) then
           call out_NEWLINE(uTeig)
        end if

        if ( CalcIeig ) then
           call out_NEWLINE(uIeig)
        end if

        if ( CalcCOOP ) then
           call out_NEWLINE(uC)
           call out_NEWLINE(uCL)
           call out_NEWLINE(uCR)
        end if

        if ( CalcAtomPDOS ) then
           call out_NEWLINE(uTOTDOS)
           call out_NEWLINE(uORBDOS)
        end if
        
     end do l_kpt

     ! Do globalization of the Averaged arrays
     ! This is the handling of:
     !  1. TAv
     !  2. TDOSAv
     !  3. PDOSAv
     !  4: TEigAv
#ifdef MPI
     buf_send(:,1) = TAv(:)
     buf_send(:,2) = TDOSAv(:)
     buf_send(:,3) = PDOSAv(:)
     if ( NeigCh > 0 ) then
        buf_send(:,4:3+NeigCh) = TEigAv(:,:)
     end if
     call MPI_Reduce(buf_send(1,1),buf_recv(1,1),size(buf_recv), &
          MPI_Double_Precision, MPI_SUM,0,MPI_Comm_World,MPIerror)
     do iE = 1 , NEn
        call out_Trans(uTAv,real(contour(iE)%c,dp),buf_recv(iE,1), &
             buf_recv(iE,2),buf_recv(iE,3))
        if ( NeigCh > 0 ) then
           call out_TEIG(uTeigAv,real(contour(iE)%c,dp),NeigCh,buf_recv(iE,4:3+NeigCh))
        end if
     end do
#else
     do iE = 1 , NEn
        call out_Trans(uTAv,real(contour(iE)%c,dp),TAv(iE), &
             TDOSAv(iE),PDOSAv(iE))
        if ( NeigCh > 0 ) then
           call out_TEIG(uTeigAv,real(contour(iE)%c,dp),NeigCh,TEigAv)
        end if
     end do
#endif

#ifdef MPI
     rTmp = 0.0_dp
     call MPI_Reduce(Current,rTmp,1,MPI_double_precision, &
          MPI_sum,0,MPI_Comm_World,MPIerror)
     Current = rTmp
#endif

     ! Calculate the current in Ampere
     Current = Current/eV*38.73_dp*1E-06

! Create the files that are needed for the calculation
     if ( IONode ) then
        call io_close(uT)
        call io_close(uTAv)
        if ( NeigCh > 0 ) then
           call io_close(uTeig)
           call io_close(uTeigAv)
        end if

        if ( CalcIeig ) then
           call io_close(uIeig)
        end if
        
        if ( CalcCOOP ) then
           call io_close(uC)
           call io_close(uCL)
           call io_close(uCR)
        end if

        if ( CalcAtomPDOS ) then
           call io_close(uTOTDOS)
           call io_close(uORBDOS)
        end if

        write(*,*) repeat('=',62)
        write(*,*) 'Results:'
        write(*,'(1x,a,2ES16.8)') 'Voltage, Current(A) = ',VoltFDF/eV,Current
        write(*,*) repeat('=',62)
        
     end if

  end do l_spin

! ###########################################
! ######   Transmission calculation    ######
! ##########      ends here       ###########
! ###########################################
  
  if ( CalcIeig ) then
     call memory('D','D',Isoo,'tbtrans')
     deallocate(eig)
     call memory('D','Z',Isoo*Isoo,'tbtrans')
     deallocate(aux)
  end if

  call memory('D','I',nuaL_GF+nuaR_GF+2,'tbtrans')
  deallocate(lastoL,lastoR)

  call memory('D','Z',noL**2+noL_GF**2*nqL*2,'tbtrans')
  deallocate(SFEL,HAAL,SAAL)

  call memory('D','Z',noR**2+noR_GF**2*nqR*2,'tbtrans')
  deallocate(SFER,HAAR,SAAR)

#ifndef MPI
  call memory('D','Z',nou*nou*2,'tbtrans')
  deallocate(Hk,Sk)
#endif

#ifdef MPI 
  call memory('D','Z',nou*nou*2,'tbtrans')
  deallocate(Hk_recv,Sk_recv)

  call memory('D','Z',nou*nou*2,'tbtrans')
  deallocate(Hk_Node,Sk_Node)
#endif

  if ( CalcIeig ) then
     call memory('D','Z',Isoo*Isoo*2,'tbtrans')
     deallocate(Hk_iso,Sk_iso)
  end if

  call memory('D','Z',noD*noD,'tbtrans')
  deallocate(Sk_D)

  call memory('D','Z',noD*noD*2,'tbtrans')
  deallocate(GF,GFRGF)

  call memory('D','Z',noD*noD,'tbtrans')
  deallocate(tt)
  call memory('D','D',NeigCh,'tbtrans')
  deallocate(TEig)

  call memory('D','D',PNEn+PNEn*NeigCh,'tbtrans')
  deallocate(TAv,TEigAv)

  call memory('D','D',PNEn*2,'tbtrans')
  deallocate(TDOSAv,PDOSAv)


! Stop time counter
  call timer( 'tbtrans', 2 )
  call timer( 'all', 3 )

! Print allocation report
  call alloc_report( printNow=.true. )

  if (IONode) then
     call timestamp("End of run")
     call wallclock("End of run")
  end if

#ifdef MPI
  call MPI_Finalize(MPIerror)
#endif

contains

  pure function fermi(e,ef,T)
    real(dp), intent(in) :: e,ef,T
    real(dp) :: tmp,fermi

    if ( T == 0.0_dp ) then
       if ( e>ef ) fermi = 0.0_dp
       if ( e<ef ) fermi = 1.0_dp
       if ( e.eq.ef ) fermi = 0.5_dp
    else
       tmp = (e - ef) / T
       if (abs(tmp).lt.40.0_dp) then
          fermi = 1.0_dp/(1.+exp(tmp))
       else
          if (tmp.lt.-40.0_dp) fermi = 1.0_dp
          if (tmp.gt. 40.0_dp) fermi = 0.0_dp
       endif
    endif

  end function fermi

end program tbtrans
