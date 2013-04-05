module m_ts_mem

  use precision, only : dp
  
  implicit none
  
  public :: transiesta_mem
  
  private
  
contains
  
! ##################################################################
! ##                                                              ##       
! ##                       "TRANSIESTA"                           ##
! ##                                                              ##       
! ##          Non-equilibrium Density Matrix Subroutine           ##
! ##                   to be called from SIESTA                   ##
! ##                                                              ## 
! ## Originally:                By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##               Kurt Stokbro, ks@mic.dtu.dk                    ## 
! ##               Mikroelektronik Centret (MIC)                  ##
! ##           Technical University of Denmark (DTU)              ##
! ##                                                              ##
! ## Currently:                 By                                ##
! ##           Nick Papior Andersen, nickpapior@gmail.com         ##
! ##           Technical University of Denmark (DTU)              ##
! ##                                                              ##
! ## This code has been fully re-created to conform with the      ##
! ## sparsity patterns in SIESTA. Thus the memory requirements    ##
! ## has been greatly reduced.                                    ##
! ##                                                              ##
! ##################################################################
!
!
! Tight-binding density-matrix/transport program for the SIESTA
! package.
! Copyright by Mads Brandbyge, 1999, 2000, 2001, 2002.
! Copyright by Nick Papior Andersen, 2013.
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the authors.
!

  subroutine transiesta_mem(nspin, &
       Gamma, sp_dist, sparse_pattern, &
       ucell, no_u, na_u, lasto, xa, n_nzs, &
       xij, Hs, Ss, DM, EDM, Ef, &
       TSiscf, Qtot)

    use units, only : Pi
    use parallel, only : Node, Nodes, IONode
#ifdef MPI
    use mpi_siesta
#endif

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpArr1D
    use class_zSpArr1D

    use m_ts_kpoints

    use m_ts_options, only : IsVolt, UseBulk, UpdateDMCR
    use m_ts_options, only : VoltL, VoltR
    use m_ts_options, only : NRepA1L, NRepA2L
    use m_ts_options, only : NRepA1R, NRepA2R
    use m_ts_options, only : GFFileL, GFFileR
    use m_ts_options, only : na_BufL => NBufAtL
    use m_ts_options, only : na_BufR => NBufAtR
    use m_ts_options, only : na_L_HS => NUsedAtomsL
    use m_ts_options, only : no_L_HS => NUsedOrbsL
    use m_ts_options, only : na_R_HS => NUsedAtomsR
    use m_ts_options, only : no_R_HS => NUsedOrbsR


    use m_ts_mem_sparsity, only : ts_sp_uc
    use m_ts_mem_sparsity, only : tsup_sp_uc
    use m_ts_mem_sparsity, only : GF_INV_EQUI_PART
    use m_ts_mem_scat

    use m_ts_contour,only : PNEn, NEn, contour
    use m_ts_cctype, only : CC_PART_EQUI
    use m_ts_cctype, only : CC_PART_LEFT_EQUI
    use m_ts_cctype, only : CC_PART_RIGHT_EQUI
    use m_ts_cctype, only : CC_PART_NON_EQUI
    use m_ts_cctype, only : CC_PART_TRANSPORT

    use m_ts_gf, only : read_Green

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in) :: nspin
    logical, intent(in) :: Gamma
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: no_u,na_u
    integer, intent(in) :: lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: n_nzs
    real(dp), intent(in) :: xij(3,n_nzs)
    real(dp), intent(in) :: Hs(n_nzs,nspin), Ss(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Ef, Qtot
    integer, intent(in) :: TSiscf

! ******************* Regional sizes *************************
! * Buffer regions
    integer :: no_BufL, no_BufR
! * Electrode regions
    integer :: na_L, no_L, na_R, no_R
    integer, allocatable :: lasto_L(:), lasto_R(:)
! * Computational region..
    integer :: no_u_TS, no_C_L, no_C_R
! * In case of using part of the GF inversion
    integer :: no_u_C, no_GF_offset
! ************************************************************

! ******************** IO descriptors ************************
    integer :: uGFL, uGFR
! ************************************************************

! ****************** Electrode variables *********************
    integer :: nqL, nqR, nkparL, nkparR
    real(dp), allocatable :: qLb(:,:), qRb(:,:)
    real(dp), allocatable :: wqL(:), wqR(:)
    real(dp), allocatable :: kparL(:,:), kparR(:,:)
    real(dp), allocatable :: wkparL(:), wkparR(:)
    complex(dp), allocatable :: HAAL(:,:,:), SAAL(:,:,:), GAAL(:,:,:)
    complex(dp), allocatable :: HAAR(:,:,:), SAAR(:,:,:), GAAR(:,:,:)
    complex(dp), allocatable :: SigmaL(:,:), SigmaR(:,:)
    real(dp), allocatable :: GammaL(:,:), GammaR(:,:)
! ************************************************************

! ******************* Computational arrays *******************
    integer :: ndwork, nzwork
    real(dp),    allocatable :: dwork(:)
    complex(dp), allocatable :: zwork(:), GF(:)
    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(dSpArr1D) :: spH, spS
    type(zSpArr1D) :: spzH, spzS
    ! The different sparse matrices...
    type(dSpArr1D) :: spDM, spEDM, spDMR, spEDMR, spDMneqL, spDMneqR
    type(zSpArr1D) :: spzDM, spzEDM, spzDMR, spzEDMR, spzDMneqL, spzDMneqR
    ! Pointers for updating the density matrices
    real(dp), pointer :: dDM(:), dEDM(:)
    complex(dp), pointer :: zDM(:), zEDM(:)
! ************************************************************

! ******************* Computational variables ****************
    complex(dp) :: Z, W, ZW
    real(dp) :: k(3)
! ************************************************************

! ******************** Loop variables ************************
    integer :: ispin, ikpt, iPE, iE, NEReqs, up_nzs, ia, ia_E
! ************************************************************

! ******************* Miscalleneous variables ****************
    integer :: ierr
! ************************************************************


#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    call timer('TS_calc',1)
    
    ! Calculate the number of used atoms in left/right
    na_L = na_L_HS * NRepA1L * NRepA2L
    no_L = no_L_HS * NRepA1L * NRepA2L
    na_R = na_R_HS * NRepA1R * NRepA2R
    no_R = no_R_HS * NRepA1R * NRepA2R

    ! Calculate the number of orbitals not used (i.e. those 
    ! in the buffer regions)
    ! Left has the first atoms
    no_BufL = lasto(na_BufL)
    ! Right has the last atoms
    no_BufR = lasto(na_u) - lasto(na_u - na_BufR)

    ! Create the lasto pointers for the electrode expansions...
    allocate(lasto_L(0:na_L_HS))
    lasto_L(0) = 0
    ia_E = 0
    do ia = na_BufL + 1 , na_BufL + na_L, NRepA1L * NRepA2L
       ia_E = ia_E + 1
       lasto_L(ia_E) = lasto_L(ia_E-1) + lasto(ia) - lasto(ia-1)
    end do
    allocate(lasto_R(0:na_R_HS))
    lasto_R(0) = 0
    ia_E = 0
    do ia = na_u - na_R - na_BufR + 1 , na_u - na_BufR , NRepA1R * NRepA2R
       ia_E = ia_E + 1
       lasto_R(ia_E) = lasto_R(ia_E-1) + lasto(ia) - lasto(ia-1)
    end do
    ! Add to memory-management...
    call memory('A','I',na_R_HS+na_L_HS+2,'transiesta')

    ! Number of orbitals in TranSIESTA
    no_u_TS = no_u - no_BufL - no_BufR

    ! The SIESTA equivalent orbital which
    ! -- starts in the left central region
    no_C_L = 1 + no_BufL + no_L
    ! -- ends in the right central region
    no_C_R = no_u - no_R - no_BufL


    ! Number of elements that are transiesta updated
    up_nzs = nnzs(tsup_sp_uc)

    ! Note that:
    ! no_BufL + no_L + no_C + no_R + no_BufR == no_u
    ! no_BufL +      no_C_TS       + no_BufR == no_u
    ! no_BufL + no_L + (no_C_R - no_C_L + 1) + no_R + no_BufR == no_u

    ! Do a crude check of the sizes
    ! if the transiesta region is equal of size to or smaller 
    ! than the size of the combined electrodes, then the system
    ! is VERY WRONG...
    if ( no_u_TS <= no_L + no_R ) then
       call die("The contact region size is &
            &smaller than the electrode size. &
            &What have you done? Please correct this insanity...")
    end if


    ! Open GF files...
    if ( IONode ) then
       call io_assign(uGFL)
       open(file=GFFileL,unit=uGFL,form='unformatted')
       call io_assign(uGFR)
       open(file=GFFileR,unit=uGFR,form='unformatted')
    end if

! Read-in header of Green's functions
! Prepare for the calculation
! We read in the k-points that the electrode was generated with.
! Furthermore we read in the expansion q-points
! They are communicated in the routine

! Read in the headers of the surface-Green's function files...
! Left
    call read_Green(uGFL,TSiscf==1,VoltL,ts_nkpnt,NEn,na_L_HS,  &
         NRepA1L,NRepA2L,.false.,no_L_HS,nspin, &
         nkparL,kparL,wkparL, &
         nqL,wqL,qLb)

! Right
    call read_Green(uGFR,TSiscf==1,VoltR,ts_nkpnt,NEn,na_R_HS, &
         NRepA1R,NRepA2R,.false.,no_R_HS,nspin,  &
         nkparR,kparR,wkparR, &
         nqR,wqR,qRb)


    ! We do need the full GF AND a single work array to handle the
    ! left-hand side of the inversion...
    ! We will provide all work arrays as single dimension arrays.
    ! This will make interfaces more stringent and allow for
    ! re-use in several other places.
    ! However, this comes at the cost of programmer book-keeping.
    ! It should be expected that the work arrays return GARBAGE
    ! on ALL routines, i.e. they are not used for anything other
    ! than, yes, work.

    ! We need to allocate the important arrays first...
    ! This will (most likely) ensure them to be placed on the 
    ! fast memory slots, whereas allocating in the end will
    ! mean swap-space... We need speed for the inversion....
    if ( ts_Gamma_SCF ) then
       ndwork = max(nnzs(ts_sp_uc),nnzs(tsup_sp_uc))
       allocate(dwork(ndwork),stat=ierr)
       if (ierr/=0) call die('Could not allocate space for dwork')
       call memory('A','D',ndwork,'transiesta')
    end if
    ! The zwork is needed to construct the LHS for solving: G^{-1} G = I
    ! Hence, we will minimum require this...
    nzwork = no_u_TS**2
    allocate(zwork(nzwork),stat=ierr)
    if (ierr/=0) call die('Could not allocate space for zwork')
    call memory('A','Z',nzwork,'transiesta')
    if ( GF_INV_EQUI_PART .and. .not. IsVolt ) then
       allocate(GF(no_u_TS*(no_u_TS-no_R-no_L)),stat=ierr)
       if (ierr/=0) call die('Could not allocate space for GFpart')
       call memory('A','Z',no_u_TS*(no_u_TS-no_R-no_L),'transiesta')
    else
       allocate(GF(no_u_TS**2),stat=ierr)
       if (ierr/=0) call die('Could not allocate space for GF')
       call memory('A','Z',no_u_TS**2,'transiesta')
    end if

    ! Allocate the left-right electrode quantities that we need
    allocate(HAAL(no_L_HS,no_L_HS,NRepA1L*NRepA2L))
    allocate(SAAL(no_L_HS,no_L_HS,NRepA1L*NRepA2L))
    allocate(GAAL(no_L_HS,no_L_HS,NRepA1L*NRepA2L))
    ispin = no_L_HS**2*NRepA1L*NRepA2L*3
    allocate(SigmaL(no_L,no_L))
    ispin = ispin + no_L**2
    allocate(HAAR(no_R_HS,no_R_HS,NRepA1R*NRepA2R))
    allocate(SAAR(no_R_HS,no_R_HS,NRepA1R*NRepA2R))
    allocate(GAAR(no_R_HS,no_R_HS,NRepA1R*NRepA2R))
    ispin = ispin + no_R_HS**2*NRepA1R*NRepA2R*3
    allocate(SigmaR(no_R,no_R))
    ispin = ispin + no_R**2
    call memory('A','Z',ispin,'transiesta')

    if ( IsVolt ) then
       ! We only need Gamma's with voltages
       allocate(GammaL(no_L,no_L),GammaR(no_R,no_R))
       call memory('A','D',no_L**2+no_R**2,'transiesta')
    end if

    ispin = 0

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
    ! Notice that we DO need it to be the SIESTA size.
#ifdef MPI
    call newDistribution(no_u,MPI_COMM_WORLD,fdist,name='TS-fake dist')
#else
    call newDistribution(no_u,-1            ,fdist,name='TS-fake dist')
#endif

    if ( ts_Gamma_SCF ) then
       call newdSpArr1D(tsup_sp_uc,fdist,spDM,name='Transiesta spDM')
       call newdSpArr1D(tsup_sp_uc,fdist,spEDM,name='Transiesta spEDM')
       if ( IsVolt ) then
          call newdSpArr1D(tsup_sp_uc,fdist,spDMR,name='Transiesta spDM-R')
          call newdSpArr1D(tsup_sp_uc,fdist,spDMneqL,name='Transiesta spDMneq-L')
          call newdSpArr1D(tsup_sp_uc,fdist,spDMneqR,name='Transiesta spDMneq-R')
          call newdSpArr1D(tsup_sp_uc,fdist,spEDMR,name='Transiesta spEDM-R')
       end if

       ! The Hamiltonian and overlap matrices (in Gamma calculations
       ! we will not have any phases, hence, it makes no sense to
       ! have the arrays in complex)
       call newdSpArr1D(ts_sp_uc,fdist,spH,name='Transiesta spH')
       call newdSpArr1D(ts_sp_uc,fdist,spS,name='Transiesta spS')

    else
       call newzSpArr1D(tsup_sp_uc,fdist,spzDM,name='Transiesta spzDM')
       call newzSpArr1D(tsup_sp_uc,fdist,spzEDM,name='Transiesta spzEDM')
       if ( IsVolt ) then
          call newzSpArr1D(tsup_sp_uc,fdist,spzDMR,name='Transiesta spzDM-R')
          call newzSpArr1D(tsup_sp_uc,fdist,spzDMneqL,name='Transiesta spzDMneq-L')
          call newzSpArr1D(tsup_sp_uc,fdist,spzDMneqR,name='Transiesta spzDMneq-R')
          call newzSpArr1D(tsup_sp_uc,fdist,spzEDMR,name='Transiesta spzEDM-R')
       end if

       ! The Hamiltonian and overlap matrices
       call newzSpArr1D(ts_sp_uc,fdist,spzH,name='Transiesta spzH')
       call newzSpArr1D(ts_sp_uc,fdist,spzS,name='Transiesta spzS')

    end if
    ! We will not write out all created sparsity patterns, it provides
    ! no purpose... other than displaying how much memory this reduces :)

    SPIN: do ispin = 1 , nspin
       
       ! This is going to get messy...
       ! we do not want to create, yet another sparsity pattern
       ! Hence we need to do the SAME checks, again and again and again....
       ! However, the extra computation should be negligible to the gain.
       
       call init_DM(sp_dist,sparse_pattern, &
            n_nzs, DM(:,ispin), EDM(:,ispin), &
            tsup_sp_uc)


! we wish to loop over the large k-points... 
!     other sub calls that kpoint is the correct array...
    KPOINT: DO ikpt = 1 , ts_nkpnt
       
       k(:) = ts_kpoint(:,ikpt)
       
       ! We initialize the updated region of DM arrays
       if ( ts_Gamma_SCF ) then
          call init_dSpArr1D(spDM)
          call init_dSpArr1D(spEDM)
          if ( IsVolt ) then
             call init_dSpArr1D(spDMR)
             call init_dSpArr1D(spDMneqL)
             call init_dSpArr1D(spDMneqR)
             call init_dSpArr1D(spEDMR)
          end if
       else
          call init_zSpArr1D(spzDM)
          call init_zSpArr1D(spzEDM)
          if ( IsVolt ) then
             call init_zSpArr1D(spzDMR)
             call init_zSpArr1D(spzDMneqL)
             call init_zSpArr1D(spzDMneqR)
             call init_zSpArr1D(spzEDMR)
          end if
       end if
       ! All of the above referenced arrays are now ZERO

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',1)
#endif

       ! Work-arrays are for MPI distribution...
       if ( ts_Gamma_SCF ) then
          call create_HS_Gamma(sp_dist,sparse_pattern, &
               Ef, &
               no_BufL, no_BufR, & ! cut-out region
               no_C_L, no_C_R, &   ! central region (junction)
               no_u, &             ! SIESTA size
               n_nzs, Hs(:,ispin), Ss, &
               spH, spS, &
               ndwork, dwork)
       else
          call create_HS_kpt(sp_dist,sparse_pattern, &
               Ef, &
               no_BufL, no_BufR, & ! cut-out region
               no_C_L, no_C_R, &   ! central region (junction)
               no_u, &             ! SIESTA size
               n_nzs, Hs(:,ispin), Ss, &
               xij, &
               spzH, spzS, k, &
               nzwork, zwork)
       end if

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',2)
#endif

       ! Energy point loop (note that this will loop in steps of the 
       ! nodes (PNEn takes into account the extra "filling" for the
       ! last step)
       EPOINTS: do iPE = Node + 1 , PNEn , Nodes

          ! obtain a valid energy point (truncate at NEn)
          iE = min(iPE,NEn)

          ! save the current weight of the point
          ! This is where we include the factor-of-two for spin and
          ! and the (1/Pi) from DM = Im[G]/Pi
          ! Furthermore we include the weight of the k-point
          W = 1._dp/Pi*contour(iE)%w * ts_kweight(ikpt)
          if ( nspin == 1 ) W = W * 2._dp

          ! save the contour energy point
          Z = contour(iE)%c
          ! Save Z*W, used for E-arrays
          ZW = Z*W

          ! the number of points we wish to read in this segment
          NEReqs = min(Nodes, NEn-(iPe-1-Node))

          ! TODO Move reading of the energy points
          ! directly into the subroutines which need them
          ! In this way we can save both GAA, Sigma AND Gamma arrays!!!!
          ! However, this will probably come at the expense 
          ! of doing the same "repetition" expansion twice, we can live with
          ! that!

#ifdef TRANSIESTA_TIMING
          call timer('TS_READ',1)
#endif

          ! Read in the left electrode
          call read_next_GS(uGFL, NEReqs, &
               ikpt,no_L_HS,nqL, HAAL, SAAL, &
               GAAL, Z, nzwork, zwork)

          ! Read in the right electrode
          call read_next_GS(uGFR, NEReqs, &
               ikpt,no_R_HS,nqR, HAAR, SAAR, &
               GAAR, Z, nzwork, zwork)

#ifdef TRANSIESTA_TIMING
          call timer('TS_READ',2)
#endif

          ! We only need to do a last communication within 
          ! the above reads. Hence we can quit the energy point loop now!
          if ( iPE > NEn ) cycle

#ifdef TRANSIESTA_DEBUG
          write(*,*)'Before calculation of the GF'
#endif

          select case ( contour(iE)%part )
          case ( CC_PART_EQUI , CC_PART_LEFT_EQUI, CC_PART_RIGHT_EQUI ) 

#ifdef TRANSIESTA_TIMING
             call timer('TS_EXPAND',1)
#endif

             ! for these contour parts we do not require to calculate
             ! Gamma's.
             ! Hence we can perform the calculation without 
             ! calculating them.

             ! Calculate the left-right Sigma
             if ( UseBulk ) then

                call UC_expansion_Sigma_Bulk(no_L_HS, no_L, NRepA1L, NRepA2L, &
                     na_L_HS, lasto_L, nqL, qLb, wqL, HAAL, SAAL, GAAL, SigmaL, &
                     nzwork,zwork)
                
                call UC_expansion_Sigma_Bulk(no_R_HS, no_R, NRepA1R, NRepA2R, &
                     na_R_HS, lasto_R, nqR, qRb, wqR, HAAR, SAAR, GAAR, SigmaR, &
                     nzwork,zwork)

             else
                call UC_expansion_Sigma(Z,no_L_HS, no_L,NRepA1L, NRepA2L, &
                     na_L_HS, lasto_L, nqL, qLb, wqL, HAAL, SAAL, GAAL, SigmaL, &
                     nzwork,zwork)

                call UC_expansion_Sigma(Z,no_R_HS, no_R, NRepA1R, NRepA2R, &
                     na_R_HS, lasto_R, nqR, qRb, wqR, HAAR, SAAR, GAAR, SigmaR, &
                     nzwork,zwork)
             end if

#ifdef TRANSIESTA_TIMING
             call timer('TS_EXPAND',2)
             call timer('TS_PREPG',1)
#endif

             if ( ts_Gamma_SCF ) then
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_D(spH , spS,Z,no_BufL,no_u_TS,zwork)
             else
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_Z(spzH,spzS,Z,no_BufL,no_u_TS,zwork)
             end if

#ifdef TRANSIESTA_TIMING
             call timer('TS_PREPG',2)
#endif

             if ( GF_INV_EQUI_PART ) then
                ! Calculate the GF
                call calc_GF_Part(no_u_TS,no_L,no_R, &
                     SigmaL, SigmaR, zwork, GF,ierr)
                no_GF_offset = no_L
               ! The size of the central region (without left-right electrodes)
                no_u_C = no_u_TS - no_R - no_L
             else
                ! Calculate the full GF
                call calc_GF(UseBulk, &
                     no_u_TS,no_L,no_R, &
                     SigmaL, SigmaR, zwork, GF,ierr)
                no_GF_offset = 0
                ! The size of the central region (with left-right electrodes)
                no_u_C = no_u_TS
             end if

             if ( contour(iE)%part == CC_PART_RIGHT_EQUI ) then
                ! We have the right equilibrium contour
                if ( ts_Gamma_SCF ) then
                   call add_DM_dE_D(spDMR , spEDMR , no_u_TS,no_u_C, &
                        GF, no_BufL, no_GF_offset, W, ZW)
                else
                   call add_DM_dE_Z(spzDMR, spzEDMR, no_u_TS,no_u_C, &
                        GF, no_BufL, no_GF_offset, W, ZW)
                end if
             else
                ! We have the left- or the equilibrium contour...
                if ( ts_Gamma_SCF ) then
                   call add_DM_dE_D(spDM ,  spEDM, no_u_TS, no_u_C, &
                        GF, no_BufL, no_GF_offset, W, ZW)
                else
                   call add_DM_dE_Z(spzDM, spzEDM, no_u_TS, no_u_C, &
                        GF, no_BufL, no_GF_offset, W, ZW)
                end if
             end if

          case ( CC_PART_NON_EQUI )
             
#ifdef TRANSIESTA_TIMING
             call timer('TS_EXPAND',1)
#endif

             ! TODO, we can actually do without the Gamma-calculations
             ! in this step. After having calculated GF, we could
             ! calculate the Gamma's and save them in Sigma
             ! This would remove no_L**2 + no_R**2
             ! and should be quite easy. Do this after full
             ! implementation... (note that it will require more computation
             ! in the case of repetition.

             ! Do the left electrode
             call UC_expansion_Sigma_Gamma(UseBulk,Z,no_L_HS,no_L, &
                  NRepA1L, NRepA2L, &
                  na_L_HS,lasto_L,nqL,qLb,wqL, &
                  HAAL, SAAL, GAAL, &
                  SigmaL, GammaL, & 
                  nzwork,zwork)

             ! Do the right electrode
             call UC_expansion_Sigma_Gamma(UseBulk,Z,no_R_HS,no_R, &
                  NRepA1R, NRepA2R, &
                  na_R_HS,lasto_R,nqR,qRb,wqR, &
                  HAAR, SAAR, GAAR, &
                  SigmaR, GammaR, & 
                  nzwork,zwork)

#ifdef TRANSIESTA_TIMING
             call timer('TS_EXPAND',2)
             call timer('TS_PREPG',1)
#endif

             if ( ts_Gamma_SCF ) then
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_D(spH ,spS ,Z,no_BufL,no_u_TS,zwork)
             else
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_Z(spzH,spzS,Z,no_BufL,no_u_TS,zwork)
             end if

#ifdef TRANSIESTA_TIMING
             call timer('TS_PREPG',2)
#endif

             ! Calculate the Greens function
             call calc_GF(UseBulk, &
                     no_u_TS,no_L,no_R, &
                     SigmaL, SigmaR, zwork, GF,ierr)

             ! We calculate the right thing.
             call GF_Gamma_GF(.FALSE., no_u_TS, no_R, GF, &
                  GammaR, zwork, no_u_TS, SigmaR) ! SigmaR is a "work" array

             ! Note that we use '-' here
             if ( ts_Gamma_SCF ) then
                call add_DM_dE_D(spDMneqR,spEDM, no_u_TS, no_u_TS, &
                     zwork, no_BufL, 0, -W, -ZW)
             else
                call add_DM_dE_Z(spzDMneqR,spzEDM, no_u_TS, no_u_TS, &
                     zwork, no_BufL, 0, -W, -ZW)
             end if

             ! We calculate the left thing.
             call GF_Gamma_GF(.TRUE., no_u_TS, no_L, GF, &
                  GammaL, zwork, no_u_TS, SigmaL) ! SigmaL is a "work" array

             ! Note that we use '-' here
             if ( ts_Gamma_SCF ) then
                call add_DM_dE_D(spDMneqL,spEDMR, no_u_TS, no_u_TS, &
                     zwork, no_BufL, 0, W, -ZW)
             else
                call add_DM_dE_Z(spzDMneqL,spzEDMR, no_u_TS, no_u_TS, &
                     zwork, no_BufL, 0, W, -ZW)
             end if

          !case ( CC_PART_TRANSPORT )
          ! For now this is commented out as it isn't implemented
          !   < we should do something special here... >

          case default
             call die('transiesta: ERROR in contour setup')
          end select

#ifdef TRANSIESTA_DEBUG
          write(*,*) 'Completed energy-point...'
#endif

       end do EPOINTS

#ifdef MPI
! Global reduction of density matrices
! this completes the energy integration of each DM
      call timer("TS_comm",1)
      
      if ( ts_Gamma_SCF ) then
         call AllReduce_dSpArr1D(spDM ,ndwork,dwork)
         call AllReduce_dSpArr1D(spEDM,ndwork,dwork)
         if ( IsVolt ) then
            call AllReduce_dSpArr1D(spDMR   ,ndwork,dwork)
            call AllReduce_dSpArr1D(spDMneqL,ndwork,dwork)
            call AllReduce_dSpArr1D(spDMneqR,ndwork,dwork)
            call AllReduce_dSpArr1D(spEDMR  ,ndwork,dwork)
         end if
      else
         call AllReduce_zSpArr1D(spzDM ,nzwork,zwork)
         call AllReduce_zSpArr1D(spzEDM,nzwork,zwork)
         if ( IsVolt ) then
            call AllReduce_zSpArr1D(spzDMR   ,nzwork,zwork)
            call AllReduce_zSpArr1D(spzDMneqL,nzwork,zwork)
            call AllReduce_zSpArr1D(spzDMneqR,nzwork,zwork)
            call AllReduce_zSpArr1D(spzEDMR  ,nzwork,zwork)
         end if
      end if

      call timer("TS_comm",2)
#endif

#ifdef TRANSIESTA_DEBUG
      write(*,*)'Completed energy integration'
#endif

#ifdef TRANSIESTA_TIMING
      call timer('TS_UPDM',1)
#endif

      if ( IsVolt ) then
         if ( ts_Gamma_SCF ) then
            call weightDM(no_C_L, no_C_R, &
                 spDM,   spDMR,  spDMneqL,  spDMneqR, &
                 spEDM,  spEDMR )
         else
            call weightDMC(no_C_L, no_C_R, &
                 spzDM, spzDMR, spzDMneqL, spzDMneqR, &
                 spzEDM, spzEDMR)
         end if
      end if

      ! The original Hamiltonian from SIESTA was shifted Ef: 
      ! Hence we need to shift EDM 
      ! TODO consider moving this to the corresponding 
      ! update region (i.e. copy of the values to the corresponding
      ! local points...
      if ( ts_Gamma_SCF ) then
         ia   =  nnzs(spDM)
         dDM  => val(spDM)
         dEDM => val(spEDM)
         call daxpy(ia,Ef,dDM,1,dEDM,1)

         ! Directly save to the correct DM
         call update_DM(sp_dist,sparse_pattern, n_nzs, &
              DM(:,ispin), EDM(:,ispin), spDM, spEDM)

      else
         Z = dcmplx(Ef,0._dp)
         ia   =  nnzs(spzDM)
         zDM  => val(spzDM)
         zEDM => val(spzEDM)
         call zaxpy(ia,Z,zDM,1,zEDM,1)

         ! Directly save to the correct DM
         call update_zDM(sp_dist,sparse_pattern, n_nzs, &
              DM(:,ispin), EDM(:,ispin), xij, spzDM, spzEDM, k)
      end if

#ifdef TRANSIESTA_TIMING
      call timer('TS_UPDM',2)
#endif

   end do KPOINT

   ! We don't need to do anything here..
   end do SPIN

#ifdef TRANSIESTA_DEBUG
   write(*,*) 'Completed TRANSIESTA SCF'
#endif

#ifdef TRANSIESTA_TIMING
   call timer('TS_HS',3)
   call timer('TS_READ',3)
   call timer('TS_EXPAND',3)
   call timer('TS_PREPG',3)
   call timer('TS_UPDM',3)
   call timer('GFGGF',3)
#endif
   
!***********************
!     Close Files
!***********************
   if ( IONode ) then
      call io_close(uGFL)
      call io_close(uGFR)
   end if

!***********************
! CLEAN UP, there is alot to do... :)
!***********************
    if ( ts_Gamma_SCF ) then
       call delete(spDM)
       call delete(spEDM)
       if ( IsVolt ) then
          call delete(spDMR)
          call delete(spDMneqL)
          call delete(spDMneqR)
          call delete(spEDMR)
       end if

       ! The Hamiltonian and overlap matrices (in Gamma calculations
       ! we will not have any phases, hence, it makes no sense to
       ! have the arrays in complex)
       call delete(spH)
       call delete(spS)

    else
       call delete(spzDM)
       call delete(spzEDM)
       if ( IsVolt ) then
          call delete(spzDMR)
          call delete(spzDMneqL)
          call delete(spzDMneqR)
          call delete(spzEDMR)
       end if

       ! The Hamiltonian and overlap matrices
       call delete(spzH)
       call delete(spzS)

    end if

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    deallocate(qLb,wqL,qRb,wqR)
    deallocate(kparL,wkparL,kparR,wkparR)
    call memory('D','D',4*(nqL+nqR+nkparL+nkparR),'transiesta')

    if ( ts_Gamma_SCF ) then
       call memory('D','D',ndwork,'transiesta')
       deallocate(dwork)
    end if
    call memory('D','Z',nzwork,'transiesta')
    deallocate(zwork)
    call memory('D','Z',size(GF),'transiesta')
    deallocate(GF)

    deallocate(lasto_L,lasto_R)
    call memory('D','I',na_R_HS+na_L_HS+2,'transiesta')

    call memory('D','Z',size(HAAL)*3+size(HAAR)*3,'transiesta')
    deallocate(HAAL,SAAL,GAAL)
    deallocate(HAAR,SAAR,GAAR)
    call memory('D','Z',size(SigmaL)+size(SigmaR),'transiesta')
    deallocate(SigmaL,SigmaR)
    
    if ( IsVolt ) then
       ! These where only allocated on voltage calculations
       call memory('D','D',size(GammaL)+size(GammaR),'transiesta')
       deallocate(GammaL,GammaR)
    end if

    ! I would like to add something here which enables 
    ! 'Transiesta' post-process

    call timer('TS_calc',2)

#ifdef TRANSIESTA_DEBUG
   call write_debug( 'POS transiesta mem' )
#endif

 end subroutine transiesta_mem

   ! Helper routine to create and distribute the sparse 
   ! k-point Hamiltonian.
   subroutine create_HS_kpt(dit,sp, &
        Ef, &
        no_BufL, no_BufR, &
        no_C_L, no_C_R, no_u, &
        maxnh, H, S, xij, SpArrH, SpArrS, k, &
        nwork, work)

     use intrinsic_missing, only : SFIND
     use geom_helper, only : UCORB
     use class_OrbitalDistribution
     use class_Sparsity
     use class_zSpArr1D
     use parallel, only : Node

! *********************
! * INPUT variables   *
! *********************
     ! the distribution that the H and S live under
     type(OrbitalDistribution), intent(inout) :: dit
     ! The (local) sparsity pattern that H, S and xij lives by
     type(Sparsity), intent(inout) :: sp
     ! Fermi-level
     real(dp), intent(in) :: Ef
     ! The number of orbitals we wish to cut-off at both ends
     integer, intent(in) :: no_BufL, no_BufR
     integer, intent(in) :: no_C_L, no_C_R, no_u
     ! The number of elements in the sparse arrays
     integer, intent(in) :: maxnh
     ! The hamiltonian and overlap sparse matrices 
     real(dp), intent(in) :: H(maxnh),S(maxnh)
     ! The orbital distance array
     real(dp), intent(in) :: xij(3,maxnh)
     ! The arrays we will save in...
     type(zSpArr1D), intent(inout) :: SpArrH, SpArrS
     ! The k-point we will create
     real(dp), intent(in) :: k(3)
     ! we pass a work array
     integer, intent(in) :: nwork
     ! work-array
     complex(dp), intent(in out) :: work(nwork)

! *********************
! * LOCAL variables   *
! *********************
     complex(dp), parameter :: z_one = dcmplx(0._dp,1._dp)
     ! Create loop-variables for doing stuff
     integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
     integer, pointer :: k_ncol(:), k_ptr(:), k_col(:)
     complex(dp), pointer :: zH(:), zS(:)
     real(dp) :: ph
     type(Sparsity), pointer :: sp_k
     integer :: no_l, lio, io, ind, jo, jg, ind_k, kn
     integer :: no_max
     
     ! obtain the local number of rows and the global...
     no_l = nrows(sp)
     if ( no_u /= nrows_g(sp) ) then
        call die('Creating the k-&point matrix in &
             &transiesta went wrong. Please TODO...')
     end if

     ! Create all the local sparsity super-cell
     l_ncol => n_col   (sp)
     l_ptr  => list_ptr(sp)
     l_col  => list_col(sp)

     ! obtain the full sparsity unit-cell
     sp_k   => spar    (SpArrH)
     k_ncol => n_col   (sp_k)
     k_ptr  => list_ptr(sp_k)
     k_col  => list_col(sp_k)

     ! The boundary at the right buffer
     no_max = no_u - no_BufR
     
     call init_zSpArr1D(SpArrH)
     call init_zSpArr1D(SpArrS)
     ! obtain the value arrays...
     zH => val(SpArrH)
     zS => val(SpArrS)

     do lio = 1 , no_l

        ! obtain the global index of the orbital.
        io = index_local_to_global(dit,lio,Node)
        kn = k_ncol(io)
        ! if there is no contribution in this row
        if ( kn == 0 ) cycle

        ! The io orbitals are in the range [1;no_u_TS]
        ! This should be redundant as it is catched by kn==0
        if ( io <= no_BufL .or. no_max < io ) cycle

        ! Loop number of entries in the row... (index frame)
        do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

           ! as the local sparsity pattern is a super-cell pattern,
           ! we need to check the unit-cell orbital
           ! The unit-cell column index
           jo = UCORB(l_col(ind),no_u)

           ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
           if ( jo <= no_BufL .or. no_max < jo ) cycle

           ! Do a check whether we have connections
           ! across the junction...
           ! This is the same as removing LEFT-RIGHT states..
           if ( io < no_C_L .and. no_C_R < jo ) cycle
           if ( jo < no_C_L .and. no_C_R < io ) cycle
           
           ! find the equivalent position in the sparsity pattern
           ! of the full unit cell
           ind_k = k_ptr(io)

           ! Notice that SFIND REQUIRES that the sparsity pattern
           ! is SORTED!
           ! Thus it will only work for UC sparsity patterns.
           ind_k = ind_k + SFIND(k_col(ind_k+1:ind_k+kn),jo)

!           if ( k_ptr(io) == ind_k ) then
!              write(*,'(a,10000(tr1,i0))') &
!                   'Something should be checked...',jo,k_col(ind_k+1:ind_k+kn)
!           end if
!              do jg = 1 , k_ncol(io)
!                 ind_k = k_ptr(io)+jg
!                 if ( k_col(ind_k) == jo ) then
!                    ph = k(1) * xij(1,ind) + &
!                         k(2) * xij(2,ind) + &
!                         k(3) * xij(3,ind)
!                    zH(ind_k) = zH(ind_k) + H(ind) * cdexp(z_one * ph)
!                    zS(ind_k) = zS(ind_k) + S(ind) * cdexp(z_one * ph)
!                 end if
!              end do
!           else
              ph = k(1) * xij(1,ind) + &
                   k(2) * xij(2,ind) + &
                   k(3) * xij(3,ind)

              zH(ind_k) = zH(ind_k) + H(ind) * cdexp(z_one * ph)

              zS(ind_k) = zS(ind_k) + S(ind) * cdexp(z_one * ph)

!           end if

        end do

     end do

#ifdef MPI
     ! Note that zH => val(SpArrH)
     ! Note that zS => val(SpArrS)
     call AllReduce_zSpArr1D(SpArrH,nwork,work)
     call AllReduce_zSpArr1D(SpArrS,nwork,work)
#endif

     ! We symmetrize AND shift
     call symmetrize_HS_kpt(Ef,SpArrH,SpArrS)
     
     ! It could be argued that MPI reduction provides
     ! numeric fluctuations.
     ! However, the block-cyclic distribution ensures that
     ! there are no two elements accessed by two or more processors.
     ! This makes all non-local elements ZERO, and there should not
     ! be flucuations on adding ZEROS as they are *only* dependent
     ! on order of summation.

   end subroutine create_HS_kpt

   ! Helper routine to create and distribute the sparse 
   ! k-point Hamiltonian.
   subroutine create_HS_Gamma(dit,sp, &
        Ef, &
        no_BufL, no_BufR, &
        no_C_L, no_C_R, no_u, &
        maxnh, H, S, SpArrH, SpArrS, &
        nwork, work)

     use intrinsic_missing, only: SFIND, UCORB => MODP
     use class_OrbitalDistribution
     use class_Sparsity
     use class_dSpArr1D
     use parallel, only : Node
! *********************
! * INPUT variables   *
! *********************
     ! the distribution that the H and S live under
     type(OrbitalDistribution), intent(inout) :: dit
     ! The (local) sparsity pattern that H, S and xij lives by
     type(Sparsity), intent(inout) :: sp
     ! Fermi-level
     real(dp), intent(in) :: Ef
     ! The number of orbitals we wish to cut-off at both ends
     integer, intent(in) :: no_BufL, no_BufR
     integer, intent(in) :: no_C_L, no_C_R, no_u
     ! The number of elements in the sparse arrays
     integer, intent(in) :: maxnh
     ! The hamiltonian and overlap sparse matrices 
     real(dp), intent(in) :: H(maxnh),S(maxnh)
     ! The arrays we will save in... these are the entire TS-region sparsity
     type(dSpArr1D), intent(inout) :: SpArrH, SpArrS
     ! we pass a work array
     integer, intent(in) :: nwork
     ! work-array
     real(dp), intent(in out) :: work(nwork)

! *********************
! * LOCAL variables   *
! *********************
     ! Create loop-variables for doing stuff
     integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
     integer, pointer  :: k_ncol(:), k_ptr(:), k_col(:)
     real(dp), pointer :: dH(:), dS(:)
     type(Sparsity), pointer :: sp_G
     integer :: no_l, lio, io, ind, jo, jg, ind_k
     
     ! obtain the local number of rows and the global...
     no_l = nrows(sp)
     if ( no_u /= nrows_g(sp) ) then
        call die('Creating the k-&point matrix in &
             &transiesta went wrong. Please TODO...')
     end if

     ! Create all the local sparsity super-cell
     l_ncol => n_col   (sp)
     l_ptr  => list_ptr(sp)
     l_col  => list_col(sp)

     ! obtain the full sparsity unit-cell
     sp_G   => spar(SpArrH)
     k_ncol => n_col   (sp_G)
     k_ptr  => list_ptr(sp_G)
     k_col  => list_col(sp_G)
     
     ! initialize to 0
     call init_dSpArr1D(SpArrH)
     call init_dSpArr1D(SpArrS)
     ! obtain the value arrays...
     dH => val(SpArrH)
     dS => val(SpArrS)

     do lio = 1 , no_l

        ! obtain the global index of the orbital.
        io = index_local_to_global(dit,lio,Node)
        ! if there is no contribution in this row
        if ( k_ncol(io) == 0 ) cycle

        ! The io orbitals are in the range [1;no_u]
        if ( io <= no_BufL )       cycle
        if ( no_u - no_BufR < io ) cycle

        ! Loop number of entries in the row... (in the index frame)
        do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

           ! as the local sparsity pattern is a super-cell pattern,
           ! we need to check the unit-cell orbital
           ! The unit-cell column index
           jo = UCORB(l_col(ind),no_u)

           ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
           if ( jo <= no_BufL .or. no_u - no_BufR < jo ) cycle

           ! Do a check whether we have connections
           ! across the junction...
           ! This is the same as removing LEFT-RIGHT states..
           if ( io < no_C_L .and. no_C_R < jo ) cycle
           if ( jo < no_C_L .and. no_C_R < io ) cycle
           
           ! find the equivalent position in the sparsity pattern
           ! of the full unit cell
           ind_k = k_ptr(io)
           ind_k = ind_k + SFIND(k_col(ind_k+1:ind_k+k_ncol(io)),jo)

           ! Todo, as this is a Gamma-calculation
           ! we probably should NOT do 'dH = dH + H'
           ! rather 'dH = H'

!           if ( ind_k == k_ptr(io) ) then
!              write(*,*) 'Something should be checked:'
!           end if
!              write(*,*) ' 1. Is there a central region orbital having a &
!                   &z-component to the following cell?'
!              write(*,*) ' 2. Try and increase the z-direction of your &
!                   &unit-cell.'
!              
!              do jg = 1 , k_ncol(io)
!                 ind_k = k_ptr(io)+jg
!                 if ( k_col(ind_k) == jo ) then
!                    dH(ind_k) = dH(ind_k) + H(ind)
!                    dS(ind_k) = dS(ind_k) + S(ind)
!                 end if
!              end do
!           else
              dH(ind_k) = dH(ind_k) + H(ind)
              dS(ind_k) = dS(ind_k) + S(ind)
!           end if

        end do

     end do
     

#ifdef MPI
     ! Note that dH => val(SpArrH)
     ! Note that dS => val(SpArrS)
     call AllReduce_dSpArr1D(SpArrH,nwork,work)
     call AllReduce_dSpArr1D(SpArrS,nwork,work)
#endif

     ! We need to do symmetrization AFTER reduction as we need the full
     ! Hamiltonian before we can do anything
     call symmetrize_HS_Gamma(Ef,SpArrH,SpArrS)

     ! It could be argued that MPI reduction provides
     ! numeric fluctuations.
     ! However, the block-cyclic distribution ensures that
     ! there are no two elements accessed by two or more processors.
     ! This makes all non-local elements ZERO, and there should not
     ! be flucuations on adding ZEROS as they are *only* dependent
     ! on order of summation.

   end subroutine create_HS_Gamma

   subroutine symmetrize_HS_Gamma(Ef,SpArrH, SpArrS)
     use intrinsic_missing, only: SFIND, UCORB => MODP
     use class_Sparsity
     use class_dSpArr1D
     use parallel, only : Node
! *********************
! * INPUT variables   *
! *********************
     real(dp), intent(in) :: Ef
     ! The arrays we will save in... these are the entire TS-region sparsity
     type(dSpArr1D), intent(inout) :: SpArrH, SpArrS

! *********************
! * LOCAL variables   *
! *********************
     ! Create loop-variables for doing stuff
     type(Sparsity), pointer :: s
     integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
     real(dp), pointer :: dH(:), dS(:)
     integer :: nr, io, ind, jo, rin, rind

     s    => spar(SpArrH)
     nr   = nrows_g(s)
     l_ncol => n_col   (s)
     l_ptr  => list_ptr(s)
     l_col  => list_col(s)

     dH     => val(SpArrH)
     dS     => val(SpArrS)

     ! This loop is across the local rows...
     do io = 1 , nr

        ! Quickly go past the empty regions... (we have nothing to update)
        if ( l_ncol(io) == 0 ) cycle

        ! Now we loop across the update region
        ! This one must *per definition* have less elements.
        ! Hence, we can exploit this, and find equivalent
        ! super-cell orbitals.
        do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

           jo = l_col(ind)

           ! As we symmetrize we do not need
           ! to cycle all points through two times...
           if ( jo < io ) cycle

           ! We will find the Hermitian part:
           ! The fact that we have a SYMMETRIC
           ! update region makes this *tricky* part easy...
           rin  = l_ptr(jo)
           ! TODO, this REQUIRES that l_col(:) is sorted
           rind = rin + SFIND(l_col(rin+1:rin+l_ncol(jo)),io)
           ! We do a check, just to be sure...
           if ( rind == rin ) then
              call die('ERROR symmetrization orbital does not &
                   &exist.')
           end if

           ! Symmetrize
           dS(ind)  = 0.5_dp * ( dS(ind) + dS(rind) )
           dH(ind)  = 0.5_dp * ( dH(ind) + dH(rind) ) &
                - Ef * dS(ind)

           ! we have a real Matrix (so imaginary part is zero)
           dH(rind) = dH(ind)
           dS(rind) = dS(ind)
                      
        end do
     end do

   end subroutine symmetrize_HS_Gamma

   subroutine symmetrize_HS_kpt(Ef,SpArrH, SpArrS)
     use intrinsic_missing, only: SFIND, UCORB => MODP
     use class_Sparsity
     use class_zSpArr1D
     use parallel, only : Node
! *********************
! * INPUT variables   *
! *********************
     real(dp), intent(in) :: Ef
     ! The arrays we will save in... these are the entire TS-region sparsity
     type(zSpArr1D), intent(inout) :: SpArrH, SpArrS

! *********************
! * LOCAL variables   *
! *********************
     ! Create loop-variables for doing stuff
     type(Sparsity), pointer :: s
     integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
     complex(dp), pointer :: zH(:), zS(:)
     integer :: nr, io, ind, jo, rin, rind

     s    => spar(SpArrH)
     nr   = nrows_g(s)
     l_ncol => n_col   (s)
     l_ptr  => list_ptr(s)
     l_col  => list_col(s)

     zH     => val(SpArrH)
     zS     => val(SpArrS)

     ! This loop is across the local rows...
     do io = 1 , nr

        ! Quickly go past the empty regions... (we have nothing to update)
        if ( l_ncol(io) == 0 ) cycle

        ! Now we loop across the update region
        ! This one must *per definition* have less elements.
        ! Hence, we can exploit this, and find equivalent
        ! super-cell orbitals.
        do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

           jo = l_col(ind)

           ! As we symmetrize we do not need
           ! to cycle all points through two times...
           if ( jo < io ) cycle

           ! We will find the Hermitian part:
           ! The fact that we have a SYMMETRIC
           ! update region makes this *tricky* part easy...
           rin  = l_ptr(jo)
           ! TODO, this REQUIRES that l_col(:) is sorted
           rind = rin + SFIND(l_col(rin+1:rin+l_ncol(jo)),io)
           ! We do a check, just to be sure...
           if ( rind == rin ) then
              call die('ERROR symmetrization orbital does not &
                   &exist.')
           end if

           ! Symmetrize (notice that this is *transposed*)
           ! See prep_GF
           zS(rind)  = 0.5_dp * ( zS(ind) + dconjg(zS(rind)) )
           zH(rind)  = 0.5_dp * ( zH(ind) + dconjg(zH(rind)) ) &
                - Ef * zS(rind)
           
           zS(ind) = dconjg(zS(rind))
           zH(ind) = dconjg(zH(rind))

           if ( ind == rind ) then
              ! This is the diagonal matrix elements
              zS(ind) = dreal(zS(ind))
              zH(ind) = dreal(zH(ind))
           end if
                      
        end do
     end do

   end subroutine symmetrize_HS_kpt


! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************

   subroutine init_dSpArr1D(sp_arr)
     use class_dSpArr1D
     type(dSpArr1D), intent(inout) :: sp_arr
     real(dp), pointer :: arr(:)
     arr    => val(sp_arr)
     arr(:) =  0._dp
   end subroutine init_dSpArr1D

   subroutine init_zSpArr1D(sp_arr)
     use class_zSpArr1D
     type(zSpArr1D), intent(inout) :: sp_arr
     complex(dp), pointer :: arr(:)
     arr    => val(sp_arr)
     arr(:) =  dcmplx(0._dp,0._dp)
   end subroutine init_zSpArr1D

#ifdef MPI
   subroutine AllReduce_zSpArr1D(sp_arr,nwork,work)
     use mpi_siesta
     use class_zSpArr1D
     type(zSpArr1D), intent(inout) :: sp_arr
     integer, intent(in) :: nwork
     complex(dp), intent(inout) :: work(nwork)
     complex(dp), pointer :: arr(:)
     integer :: s, MPIerror
     s = nnzs(sp_arr)
     ! This should never happen, exactly due to the sparsity
     if ( s > nwork ) call die('Sparsity seems larger than &
          &work arrays, Transiesta????')
     arr => val(sp_arr)
     call MPI_AllReduce(arr,work,s, &
          MPI_Double_Complex, MPI_Sum, MPI_Comm_World, MPIerror)
     arr(:) = work(1:s)
   end subroutine AllReduce_zSpArr1D

   subroutine AllReduce_dSpArr1D(sp_arr,nwork,work)
     use mpi_siesta
     use class_dSpArr1D
     type(dSpArr1D), intent(inout) :: sp_arr
     integer, intent(in) :: nwork
     real(dp), intent(inout) :: work(nwork)
     real(dp), pointer :: arr(:)
     integer :: s, MPIerror
     s = nnzs(sp_arr)
     ! This should never happen, exactly due to the sparsity
     if ( s > nwork ) call die('Sparsity seems larger than &
          &work arrays, Transiesta????')
     arr => val(sp_arr)
     call MPI_AllReduce(arr,work,s, &
          MPI_Double_Precision, MPI_Sum, MPI_Comm_World, MPIerror)
     arr(:) = work(1:s)
   end subroutine AllReduce_dSpArr1D
#endif


! Update DM
! These routines are supplied for easy update of the update region
! sparsity patterns
! Note that these routines implement the usual rho(Z) \propto - GF
   subroutine add_DM_dE_Z(DM,EDM,no1,no2,GF,no_BufL,GF_offset,DMfact,EDMfact)
     use class_zSpArr1D
     use class_Sparsity
     ! The DM and EDM equivalent matrices
     type(zSpArr1D), intent(inout) :: DM,EDM
     ! The size of GF
     integer, intent(in) :: no1,no2
     ! The Green's function
     complex(dp), intent(in) :: GF(no1,no2)
     ! The number of buffer atoms (needed for the offset in the sparsity
     ! patterns), and the offset in the GF
     integer, intent(in) :: no_BufL, GF_offset
     ! Complex numbers that are used in the factor of GF
     complex(dp), intent(in) :: DMfact, EDMfact

     ! Arrays needed for looping the sparsity
     type(Sparsity), pointer :: s
     integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
     complex(dp), pointer :: zD(:), zE(:)
     integer :: io, ind, nr, iu, ju

     s      => spar(DM)
     l_ncol => n_col   (s)
     l_ptr  => list_ptr(s)
     l_col  => list_col(s)
     zD     => val(DM)
     zE     => val(EDM)

     ! Number of orbitals in the SIESTA unit-cell
     ! Remember that this is a sparsity pattern which contains
     ! a subset of the SIESTA pattern.
     nr = nrows(s)
     
     if ( no1 < no2 ) call die('The GF format is not as &
          &expected.')
     
     do io = 1 , nr
        ! Quickly go past the buffer atoms...
        if ( l_ncol(io) == 0 ) cycle

        ! The update region equivalent GF part
        iu = io - no_BufL
        
        do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

           ! We need to subtract the offset of
           ! the Green's function
           ! Any offset will ONLY be in the column
           ! index. See explanation in the 
           ! mem_sparsity module...
           ju = l_col(ind) - no_BufL - GF_offset

           zD(ind) = zD(ind) - GF(iu,ju) * DMfact
           zE(ind) = zE(ind) - GF(iu,ju) * EDMfact
        end do
     end do

   end subroutine add_DM_dE_Z

   subroutine add_DM_dE_D(DM,EDM,no1,no2,GF,no_BufL,GF_offset,DMfact,EDMfact)
     use class_dSpArr1D
     use class_Sparsity
     ! The DM and EDM equivalent matrices
     type(dSpArr1D), intent(inout) :: DM,EDM
     ! The size of GF
     integer, intent(in) :: no1,no2
     ! The Green's function
     complex(dp), intent(in) :: GF(no1,no2)
     ! The number of buffer atoms (needed for the offset in the sparsity
     ! patterns), and the offset in the GF
     integer, intent(in) :: no_BufL, GF_offset
     ! Complex numbers that are used in the factor of GF
     complex(dp), intent(in) :: DMfact, EDMfact

     ! Arrays needed for looping the sparsity
     type(Sparsity), pointer :: s
     integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
     real(dp), pointer :: dD(:), dE(:)
     integer :: io, ind, nr, iu, ju

     s      => spar(DM)
     l_ncol => n_col   (s)
     l_ptr  => list_ptr(s)
     l_col  => list_col(s)
     dD     => val(DM)
     dE     => val(EDM)

     ! Number of orbitals in the SIESTA unit-cell
     ! Remember that this is a sparsity pattern which contains
     ! a subset of the SIESTA pattern.
     nr = nrows(s)
     
     if ( no1 < no2 ) call die('The GF format is not as &
          &expected.')
     
     do io = 1 , nr
        ! Quickly go past the buffer atoms...
        if ( l_ncol(io) == 0 ) cycle

        ! The update region equivalent GF part
        iu = io - no_BufL
        
        do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

           ! We need to subtract the offset of
           ! the Green's function
           ! Any offset will ONLY be in the column
           ! index. See explanation in the 
           ! mem_sparsity module...
           ju = l_col(ind) - no_BufL - GF_offset
     
           dD(ind) = dD(ind) - dimag( GF(iu,ju) * DMfact  )
           dE(ind) = dE(ind) - dimag( GF(iu,ju) * EDMfact )

        end do
     end do

   end subroutine add_DM_dE_D

   subroutine update_DM(dit,sp,maxn,DM,EDM, spDM, spEDM)
          ! The DM and EDM equivalent matrices
     use class_OrbitalDistribution
     use class_Sparsity
     use class_dSpArr1D
     use intrinsic_missing, only : UCORB => MODP
     use parallel, only : Node
     type(OrbitalDistribution), intent(inout) :: dit
     type(Sparsity), intent(inout) :: sp
     ! Size of the sparsity arrays
     integer, intent(in) :: maxn
     ! Sparse DM-arrays (local)
     real(dp), intent(inout) :: DM(maxn), EDM(maxn)
     ! Updated sparsity arrays (they contain the current integration)
     type(dSpArr1D), intent(inout) :: spDM, spEDM

     ! Arrays needed for looping the sparsity
     type(Sparsity), pointer :: s
     integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
     integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
     real(dp), pointer :: dD(:), dE(:)
     integer :: lnr, lio, lind, io, ind, nr, ljo

     l_ncol => n_col   (sp)
     l_ptr  => list_ptr(sp)
     l_col  => list_col(sp)
     s        => spar(spDM)
     lup_ncol => n_col   (s)
     lup_ptr  => list_ptr(s)
     lup_col  => list_col(s)
     dD     => val(spDM)
     dE     => val(spEDM)
     
     ! Number of orbitals in the SIESTA unit-cell
     ! Remember that this is a sparsity pattern which contains
     ! a subset of the SIESTA pattern.
     lnr = nrows(sp)
     nr  = nrows_g(sp)
     
     if ( nr /= nrows(s) ) call die('The sparsity format is not as &
          &expected.')
     
     ! This loop is across the local rows...
     do lio = 1 , lnr

        ! obtain the global index of the local orbital.
        io = index_local_to_global(dit,lio,Node)

        ! Quickly go past the empty regions... (we have nothing to update)
        if ( lup_ncol(io) == 0 ) cycle

        ! Do a loop in the local sparsity pattern...
        ! The local sparsity pattern is more "spread", hence
        ! we do fewer operations by having this as an outer loop
        do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

           ljo = UCORB(l_col(lind),nr)

           ! Now we loop across the update region
           ! This one must *per definition* have less elements.
           ! Hence, we can exploit this, and find equivalent
           ! super-cell orbitals.
           ! Ok, this is Gamma (but to be consistent)
           do ind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)

              if ( ljo /= lup_col(ind) ) cycle

              ! Probably we dont need to "add"
              ! We only have one k-point...
              DM(lind)  = DM(lind)  + dD(ind)
              EDM(lind) = EDM(lind) + dE(ind)

           end do
        end do
     end do

   end subroutine update_DM

   subroutine init_DM(dit,sp,maxn,DM,EDM, up_sp)
          ! The DM and EDM equivalent matrices
     use class_OrbitalDistribution
     use class_Sparsity
     use intrinsic_missing, only : UCORB => MODP
     use parallel, only : Node
     type(OrbitalDistribution), intent(inout) :: dit
     type(Sparsity), intent(inout) :: sp
     ! Size of the sparsity arrays
     integer, intent(in) :: maxn
     ! Sparse DM-arrays (local)
     real(dp), intent(inout) :: DM(maxn), EDM(maxn)
     ! The updated sparsity arrays...
     type(Sparsity), intent(inout) :: up_sp

     integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
     integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
     integer :: lnr, lio, lind, io, jo, ind, nr

     l_ncol => n_col   (sp)
     l_ptr  => list_ptr(sp)
     l_col  => list_col(sp)
     lup_ncol => n_col   (up_sp)
     lup_ptr  => list_ptr(up_sp)
     lup_col  => list_col(up_sp)
     
     ! Number of orbitals in the SIESTA unit-cell
     ! Remember that this is a sparsity pattern which contains
     ! a subset of the SIESTA pattern.
     lnr = nrows(sp)
     nr  = nrows_g(sp)
     
     if ( nr /= nrows(up_sp) ) call die('The sparsity format is not as &
          &expected.')
     
     ! This loop is across the local rows...
     do lio = 1 , lnr

        ! obtain the global index of the local orbital.
        io = index_local_to_global(dit,lio,Node)

        ! Quickly go past the empty regions... (we have nothing to update)
        if ( lup_ncol(io) == 0 ) cycle

        ! Now we loop across the update region
        ! This one must *per definition* have less elements.
        ! Hence, we can exploit this, and find equivalent
        ! super-cell orbitals.
        ! Ok, this is Gamma (but to be consistent)
        do ind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)

           jo = lup_col(ind)

           ! Do a loop in the local sparsity pattern...
           do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

              ! We know that the update region is in 
              ! UC-format. Hence we can compare directly, via
              ! the orbital index in the unit-cell.
              if ( UCORB(l_col(lind),nr) == jo ) then
                 DM (lind) = 0._dp
                 EDM(lind) = 0._dp
              end if

           end do
        end do
     end do

   end subroutine init_DM

   subroutine update_zDM(dit,sp,n_nzs,DM,EDM,xij,spDM, spEDM,k)
     use class_OrbitalDistribution
     use class_Sparsity
     use class_zSpArr1D

     use intrinsic_missing, only : SFIND, UCORB => MODP
     use parallel, only : Node
     type(OrbitalDistribution), intent(inout) :: dit
     type(Sparsity), intent(inout) :: sp
     ! Size of the sparsity arrays
     integer, intent(in) :: n_nzs
     ! Sparse DM-arrays (local)
     real(dp), intent(inout) :: DM(n_nzs), EDM(n_nzs)
     ! The orbital distances
     real(dp), intent(in) :: xij(3,n_nzs)
     ! Updated sparsity arrays (they contain the current integration)
     type(zSpArr1D), intent(inout) :: spDM, spEDM
     ! The k-point...
     real(dp), intent(in) :: k(3)

     ! Arrays needed for looping the sparsity
     type(Sparsity), pointer :: s
     integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
     integer, pointer :: lup_ncol(:), lup_ptr(:), lup_col(:)
     complex(dp), pointer :: zD(:), zE(:)
     complex(dp) :: ph_m, ph_p, kx
     integer :: lio, io, jo, ind, nr, ljo
     integer :: lnr, lind, rin, rind

     l_ncol => n_col   (sp)
     l_ptr  => list_ptr(sp)
     l_col  => list_col(sp)

     s        => spar(spDM)
     lup_ncol => n_col   (s)
     lup_ptr  => list_ptr(s)
     lup_col  => list_col(s)
     zD     => val(spDM)
     zE     => val(spEDM)
     
     ! Number of orbitals in the SIESTA unit-cell
     ! Remember that this is a sparsity pattern which contains
     ! a subset of the SIESTA pattern.
     lnr = nrows(sp)
     nr  = nrows_g(sp)
     
     if ( nr /= nrows(s) ) call die('The sparsity format is not as &
          &expected.')
     
     ! This loop is across the local rows...
     do lio = 1 , lnr

        ! obtain the global index of the local orbital.
        io = index_local_to_global(dit,lio,Node)

        ! Quickly go past the empty regions... (we have nothing to update)
        if ( lup_ncol(io) == 0 ) cycle

        ! Do a loop in the local sparsity pattern...
        ! The local sparsity pattern is more "spread", hence
        ! we do fewer operations by having this as an outer loop
        do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

           ljo = UCORB(l_col(lind),nr)
           
           ! Now we loop across the update region
           ! This one must *per definition* have less elements.
           ! Hence, we can exploit this, and find equivalent
           ! super-cell orbitals.
           do ind = lup_ptr(io) + 1 , lup_ptr(io) + lup_ncol(io)
              
              jo = lup_col(ind)

              ! We know that the update region is in 
              ! UC-format. Hence we can compare directly, via
              ! the orbital index in the unit-cell.
              if ( ljo /= jo ) cycle

              kx = k(1) * xij(1,lind) + &
                   k(2) * xij(2,lind) + &
                   k(3) * xij(3,lind)
              
              ! The fact that we have a SYMMETRIC
              ! update region makes this *tricky* part easy...
              rin  = lup_ptr(jo)
              ! TODO, this REQUIRES that lup_col(:) is sorted
              rind = rin+SFIND(lup_col(rin+1:rin+lup_ncol(jo)),io)
              ! We do a check, just to be sure...
              if ( rind == rin ) then
                 call die('ERROR: symmetrization points does not exist')
              end if
              
              ph_p = cdexp(dcmplx(0._dp,+1._dp)*kx)
              ph_m = cdexp(dcmplx(0._dp,-1._dp)*kx)

              DM(lind)  = DM(lind)  + 0.5_dp * dimag( &
                   ph_p*zD(rind) + ph_m*zD(ind) )

              EDM(lind) = EDM(lind) + 0.5_dp * dimag( &
                   ph_p*zE(rind) + ph_m*zE(ind) )

           end do
        end do
     end do

   end subroutine update_zDM

   ! creation of the GF^{-1}.
   ! this routine will ONLY insert the zS-H terms in the GF 
   subroutine prepare_GF_inv_D(spH,spS, ZE, no_BufL,no_u,GFinv)
     use class_dSpArr1D
     use class_Sparsity

     ! The Hamiltonian and overlap sparse matrices
     type(dSpArr1D), intent(inout) :: spH, spS
     ! the current energy point
     complex(dp), intent(in) :: ZE
     ! Remark that we need the left buffer orbitals
     ! to calculate the actual orbital of the sparse matrices...
     integer, intent(in) :: no_BufL,no_u
     complex(dp), intent(out) :: GFinv(no_u**2)

     ! Local variables
     type(Sparsity), pointer :: s
     integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
     real(dp), pointer :: dH(:), dS(:)
     integer :: io, iu, ind, ioff

     s      => spar    (spH)
     l_ncol => n_col   (s)
     l_ptr  => list_ptr(s)
     l_col  => list_col(s)

     ioff = no_BufL + 1

     dH     => val(spH)
     dS     => val(spS)

     ! initialize 
     GFinv(1:no_u**2) = dcmplx(0._dp,0._dp)

     ! We will only loop in the central region
     ! We have constructed the sparse array to only contain
     ! values in this part...
     do io = no_BufL + 1, no_BufL + no_u

        iu = (io - ioff) * no_u - no_BufL

        do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 

           !ju = l_col(ind) ! the '- no_BufL' is moved outside the loop

           GFinv(l_col(ind)+iu) = ZE * dS(ind) - dH(ind)
        end do
     end do

   end subroutine prepare_GF_inv_D

   ! creation of the GF^{-1}.
   ! this routine will ONLY insert the zS-H terms in the GF 
   subroutine prepare_GF_inv_Z(spH,spS,ZE, no_BufL,no_u,GFinv)
     use class_zSpArr1D
     use class_Sparsity
     ! The Hamiltonian and overlap sparse matrices
     type(zSpArr1D), intent(inout) :: spH, spS
     ! the current energy point
     complex(dp), intent(in) :: ZE
     ! Remark that we need the left buffer orbitals
     ! to calculate the actual orbital of the sparse matrices...
     integer, intent(in) :: no_BufL,no_u
     complex(dp), intent(out) :: GFinv(no_u**2)

     ! Local variables
     type(Sparsity), pointer :: s
     integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
     complex(dp), pointer :: zH(:), zS(:)
     integer :: io, iu,ind, ioff
     
     s      => spar    (spH)
     l_ncol => n_col   (s)
     l_ptr  => list_ptr(s)
     l_col  => list_col(s)
     
     ! Offset
     ioff = no_BufL + 1

     zH     => val(spH)
     zS     => val(spS)

     ! Initialize
     GFinv(1:no_u**2) = dcmplx(0._dp,0._dp)

     ! We will only loop in the central region
     ! We have constructed the sparse array to only contain
     ! values in this part...
     do io = no_BufL + 1, no_BufL + no_u

        iu = (io-ioff) * no_u - no_BufL

        do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 

           !ju = l_col(ind) ! The '- no_BufL' is moved outside the loop

           ! Notice that we transpose S and H back here
           ! See symmetrize_HS_k
           GFinv(l_col(ind)+iu) = ZE * zS(ind) - zH(ind)
        end do
     end do
     
   end subroutine prepare_GF_inv_Z
 
end module m_ts_mem
