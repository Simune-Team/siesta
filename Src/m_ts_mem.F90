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
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.
! * It has been heavily inspired by the original authors of the 
!   Transiesta code (hence the references here are still remaining) *

module m_ts_mem

  use precision, only : dp

  use m_ts_sparse_helper, only : create_HS_kpt
  use m_ts_sparse_helper, only : create_HS_Gamma
  use m_ts_sparse_helper, only : symmetrize_HS_kpt
  use m_ts_sparse_helper, only : symmetrize_HS_Gamma
#ifdef MPI
  use m_ts_sparse_helper, only : AllReduce_dSpData1D
  use m_ts_sparse_helper, only : AllReduce_zSpData1D
#endif
  use m_ts_sparse_helper, only : init_DM
  use m_ts_sparse_helper, only : update_DM
  use m_ts_sparse_helper, only : update_zDM

  use m_ts_sparse_helper, only : ts_print_charges
  use m_ts_sparse_helper, only : TS_INFO_SCF
  
  use m_ts_sparse_helper, only : weightDM
  use m_ts_sparse_helper, only : weightDMC
  
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
    use class_dSpData1D
    use class_zSpData1D

    use m_ts_kpoints

    use m_ts_options, only : IsVolt, UseBulk, UpdateDMCR
    use m_ts_options, only : VoltL, VoltR
    use m_ts_electype
    use m_ts_options, only : ElLeft, ElRight
    use m_ts_options, only : GFFileL, GFFileR
    use m_ts_options, only : na_BufL, no_BufL
    use m_ts_options, only : na_BufR, no_BufR

    use m_ts_sparse, only : ts_sp_uc
    use m_ts_sparse, only : tsup_sp_uc
    
    ! Self-energy retrival and expansion
    use m_ts_elec_se
    ! Gf calculation
    use m_ts_mem_scat

    use m_ts_method, only : GF_INV_EQUI_PART

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
    integer, intent(in)  :: no_u,na_u
    integer, intent(in)  :: lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: xij(3,n_nzs)
    real(dp), intent(in) :: Hs(n_nzs,nspin), Ss(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Ef, Qtot
    integer, intent(in)  :: TSiscf

! ******************* Regional sizes *************************
! * Electrode regions
    integer :: na_L_HS, na_R_HS
    integer :: no_L_HS, no_R_HS
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
    complex(dp), allocatable :: HAAL(:,:,:), SAAL(:,:,:)
    complex(dp), allocatable :: HAAR(:,:,:), SAAR(:,:,:)
    complex(dp), pointer :: GAAL(:,:), GAAR(:,:)
    complex(dp), allocatable :: SigmaL(:,:), SigmaR(:,:)
    complex(dp), target, allocatable :: GammaLT(:,:), GammaRT(:,:)
!    real(dp), allocatable :: GammaLT(:,:), GammaRT(:,:)
! ************************************************************

! ******************* Computational arrays *******************
    integer :: ndwork, nzwork
    real(dp),    allocatable :: dwork(:)
    complex(dp), allocatable :: zwork(:), GF(:)

    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData1D) :: spH, spS
    type(zSpData1D) :: spzH, spzS
    ! The different sparse matrices...
    type(dSpData1D) :: spDM, spEDM, spDMR, spEDMR, spDMneqL, spDMneqR
    type(zSpData1D) :: spzDM, spzEDM, spzDMR, spzEDMR, spzDMneqL, spzDMneqR
    ! Pointers for updating the density matrices
    real(dp),    pointer :: dDM(:), dEDM(:)
    complex(dp), pointer :: zDM(:), zEDM(:)
! ************************************************************

! ******************* Computational variables ****************
    complex(dp) :: Z, W, ZW
    real(dp)    :: k(3)
    complex(dp), parameter :: zmi = dcmplx(0._dp,-1._dp)
! ************************************************************

! ******************** Loop variables ************************
    integer :: ispin, ikpt, iPE, iE, NEReqs, up_nzs, ia, ia_E
    integer :: ind
! ************************************************************

! ******************* Miscalleneous variables ****************
    integer :: ierr
! ************************************************************


#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    call timer('TS_calc',1)
    
    ! Calculate the number of used atoms in left/right
    na_L_HS = UsedAtoms(ElLeft)
    na_R_HS = UsedAtoms(ElRight)
    no_L_HS = UsedOrbs(ElLeft)
    no_R_HS = UsedOrbs(ElRight)
    na_L = TotUsedAtoms(ElLeft)
    no_L = TotUsedOrbs(ElLeft)
    na_R = TotUsedAtoms(ElRight)
    no_R = TotUsedOrbs(ElRight)

    ! Create the lasto pointers for the electrode expansions...
    allocate(lasto_L(0:na_L_HS))
    lasto_L(0) = 0
    ia_E = 0
    do ia = na_BufL + 1 , na_BufL + na_L, Rep(ElLeft)
       ia_E = ia_E + 1
       lasto_L(ia_E) = lasto_L(ia_E-1) + lasto(ia) - lasto(ia-1)
    end do
    allocate(lasto_R(0:na_R_HS))
    lasto_R(0) = 0
    ia_E = 0
    do ia = na_u - na_R - na_BufR + 1 , na_u - na_BufR , Rep(ElRight)
       ia_E = ia_E + 1
       lasto_R(ia_E) = lasto_R(ia_E-1) + lasto(ia) - lasto(ia-1)
    end do
    ! Add to memory-management...
    call memory('A','I',na_L_HS+na_R_HS+2,'transiesta')

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
    call read_Green(uGFL,TSiscf==1,VoltL,ts_nkpnt,NEn,  &
         ElLeft,.false.,nspin, &
         nkparL,kparL,wkparL, &
         nqL,wqL,qLb)

! Right
    call read_Green(uGFR,TSiscf==1,VoltR,ts_nkpnt,NEn, &
         ElRight,.false.,nspin,  &
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

    if ( GF_INV_EQUI_PART ) then
       ispin = no_u_TS*(no_u_TS-no_R-no_L)
       if ( IsVolt ) then
          ispin = no_u_TS*max(no_u_TS-no_R-no_L,no_L+no_R)
       end if
       allocate(GF(ispin),stat=ierr)
       if (ierr/=0) call die('Could not allocate space for GFpart')
       call memory('A','Z',ispin,'transiesta')
    else
       ! Maybe we could introduce a blocking inversion, so that we can limit
       ! memory size... (perhaps to much work as we need a new solver)
       allocate(GF(no_u_TS**2),stat=ierr)
       if (ierr/=0) call die('Could not allocate space for GF')
       call memory('A','Z',no_u_TS**2,'transiesta')
    end if

    ! Allocate the left-right electrode quantities that we need
    allocate(HAAL(no_L_HS,no_L_HS,Rep(ElLeft)))
    allocate(SAAL(no_L_HS,no_L_HS,Rep(ElLeft)))
    ispin = no_L_HS**2*Rep(ElLeft)*2
    allocate(SigmaL(no_L,no_L))
    ispin = ispin + no_L**2
    allocate(HAAR(no_R_HS,no_R_HS,Rep(ElRight)))
    allocate(SAAR(no_R_HS,no_R_HS,Rep(ElRight)))
    ispin = ispin + no_R_HS**2*Rep(ElRight)*2
    allocate(SigmaR(no_R,no_R))
    ispin = ispin + no_R**2
    call memory('A','Z',ispin,'transiesta')

    if ( IsVolt ) then
       ! We need Gamma's with voltages (now they are both GAA and GammaT)
       allocate(GammaLT(no_L,no_L),GammaRT(no_R,no_R))
       call memory('A','Z',no_L**2+no_R**2,'transiesta')
    else
       ! Now they are actually only the GAA arrays...
       allocate(GammaLT(no_L_HS,no_L), GammaRT(no_R_HS,no_R))
       call memory('A','Z',no_L_HS*no_L+no_R_HS*no_R,'transiesta')
    end if
    ! This seems stupid, however, we never use the GAAL and
    ! GammaL at the same time. Hence it will be safe
    ! to have them point to the same array.
    ! When the UC_expansion_Sigma_GammaT is called:
    ! first the GAA is "emptied" of information and then
    ! Gamma is filled.
    GAAL => GammaLT
    GAAR => GammaRT

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
       call newdSpData1D(tsup_sp_uc,fdist,spDM,name='TS spDM')
       call newdSpData1D(tsup_sp_uc,fdist,spEDM,name='TS spEDM')
       if ( IsVolt ) then
          call newdSpData1D(tsup_sp_uc,fdist,spDMR,name='TS spDM-R')
          call newdSpData1D(tsup_sp_uc,fdist,spDMneqL,name='TS spDMneq-L')
          call newdSpData1D(tsup_sp_uc,fdist,spDMneqR,name='TS spDMneq-R')
          call newdSpData1D(tsup_sp_uc,fdist,spEDMR,name='TS spEDM-R')
       end if

       ! The Hamiltonian and overlap matrices (in Gamma calculations
       ! we will not have any phases, hence, it makes no sense to
       ! have the arrays in complex)
       call newdSpData1D(ts_sp_uc,fdist,spH,name='TS spH')
       call newdSpData1D(ts_sp_uc,fdist,spS,name='TS spS')

    else
       call newzSpData1D(tsup_sp_uc,fdist,spzDM,name='TS spzDM')
       call newzSpData1D(tsup_sp_uc,fdist,spzEDM,name='TS spzEDM')
       if ( IsVolt ) then
          call newzSpData1D(tsup_sp_uc,fdist,spzDMR,name='TS spzDM-R')
          call newzSpData1D(tsup_sp_uc,fdist,spzDMneqL,name='TS spzDMneq-L')
          call newzSpData1D(tsup_sp_uc,fdist,spzDMneqR,name='TS spzDMneq-R')
          call newzSpData1D(tsup_sp_uc,fdist,spzEDMR,name='TS spzEDM-R')
       end if

       ! The Hamiltonian and overlap matrices
       call newzSpData1D(ts_sp_uc,fdist,spzH,name='TS spzH')
       call newzSpData1D(ts_sp_uc,fdist,spzS,name='TS spzS')

    end if
    ! We will not write out all created sparsity patterns, it provides
    ! no purpose... other than displaying how much memory this reduces :)

    ! We just check that the work for reducing the matrices can be made
    if ( nzwork < nnzs(ts_sp_uc) ) call die('The memory for transiesta cannot &
         &sustain the implementation, contact the developers.')

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
          call init_val(spDM)
          call init_val(spEDM)
          if ( IsVolt ) then
             call init_val(spDMR)
             call init_val(spDMneqL)
             call init_val(spDMneqR)
             call init_val(spEDMR)
          end if
       else
          call init_val(spzDM)
          call init_val(spzEDM)
          if ( IsVolt ) then
             call init_val(spzDMR)
             call init_val(spzDMneqL)
             call init_val(spzDMneqR)
             call init_val(spzEDMR)
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

             ! for these contour parts we do not require to calculate
             ! Gamma's.
             ! Hence we can perform the calculation without 
             ! calculating them.

             ! Calculate the left-right Sigma
             if ( UseBulk ) then

                call UC_expansion_Sigma_Bulk(no_L_HS, no_L, &
                     RepA1(ElLeft), RepA2(ElLeft), &
                     na_L_HS, lasto_L, nqL, qLb, wqL, HAAL, SAAL, GAAL, SigmaL, &
                     nzwork,zwork)
                
                call UC_expansion_Sigma_Bulk(no_R_HS, no_R, &
                     RepA1(ElRight), RepA2(ElRight), &
                     na_R_HS, lasto_R, nqR, qRb, wqR, HAAR, SAAR, GAAR, SigmaR, &
                     nzwork,zwork)

             else
                call UC_expansion_Sigma(Z,no_L_HS, no_L, &
                     RepA1(ElLeft), RepA2(ElLeft), &
                     na_L_HS, lasto_L, nqL, qLb, wqL, HAAL, SAAL, GAAL, SigmaL, &
                     nzwork,zwork)

                call UC_expansion_Sigma(Z,no_R_HS, no_R, &
                     RepA1(ElRight), RepA2(ElRight), &
                     na_R_HS, lasto_R, nqR, qRb, wqR, HAAR, SAAR, GAAR, SigmaR, &
                     nzwork,zwork)
             end if

#ifdef TRANSIESTA_TIMING
             call timer('TS_PREPG',1)
#endif

             if ( ts_Gamma_SCF ) then
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_D(UseBulk,spH , spS,Z,no_BufL, &
                     no_u_TS,zwork, &
                     no_L, SigmaL, no_R, SigmaR)
             else
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_Z(UseBulk,spzH,spzS,Z,no_BufL, &
                     no_u_TS,zwork, &
                     no_L, SigmaL, no_R, SigmaR)
             end if

#ifdef TRANSIESTA_TIMING
             call timer('TS_PREPG',2)
#endif

             if ( GF_INV_EQUI_PART ) then
                ! Calculate the GF
                call calc_GF_Part(no_u_TS,no_L,no_R, zwork, GF,ierr)
                no_GF_offset = no_L
                ! The size of the central region (without left-right electrodes)
                no_u_C = no_u_TS - no_R - no_L
             else
                ! Calculate the full GF
                call calc_GF(no_u_TS, zwork, GF,ierr)
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

             ! The non-equilibrium integration points have the density
             ! in the real part of the Gf.Gamma.Gf^\dagger
             ! Hence we simply multiply W by -i to move the density
             ! to the same scheme i.e. \rho = - Im(Gf.Gamma.Gf^\dagger)
             ! Hence the signs on "add_DM" are the same as the integration
             ! signs
             W  = zmi * W
             ZW = Z * W
             
             ! Do the left electrode
             call UC_expansion_Sigma_GammaT(UseBulk,Z,no_L_HS,no_L, &
                  RepA1(ElLeft), RepA2(ElLeft), &
                  na_L_HS,lasto_L,nqL,qLb,wqL, &
                  HAAL, SAAL, GAAL, &
                  SigmaL, GammaLT, & 
                  nzwork, zwork)

             ! Do the right electrode
             call UC_expansion_Sigma_GammaT(UseBulk,Z,no_R_HS,no_R, &
                  RepA1(ElRight), RepA2(ElRight), &
                  na_R_HS,lasto_R,nqR,qRb,wqR, &
                  HAAR, SAAR, GAAR, &
                  SigmaR, GammaRT, & 
                  nzwork, zwork)

#ifdef TRANSIESTA_TIMING
             call timer('TS_PREPG',1)
#endif

             if ( ts_Gamma_SCF ) then
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_D(UseBulk,spH , spS,Z,no_BufL, &
                     no_u_TS,zwork, &
                     no_L, SigmaL, no_R, SigmaR)
             else
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_Z(UseBulk,spzH,spzS,Z,no_BufL, &
                     no_u_TS,zwork, &
                     no_L, SigmaL, no_R, SigmaR)
             end if

#ifdef TRANSIESTA_TIMING
             call timer('TS_PREPG',2)
#endif

             ! Calculate the Greens function
             call calc_GF_Bias(no_u_TS,no_L,no_R, zwork, GF,ierr)

             ! We calculate the right contribution
             call GF_Gamma_GF(no_L+1,no_u_TS,no_L+no_R, &
                  no_R, Gf, &
                  GammaRT, zwork, no_R*no_R, SigmaR) ! SigmaR is a "work" array
             ! zwork is now GFGGF

             ! Note that we use '--' here
             ! The integral is defined with:
             !   DMneqL = ... [n_F(\epsilon - \mu_L) - n_F(\epsilon - \mu_R)]
             ! and the weight is defined with:
             !   n_F(\epsilon - \mu_L) - n_F(\epsilon - \mu_R)
             ! hence the weight function on the right has to be the negative.
             if ( ts_Gamma_SCF ) then
                call add_DM_dE_D(spDMneqR,spEDM, no_u_TS, no_u_TS, &
                     zwork, no_BufL, 0, -W, -ZW)
             else
                call add_DM_dE_Z(spzDMneqR,spzEDM, no_u_TS, no_u_TS, &
                     zwork, no_BufL, 0, -W, -ZW)
             end if

             ! We calculate the left contribution
             call GF_Gamma_GF(1,no_u_TS, no_L+no_R, &
                  no_L, Gf, &
                  GammaLT, zwork, no_L*no_L, SigmaL) ! SigmaL is a "work" array
             ! zwork is now GF.G.GF

             ! Note that we use '++' here
             ! See sign convention above...
             if ( ts_Gamma_SCF ) then
                call add_DM_dE_D(spDMneqL,spEDMR, no_u_TS, no_u_TS, &
                     zwork, no_BufL, 0, W, ZW)
             else
                call add_DM_dE_Z(spzDMneqL,spzEDMR, no_u_TS, no_u_TS, &
                     zwork, no_BufL, 0, W, ZW)
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
      ! We insert a barrier to better see the actual
      ! communication times
      call MPI_Barrier(MPI_Comm_World,ind)

! Global reduction of density matrices
! this completes the energy integration of each DM

      call timer("TS_comm",1)

      if ( ts_Gamma_SCF ) then
         call AllReduce_dSpData1D(spDM , ndwork,dwork)
         call AllReduce_dSpData1D(spEDM, ndwork,dwork)
         if ( IsVolt ) then
            call AllReduce_dSpData1D(spDMR   , ndwork,dwork)
            call AllReduce_dSpData1D(spDMneqL, ndwork,dwork)
            call AllReduce_dSpData1D(spDMneqR, ndwork,dwork)
            call AllReduce_dSpData1D(spEDMR  , ndwork,dwork)
         end if
      else
         call AllReduce_zSpData1D(spzDM , nzwork,zwork)
         call AllReduce_zSpData1D(spzEDM, nzwork,zwork)
         if ( IsVolt ) then
            call AllReduce_zSpData1D(spzDMR   , nzwork,zwork)
            call AllReduce_zSpData1D(spzDMneqL, nzwork,zwork)
            call AllReduce_zSpData1D(spzDMneqR, nzwork,zwork)
            call AllReduce_zSpData1D(spzEDMR  , nzwork,zwork)
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
         call timer('ts_weight',1)
         if ( ts_Gamma_SCF ) then
            call weightDM(no_C_L, no_C_R, &
                 spDM,   spDMR,  spDMneqL,  spDMneqR, &
                 spEDM,  spEDMR )
         else
            call weightDMC(no_C_L, no_C_R, &
                 spzDM, spzDMR, spzDMneqL, spzDMneqR, &
                 spzEDM, spzEDMR)
         end if
         call timer('ts_weight',2)
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

    call memory('D','Z',size(zwork)+size(GF),'transiesta')
    deallocate(zwork,GF)

    deallocate(lasto_L,lasto_R)
    call memory('D','I',na_R_HS+na_L_HS+2,'transiesta')

    call memory('D','Z',size(HAAL)*2+size(HAAR)*2,'transiesta')
    deallocate(HAAL,SAAL)
    deallocate(HAAR,SAAR)
    call memory('D','Z',size(SigmaL)+size(SigmaR),'transiesta')
    deallocate(SigmaL,SigmaR)
    
    ! These are allocated instead of the GAA[LR] arrays.
    ! Hence they are used in both non-bias and bias calculations.
    call memory('D','Z',size(GammaLT)+size(GammaRT),'transiesta')
    deallocate(GammaLT,GammaRT)
    
    ! I would like to add something here which enables 
    ! 'Transiesta' post-process

    call timer('TS_calc',2)

    call ts_print_charges(sp_dist, sparse_pattern, na_u, lasto, &
         nspin, n_nzs, DM, Ss, TS_INFO_SCF)

#ifdef TRANSIESTA_DEBUG
   call write_debug( 'POS transiesta mem' )
#endif

 end subroutine transiesta_mem


! Update DM
! These routines are supplied for easy update of the update region
! sparsity patterns
! Note that these routines implement the usual rho(Z) \propto - GF
   subroutine add_DM_dE_Z(DM,EDM,no1,no2,GF,no_BufL,GF_offset,DMfact,EDMfact)
     use class_zSpData1D
     use class_Sparsity
     ! The DM and EDM equivalent matrices
     type(zSpData1D), intent(inout) :: DM,EDM
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
     integer :: io, ind, nr
     integer :: iu, ju
     
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
     use class_dSpData1D
     use class_Sparsity
     ! The DM and EDM equivalent matrices
     type(dSpData1D), intent(inout) :: DM,EDM
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
     integer :: io, ind, nr
     integer :: iu, ju

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
     
     do io = 1 , nr ! TODO introduce reduced loop
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


   ! creation of the GF^{-1}.
   ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
   subroutine prepare_GF_inv_D(UseBulk,spH,spS, ZE, no_BufL,no_u,GFinv, &
        no_L, SigmaL, no_R, SigmaR)
     use class_dSpData1D
     use class_Sparsity
     logical, intent(in) :: UseBulk
     ! The Hamiltonian and overlap sparse matrices
     type(dSpData1D), intent(inout) :: spH, spS
     ! the current energy point
     complex(dp), intent(in) :: ZE
     ! Remark that we need the left buffer orbitals
     ! to calculate the actual orbital of the sparse matrices...
     integer, intent(in) :: no_BufL,no_u
     complex(dp), intent(out) :: GFinv(no_u**2)
     integer, intent(in) :: no_L, no_R
     complex(dp), intent(in) :: SigmaL(no_L**2)
     complex(dp), intent(in) :: SigmaR(no_R**2)

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

           ! Notice that we transpose S and H back here
           ! See symmetrize_HS_Gamma (H is hermitian)
           GFinv(l_col(ind)+iu) = ZE * dS(ind) - dH(ind)
        end do
     end do

     call insert_Self_Energies(UseBulk,no_u,Gfinv, &
          no_L, SigmaL, no_R, SigmaR)

   end subroutine prepare_GF_inv_D

   ! creation of the GF^{-1}.
   ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
   subroutine prepare_GF_inv_Z(UseBulk, spH,spS,ZE, no_BufL,no_u,GFinv, &
        no_L, SigmaL, no_R, SigmaR)
     use class_zSpData1D
     use class_Sparsity
     logical, intent(in) :: UseBulk
     ! The Hamiltonian and overlap sparse matrices
     type(zSpData1D), intent(inout) :: spH, spS
     ! the current energy point
     complex(dp), intent(in) :: ZE
     ! Remark that we need the left buffer orbitals
     ! to calculate the actual orbital of the sparse matrices...
     integer, intent(in) :: no_BufL,no_u
     complex(dp), intent(out) :: GFinv(no_u**2)
     integer, intent(in) :: no_L, no_R
     complex(dp), intent(in) :: SigmaL(no_L**2)
     complex(dp), intent(in) :: SigmaR(no_R**2)

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
           ! See symmetrize_HS_kpt
           GFinv(l_col(ind)+iu) = ZE * zS(ind) - zH(ind)
        end do
     end do

     call insert_Self_Energies(UseBulk,no_u,Gfinv, &
          no_L, SigmaL, no_R, SigmaR)

   end subroutine prepare_GF_inv_Z
   
   subroutine insert_Self_Energies(UseBulk,no_u,Gfinv, &
        no_L, SigmaL, no_R, SigmaR)
     logical, intent(in) :: UseBulk
     integer, intent(in) :: no_u
     complex(dp), intent(in out) :: GFinv(no_u,no_u)
     integer, intent(in) :: no_L, no_R
     complex(dp), intent(in) :: SigmaL(no_L,no_L)
     complex(dp), intent(in) :: SigmaR(no_R,no_R)
     
     integer :: i, j, jj, off

     off = no_u - no_R
     if ( UseBulk ) then
        do j = 1 , no_L
           do i = 1 , no_L
              Gfinv(i,j) = SigmaL(i,j)
           end do
        end do
        do j = 1 , no_R
           do i = 1 , no_R
              Gfinv(off+i,off+j) = SigmaR(i,j)
           end do
        end do
     else
        do j = 1 , no_L
           do i = 1 , no_L
              Gfinv(i,j) = Gfinv(i,j) - SigmaL(i,j)
           end do
        end do
        do j = 1 , no_R
           jj = off + j
           do i = 1 , no_R
              Gfinv(off+i,jj) = Gfinv(off+i,jj) - SigmaR(i,j) 
           end do
        end do
     end if

   end subroutine insert_Self_Energies

end module m_ts_mem
