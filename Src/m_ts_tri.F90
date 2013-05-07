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

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.

module m_ts_tri

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
  
  implicit none

  
  public :: transiesta_tri

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

  subroutine transiesta_tri(nspin, &
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
    use class_zTriMat3

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

    use m_ts_method, only : GF_INV_EQUI_PART

    use m_ts_contour,only : PNEn, NEn, contour
    use m_ts_cctype, only : CC_PART_EQUI
    use m_ts_cctype, only : CC_PART_LEFT_EQUI
    use m_ts_cctype, only : CC_PART_RIGHT_EQUI
    use m_ts_cctype, only : CC_PART_NON_EQUI
    use m_ts_cctype, only : CC_PART_TRANSPORT

    use m_ts_gf, only : read_Green

    use m_ts_tri_scat

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
    complex(dp), allocatable :: HAAL(:,:,:), SAAL(:,:,:)
    complex(dp), allocatable :: HAAR(:,:,:), SAAR(:,:,:)
    complex(dp), pointer :: GAAL(:,:), GAAR(:,:)
    complex(dp), allocatable :: SigmaL(:,:), SigmaR(:,:)
    complex(dp), target, allocatable :: GammaLT(:,:), GammaRT(:,:)
! ************************************************************

! ******************* Computational arrays *******************
    integer :: ndwork, nzwork
    real(dp),    allocatable :: dwork(:)
    complex(dp), pointer :: GF22(:)
    complex(dp), pointer :: zwork(:)
    type(zTriMat3) :: zwork_tri, GF_tri
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
    integer :: i, j, ind
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

    ! We must ensure that the Sigma[LR] have enough space to hold
    ! one line of the full matrix (we use them as work arrays
    ! when calculating the Gamma[LR]
    if ( no_L ** 2 < no_u_TS .or. no_R ** 2 < no_u_TS ) then
       call die('The current implementation requires that the &
            &square of the orbitals in the electrodes are larger &
            &than the dimension of the problem.')
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
    no_u_C = no_u_TS - no_L - no_R
    nzwork =          (no_L + no_u_C)        * no_L
    nzwork = nzwork + (no_L + no_u_C + no_R) * no_u_C
    nzwork = nzwork + (       no_u_C + no_R) * no_R
    call newzTriMat3(zwork_tri,no_L,no_u_C,no_R,'GFinv')

    ! Save the work-space
    ! Now the programmer should "keep a straight tongue"
    ! The zwork points to the array in the zwork_tri
    ! tri-diagonal array. This means that there are two
    ! arrays that point to the same.
    ! Generally the zwork need only to retain the value in
    ! one call!
    zwork => val(zwork_tri)

    call newzTriMat3(GF_tri,no_L,no_u_C,no_R,'GF')
    if ( GF_INV_EQUI_PART ) then
       Gf22 => val22(Gf_tri)
    end if

    ! Allocate the left-right electrode quantities that we need
    allocate(HAAL(no_L_HS,no_L_HS,NRepA1L*NRepA2L))
    allocate(SAAL(no_L_HS,no_L_HS,NRepA1L*NRepA2L))
    ispin = no_L_HS**2*NRepA1L*NRepA2L*2
    allocate(SigmaL(no_L,no_L))
    ispin = ispin + no_L**2
    allocate(HAAR(no_R_HS,no_R_HS,NRepA1R*NRepA2R))
    allocate(SAAR(no_R_HS,no_R_HS,NRepA1R*NRepA2R))
    ispin = ispin + no_R_HS**2*NRepA1R*NRepA2R*2
    allocate(SigmaR(no_R,no_R))
    ispin = ispin + no_R**2
    call memory('A','Z',ispin,'transiesta')

    if ( IsVolt ) then
       ! We only need Gamma's with voltages
       allocate(GammaLT(no_L,no_L),GammaRT(no_R,no_R))
       call memory('A','Z',no_L**2+no_R**2,'transiesta')
    else
       allocate(GammaLT(no_L_HS,no_L),GammaRT(no_R_HS,no_R))
       call memory('A','Z',no_L_HS*no_L+no_R_HS*no_R,'transiesta')
    end if
    ! This seems stupid, however, we never use the GAAL and
    ! GammaL at the same time. Hence it will be safe
    ! to have them point to the same array.
    ! When the UC_expansion_Sigma_GammaT is called
    ! first the GAA is "emptied" of information" and then
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
             call timer('TS_PREPG',1)
#endif

             if ( ts_Gamma_SCF ) then
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_D(UseBulk, spH , spS,Z,no_BufL, &
                     no_u_TS,zwork_tri, &
                     no_L, SigmaL, no_R, SigmaR)
             else
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_Z(UseBulk, spzH,spzS,Z,no_BufL, &
                     no_u_TS,zwork_tri, &
                     no_L, SigmaL, no_R, SigmaR)
             end if

#ifdef TRANSIESTA_TIMING
             call timer('TS_PREPG',2)
#endif

             if ( GF_INV_EQUI_PART ) then
                ! Calculate the GF22 (note that GF22 points to the 
                ! tri-diag array...
                call calc_GF_Part(no_u_TS, no_L,no_R,zwork_tri, GF22, ierr)
             else
                ! Calculate the full GF
                call calc_GF(.false., no_u_TS, zwork_tri, GF_tri, ierr)
             end if

             if ( contour(iE)%part == CC_PART_RIGHT_EQUI ) then
                ! We have the right equilibrium contour
                if ( ts_Gamma_SCF ) then
                   call add_DM_dE_D(spDMR , spEDMR, &
                        GF_tri, no_BufL, W, ZW)
                else
                   call add_DM_dE_Z(spzDMR, spzEDMR, &
                        GF_tri, no_BufL, W, ZW)
                end if
             else
                ! We have the left- or the equilibrium contour...
                if ( ts_Gamma_SCF ) then
                   call add_DM_dE_D(spDM ,  spEDM, &
                        GF_tri, no_BufL, W, ZW)
                else
                   call add_DM_dE_Z(spzDM, spzEDM, &
                        GF_tri, no_BufL, W, ZW)
                end if
             end if

          case ( CC_PART_NON_EQUI )
             
             ! The non-equilibrium integration points have the density
             ! in the real part of the Gf.Gamma.Gf^\dagger
             ! Hence we simply multiply W by -i to move the density
             ! to the same scheme i.e. \rho = - Im(Gf.Gamma.Gf^\dagger)
             W  = zmi * W
             ZW = Z * W

             ! Do the left electrode
             call UC_expansion_Sigma_GammaT(UseBulk,Z,no_L_HS,no_L, &
                  NRepA1L, NRepA2L, &
                  na_L_HS,lasto_L,nqL,qLb,wqL, &
                  HAAL, SAAL, GAAL, &
                  SigmaL, GammaLT, & 
                  nzwork, zwork)

             ! Do the right electrode
             call UC_expansion_Sigma_GammaT(UseBulk,Z,no_R_HS,no_R, &
                  NRepA1R, NRepA2R, &
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
                call prepare_GF_inv_D(UseBulk, spH , spS,Z,no_BufL, &
                     no_u_TS,zwork_tri, &
                     no_L, SigmaL, no_R, SigmaR)
             else
                ! Notice that we now actually need to retain the values
                ! in zwork...!!!
                call prepare_GF_inv_Z(UseBulk, spzH,spzS,Z,no_BufL, &
                     no_u_TS,zwork_tri, &
                     no_L, SigmaL, no_R, SigmaR)
             end if

#ifdef TRANSIESTA_TIMING
             call timer('TS_PREPG',2)
#endif

             ! Calculate the Greens function
             call calc_GF_Bias(UseBulk, &
                     no_u_TS,zwork_tri,GF_tri, &
                     no_L, SigmaL, & ! These are work-arrays here...
                     no_R, SigmaR)

             ! We calculate the right thing.
             call GF_Gamma_GF_Right_All(no_R, Gf_tri, GammaRT, zwork_tri)
             ! work is now GFGGF

             ! Note that we use '--' here
             if ( ts_Gamma_SCF ) then
                call add_DM_dE_D(spDMneqR,spEDM, &
                     zwork_tri, no_BufL, -W, -ZW)
             else
                call add_DM_dE_Z(spzDMneqR,spzEDM, &
                     zwork_tri, no_BufL, -W, -ZW)
             end if

             ! We calculate the left thing.
             call GF_Gamma_GF_Left_All(no_L, Gf_tri,GammaLT, zwork_tri)
             ! work is now GFGGF

             ! Note that we use '++' here
             if ( ts_Gamma_SCF ) then
                call add_DM_dE_D(spDMneqL,spEDMR, &
                     zwork_tri, no_BufL, +W, +ZW)
             else
                call add_DM_dE_Z(spzDMneqL,spzEDMR, &
                     zwork_tri, no_BufL, +W, +ZW)
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
         ind = nnzs(spDM)
         call AllReduce_dSpData1D(spDM ,ind, ndwork,dwork)
         call AllReduce_dSpData1D(spEDM,ind, ndwork,dwork)
         if ( IsVolt ) then
            call AllReduce_dSpData1D(spDMR   ,ind, ndwork,dwork)
            call AllReduce_dSpData1D(spDMneqL,ind, ndwork,dwork)
            call AllReduce_dSpData1D(spDMneqR,ind, ndwork,dwork)
            call AllReduce_dSpData1D(spEDMR  ,ind, ndwork,dwork)
         end if
      else
         ind = nnzs(spzDM)
         call AllReduce_zSpData1D(spzDM ,ind, nzwork,zwork)
         call AllReduce_zSpData1D(spzEDM,ind, nzwork,zwork)
         if ( IsVolt ) then
            call AllReduce_zSpData1D(spzDMR   ,ind, nzwork,zwork)
            call AllReduce_zSpData1D(spzDMneqL,ind, nzwork,zwork)
            call AllReduce_zSpData1D(spzDMneqR,ind, nzwork,zwork)
            call AllReduce_zSpData1D(spzEDMR  ,ind, nzwork,zwork)
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

    call delete(zwork_tri)
    call delete(GF_tri)

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

  end subroutine transiesta_tri

! Update DM
! These routines are supplied for easy update of the update region
! sparsity patterns
! Note that these routines implement the usual rho(Z) \propto - GF
  subroutine add_DM_dE_Z(DM,EDM,GF_tri,no_BufL,DMfact,EDMfact)
    use class_zSpData1D
    use class_Sparsity
    use class_zTriMat3
    ! The DM and EDM equivalent matrices
    type(zSpData1D), intent(inout) :: DM,EDM
    ! The Green's function
    type(zTriMat3), intent(inout) :: GF_tri
    ! The number of buffer atoms (needed for the offset in the sparsity
    ! patterns)
    integer, intent(in) :: no_BufL
    ! Complex numbers that are used in the factor of GF
    complex(dp), intent(in) :: DMfact, EDMfact

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: zD(:), zE(:), Gf(:)
    integer :: io, ind, nr, iu, ju, idx

    s      => spar(DM)
    l_ncol => n_col   (s)
    l_ptr  => list_ptr(s)
    l_col  => list_col(s)
    zD     => val(DM)
    zE     => val(EDM)
    Gf     => val(Gf_tri)

    ! Number of orbitals in the SIESTA unit-cell
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    nr = nrows(s)
    
    do io = 1 , nr
       ! Quickly go past the buffer atoms... (the right side)
       if ( l_ncol(io) == 0 ) cycle

       ! The update region equivalent GF part
       iu = io - no_BufL
       
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          idx = index(Gf_tri,iu,l_col(ind) - no_BufL)

          zD(ind) = zD(ind) - GF(idx) * DMfact
          zE(ind) = zE(ind) - GF(idx) * EDMfact
       end do
    end do

  end subroutine add_DM_dE_Z

  subroutine add_DM_dE_D(DM,EDM,GF_tri,no_BufL,DMfact,EDMfact)
    use class_dSpData1D
    use class_Sparsity
    use class_zTriMat3
    ! The DM and EDM equivalent matrices
    type(dSpData1D), intent(inout) :: DM,EDM
    ! The Green's function
    type(zTriMat3), intent(inout) :: GF_tri
    ! The number of buffer atoms (needed for the offset in the sparsity
    ! patterns)
    integer, intent(in) :: no_BufL
    ! Complex numbers that are used in the factor of GF
    complex(dp), intent(in) :: DMfact, EDMfact

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: dD(:), dE(:)
    complex(dp), pointer :: Gf(:)
    integer :: io, ind, nr, iu, idx

    s      => spar(DM)
    l_ncol => n_col   (s)
    l_ptr  => list_ptr(s)
    l_col  => list_col(s)
    dD     => val(DM)
    dE     => val(EDM)
    Gf     => val(Gf_tri)

    ! Notice that we do not need to do any transposing here...
    ! The tri-diagonal calculation of GF_Gamma_GF will always be correct

    ! Number of orbitals in the SIESTA unit-cell
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    nr = nrows(s)
    
    do io = 1 , nr !TODO introduce reduced loop
       ! Quickly go past the buffer atoms... (in the right side)
       if ( l_ncol(io) == 0 ) cycle

       ! The update region equivalent GF part
       iu = io - no_BufL
       
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          idx = index(Gf_tri,iu,l_col(ind) - no_BufL)
    
          dD(ind) = dD(ind) - dimag( GF(idx) * DMfact  )
          dE(ind) = dE(ind) - dimag( GF(idx) * EDMfact )

       end do
    end do

  end subroutine add_DM_dE_D


  ! creation of the GF^{-1}.
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_GF_inv_D(UseBulk, spH,spS, ZE, no_BufL,no_u,GFinv_tri, &
       no_L, SigmaL, no_R, SigmaR)
    use class_dSpData1D
    use class_Sparsity
    use class_zTriMat3

    logical, intent(in) :: UseBulk
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData1D), intent(inout) :: spH, spS
    ! the current energy point
    complex(dp), intent(in) :: ZE
    ! Remark that we need the left buffer orbitals
    ! to calculate the actual orbital of the sparse matrices...
    integer, intent(in) :: no_BufL,no_u
    type(zTriMat3), intent(inout) :: GFinv_tri
    integer, intent(in) :: no_L, no_R
    complex(dp), intent(in) :: SigmaL(no_L**2), SigmaR(no_R**2)

    ! Local variables
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: dH(:), dS(:)
    complex(dp), pointer :: Gfinv(:)
    integer :: io, iu, ind, idx

    s      => spar    (spH)
    l_ncol => n_col   (s)
    l_ptr  => list_ptr(s)
    l_col  => list_col(s)

    dH     => val(spH)
    dS     => val(spS)
    Gfinv  => val(Gfinv_tri)

    ! initialize 
    GFinv(:) = dcmplx(0._dp,0._dp)

    ! We will only loop in the central region
    ! We have constructed the sparse array to only contain
    ! values in this part...
    do io = no_BufL + 1, no_BufL + no_u

       iu = io - no_BufL

       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          ! we could transpose... but...
          idx = index(Gfinv_tri, iu, l_col(ind)-no_BufL)

          GFinv(idx) = ZE * dS(ind) - dH(ind)
       end do
    end do

    call insert_Self_Energies(UseBulk, Gfinv_tri, &
         no_L, SigmaL, no_R, SigmaR)

  end subroutine prepare_GF_inv_D

  ! creation of the GF^{-1}.
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_GF_inv_Z(UseBulk,spH,spS,ZE, no_BufL,no_u,GFinv_tri, &
       no_L, SigmaL, no_R, SigmaR)
    use class_zSpData1D
    use class_Sparsity
    use class_zTriMat3

    logical, intent(in) :: UseBulk
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D), intent(inout) :: spH, spS
    ! the current energy point
    complex(dp), intent(in) :: ZE
    ! Remark that we need the left buffer orbitals
    ! to calculate the actual orbital of the sparse matrices...
    integer, intent(in) :: no_BufL,no_u
    type(zTriMat3), intent(inout) :: GFinv_tri
    integer, intent(in) :: no_L, no_R
    complex(dp), intent(in) :: SigmaL(no_L**2), SigmaR(no_R**2)

    ! Local variables
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: zH(:), zS(:), Gfinv(:)
    integer :: io, iu,ind, idx
    
    s      => spar    (spH)
    l_ncol => n_col   (s)
    l_ptr  => list_ptr(s)
    l_col  => list_col(s)
    
    zH     => val(spH)
    zS     => val(spS)
    Gfinv  => val(Gfinv_tri)

    ! Initialize
    GFinv(:) = dcmplx(0._dp,0._dp)

    ! We will only loop in the central region
    ! We have constructed the sparse array to only contain
    ! values in this part...
    do io = no_BufL + 1, no_BufL + no_u

       iu = io - no_BufL

       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 

          !ju = l_col(ind) ! The '- no_BufL' is moved outside the loop
          idx = index(Gfinv_tri,l_col(ind)-no_BufL,iu)

          ! Notice that we transpose back here...
          ! See symmetrize_HS_kpt
          GFinv(idx) = ZE * zS(ind) - zH(ind)
       end do
    end do

    call insert_Self_Energies(UseBulk, Gfinv_tri, &
         no_L, SigmaL, no_R, SigmaR)
    
  end subroutine prepare_GF_inv_Z

  subroutine insert_Self_Energies(UseBulk,Gfinv_tri, &
       no_L, SigmaL, no_R, SigmaR)
    use class_zTriMat3

    logical, intent(in) :: UseBulk
    type(zTriMat3), intent(inout) :: GFinv_tri
    integer, intent(in) :: no_L, no_R
    complex(dp), intent(in) :: SigmaL(no_L**2)
    complex(dp), intent(in) :: SigmaR(no_R**2)

    complex(dp), pointer :: Gfpart(:)
    integer :: i

    Gfpart => val11(GFinv_tri)
    
    if ( UseBulk ) then
       Gfpart(:) = SigmaL(:)
    else
       do i = 1 , no_L**2
          Gfpart(i) = Gfpart(i) - SigmaL(i)
       end do
    end if

    Gfpart => val33(GFinv_tri)
    
    if ( UseBulk ) then
       Gfpart(:) = SigmaR(:)
    else
       do i = 1 , no_R**2
          Gfpart(i) = Gfpart(i) - SigmaR(i)
       end do
    end if

  end subroutine insert_Self_Energies


  subroutine test_tri(tri)
    use class_zTriMat3
    type(zTriMat3), intent(inout) :: tri
    integer :: i,j,ind,nL,nC,nR
    nL =nrows_g_left(tri)
    nC =nrows_g_center(tri)
    nR =nrows_g_right(tri)
    ind = 0
    do j = 1 , nL
       do i = 1 , nL
          ind = ind + 1
          if ( ind /= index(tri,i,j) ) then
             print *,'got',index(tri,i,j)
             print *,'expected',ind
             call die('wrong 11')
          end if
       end do
    end do
    do j = 1 , nL
       do i = nL+1 , nL+nC
          ind = ind + 1
          if ( ind /= index(tri,i,j) ) then
             print *,'got',index(tri,i,j)
             print *,'expected',ind
             call die('wrong 21')
          end if
       end do
    end do
    do j = nL+1 , nL+nC
       do i = 1 , nL
          ind = ind + 1
          if ( ind /= index(tri,i,j) ) then
             print *,'got',index(tri,i,j)
             print *,'expected',ind
             call die('wrong 12')
          end if
       end do
    end do
    do j = nL+1 , nL+nC
       do i = nL+1 , nL+nC
          ind = ind + 1
          if ( ind /= index(tri,i,j) ) then
             print *,'got',index(tri,i,j)
             print *,'expected',ind
             call die('wrong 22')
          end if
       end do
    end do
    do j = nL+1 , nL+nC
       do i = nL+nC+1 , nL+nC+nR
          ind = ind + 1
          if ( ind /= index(tri,i,j) ) then
             print *,'got',index(tri,i,j)
             print *,'expected',ind
             call die('wrong 32')
          end if
       end do
    end do
    do j = nL+nC+1 , nL+nC+nR
       do i = nL+1 , nL+nC
          ind = ind + 1
          if ( ind /= index(tri,i,j) ) then
             print *,'got',index(tri,i,j)
             print *,'expected',ind
             call die('wrong 23')
          end if
       end do
    end do
    do j = nL+nC+1 , nL+nC+nR
       do i = nL+nC+1 , nL+nC+nR
          ind = ind + 1
          if ( ind /= index(tri,i,j) ) then
             print *,'got',index(tri,i,j)
             print *,'expected',ind
             call die('wrong 33')
          end if
       end do
    end do
    print *,'successfull'
  end subroutine test_tri

  subroutine zTriMat3_transpose(m_tri,mT_tri)
    use class_zTriMat3
    
    ! The input matrix
    type(zTriMat3), intent(inout) :: m_tri
    ! The returned transposed matrix...
    type(zTriMat3), intent(inout) :: mT_tri

    ! Local variables
    complex(dp), pointer :: m(:)
    complex(dp), pointer :: mT(:)
    integer :: nL, nC, nR
    integer :: i, j, ii

    nL = nrows_g_left  (m_tri)
    nC = nrows_g_center(m_tri)
    nR = nrows_g_right (m_tri)

    ! Transpose m11
    m  => val11(m_tri)
    mT => val11(mT_tri)
    call t(nL,nL,m,mT)
    ! Transpose m12
    m  => val12(m_tri)
    mT => val21(mT_tri)
    call t(nL,nC,m,mT)
    ! Transpose m21
    m  => val21(m_tri)
    mT => val12(mT_tri)
    call t(nC,nL,m,mT)
    ! Transpose m22
    m  => val22(m_tri)
    mT => val22(mT_tri)
    call t(nC,nC,m,mT)
    ! Transpose m23
    m  => val23(m_tri)
    mT => val32(mT_tri)
    call t(nC,nR,m,mT)
    ! Transpose m32
    m  => val32(m_tri)
    mT => val23(mT_tri)
    call t(nR,nC,m,mT)
    ! Transpose m33
    m  => val33(m_tri)
    mT => val33(mT_tri)
    call t(nR,nR,m,mT)
    ! Done with transposing...

  contains
    subroutine t(N1,N2,m,mt)
      integer, intent(in) :: N1,N2
      complex(dp), intent(in) :: m(n1,n2)
      complex(dp), intent(out) :: mT(n2,n1)
      integer :: i,j
      
      do j = 1 , N2
         do i = 1 , N1
            mT(j,i) = m(i,j)
         end do
      end do

    end subroutine t

  end subroutine zTriMat3_transpose

end module m_ts_tri
