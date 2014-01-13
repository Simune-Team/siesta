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

module m_ts_full

  use precision, only : dp

  use m_ts_sparse_helper, only : create_HS_kpt
  use m_ts_sparse_helper, only : create_HS_Gamma
  use m_ts_sparse_helper, only : symmetrize_HS_kpt
  use m_ts_sparse_helper, only : symmetrize_HS_Gamma

  use m_ts_sparse_helper, only : d_DM_EDM_Reduce_Shift
  use m_ts_sparse_helper, only : z_DM_EDM_Reduce_Shift

  use m_ts_dm_update, only : init_DM
  use m_ts_dm_update, only : update_DM
  use m_ts_dm_update, only : update_zDM
  use m_ts_dm_update, only : add_Gamma_DM
  use m_ts_dm_update, only : add_k_DM
  
  use m_ts_weight, only : TS_W_METHOD
  use m_ts_weight, only : TS_W_CORRELATED
  use m_ts_weight, only : TS_W_UNCORRELATED
  use m_ts_weight, only : TS_W_K_UNCORRELATED
  use m_ts_weight, only : weight_DM
  
  implicit none
  
  public :: transiesta_full
  
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
  
  subroutine transiesta_full(nspin, &
       Gamma, sp_dist, sparse_pattern, &
       ucell, no_u, na_u, lasto, xa, n_nzs, &
       xij, Hs, Ss, DM, EDM, Ef, &
       TSiscf, Qtot)

    use units, only : Pi
    use parallel, only : Node, Nodes, IONode, operator(.PARCOUNT.)
#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use class_zSpData1D

    use m_ts_kpoints

    use m_ts_electype

    use m_ts_options, only : N_Elec, Elecs
    use m_ts_options, only : N_mu, mus
    use m_ts_options, only : na_BufL, no_BufL
    use m_ts_options, only : na_BufR, no_BufR

    use m_ts_options, only : IsVolt

    use m_ts_sparse, only : ts_sp_uc
    use m_ts_sparse, only : tsup_sp_uc
    use m_ts_sparse, only : ltsup_sp_sc
    use m_ts_sparse, only : ltsup_sc_pnt

    ! Self-energy retrival and expansion
    use m_ts_elec_se
    ! Gf calculation
    use m_ts_full_scat

    use m_ts_cctype
    use m_ts_contour, only : nextE
    
    use m_ts_gf, only : read_Green

    use m_ts_cctype

    use m_iterator

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

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
    integer,  intent(in) :: TSiscf

! ******************** IO descriptors ************************
    integer, allocatable :: uGF(:)
! ************************************************************

! ****************** Electrode variables *********************
    integer, allocatable :: nq(:)
    complex(dp), pointer :: GFGGF_work(:) => null()
! ************************************************************

! ******************* Computational arrays *******************
    integer :: ndwork, nzwork
    real(dp),    allocatable :: dwork(:)
    complex(dp), allocatable, target :: zwork(:), GF(:)

    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData1D) :: spH, spS
    type(zSpData1D) :: spzH, spzS
    ! The different sparse matrices that will surmount to the integral
    ! Tohese two lines are in local update sparsity pattern
    type(dSpData2D) ::  spDM,  spDMneq
    type(dSpData2D) :: spEDM
    ! The different sparse matrices... (these two lines are in global update sparsity pattern)
    type(dSpData1D) ::  spDMu,  spEDMu
    type(zSpData1D) :: spzDMu, spzEDMu
    ! Pointers for updating the density matrices
    real(dp),    pointer :: dDM(:), dEDM(:)
    complex(dp), pointer :: zDM(:), zEDM(:), tmp(:)
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c)  :: cE
    complex(dp) :: Z, W, ZW
    real(dp)    :: k(3)
    complex(dp), parameter :: zmi = dcmplx(0._dp,-1._dp)
! ************************************************************

! ******************** Loop variables ************************
    type(itt2) :: SpKp
    integer, pointer :: ispin, ikpt
    integer :: ispin, ikpt, iPE, iE, NEReqs, up_nzs, ia, ia_E
    integer :: ind, El
#ifdef TRANSIESTA_DEBUG
    integer :: iu_GF, iu_GFinv
    integer :: iu_SL, iu_SR
#endif
! ************************************************************

! ******************* Miscalleneous variables ****************
    integer :: ierr
! ************************************************************

#ifdef MPI
! ******************* MPI-related variables   ****************
    type(ts_c) :: nE
    integer :: MPIerror
    logical :: tmp
    logical, pointer :: lhas_elec(:) => null()
    integer, pointer :: has_elec(:,:) => null()
! ************************************************************
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    call timer('TS_calc',1)

    ! Number of orbitals in TranSIESTA
    no_u_TS = no_u - no_BufL - no_BufR

    ! Number of elements that are transiesta updated
    up_nzs = nnzs(tsup_sp_uc)

    ! Open GF files...
    ! Read-in header of Green's functions
    ! Prepare for the calculation
    ! We read in the k-points that the electrode was generated with.
    ! Furthermore we read in the expansion q-points
    ! They are communicated in the routine

    allocate(uGF(N_Elec))
    allocate(nq(N_Elec))
    do iEl = 1 , N_Elec
       if ( IONode ) then
          call io_assign(uGF(iEl))
          open(file=GFFile(Elecs(iEl)),unit=uGF(iEl),form='unformatted')
       end if
       call read_Green(uGF(iEl),Elecs(iEl), ts_nkpnt, NEn, .false. )
       nq(iEl) = Rep(Elecs(iEl))
    end do

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

    ! We only need a partial size of the Green's function
    ind = no_u_TS * (no_u_TS-sum(TotUsedOrbs(Elecs)))
    do iEl = 1 , N_Elec
       if ( Elecs(iEl)%DM_CrossTerms ) then
          ind = ind + no_u_TS * TotUsedOrbs(Elecs(iEl))
       end if
    end do
    ! when bias is needed we need the entire GF column
    ! for all the electrodes (at least some of the contour points needs this)
    if ( IsVolt ) then
       ind = max(ind,no_u_TS * sum(TotUsedOrbs(Elecs)))
    end if
    allocate(GF(ind),stat=ierr)
    if (ierr/=0) call die('Could not allocate space for GFpart')
    call memory('A','Z',ind,'transiesta')

    ! Allocate the electrode quantities
    ind = 0
    do iEl = 1 , N_Elec
       nullify(Elecs(iEl)%HA,Elecs(iEl)%SA,Elecs(iEl)%Gamma)

       ! We allocate for once as much space as needed,

       ! Allocate the non-repeated hamiltonian and overlaps...
       ia = UsedOrbs(Elecs(iEl))
       call re_alloc(Elecs(iEl)%HA,1,ia,1,ia,1,Rep(Elecs(iEl)),routine='transiesta')
       call re_alloc(Elecs(iEl)%SA,1,ia,1,ia,1,Rep(Elecs(iEl)),routine='transiesta')

       ! This seems stupid, however, we never use the Sigma and
       ! GF at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called
       ! first the Sigma is assigned and then 
       ! it is required that prepare_GF_inv is called
       ! immediately (which it is)
       ! Hence the GF must NOT be used in between these two calls!
       ia = TotUsedOrbs(Elecs(iEl))
       Elecs(iEl)%Sigma => GF(ind+1:ind+ia**2)
       ind = ind + ia ** 2

       if ( IsVolt ) then
          ! We need Gamma's with voltages (now they are both GAA and GammaT)
          ind = ia
       else
          ! This is only for having space for GA
          ind = UsedOrbs(Elecs(iEl))
       end if
       call re_alloc(Elecs(iEl)%Gamma,1,ind,1,ia,routine='transiesta')

       ! This seems stupid, however, we never use the GAAL and
       ! GammaT at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called:
       ! first the GAA is "emptied" of information and then
       ! Gamma is filled.
       Elecs(iEl)%GA => Elecs(iEl)%Gamma

    end do

    if ( IsVolt ) then
       ! we need only allocate one work-array for
       ! Gf.G.Gf^\dagger
       call re_alloc(GFGGF_work,1,maxval(TotUsedOrbs(Elecs))**2,routine='transiesta')
    end if

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
    ! Notice that we DO need it to be the SIESTA size.
#ifdef MPI
    call newDistribution(no_u,MPI_COMM_WORLD,fdist,name='TS-fake dist')
#else
    call newDistribution(no_u,-1            ,fdist,name='TS-fake dist')
#endif

    ! The Hamiltonian and overlap matrices (in Gamma calculations
    ! we will not have any phases, hence, it makes no sense to
    ! have the arrays in complex)
    if ( ts_Gamma_SCF ) then
       call newdSpData1D(ts_sp_uc,fdist,spH,name='TS spH')
       call newdSpData1D(ts_sp_uc,fdist,spS,name='TS spS')

    else
       call newzSpData1D(ts_sp_uc,fdist,spzH,name='TS spH')
       call newzSpData1D(ts_sp_uc,fdist,spzS,name='TS spS')

    end if

    if ( ts_Gamma_SCF ) then
       ! The (temporary) update arrays
       call newdSpData1D(tsup_sp_uc,fdist,spDMu,name='TS up DM')
       if ( .not. IsVolt .and. Calc_Forces ) &
            call newdSpData1D(tsup_sp_uc,fdist,spEDMu,name='TS up EDM')
       
    else
       call newzSpData1D(tsup_sp_uc,fdist,spzDMu,name='TS up DM')
       if ( .not. IsVolt .and. Calc_Forces ) &
            call newzSpData1D(tsup_sp_uc,fdist,spzEDMu,name='TS up EDM')
       
    end if

    if ( IsVolt ) then
       ! If we have a bias calculation we need additional arrays.
       ! If not bias we don't need the update arrays (we already have
       ! all information in tsup_sp_uc (spDMu))

       ! Allocate space for update arrays
       call newdSpData2D(ltsup_sp_sc,N_mu,sp_dist, &
            spDM,name='TS spDM')
       if ( Calc_Forces ) then
          call newdSpData2D(ltsup_sp_sc,N_mu,sp_dist, &
               spEDM,name='TS spEDM')
       end if

       ! The density matrix arrays of the non-equilibrium part
       call newdSpData1D(ltsup_sp_sc,N_nEq_id,sp_dist, &
            spDMneq,name='TS spDMneq')
       
    end if

    ! We just check that the work for reducing the matrices can be made
    if ( nzwork < nnzs(ts_sp_uc) ) call die('The memory for transiesta cannot &
         &sustain the implementation, contact the developers.')

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Writing GF^-1s (1000)'
    if(IONode)write(*,*)'Writing GFs (2000)'
    iu_GFinv = 1000 + Node
    iu_GF = 2000 + Node
    if(IONode)write(*,*)'Writing SigmaLs (3000)'
    if(IONode)write(*,*)'Writing SigmaRs (4000)'
    iu_SL = 3000 + Node
    iu_SR = 4000 + Node
#endif

    ! start the itterators
    call itt_init(SpKp,end1=nspin,end2=ts_nkpnt)
    ! point to the index iterators
    call itt_attach(SpKp,cur1=ispin,cur2=ikpt)

    do while ( .not. itt_next(SpKp) )

       ! This is going to get messy...
       ! we do not want to create, yet another sparsity pattern
       ! Hence we need to do the SAME checks, again and again and again....
       ! However, the extra computation should be negligible to the gain.

       if ( itt_stepped(SpKp,1) ) then
          < fix usage of calc_forces, dont reset EDM if .not. Calc_forces >
          call init_DM(sp_dist,sparse_pattern, &
               n_nzs, DM(:,ispin), EDM(:,ispin), &
               tsup_sp_uc)

          if ( IsVolt .and. TS_W_METHOD /= TS_W_K_UNCORRELATED ) then
             call init_val(spDM)
             call init_val(spEDM)
             call init_val(spDMneq)
          end if

       end if

       k(:) = ts_kpoint(:,ikpt)
       ! get the weight of the k-point
       kw = 1._dp / Pi * ts_kweight(ikpt)
       if ( nspin == 1 ) kw = kw * 2._dp

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',1)
#endif

       ! Work-arrays are for MPI distribution...
       if ( ts_Gamma_SCF ) then
          call create_HS_Gamma(sp_dist,sparse_pattern, &
               Ef, &
               no_BufL, no_BufR, & ! cut-out region
               N_Elec, Elecs, no_u, & ! electrodes, SIESTA size
               n_nzs, Hs(:,ispin), Ss, &
               spH, spS, &
               ndwork, dwork)
       else
          call create_HS_kpt(sp_dist,sparse_pattern, &
               Ef, &
               no_BufL, no_BufR, & ! cut-out region
               N_Elec, Elecs, no_u, & ! electrodes, SIESTA size
               n_nzs, Hs(:,ispin), Ss, &
               xij, &
               spzH, spzS, k, &
               nzwork, zwork)
       end if

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',2)
#endif

       ! loop across all energy-points
       iE = 1
       cE = nextE(iE+Node,steps=Nodes)
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          ! read in the next surface Green's function
          call read_next_GS(ikpt, cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = cE%idx(1) /= CONTOUR_EQ )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(no_BufL, no_u_TS, zwork, &
               N_Elec, Elecs, &
               spH =spH , spS =spS, &
               spzH=spzH, spzS=spzS )

          ! *******************
          ! * calc GF         *
          ! *******************
          select case ( cE%idx(1) )
          case ( CONTOUR_EQ )
             
             if ( all(Elecs(:)%DM_CrossTerms) ) then
                call calc_GF(cE,no_u_TS, zwork, GF, ierr)
                
             else
                call calc_GF_part(cE,no_BufL, no_u_TS, &
                     N_Elec, Elecs, &
                     zwork, GF, ierr)
             end if
             
          case ( CONTOUR_NEQ , CONTOUR_NEQ_TAIL )
             
             ! ensure that we do not do too many calculations
             ! capture
             call calc_GF_Bias(cE,no_BufL, no_u_TS, &
                  N_Elec, Elecs, &
                  zwork, GF, ierr)
             
          case default
             call die('Error in contour index')
             
          end select

          ! ** At this point we have calculated the Green's function

          ! ****************
          ! * save GF      *
          ! ****************
          if ( .not. IsVolt ) then
             if ( cE%idx(1) /= CONTOUR_EQ ) &
                  call die('Error in algorithm')

             ! TODO add the weight in case of the equilibrium thing
             < get weight of energy-point >
             if ( ts_Gamma_SCF ) then
                call add_DM_dE_D(spDMu, no_u_TS, no, &
                     GF, no_BufL, N_Elec, Elecs, cE, &
                     EDM = spEDMu)
             else
                call add_DM_dE_Z(spzDMu, no_u_TS, no, &
                     GF, no_BufL, N_Elec, Elecs, cE, &
                     EDM = spzEDMu)
             end if
             
          else if ( cE%idx(1) == CONTOUR_EQ ) then

             ! TODO move into add
             !if ( .not. cE%fake ) then
             ! We need to transfer the result over to the
             ! update array, and ready it for distribution
             call add_DM_dE_Z(spzDMu, no_u_TS, no, &
                  GF, no_BufL, N_Elec, Elecs, cE)

#ifdef MPI
             do i = 0 , Nodes
                ! retrieve the contour point from the node
                nE = nextE(iE+i,steps=Nodes)
                if ( nE%fake ) cycle

                call MPI_ScatterV(zDM,d_dist(:,3),d_dist(:,4), &
                     MPI_Double_Complex, Gf(1), d_dist(Node,3), &
                     MPI_Double_Complex, i, MPI_Comm_World, MPIerror)
                
                do imu = 1 , N_mu
                   ! we just need to take the first electrode
                   if ( Elec_hasE(Elecs(mus(imu)%el(1)),nE) ) then
                      call add_to_dm(idx=imu)
                   end if
                end do
             end do
#else
             do imu = 1 , N_mu
                ! we just need to take the first electrode
                if ( Elec_hasE(Elecs(mus(imu)%el(1)),cE) ) then
                   call add_to_dm(idx=imu)
                end if
             end do
#endif
          else
             ! non-equilibrium contour
             
             off = 0
             do iEl = 1 , N_Elec
                
                if ( .not. cE%fake ) then
                   if ( Elec_hasE(Elecs(iEl),cE) ) then

                   ! offset and number of orbitals
                   no = TotUsedOrbs(Elecs(iEl))
                   ! Calculate the triple-product
                   call GF_Gamma_GF(no_BufL, Elecs(iEl), no_u_TS, no, &
                        Gf(no_u_TS*off+1), zwork, size(GFGGF_work), GFGGF_work)
                   
                   ! we only calculate the Gf-column if the electrode
                   ! requires it
                   off = off + no
                end if

                ! We need to transfer the result over to the
                ! update array, and ready it for distribution
                call add_DM_dE_Z(spzDMu, no_u_TS, no_u_TS, &
                     zwork, no_BufL, N_Elec, Elecs, cE)

                end if

#ifdef MPI
                do i = 0 , Nodes

                   ! the nodes energy-point
                   nE = nextE(iE+i,steps=Nodes)
                   if ( nE%fake ) cycle

                   if ( Elec_hasE(Elecs(iEl),nE) ) then
                      < todo check that GFGGF_work is large enough >
                      call MPI_ScatterV(zDM,d_dist(:,3),d_dist(:,4), &
                           MPI_Double_Complex, GFGGF_work(1), d_dist(Node,3), &
                           MPI_Double_Complex, i, MPI_Comm_World, MPIerror)
                      
                      ! we have the contour, add to our system
                      < check which segment has this contour point and electrode >
                      
                   end if
                end do
#else
                < add to correct DM array >
#endif
             end do
                   

                ! check whether we should calculate
                ! Gf.Gamma.Gf for this electrode
                calc = .false.
                do j = 1 , N_nEq_segs 
                   if ( segment_has_El(nEq_segs(j),Elecs(i)) ) then
                      if ( cE%idx(1) == CONTOUR_NEQ .and. &
                           segment_has_c(nEq_segs(j),cE%idx(2)) ) then
                         calc = .true.
                         exit
                      else if ( cE%idx(1) == CONTOUR_NEQ_TAIL .and. &
                           segment_has_tail_c(nEq_segs(j),cE%idx(2)) ) then
                         calc = .true.
                         exit
                      end if
                   end if
                end do
                if ( .not. calc ) cycle

             end do

          end if

          ! step energy-point
          iE = iE + Nodes
          cE = next(iE+Node,steps=Nodes)
       end do

       ! reduce and shift to fermi-level
       call timer("TS_comm",1)
       if ( ts_Gamma_SCF ) then
          call d_DM_EDM_Reduce_Shift(Ef,spDMu, spEDMu, ndwork, dwork)
       else
          call z_DM_EDM_Reduce_Shift(Ef,spzDMu, spzEDMu, nzwork, zwork)
       end if
       call timer("TS_comm",2)

       if ( .not. IsVolt ) then
          if ( ts_Gamma_SCF ) then
             ! Directly save to the correct DM
             call update_DM(sp_dist,sparse_pattern, n_nzs, &
                  DM(:,ispin), EDM(:,ispin), spDMu, spEDMu, UpGlobal=.true.)
          else
             ! Directly save to the correct DM
             call update_zDM(sp_dist,sparse_pattern, n_nzs, &
                  DM(:,ispin), EDM(:,ispin), xij, spzDMu, spzEDMu, k)
          end if
       end if

       if ( ts_Gamma_SCF ) then
          ! Directly save to the correct DM
          call add_Gamma_DM(sp_dist,spDML, spEDML, spDMu, spEDMu)
       else
          ! Directly save to the correct DM
          call add_k_DM(sp_dist,spDML, spEDML, spzDMu, spzEDMu, &
               k, ltsup_sc_pnt, n_nzs, xij , non_Eq = .false. )
       end if

#ifdef TRANSIESTA_DEBUG
       call timer('TS_calc',2)
       if ( IONode ) then
          do iEl = 1 , N_Elec
          call io_close(uGF(i))
       end if
       return
#endif


       if ( ts_Gamma_SCF ) then
          ! Directly save to the correct DM
          ! Notice that we here save EDM to the correct EDM see weight_DM
          call add_Gamma_DM(sp_dist,spDMneqL, spEDMR, spDMu, spEDMu)
          call add_Gamma_DM(sp_dist,spDMneqR, spEDML, spDMuR, spEDMuR)
       else
          ! Here we have a couple of things to do
          if ( TS_W_METHOD == TS_W_CORRELATED ) then

             call add_k_DM(sp_dist,spDMneqL, spEDMR, spzDMu, spzEDMu, &
                  k, ltsup_sc_pnt, n_nzs, xij , non_Eq = .true.)
             call add_k_DM(sp_dist,spDMneqR, spEDML, spzDMuR, spzEDMuR, &
                  k, ltsup_sc_pnt, n_nzs, xij , non_Eq = .true.)

          else if ( TS_W_METHOD == TS_W_UNCORRELATED ) then

             call add_k_DM(sp_dist,spDMR, spEDMR, spzDMu, spzEDMu, &
                  k, ltsup_sc_pnt, n_nzs, xij , non_Eq = .true., spW=spDMneqL)
             call add_k_DM(sp_dist,spDML, spEDML, spzDMuR, spzEDMuR, &
                  k, ltsup_sc_pnt, n_nzs, xij , non_Eq = .true., spW=spDMneqR)

          else if ( TS_W_METHOD == TS_W_K_UNCORRELATED ) then

             call add_k_DM(sp_dist,spDMneqL, spEDMR, spzDMu, spzEDMu, &
                  k, ltsup_sc_pnt, n_nzs, xij , non_Eq = .true.)
             call add_k_DM(sp_dist,spDMneqR, spEDML, spzDMuR, spzEDMuR, &
                  k, ltsup_sc_pnt, n_nzs, xij , non_Eq = .true.)

          end if
       end if
       
       if ( IsVolt .and. TS_W_METHOD == TS_W_K_UNCORRELATED ) then
          call weight_DM( spDML, spDMR, spDMneqL, spDMneqR, &
               spEDML, spEDMR, nonEq_IsWeight = .false.)
          ! Directly save to the correct DM
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM(:,ispin), EDM(:,ispin), spDML, spEDML, ipnt=ltsup_sc_pnt)
       end if

       if ( itt_last(SpKp,2) ) then
          if ( IsVolt .and. TS_W_METHOD /= TS_W_K_UNCORRELATED ) then
             call weight_DM( spDML, spDMR, spDMneqL, spDMneqR, &
                  spEDML, spEDMR, nonEq_IsWeight = (TS_W_METHOD == TS_W_UNCORRELATED) )
             
             ! Directly save to the correct DM
             call update_DM(sp_dist, sparse_pattern, n_nzs, &
                  DM(:,ispin), EDM(:,ispin), spDML, spEDML, ipnt=ltsup_sc_pnt)
             
          end if
          
       end if

    ! We don't need to do anything here..
    end do
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
       do iEl = 1 , N_Elec
          call io_close(uGF(i))
       end if
    end if

!***********************
! CLEAN UP, there is alot to do... :)
!***********************
    if ( ts_Gamma_SCF ) then
       call delete(spDMu)
       call delete(spEDMu)
       if ( IsVolt ) then
          call delete(spDMuR)
          call delete(spEDMuR)
       end if

       ! The Hamiltonian and overlap matrices (in Gamma calculations
       ! we will not have any phases, hence, it makes no sense to
       ! have the arrays in complex)
       call delete(spH)
       call delete(spS)

    else
       call delete(spzDMu)
       call delete(spzEDMu)
       if ( IsVolt ) then
          call delete(spzDMuR)
          call delete(spzEDMuR)
       end if
       
       ! The Hamiltonian and overlap matrices
       call delete(spzH)
       call delete(spzS)
       
    end if
    
    if ( IsVolt ) then
       call delete(spDML)
       call delete(spEDML)
       call delete(spDMR)
       call delete(spDMneqL)
       call delete(spDMneqR)
       call delete(spEDMR)
    end if

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    if ( ts_Gamma_SCF ) then
       call memory('D','D',ndwork,'transiesta')
       deallocate(dwork)
    end if

    call memory('D','Z',size(zwork)+size(GF),'transiesta')
    deallocate(zwork,GF)

    do iEl = 1 , N_Elec
       call de_alloc(Elecs(iEl)%HA,routine='transiesta')
       call de_alloc(Elecs(iEl)%SA,routine='transiesta')
       call de_alloc(Elecs(iEl)%Gamma,routine='transiesta')
    end do

    ! In case of voltage calculations we need a work-array for
    ! handling the GF.Gamma.Gf^\dagger multiplication
    if ( IsVolt ) then
       call de_alloc(GFGGF_work, routine='transiesta')
    end if
   
    ! I would like to add something here which enables 
    ! 'Transiesta' post-process


    call timer('TS_calc',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS transiesta mem' )
#endif

  contains

    subroutine Equilibrium_Density(Z,W,ZW)
      complex(dp), intent(in) :: Z,W,ZW

      ! for these contour parts we do not require to calculate
      ! Gamma's.
      ! Hence we can perform the calculation without 
      ! calculating them.
      call UC_expansion(.false.,Z,no_L_HS,no_L, &
           Elecs(1), &
           na_L_HS,lasto_L,nqL, &
           HAAL, SAAL, GAAL, &
           SigmaL, GammaLT, & 
           nzwork, zwork)

      call UC_expansion(.false.,Z,no_R_HS,no_R, &
           Elecs(2), &
           na_R_HS,lasto_R,nqR, &
           HAAR, SAAR, GAAR, &
           SigmaR, GammaRT, & 
           nzwork, zwork)

#ifdef TRANSIESTA_DEBUG
      call write_Full(iu_SL,no_L,SigmaL)
      call write_Full(iu_SR,no_R,SigmaR)
#endif

      call prepare_GF_inv(UseBulk, Z, no_BufL, &           no_u_TS, zwork, &
           no_L, SigmaL, no_R, SigmaR, &
           spH =spH , spS =spS, &
           spzH=spzH, spzS=spzS )

         
      if ( UpdateDMCR ) then
         ! Only calculate the middle part of the Gf
         call calc_GF_Part(no_u_TS,no_L,no_R, zwork, GF,ierr)
         no_GF_offset = no_L
         ! The size of the central region (without left-right electrodes)
         no_u_C = no_u_TS - no_R - no_L
      else

#ifdef TRANSIESTA_DEBUG
         call write_Full(iu_GFinv,no_u_TS,zwork)
#endif

         ! Calculate the full GF
         call calc_GF(no_u_TS, zwork, GF,ierr)
         no_GF_offset = 0
         ! The size of the central region (with left-right electrodes)
         no_u_C = no_u_TS

#ifdef TRANSIESTA_DEBUG
         ! currently we will only write out the equilibrium GF
         call write_Full(iu_GF,no_u_TS,GF)
#endif

      end if

      if ( ts_Gamma_SCF ) then
         call add_DM_dE_D(spDMu ,  spEDMu, no_u_TS, no_u_C, &
              GF, no_BufL, no_GF_offset, W, ZW)
      else
         call add_DM_dE_Z(spzDMu, spzEDMu, no_u_TS, no_u_C, &
              GF, no_BufL, no_GF_offset, W, ZW)
      end if
      
    end subroutine Equilibrium_Density
    
    subroutine init_update_regions(BiasContour)
      logical, intent(in) :: BiasContour
      
      call init_val(spDMu)
      call init_val(spEDMu)
      call init_val(spzDMu)
      call init_val(spzEDMu)
      if ( BiasContour ) then
         call init_val(spDMuR)
         call init_val(spEDMuR)
         call init_val(spzDMuR)
         call init_val(spzEDMuR)
      end if

    end subroutine init_update_regions
    
    subroutine non_Equilibrium_Density(Z,i_W,i_ZW)
      complex(dp), intent(in) :: Z,i_W,i_ZW
      complex(dp) :: W,ZW

      ! The non-equilibrium integration points have the density
      ! in the real part of the Gf.Gamma.Gf^\dagger
      ! Hence we simply multiply W by -i to move the density
      ! to the same scheme i.e. \rho = - Im(Gf.Gamma.Gf^\dagger)
      W  = zmi * i_W
      ZW = Z * W

      call UC_expansion(.true.,Z,no_L_HS,no_L, &
           Elecs(1), &
           na_L_HS,lasto_L,nqL, &
           HAAL, SAAL, GAAL, &
           SigmaL, GammaLT, & 
           nzwork, zwork)

      call UC_expansion(.true.,Z,no_R_HS,no_R, &
           Elecs(2), &
           na_R_HS,lasto_R,nqR, &
           HAAR, SAAR, GAAR, &
           SigmaR, GammaRT, & 
           nzwork, zwork)

      call prepare_GF_inv(UseBulk, Z, no_BufL, &
           no_u_TS,zwork, &
           no_L, SigmaL, no_R, SigmaR, &
           spH =spH , spS =spS, &
           spzH=spzH, spzS=spzS )


      ! Calculate the Greens function
      call calc_GF_Bias(no_u_TS,no_L,no_R, zwork, GF,ierr)

#ifdef TRANSIESTA_DEBUG
      ind = 0 
      do ia_E = 1 , no_L
         do ia = 1 , no_u_TS
            ind = ind + 1
            call out_write(10000+iu_GF,ia,ia_E,GF(ind))
         end do
      end do
      do ia_E = no_u_TS - no_R + 1 , no_u_TS
         do ia = 1 , no_u_TS
            ind = ind + 1
            call out_write(10000+iu_GF,ia,ia_E,GF(ind))
         end do
      end do
#endif

      ! We calculate the right contribution
      call GF_Gamma_GF(UpdateDMCR,no_L+1,no_u_TS,no_L+no_R, &
           no_R, Gf, &
           GammaRT, zwork, no_R*no_R, GFGGF_work)
      ! zwork is now GF.G.GF

#ifdef TRANSIESTA_DEBUG
      call write_full(iu_GF,no_u_TS,zwork)
#endif

      ! Note that we use '--' here
      if ( ts_Gamma_SCF ) then
         call add_DM_dE_D(spDMuR,spEDMuR, no_u_TS, no_u_TS, &
              zwork, no_BufL, 0, -W, -ZW)
      else
         call add_DM_dE_Z(spzDMuR,spzEDMuR, no_u_TS, no_u_TS, &
              zwork, no_BufL, 0, -W, -ZW)
      end if
         
      ! We calculate the left contribution
      call GF_Gamma_GF(UpdateDMCR,1,no_u_TS, no_L+no_R, &
           no_L, Gf, &
           GammaLT, zwork, no_L*no_L, GFGGF_work)
      ! zwork is now GF.G.GF

#ifdef TRANSIESTA_DEBUG
      call write_full(iu_GF,no_u_TS,zwork)
#endif

      ! Note that we use '++' here
      if ( ts_Gamma_SCF ) then
         call add_DM_dE_D(spDMu,spEDMu, no_u_TS, no_u_TS, &
              zwork, no_BufL, 0, +W, +ZW)
      else
         call add_DM_dE_Z(spzDMu,spzEDMu, no_u_TS, no_u_TS, &
              zwork, no_BufL, 0, +W, +ZW)
      end if
      
    end subroutine non_Equilibrium_Density
    
  end subroutine transiesta_full


  ! Update DM
  ! These routines are supplied for easy update of the update region
  ! sparsity patterns
  ! Note that these routines implement the usual rho(Z) \propto - GF
  subroutine add_DM_dE_Z(DM,no1,no2,GF,no_BufL,N_Elec,Elecs, cE, &
       EDM)
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
    integer, intent(in) :: no_BufL
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
  subroutine prepare_invGF(cE, no_BufL,no_u,GFinv, &
       N_Elec, Elecs, spH, spS, spzH, spzS)
    use class_dSpData1D
    use class_zSpData1D
    use class_Sparsity
    ! the current energy point
    type(ts_c), intent(in) :: cE
    ! Remark that we need the left buffer orbitals
    ! to calculate the actual orbital of the sparse matrices...
    integer, intent(in) :: no_BufL, no_u
    complex(dp), intent(out) :: GFinv(no_u**2)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData1D), intent(inout), optional :: spH,  spS
    type(zSpData1D), intent(inout), optional :: spzH, spzS

    ! Local variables
    complex(dp) :: Z
    type(Sparsity), pointer :: s
    logical :: Is_Gamma
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: dH(:), dS(:)
    complex(dp), pointer :: zH(:), zS(:)
    integer :: io, iu,ind, ioff

    if ( cE%fake ) return

    Z = cE%e
    ! Determine whether we have a Gamma or k-point Hamiltonian
    if ( initialized(spH) .eqv. initialized(spzH) ) then
       call die('Transiesta error, not &
            &two initialized arrays are allowed, check with the &
            &developers.')
    end if
    
    if ( initialized(spH) ) then
       Is_Gamma = .true.
       s  => spar(spH)
       dH => val (spH)
       dS => val (spS)
    else
       Is_Gamma = .false.
       s  => spar(spzH)
       zH => val (spzH)
       zS => val (spzS)
    end if

    l_ncol => n_col   (s)
    l_ptr  => list_ptr(s)
    l_col  => list_col(s)
     
    ! Offset
    ioff = no_BufL + 1
    
    ! Initialize
    GFinv(1:no_u**2) = dcmplx(0._dp,0._dp)

    if ( Is_Gamma ) then
       ! We will only loop in the central region
       ! We have constructed the sparse array to only contain
       ! values in this part...
       do io = ioff, no_BufL + no_u
       
          iu = (io - ioff) * no_u - no_BufL
          
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 
             
             !ju = l_col(ind) ! the '- no_BufL' is moved outside the loop
          
             ! Notice that we transpose S and H back here
             ! See symmetrize_HS_Gamma (H is hermitian)
             GFinv(l_col(ind)+iu) = Z * dS(ind) - dH(ind)
          end do
       end do
    else
    
       ! We will only loop in the central region
       ! We have constructed the sparse array to only contain
       ! values in this part...
       do io = ioff, no_BufL + no_u
          
          iu = (io - ioff) * no_u - no_BufL
          
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 
             
             !ju = l_col(ind) ! The '- no_BufL' is moved outside the loop
             
             ! Notice that we transpose S and H back here
             ! See symmetrize_HS_kpt
             GFinv(l_col(ind)+iu) = Z * zS(ind) - zH(ind)
          end do
       end do

    end if

    do io = 1 , N_Elec
       call insert_Self_Energies(no_BufL, no_u, Gfinv, Elecs(io))
    end do

  end subroutine prepare_invGF
   
  subroutine insert_Self_Energies(no_BufL, no_u, Gfinv, El)
    integer, intent(in) :: no_BufL, no_u
    complex(dp), intent(in out) :: GFinv(no_u,no_u)
    type(Elec), intent(in) :: El

    integer :: i, j, ii, jj, off, no

    no = TotUsedOrbs(El)
    off = El%idx_no - no_BufL - 1

    if ( El%Bulk ) then
       do j = 1 , no
          do i = 1 , no
             Gfinv(off+i,off+j) = El%Sigma(i,j)
          end do
       end do
    else
       do j = 1 , no
          jj = off + j
          do i = 1 , no
             ii = off + i
             Gfinv(ii,jj) = Gfinv(ii,jj) - El%Sigma(i,j) 
          end do
       end do
    end if

  end subroutine insert_Self_Energies

end module m_ts_full
