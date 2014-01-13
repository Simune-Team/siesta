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

module m_ts_fullg

  use precision, only : dp

  use m_ts_sparse_helper

  use m_ts_dm_update, only : init_DM
  use m_ts_dm_update, only : update_DM
  use m_ts_dm_update, only : add_Gamma_DM
  
  use m_ts_weight, only : weight_DM
  
  implicit none
  
  public :: ts_fullg
  
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
  
  subroutine ts_fullg(N_Elec,Elecs, &
       nq,uGF, &
       nspin, &
       sp_dist, sparse_pattern, &
       ucell, no_u, na_u, lasto, xa, n_nzs, &
       Hs, Ss, DM, EDM, Ef, kT)

    use units, only : Pi, eV
    use parallel, only : Node, Nodes, IONode
#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use m_ts_electype
    ! Self-energy retrival and expansion
    use m_ts_elec_se

    use m_ts_options, only : Calc_Forces
    use m_ts_options, only : N_mu, mus
    use m_ts_options, only : no_BufL, no_BufR

    use m_ts_options, only : IsVolt

    use m_ts_sparse, only : ts_sp_uc
    use m_ts_sparse, only : tsup_sp_uc
    use m_ts_sparse, only : ltsup_sp_sc
    use m_ts_sparse, only : ltsup_sc_pnt

    use m_ts_cctype
    use m_ts_contour,     only : has_cE
    use m_ts_contour_eq,  only : Eq_E, ID2idx, c2weight_eq
    use m_ts_contour_neq, only : nEq_E, indices2eq
    use m_ts_contour_neq, only : N_nEq_ID, c2weight_neq
    
    use m_iterator

    ! Gf calculation
    use m_ts_full_scat

use m_object_debug

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    integer, intent(in) :: nq(N_Elec), uGF(N_Elec)
    integer, intent(in) :: nspin
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in)  :: no_u, na_u
    integer, intent(in)  :: lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: Hs(n_nzs,nspin), Ss(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Ef, kT

! ****************** Electrode variables *********************
    complex(dp), pointer :: GFGGF_work(:) => null()
! ************************************************************

! ******************* Computational arrays *******************
    integer :: ndwork, nzwork
    real(dp), pointer :: dwork(:,:)
    complex(dp), allocatable, target :: zwork(:), GF(:)

    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData1D) :: spH, spS
    ! local sparsity pattern in local SC pattern
    type(dSpData2D) :: spDM, spDMneq
    type(dSpData2D) :: spEDM ! only used if calc_forces
    ! The different sparse matrices that will surmount to the integral
    ! These two lines are in global update sparsity pattern (UC)
    type(dSpData2D) ::  spuDM
    type(dSpData2D) :: spuEDM ! only used if calc_forces
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    real(dp) :: kw
    complex(dp) :: Z, W, ZW
    complex(dp), parameter :: zmi = dcmplx(0._dp,-1._dp)
! ************************************************************

! ******************** Loop variables ************************
    type(itt1) :: Sp
    integer, pointer :: ispin
    integer :: iEl, iID, up_nzs, ia, ia_E
    integer :: ind, iE, imu, io, idx
#ifdef TRANSIESTA_DEBUG
    integer :: iu_GF, iu_GFinv
    integer :: iu_SL, iu_SR
#endif
! ************************************************************

! ******************* Miscalleneous variables ****************
    integer :: ierr, no_u_TS, off, no
! ************************************************************

#ifdef MPI
! ******************* MPI-related variables   ****************
    integer :: MPIerror
! ************************************************************
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    ! Number of orbitals in TranSIESTA
    no_u_TS = no_u - no_BufL - no_BufR

    ! Number of elements that are transiesta updated
    up_nzs = nnzs(tsup_sp_uc)

    ! We do need the full GF AND a single work array to handle the
    ! left-hand side of the inversion...
    ! We will provide all work arrays as single dimension arrays.
    ! This will make interfaces more stringent and allow for
    ! re-use in several other places.
    ! However, this comes at the cost of programmer book-keeping.
    ! It should be expected that the work arrays return GARBAGE
    ! on ALL routines, i.e. they are not used for anything other
    ! than, yes, work.

    ! The zwork is needed to construct the LHS for solving: G^{-1} G = I
    ! Hence, we will minimum require the full matrix...
    nzwork = no_u_TS ** 2
    allocate(zwork(nzwork),stat=ierr)
    if (ierr/=0) call die('Could not allocate space for zwork')
    call memory('A','Z',nzwork,'transiesta')

    ! We only need a partial size of the Green's function
    no = no_u_TS
    do iEl = 1 , N_Elec
       if ( .not. Elecs(iEl)%DM_CrossTerms ) then
          no = no - TotUsedOrbs(Elecs(iEl))
       end if
    end do
    ! when bias is needed we need the entire GF column
    ! for all the electrodes (at least some of the contour points needs this)
    if ( IsVolt ) then
       no = max(no,sum(TotUsedOrbs(Elecs)))
    end if
    no = no * no_u_TS
    allocate(GF(no),stat=ierr)
    if (ierr/=0) call die('Could not allocate space for GFpart')
    call memory('A','Z',no,'transiesta')

    no = 0
    do iEl = 1 , N_Elec

       ! This seems stupid, however, we never use the Sigma and
       ! GF at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called
       ! first the Sigma is assigned and then 
       ! it is required that prepare_GF_inv is called
       ! immediately (which it is)
       ! Hence the GF must NOT be used in between these two calls!
       io = TotUsedOrbs(Elecs(iEl))
       Elecs(iEl)%Sigma => GF(no+1:no+io**2)
       no = no + io ** 2

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
    call newdSpData1D(ts_sp_uc,fdist,spH,name='TS spH')
    call newdSpData1D(ts_sp_uc,fdist,spS,name='TS spS')

    ! If we have a bias calculation we need additional arrays.
    ! If not bias we don't need the update arrays (we already have
    ! all information in tsup_sp_uc (spDMu))

    ! Allocate space for global sparsity arrays
    no = max(N_mu,N_nEq_id)
    call newdSpData2D(tsup_sp_uc,no,fdist, spuDM, name='TS spuDM')
    ! assign dwork, this will problably come at the expence of
    ! two full reductions, however, we save some memory...
    ndwork = nnzs(tsup_sp_uc)
    dwork => val(spuDM)
    if ( Calc_Forces ) then
       call newdSpData2D(tsup_sp_uc,N_mu,fdist, spuEDM, name='TS spuEDM')
    end if
    
    if ( IsVolt ) then
       ! Allocate space for update arrays, local sparsity arrays
       call newdSpData2D(ltsup_sp_sc,N_mu,    sp_dist,spDM   ,name='TS spDM')
       call newdSpData2D(ltsup_sp_sc,N_nEq_id,sp_dist,spDMneq,name='TS spDM-neq')
       if ( nnzs(ltsup_sp_sc) > ndwork ) then
          ! only update if this array is larger (should only happen in 
          ! few processor setups
          ndwork = nnzs(ltsup_sp_sc)
          dwork => val(spDMneq)
       end if
       if ( Calc_Forces ) then
          call newdSpData2D(ltsup_sp_sc,N_mu, sp_dist,spEDM  ,name='TS spEDM')
       end if
    end if

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
    call itt_init  (Sp,end=nspin)
    ! point to the index iterators
    call itt_attach(Sp,cur=ispin)

    do while ( .not. itt_step(Sp) )

       ! This is going to get messy...
       ! we do not want to create, yet another sparsity pattern
       ! Hence we need to do the SAME checks, again and again and again....
       ! However, the extra computation should be negligible to the gain.

       call init_DM(sp_dist,sparse_pattern, &
            n_nzs, DM(:,ispin), EDM(:,ispin), &
            tsup_sp_uc, Calc_Forces)

       ! Include spin factor and 1/\pi
       kw = 1._dp / Pi
       if ( nspin == 1 ) kw = kw * 2._dp

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',1)
#endif

       ! Work-arrays are for MPI distribution...
       call create_HS(sp_dist,sparse_pattern, &
            Ef, &
            no_BufL, no_BufR, & ! cut-out region
            N_Elec, Elecs, no_u, & ! electrodes, SIESTA size
            n_nzs, Hs(:,ispin), Ss, &
            spH, spS, &
            ndwork, dwork(:,1)) ! annoyingly we can't pass the full array!!!!!

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',2)
#endif

       ! initialize to zero
       ! local sparsity update patterns
       if ( IsVolt ) then
          call init_val(spDM)
          call init_val(spDMneq)
          if ( Calc_Forces ) call init_val(spEDM)
       end if

       ! initialize to zero
       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)


       ! ***************
       ! * EQUILIBRIUM *
       ! ***************
       no = no_u_TS
       do iEl = 1 , N_Elec
          if ( .not. Elecs(iEl)%DM_CrossTerms ) then
             no = no - TotUsedOrbs(Elecs(iEl))
          end if
       end do
       iE = 0
       cE = Eq_E(iE+Nodes-Node,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(1, cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .false. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(cE, no_BufL, no_u_TS, zwork, &
               N_Elec, Elecs, &
               spH=spH , spS=spS)

          ! *******************
          ! * calc GF         *
          ! *******************
          if ( all(Elecs(:)%DM_CrossTerms) ) then
             call calc_GF(cE,no_u_TS, zwork, GF, ierr)
             
          else
             call calc_GF_part(cE,no_BufL, no_u_TS, &
                  N_Elec, Elecs, &
                  zwork, GF, ierr)
          end if
          
          ! ** At this point we have calculated the Green's function

          ! ****************
          ! * save GF      *
          ! ****************
          do imu = 1 , N_mu
             if ( cE%fake ) cycle ! in case we implement an MPI communication solution...
             call ID2idx(cE,mus(imu)%ID,idx)
             if ( idx < 1 ) cycle
             
             call c2weight_eq(cE,idx, kw, W,ZW)
!print '(f7.3,tr1,i0,tr1,4(f14.9,tr1))',mus(imu)%mu/eV,imu,cE%E/eV,W/kw
             if ( Calc_Forces ) then
                call add_DM(spuDM, mus(imu)%ID, &
                     no_u_TS, no, &
                     GF, no_BufL, N_Elec, Elecs, W, &
                     EDM=spuEDM, EDMfact=ZW)
             else
                call add_DM(spuDM, mus(imu)%ID, &
                     no_u_TS, no, &
                     GF, no_BufL, N_Elec, Elecs, W)
             end if

          end do

          ! step energy-point
          iE = iE + Nodes
          cE = Eq_E(iE+Nodes-Node,step=Nodes) ! we read them backwards
       end do

#ifdef MPI
       ! We need to reduce all the arrays
       call timer('TS_comm',1)
       call my_full_G_reduce(spuDM,nzwork*2,zwork,N_mu)
       if ( Calc_Forces ) then
          call my_full_G_reduce(spuEDM,nzwork*2,zwork,N_mu)
       end if
       call timer('TS_comm',2)
#endif

       if ( .not. IsVolt ) then
          if ( Calc_Forces ) then
             call update_DM(sp_dist,sparse_pattern, n_nzs, &
                  DM(:,ispin), spuDM, Ef=Ef, &
                  EDM=EDM(:,ispin), spEDM=spuEDM, &
                  UpSpGlobal = .true.)
          else
             call update_DM(sp_dist,sparse_pattern, n_nzs, &
                  DM(:,ispin), spuDM, UpSpGlobal = .true.)
          end if

          ! The remaining code segment only deals with 
          ! bias integration...
          cycle
       end if

       ! *****************
       ! * only things with non-Equilibrium contour...
       ! *****************

       ! transfer data to local sparsity arrays
       if ( Calc_Forces ) then
          call add_Gamma_DM(spDM,   spuDM, D_dim2=N_mu, &
               spEDM=spEDM, spuEDM=spuEDM, E_dim2=N_mu)
       else
          call add_Gamma_DM(spDM,   spuDM, D_dim2=N_mu)
       end if

       if ( IsVolt ) then
          ! initialize to zero (only if we have a bias)
          call init_val(spuDM)
          if ( Calc_Forces ) call init_val(spuEDM)
       end if

       ! *******************
       ! * NON-EQUILIBRIUM *
       ! *******************
       iE = 0
       cE = nEq_E(iE+Nodes-Node,step=Nodes) ! we read them backwards
       do while ( cE%exist )

!print*,'nEq: read next GS...',Node,cE%fake
          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(1, cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .true. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(cE, no_BufL, no_u_TS, zwork, &
               N_Elec, Elecs, &
               spH =spH , spS =spS)
          
          ! *******************
          ! * calc GF         *
          ! *******************
          call calc_GF_Bias(cE, no_BufL, no_u_TS, &
               N_Elec, Elecs, &
               zwork, GF, ierr)

          ! ** At this point we have calculated the Green's function

          ! ****************
          ! * save GF      *
          ! ****************
          off = 0
          do iEl = 1 , N_Elec
             if ( cE%fake ) cycle ! in case we implement an MPI communication solution

             if ( .not. has_cE(cE,iEl=iEl) ) cycle

             ! offset and number of orbitals
             no = TotUsedOrbs(Elecs(iEl))

             call GF_Gamma_GF(no_BufL, Elecs(iEl), no_u_TS, no, &
                  Gf(no_u_TS*off+1), zwork, size(GFGGF_work), GFGGF_work)

             ! step to the next electrode position
             off = off + no
                
             do iID = 1 , N_nEq_ID
                
                if ( .not. has_cE(cE,iEl=iEl,ineq=iID) ) cycle
                
                call c2weight_neq(cE,kT,iEl,iID, kw, W,ZW)

                     call indices2eq(iID,imu,ZW)
!print '(3(tr1,i0),tr2,4(f14.9,tr1))',iEl,imu,iID,cE%E/eV,W/kw

                ! the energy-density matrix contribution gets added appropriately...
                call add_Bias_DM(spuDM, iID, &
                     no_u_TS, no_u_TS, &
                     zwork, no_BufL, N_Elec, Elecs, W, &
                     has_offset = .false.)

                if ( Calc_Forces ) then
                   call indices2eq(iID,imu,ZW)
                   call add_Bias_DM(spuEDM, imu, &
                        no_u_TS, no_u_TS, &
                        zwork, no_BufL, N_Elec, Elecs, ZW, &
                        has_offset = .false.)
                end if
                
             end do
          end do

          ! step energy-point
          iE = iE + Nodes
          cE = nEq_E(iE+Nodes-Node,step=Nodes) ! we read them backwards
       end do

#ifdef TRANSIESTA_DEBUG
          call timer('TS_calc',2)
          if ( IONode ) then
             do iEl = 1 , N_Elec
             call io_close(uGF(iEl))
          end do
       end if
       return
#endif

       if ( IsVolt ) then

#ifdef MPI
          ! We need to reduce all the arrays
          call timer('TS_comm',1)
          call my_full_G_reduce(spuDM, nzwork*2, zwork, N_nEq_id)
          if ( Calc_Forces ) then
             call my_full_G_reduce(spuEDM, nzwork*2, zwork, N_mu)
          end if
          call timer('TS_comm',2)
#endif

          ! 1. move from global UC to local SC
          ! 2. calculate the correct contribution by applying the weight
          ! 3. add the density to the real arrays
          if ( Calc_Forces ) then
             call add_Gamma_DM(spDMneq, spuDM, D_dim2=N_nEq_id, &
                  spEDM=spEDM,  spuEDM=spuEDM, E_dim2=N_mu)
             
             call weight_DM( N_Elec, N_mu, spDM, spDMneq, &
                  spEDM = spEDM, &
                  nonEq_IsWeight = .false.)

             call update_DM(sp_dist,sparse_pattern, n_nzs, &
                  DM(:,ispin), spDM, Ef=Ef, &
                  EDM=EDM(:,ispin), spEDM=spEDM, ipnt=ltsup_sc_pnt)
          else
             call add_Gamma_DM(spDMneq, spuDM, D_dim2=N_nEq_id)

             call weight_DM( N_Elec, N_mu, spDM, spDMneq, &
                  nonEq_IsWeight = .false.)

             call update_DM(sp_dist,sparse_pattern, n_nzs, &
                  DM(:,ispin), spDM, ipnt=ltsup_sc_pnt)
          end if

       end if

       ! We don't need to do anything here..
    end do ! spin

    call itt_destroy(Sp)

#ifdef TRANSIESTA_DEBUG
    write(*,*) 'Completed TRANSIESTA SCF'
#endif

!***********************
! CLEAN UP
!***********************

    call delete(spH)
    call delete(spS)

    call delete(spuDM)
    if ( Calc_Forces ) call delete(spuEDM)

    if ( IsVolt ) then
       call delete(spDM)
       call delete(spDMneq)
       if ( Calc_Forces ) call delete(spEDM)
    end if

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    call memory('D','Z',size(zwork)+size(GF),'transiesta')
    deallocate(zwork,GF)

    ! In case of voltage calculations we need a work-array for
    ! handling the GF.Gamma.Gf^\dagger multiplication
    if ( IsVolt ) then
       call de_alloc(GFGGF_work, routine='transiesta')
    end if
   
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS transiesta mem' )
#endif

  end subroutine ts_fullg

  ! Update DM
  ! These routines are supplied for easy update of the update region
  ! sparsity patterns
  ! Note that these routines implement the usual rho(Z) \propto - GF
  subroutine add_DM(DM,idx,no1,no2,GF, &
       no_BufL,N_Elec,Elecs,DMfact,EDM,EDMfact, &
       has_offset)
    use class_Sparsity
    use class_dSpData2D
    use m_ts_electype
    ! The DM and EDM equivalent matrices
    type(dSpData2D), intent(inout) :: DM
    ! the index of the partition
    integer, intent(in) :: idx
    ! The size of GF
    integer, intent(in) :: no1, no2
    ! The Green's function
    complex(dp), intent(in) :: GF(no1,no2)
    ! The number of buffer atoms (needed for the offset in the sparsity
    ! patterns), and the offset in the GF
    integer, intent(in) :: no_BufL, N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! Complex numbers that are used in the factor of GF
    complex(dp), intent(in) :: DMfact
    type(dSpData2D), intent(inout), optional :: EDM
    complex(dp), intent(in), optional :: EDMfact
    logical, intent(in), optional :: has_offset

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: D(:,:), E(:,:)
    integer :: io, ind, nr
    integer :: iu, ju
    logical :: lh_off

    s      => spar(DM)
    l_ncol => n_col   (s)
    l_ptr  => list_ptr(s)
    l_col  => list_col(s)
    D      => val(DM)
    if ( present(EDM) ) E => val(EDM)
    if ( present(EDM) .neqv. present(EDMfact) ) &
         call die('add_DM: Error in code') ! TODO DELETE

    lh_off = .true.
    if ( present(has_offset) ) lh_off = has_offset

    ! Number of orbitals in the SIESTA unit-cell
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    nr = nrows(s)

    if ( present(EDM) ) then
       if ( lh_off ) then

          do io = 1 , nr ! TODO introduce reduced loop
             ! Quickly go past the buffer atoms...
             if ( l_ncol(io) == 0 ) cycle

             ! The update region equivalent GF part
             iu = io - no_BufL
        
             do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                
                ju = l_col(ind) - no_BufL - offset(N_Elec,Elecs,l_col(ind))
                
                D(ind,idx) = D(ind,idx) - dimag( GF(iu,ju) * DMfact  )
                E(ind,idx) = E(ind,idx) - dimag( GF(iu,ju) * EDMfact )
                
             end do
          end do
     
       else
          do io = 1 , nr
             if ( l_ncol(io) == 0 ) cycle
             iu = io - no_BufL
             do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                ju = l_col(ind) - no_BufL
                D(ind,idx) = D(ind,idx) - dimag( GF(iu,ju) * DMfact  )
                E(ind,idx) = E(ind,idx) - dimag( GF(iu,ju) * EDMfact )
             end do
          end do

       end if
    else

       if ( lh_off ) then
          do io = 1 , nr ! TODO introduce reduced loop
             ! Quickly go past the buffer atoms...
             if ( l_ncol(io) == 0 ) cycle

             ! The update region equivalent GF part
             iu = io - no_BufL
             
             do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                
                ju = l_col(ind) - no_BufL - offset(N_Elec,Elecs,l_col(ind))
                
                D(ind,idx) = D(ind,idx) - dimag( GF(iu,ju) * DMfact )
                
             end do
          end do

       else
          do io = 1 , nr
             if ( l_ncol(io) == 0 ) cycle
             iu = io - no_BufL
             do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                ju = l_col(ind) - no_BufL
                D(ind,idx) = D(ind,idx) - dimag( GF(iu,ju) * DMfact )
             end do
          end do

       end if
    end if

  contains
    
    pure function offset(N_Elec,Elecs,io)
      integer, intent(in) :: N_Elec
      type(Elec), intent(in) :: Elecs(N_Elec)
      integer, intent(in) :: io
      integer :: offset
      ! TODO figure out if offset really works!!!
      offset = sum(TotUsedOrbs(Elecs(:)), &
           MASK=.not. Elecs(:)%DM_CrossTerms .and. Elecs(:)%idx_no <= io )
    end function offset

  end subroutine add_DM

  ! Update DM
  ! These routines are supplied for easy update of the update region
  ! sparsity patterns
  ! Note that these routines implement the usual rho(Z) \propto - GF
  subroutine add_Bias_DM(DM,idx,no1,no2,GF, &
       no_BufL,N_Elec,Elecs,DMfact,EDM,EDMfact, &
       has_offset)
    use class_Sparsity
    use class_dSpData2D
    use m_ts_electype
    ! The DM and EDM equivalent matrices
    type(dSpData2D), intent(inout) :: DM
    ! the index of the partition
    integer, intent(in) :: idx
    ! The size of GF
    integer, intent(in) :: no1, no2
    ! The Green's function
    complex(dp), intent(in) :: GF(no1,no2)
    ! The number of buffer atoms (needed for the offset in the sparsity
    ! patterns), and the offset in the GF
    integer, intent(in) :: no_BufL, N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! Complex numbers that are used in the factor of GF
    complex(dp), intent(in) :: DMfact
    type(dSpData2D), intent(inout), optional :: EDM
    complex(dp), intent(in), optional :: EDMfact
    logical, intent(in), optional :: has_offset

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: D(:,:), E(:,:)
    integer :: io, ind, nr
    integer :: iu, ju
    logical :: lh_off

    s      => spar(DM)
    l_ncol => n_col   (s)
    l_ptr  => list_ptr(s)
    l_col  => list_col(s)
    D      => val(DM)
    if ( present(EDM) ) E => val(EDM)
    if ( present(EDM) .neqv. present(EDMfact) ) &
         call die('add_DM: Error in code') ! TODO DELETE

    lh_off = .true.
    if ( present(has_offset) ) lh_off = has_offset

    ! Number of orbitals in the SIESTA unit-cell
    ! Remember that this is a sparsity pattern which contains
    ! a subset of the SIESTA pattern.
    nr = nrows(s)

    if ( present(EDM) ) then
       if ( lh_off ) then

          do io = 1 , nr ! TODO introduce reduced loop
             ! Quickly go past the buffer atoms...
             if ( l_ncol(io) == 0 ) cycle

             ! The update region equivalent GF part
             iu = io - no_BufL
        
             do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                
                ju = l_col(ind) - no_BufL - offset(N_Elec,Elecs,l_col(ind))
                
                D(ind,idx) = D(ind,idx) - GF(iu,ju) * DMfact
                E(ind,idx) = E(ind,idx) - GF(iu,ju) * EDMfact
                
             end do
          end do
     
       else
          do io = 1 , nr
             if ( l_ncol(io) == 0 ) cycle
             iu = io - no_BufL
             do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                ju = l_col(ind) - no_BufL
                D(ind,idx) = D(ind,idx) - GF(iu,ju) * DMfact
                E(ind,idx) = E(ind,idx) - GF(iu,ju) * EDMfact
             end do
          end do

       end if
    else

       if ( lh_off ) then
          do io = 1 , nr ! TODO introduce reduced loop
             ! Quickly go past the buffer atoms...
             if ( l_ncol(io) == 0 ) cycle

             ! The update region equivalent GF part
             iu = io - no_BufL
             
             do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                
                ju = l_col(ind) - no_BufL - offset(N_Elec,Elecs,l_col(ind))
                
                D(ind,idx) = D(ind,idx) - GF(iu,ju) * DMfact
                
             end do
          end do

       else
          do io = 1 , nr
             if ( l_ncol(io) == 0 ) cycle
             iu = io - no_BufL
             do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
                ju = l_col(ind) - no_BufL
                D(ind,idx) = D(ind,idx) - GF(iu,ju) * DMfact
             end do
          end do

       end if
    end if

  contains
    
    pure function offset(N_Elec,Elecs,io)
      integer, intent(in) :: N_Elec
      type(Elec), intent(in) :: Elecs(N_Elec)
      integer, intent(in) :: io
      integer :: offset
      ! TODO figure out if offset really works!!!
      offset = sum(TotUsedOrbs(Elecs(:)), &
           MASK=.not. Elecs(:)%DM_CrossTerms .and. Elecs(:)%idx_no <= io )
    end function offset

  end subroutine add_Bias_DM


  ! creation of the GF^{-1}.
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_invGF(cE, no_BufL,no_u,GFinv, &
       N_Elec, Elecs, spH, spS)
    use class_dSpData1D
    use class_Sparsity
    use m_ts_electype
    use m_ts_cctype, only : ts_c_idx
    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    ! Remark that we need the left buffer orbitals
    ! to calculate the actual orbital of the sparse matrices...
    integer, intent(in) :: no_BufL, no_u
    complex(dp), intent(out) :: GFinv(no_u**2)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData1D), intent(inout), optional :: spH,  spS

    ! Local variables
    complex(dp) :: Z
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: H(:), S(:)
    integer :: io, iu, ind, ioff

    if ( cE%fake ) return

    Z = cE%e
    
    sp => spar(spH)
    H  => val (spH)
    S  => val (spS)

    l_ncol => n_col   (sp)
    l_ptr  => list_ptr(sp)
    l_col  => list_col(sp)
     
    ! Offset
    ioff = no_BufL + 1
    
    ! Initialize
    GFinv(1:no_u**2) = dcmplx(0._dp,0._dp)

    ! We will only loop in the central region
    ! We have constructed the sparse array to only contain
    ! values in this part...
    do io = ioff, no_BufL + no_u
       
       iu = (io - ioff) * no_u - no_BufL
       
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 
             
          ! Notice that we transpose S and H back here
          ! See symmetrize_HS_Gamma (H is hermitian)
          GFinv(iu+l_col(ind)) = Z * S(ind) - H(ind)

       end do
    end do

    do io = 1 , N_Elec
       call insert_Self_Energies(no_BufL, no_u, Gfinv, Elecs(io))
    end do

  end subroutine prepare_invGF
   
  subroutine insert_Self_Energies(no_BufL, no_u, Gfinv, El)
    use m_ts_electype
    integer, intent(in) :: no_BufL, no_u
    complex(dp), intent(in out) :: GFinv(no_u,no_u)
    type(Elec), intent(in) :: El

    integer :: i, j, ii, jj, iii, off, no

    no = TotUsedOrbs(El)
    off = El%idx_no - no_BufL - 1

    if ( El%Bulk ) then
       ii = 0
       do j = 1 , no
          jj = off + j
          do i = 1 , no
             ii = ii + 1
             Gfinv(off+i,jj) = El%Sigma(ii)
          end do
       end do
    else
       iii = 0
       do j = 1 , no
          jj = off + j
          do i = 1 , no
             ii = off + i
             iii = iii + 1
             Gfinv(ii,jj) = Gfinv(ii,jj) - El%Sigma(iii) 
          end do
       end do
    end if

  end subroutine insert_Self_Energies

end module m_ts_fullg


#ifdef MPI
subroutine my_full_G_reduce(sp_arr,nwork,work,dim2_count)
  use precision, only : dp
  use class_dSpData2D
  use m_ts_sparse_helper, only : AllReduce_SpData
  type(dSpData2D), intent(inout) :: sp_arr
  integer, intent(in)     :: nwork
  real(dp), intent(inout) :: work(nwork)
  integer, intent(in) :: dim2_count
  call AllReduce_SpData(sp_arr,nwork,work,dim2_count)
end subroutine my_full_G_reduce
#endif

