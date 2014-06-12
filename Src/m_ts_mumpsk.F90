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

module m_ts_mumpsk

  use precision, only : dp

#ifdef MUMPS  
  use m_ts_sparse_helper

  use m_ts_dm_update, only : init_DM
  use m_ts_dm_update, only : update_DM, update_zDM
  use m_ts_dm_update, only : add_k_DM
  
  use m_ts_weight, only : weight_DM
  use m_ts_weight, only : TS_W_K_METHOD
  use m_ts_weight, only : TS_W_K_CORRELATED
  use m_ts_weight, only : TS_W_K_HALF_CORRELATED
  use m_ts_weight, only : TS_W_K_UNCORRELATED

  use m_ts_method, only: orb_offset, no_Buf, ts2s_orb

  use m_ts_mumps_init
  
  implicit none
  
  public :: ts_mumpsk
  
  private
  
contains

  subroutine ts_mumpsk(N_Elec,Elecs, &
       nq,uGF, &
       ucell, nspin, na_u, lasto, &
       sp_dist, sparse_pattern, &
       no_u, n_nzs, &
       Hs, Ss, xij, DM, EDM, Ef, kT)

    use units, only : Pi, eV
    use parallel, only : Node, Nodes, IONode
#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use class_zSpData1D
    use class_zSpData2D

    use m_ts_electype
    ! Self-energy retrival and expansion
    use m_ts_elec_se

    use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight

    use m_ts_options, only : Calc_Forces
    use m_ts_options, only : N_mu, mus

    use m_ts_options, only : IsVolt

    use m_ts_sparse, only : ts_sp_uc
    use m_ts_sparse, only : tsup_sp_uc
    use m_ts_sparse, only : ltsup_sp_sc
    use m_ts_sparse, only : ltsup_sc_pnt

    use m_ts_cctype
    use m_ts_contour,     only : has_cE
    use m_ts_contour_eq,  only : Eq_E, ID2idx, c2weight_eq
    use m_ts_contour_neq, only : nEq_E
    use m_ts_contour_neq, only : N_nEq_ID, c2weight_neq
    
    use m_iterator

    ! Gf calculation
    use m_ts_full_scat

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    include 'zmumps_struc.h'

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    integer, intent(in) :: nq(N_Elec), uGF(N_Elec)
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: nspin, na_u, lasto(0:na_u)
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    integer, intent(in)  :: no_u
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: Hs(n_nzs,nspin), Ss(n_nzs), xij(3,n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Ef, kT

! ****************** Electrode variables *********************
    complex(dp), pointer :: GFGGF_work(:) => null()
! ************************************************************

! ******************* Computational arrays *******************
    integer :: nzwork
    ! The solution arrays
    type(zMUMPS_STRUC) :: mum
    complex(dp), pointer :: zwork(:), GF(:)

    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D) :: spH, spS
    ! local sparsity pattern in local SC pattern
    type(dSpData2D) :: spDM, spDMneq
    type(dSpData2D) :: spEDM ! only used if calc_forces
    ! The different sparse matrices that will surmount to the integral
    ! These two lines are in global update sparsity pattern (UC)
    type(zSpData2D) ::  spuDM
    type(zSpData2D) :: spuEDM ! only used if calc_forces
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    real(dp)    :: kw, kpt(3), bkpt(3)
    complex(dp) :: W, ZW
! ************************************************************

! ******************** Loop variables ************************
    type(itt2) :: SpKp
    integer, pointer :: ispin, ikpt
    integer :: iEl, iID, up_nzs
    integer :: iE, imu, io, idx
! ************************************************************

! ******************* Miscalleneous variables ****************
    integer :: ierr, no_u_TS, off, no
! ************************************************************

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    ! Number of orbitals in TranSIESTA
    no_u_TS = no_u - no_Buf

    ! Number of elements that are transiesta updated
    up_nzs = nnzs(tsup_sp_uc)

    nullify(Gf)

    ! We do need the full GF AND a single work array to handle the
    ! left-hand side of the inversion...
    ! We will provide all work arrays as single dimension arrays.
    ! This will make interfaces more stringent and allow for
    ! re-use in several other places.
    ! However, this comes at the cost of programmer book-keeping.
    ! It should be expected that the work arrays return GARBAGE
    ! on ALL routines, i.e. they are not used for anything other
    ! than, yes, work.

#ifdef TRANSIESTA_TIMING
    call timer('TS_MUMPS_INIT',1)
#endif

    ! initialize MUMPS
    call init_MUMPS(mum,Node)

    ! Initialize for the correct size
    mum%N    = no_u_TS ! order of matrix
    mum%NRHS = no_u_TS ! RHS-vectors

    ! prepare the LHS of MUMPS-solver.
    ! This is the same no matter the contour type
    call prep_LHS(mum,N_Elec,Elecs)
    nzwork = mum%NZ
    zwork => mum%A(:)

    ! analyzation step
    call analyze_MUMPS(mum)

#ifdef TRANSIESTA_TIMING
    call timer('TS_MUMPS_INIT',2)
#endif

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
    call newzSpData1D(ts_sp_uc,fdist,spH,name='TS spH')
    call newzSpData1D(ts_sp_uc,fdist,spS,name='TS spS')

    ! If we have a bias calculation we need additional arrays.
    ! If not bias we don't need the update arrays (we already have
    ! all information in tsup_sp_uc (spDMu))

    ! Allocate space for global sparsity arrays
    no = max(N_mu,N_nEq_id)
    call newzSpData2D(tsup_sp_uc,no,fdist, spuDM, name='TS spuDM')
    if ( Calc_Forces ) then
       call newzSpData2D(tsup_sp_uc,N_mu,fdist, spuEDM, name='TS spuEDM')
    end if
    
    if ( IsVolt ) then
       ! Allocate space for update arrays, local sparsity arrays
       call newdSpData2D(ltsup_sp_sc,N_mu,    sp_dist,spDM   ,name='TS spDM')
       call newdSpData2D(ltsup_sp_sc,N_nEq_id,sp_dist,spDMneq,name='TS spDM-neq')
       if ( Calc_Forces ) then
          call newdSpData2D(ltsup_sp_sc,N_mu, sp_dist,spEDM  ,name='TS spEDM')
       end if
    end if

    ! start the itterators
    call itt_init  (SpKp,end1=nspin,end2=ts_nkpnt)
    ! point to the index iterators
    call itt_attach(SpKp,cur1=ispin,cur2=ikpt)

    do while ( .not. itt_step(SpKp) )

       if ( itt_stepped(SpKp,1) ) then ! spin has incremented
          call init_DM(sp_dist, sparse_pattern, &
               n_nzs, DM(:,ispin), EDM(:,ispin), &
               tsup_sp_uc, Calc_Forces)
       end if

       ! Include spin factor and 1/\pi
       kpt(:) = ts_kpoint(:,ikpt)
       ! create the k-point in reciprocal space
       call kpoint_convert(Ucell,kpt,bkpt,1)
       kw = 1._dp / Pi * ts_kweight(ikpt)
       if ( nspin == 1 ) kw = kw * 2._dp

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',1)
#endif

       ! Work-arrays are for MPI distribution...
       call create_HS(sp_dist,sparse_pattern, &
            Ef, &
            N_Elec, Elecs, no_u, & ! electrodes, SIESTA size
            n_nzs, Hs(:,ispin), Ss, xij, &
            spH, spS, kpt, &
            nzwork, zwork)

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',2)
#endif

#ifdef TRANSIESTA_TIMING
       call timer('TS_EQ',1)
#endif

       call prep_RHS_Eq(mum,no_u_TS,up_nzs,N_Elec,Elecs,Gf)

       ! ***************
       ! * EQUILIBRIUM *
       ! ***************
       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
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
          call read_next_GS(ispin, ikpt, kpt, &
               cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, .false., forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .false. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
#ifdef TRANSIESTA_TIMING
          call timer('TS_PREP',1)
#endif
          call prepare_invGF(cE, mum, &
               N_Elec, Elecs, &
               spH=spH , spS=spS)
#ifdef TRANSIESTA_TIMING
          call timer('TS_PREP',2)
#endif

          ! *******************
          ! * calc GF         *
          ! *******************
#ifdef TRANSIESTA_TIMING
          call timer('TS_MUMPS_SOLVE',1)
#endif
          write(mum%ICNTL(1),'(a,i0,2(a,i0))') &
               '### Solving Eq Node/iC: ',Node,'/',cE%idx(2),',',cE%idx(3)
          mum%JOB = 5
          call zMUMPS(mum)
          if ( mum%INFO(1) < 0 .or. mum%INFOG(1) < 0 ) then
             call die('MUMPS failed the Eq. inversion, check the output logs')
          end if
#ifdef TRANSIESTA_TIMING
          call timer('TS_MUMPS_SOLVE',2)
#endif
          
          ! ** At this point we have calculated the Green's function

          ! ****************
          ! * save GF      *
          ! ****************
          do imu = 1 , N_mu
             if ( cE%fake ) cycle ! in case we implement an MPI communication solution...
             call ID2idx(cE,mus(imu)%ID,idx)
             if ( idx < 1 ) cycle
             
             call c2weight_eq(cE,idx, kw, W ,ZW)
             call add_DM( spuDM, W, spuEDM, ZW, &
                  mum, &
                  N_Elec, Elecs, &
                  DMidx=mus(imu)%ID)
          end do

          ! step energy-point
          iE = iE + Nodes
          cE = Eq_E(iE+Nodes-Node,step=Nodes) ! we read them backwards
       end do

#ifdef TRANSIESTA_TIMING
       call timer('TS_EQ',2)
#endif

#ifdef MPI
       ! We need to reduce all the arrays
       call timer('TS_comm',1)
       call AllReduce_SpData(spuDM,nzwork,zwork,N_mu)
       if ( Calc_Forces ) then
          call AllReduce_SpData(spuEDM,nzwork,zwork,N_mu)
       end if
       call timer('TS_comm',2)
#endif

       if ( .not. IsVolt ) then
          call update_zDM(sp_dist,sparse_pattern, n_nzs, &
               DM(:,ispin) ,  spuDM, Ef, &
               EDM(:,ispin), spuEDM, kpt, xij)

          ! The remaining code segment only deals with 
          ! bias integration... So we skip instantly

          cycle

       end if

       ! *****************
       ! * only things with non-Equilibrium contour...
       ! *****************

       ! initialize to zero
       ! local sparsity update patterns
       ! if (tsweightmethod...)
       if ( TS_W_K_METHOD == TS_W_K_UNCORRELATED ) then
          call init_val(spDM)
          call init_val(spDMneq)
          if ( Calc_Forces ) call init_val(spEDM)
       else if ( itt_stepped(SpKp,1) ) then
          ! we only need to initialize once per spin
          call init_val(spDM)
          call init_val(spDMneq)
          if ( Calc_Forces ) call init_val(spEDM)
       end if

       ! transfer equilibrium data to local sparsity arrays
       call add_k_DM(spDM, spuDM, N_mu, &
            spEDM, spuEDM, N_mu, &
            n_nzs, xij, kpt, ipnt=ltsup_sc_pnt, non_Eq = .false. )

#ifdef TRANSIESTA_TIMING
       call timer('TS_NEQ',1)
#endif

       call prep_RHS_nEq(mum,no_u_TS,N_Elec,Elecs,Gf)

       ! *******************
       ! * NON-EQUILIBRIUM *
       ! *******************
       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
       iE = 0
       cE = nEq_E(iE+Nodes-Node,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(ispin, ikpt, kpt, &
               cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, .false., forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .true. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(cE, mum, &
               N_Elec, Elecs, &
               spH =spH , spS =spS)
          
          ! *******************
          ! * calc GF         *
          ! *******************
          write(mum%ICNTL(1),'(a,i0,2(a,i0))') &
               '### Solving nEq Node/iC: ',Node,'/',cE%idx(2),',',cE%idx(3)
          mum%JOB = 5
          call zMUMPS(mum)
          if ( mum%INFO(1) < 0 .or. mum%INFOG(1) < 0 ) then
             call die('MUMPS failed the nEq. inversion, check the output logs')
          end if

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

             ! Currently this is the *ONLY* thing
             ! that we need to correct to get MUMPS nEq to work
             !call GF_Gamma_GF(Elecs(iEl), no_u_TS, no, &
             !     Gf(no_u_TS*off+1), zwork, size(GFGGF_work), GFGGF_work)

             ! step to the next electrode position
             off = off + no
                
             do iID = 1 , N_nEq_ID
                
                if ( .not. has_cE(cE,iEl=iEl,ineq=iID) ) cycle
                
                call c2weight_neq(cE,kT,iEl,iID, kw,W,imu,ZW)

                call add_DM( spuDM, W, spuEDM, ZW, &
                     mum, &
                     N_Elec, Elecs, &
                     DMidx=iID, EDMidx=imu)
             end do
          end do

          ! step energy-point
          iE = iE + Nodes
          cE = nEq_E(iE+Nodes-Node,step=Nodes) ! we read them backwards
       end do

#ifdef TRANSIESTA_TIMING
       call timer('TS_NEQ',2)
#endif


#ifdef MPI
       ! We need to reduce all the arrays
       call timer('TS_comm',1)
       call AllReduce_SpData(spuDM, nzwork, zwork, N_nEq_id)
       if ( Calc_Forces ) then
          call AllReduce_SpData(spuEDM, nzwork, zwork, N_mu)
       end if
       call timer('TS_comm',2)
#endif

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',1)
#endif

       ! 1. move from global UC to local SC
       ! 2. calculate the correct contribution by applying the weight
       ! 3. add the density to the real arrays
       call add_k_DM(spDMneq, spuDM, N_nEq_id, &
            spEDM, spuEDM, N_mu, &
            n_nzs, xij, kpt, ipnt=ltsup_sc_pnt, non_Eq = .true. )

       if ( TS_W_K_METHOD == TS_W_K_UNCORRELATED ) then
          call weight_DM( N_Elec, Elecs, N_mu, na_u, lasto, &
               spDM, spDMneq, spEDM, &
               nonEq_IsWeight = .false.)
          
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM(:,ispin), spDM, Ef=Ef, &
               EDM=EDM(:,ispin), spEDM=spEDM, ipnt=ltsup_sc_pnt)
       else if ( TS_W_K_METHOD == TS_W_K_HALF_CORRELATED ) then
          call die('not functioning yet')
       else if ( itt_last(SpKp,2) ) then ! TS_W_K_METHOD == TS_W_K_CORRELATED
          call weight_DM( N_Elec, Elecs, N_mu, na_u, lasto, &
               spDM, spDMneq, spEDM, &
               nonEq_IsWeight = .false.)
          
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM(:,ispin), spDM, Ef=Ef, &
               EDM=EDM(:,ispin), spEDM=spEDM, ipnt=ltsup_sc_pnt)          
       end if

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',2)
#endif
       
       ! We don't need to do anything here..

    end do ! spin

    call itt_destroy(SpKp)

#ifdef TRANSIESTA_DEBUG
    write(*,*) 'Completed TRANSIESTA SCF'
#endif

!***********************
! CLEAN UP
!***********************

    call delete(spH)
    call delete(spS)

    call delete(spuDM)
    call delete(spuEDM)

    call delete(spDM)
    call delete(spDMneq)
    call delete(spEDM)

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    ! Deallocate user data
    deallocate( mum%IRN )
    deallocate( mum%JCN )
    deallocate( mum%A )
    deallocate( mum%IRHS_PTR )
    deallocate( mum%IRHS_SPARSE )
    deallocate( GF ) ! mum%RHS_SPARSE is pointing to GF
    nullify( mum%RHS_SPARSE )
    no = mum%ICNTL(1)

    ! Destroy the instance (deallocate internal data structures)
    mum%JOB = -2
    call zMUMPS(mum)
    if ( mum%INFO(1) < 0 .or. mum%INFOG(1) < 0 ) then
       call die('Cleaning of MUMPS failed.')
    end if

    ! close the file
    call io_close(no)

    ! In case of voltage calculations we need a work-array for
    ! handling the GF.Gamma.Gf^\dagger multiplication
    if ( IsVolt ) then
       call de_alloc(GFGGF_work, routine='transiesta')
    end if
   
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS transiesta mem' )
#endif

  end subroutine ts_mumpsk


  ! Update DM
  ! These routines are supplied for easy update of the update region
  ! sparsity patterns
  ! Note that these routines implement the usual rho(Z) \propto - GF
  subroutine add_DM(DM, DMfact,EDM, EDMfact, &
       mum, &
       N_Elec,Elecs, &
       DMidx, EDMidx)

    use class_Sparsity
    use class_zSpData2D
    use m_ts_electype
    use m_ts_method, only : ts2s_orb
    use intrinsic_missing, only : SFIND
    include 'zmumps_struc.h'

    ! The DM and EDM equivalent matrices
    type(zSpData2D), intent(inout) :: DM
    complex(dp), intent(in) :: DMfact
    type(zSpData2D), intent(inout) :: EDM
    complex(dp), intent(in) :: EDMfact
    ! the mumps structure
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! the index of the partition
    integer, intent(in) :: DMidx
    integer, intent(in), optional :: EDMidx

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: D(:,:), E(:,:), GF(:)
    integer :: io, jo, ind, ir, nr, Hn, ind_H
    integer :: iu, ju, i1, i2
    logical :: hasEDM

    ! Remember that this sparsity pattern HAS to be in Global UC
    s => spar(DM)
    call attach(s,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=nr)
    D => val(DM)
    hasEDM = initialized(EDM)
    if ( hasEDM ) E => val(EDM)

    i1 = DMidx
    i2 = i1
    if ( present(EDMidx) ) i2 = EDMidx

    GF => mum%RHS_SPARSE

    if ( hasEDM ) then
       
       do ir = 1 , mum%NRHS
             
          ! this is column index
          jo = ts2s_orb(ir)

          ! The update region equivalent GF part
          do ind = mum%IRHS_PTR(ir) , mum%IRHS_PTR(ir+1)-1
                
             io = ts2s_orb(mum%IRHS_SPARSE(ind))

             Hn    = l_ncol(io)
             ind_H = l_ptr(io)
             ! Requires that l_col is sorted
             ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)
             if ( ind_H == l_ptr(io) ) cycle
                
             D(ind_H,i1) = D(ind_H,i1) - GF(ind) * DMfact
             E(ind_H,i2) = E(ind_H,i2) - GF(ind) * EDMfact
             
          end do
       end do
       
    else

       do ir = 1 , mum%NRHS
             
          ! this is column index
          jo = ts2s_orb(ir)

          ! The update region equivalent GF part
          do ind = mum%IRHS_PTR(ir) , mum%IRHS_PTR(ir+1)-1
                
             io = ts2s_orb(mum%IRHS_SPARSE(ind))

             Hn    = l_ncol(io)
             ind_H = l_ptr(io)
             ! Requires that l_col is sorted
             ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)
             if ( ind_H == l_ptr(io) ) cycle
                
             D(ind_H,i1) = D(ind_H,i1) - GF(ind) * DMfact
             
          end do
       end do

    end if

  end subroutine add_DM


  ! creation of the GF^{-1}.
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_invGF(cE, mum, N_Elec, Elecs, spH, spS)

    use intrinsic_missing, only : SFIND
    use class_zSpData1D
    use class_Sparsity
    use m_ts_electype
    use m_ts_cctype, only : ts_c_idx
    use m_ts_method, only : ts2s_orb
    include 'zmumps_struc.h'

    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D), intent(inout) :: spH, spS

    ! Local variables
    complex(dp) :: Z
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: H(:), S(:), iG(:)
    integer :: io, jo, ind, Hn, ind_H

    if ( cE%fake ) return

    Z = cE%e
    
    sp => spar(spH)
    H  => val (spH)
    S  => val (spS)

    l_ncol => n_col   (sp)
    l_ptr  => list_ptr(sp)
    l_col  => list_col(sp)

    ! Initialize
    iG => mum%A
    iG = 0._dp ! possibly this is not needed...

    do ind = 1, mum%NZ

       io = ts2s_orb(mum%JCN(ind))
       jo = ts2s_orb(mum%IRN(ind))

       Hn = l_ncol(io)
       if ( Hn == 0 ) cycle

       ind_H = l_ptr(io)
       ! Requires that l_col is sorted
       ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)

       if ( ind_H == l_ptr(io) ) cycle
       
       ! Notice that we transpose S and H back here
       ! See symmetrize_HS_Gamma (H is hermitian)
       iG(ind) = Z * S(ind_H) - H(ind_H)

    end do

    do io = 1 , N_Elec
       call insert_Self_Energies(mum, Elecs(io))
    end do

  end subroutine prepare_invGF
   
#endif
end module m_ts_mumpsk
