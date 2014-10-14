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

module m_ts_trik

  use precision, only : dp

  use m_ts_sparse_helper

  use m_ts_dm_update, only : init_DM
  use m_ts_dm_update, only : update_DM, update_zDM
  use m_ts_dm_update, only : add_k_DM
  
  use m_ts_weight, only : weight_DM
  use m_ts_weight, only : TS_W_K_METHOD
  use m_ts_weight, only : TS_W_K_CORRELATED
  use m_ts_weight, only : TS_W_K_UNCORRELATED

  use m_ts_tri_init, only : N_tri_part, tri_parts

  use m_ts_method, only : orb_offset, no_Buf
  
  implicit none

  public :: ts_trik, ts_trik_Fermi

  private
  
contains

  subroutine ts_trik(N_Elec,Elecs, &
       nq, uGF, ucell, nspin, na_u, lasto, &
       sp_dist, sparse_pattern, &
       no_u, n_nzs, &
       Hs, Ss, DM, EDM, Ef, kT)

    use units, only : Pi
    use parallel, only : Node, Nodes, IONode

#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D
    use class_dSpData2D
    use class_zSpData2D
    use class_zTriMat

    use m_ts_electype
    ! Self-energy read
    use m_ts_gf
    ! Self-energy expansion
    use m_ts_elec_se

    use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight

    use m_ts_options, only : Calc_Forces
    use m_ts_options, only : N_mu, mus

    use m_ts_options, only : IsVolt

    use m_ts_sparse, only : ts_sp_uc, tsup_sp_uc
    use m_ts_sparse, only : ltsup_sp_sc, ltsup_sc_pnt
    use m_ts_sparse, only : sc_off

    use m_ts_cctype
    use m_ts_contour,     only : has_cE
    use m_ts_contour_eq,  only : Eq_E, ID2idx, c2weight_eq
    use m_ts_contour_neq, only : nEq_E
    use m_ts_contour_neq, only : N_nEq_ID, c2weight_neq

    use m_iterator
    use m_mat_invert

    use m_trimat_invert

    ! Gf calculation
    use m_ts_trimat_invert

    ! Gf.Gamma.Gf
    use m_ts_tri_scat

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

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
    real(dp), intent(in) :: Hs(n_nzs,nspin), Ss(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Ef, kT

! ******************* Computational arrays *******************
    integer :: nzwork, n_s
    complex(dp), pointer :: zwork(:)
    type(zTriMat) :: zwork_tri, GF_tri
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
    ! To figure out which parts of the tri-diagonal blocks we need
    ! to calculate
    logical, pointer :: calc_parts(:) => null()
! ************************************************************

! ****************** Electrode variables *********************
    integer :: padding, GFGGF_size ! with IsVolt we need padding and work-array
    complex(dp), pointer :: GFGGF_work(:) => null()
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    logical     :: has_El(N_Elec)
    real(dp)    :: kw, kpt(3), bkpt(3)
    complex(dp) :: W, ZW
! ************************************************************

! ******************** Loop variables ************************
    type(itt2) :: SpKp
    integer, pointer :: ispin, ikpt
    integer :: iEl, iID, ia
    integer :: iE, imu, io, idx
    integer :: no, no_u_TS
! ************************************************************

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    ! Number of supercells
    n_s = size(sc_off,dim=2)

    ! Number of orbitals in TranSIESTA
    no_u_TS = no_u - no_Buf

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
    ! Hence, we will minimum require this...
    if ( IsVolt ) then
       io = maxval(TotUsedOrbs(Elecs(:)))
       call GFGGF_needed_worksize(N_tri_part,tri_parts, &
            io, padding, GFGGF_size)
    else
       padding = 0
       GFGGF_size = 0
    end if
    call newzTriMat(zwork_tri,N_tri_part,tri_parts,'GFinv', &
         padding=padding)
    nzwork = elements(zwork_tri)

    ! Initialize the tri-diagonal inversion routine
    call init_TriMat_inversion(zwork_tri)
    if ( IsVolt ) then
       call init_BiasTriMat_inversion(zwork_tri)
    end if

    call newzTriMat(GF_tri,N_tri_part,tri_parts,'GF')

    ! initialize the matrix inversion tool
    call init_mat_inversion(maxval(tri_parts))

    ! Allocate the logical array to handle calculated
    ! entries in the block-matrix
    call re_alloc(calc_parts,1,N_tri_part)
    ! initialize to calculate all blocks
    calc_parts(:) = .true.

    ! we use the GF as a placement for the self-energies
    no = 0
    zwork => val(GF_tri)
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
       Elecs(iEl)%Sigma => zwork(no+1:no+io**2)
       no = no + io ** 2

       ! if we need the cross-terms we can not skip the blocks
       ! that are fully inside the electrode
       if ( Elecs(iEl)%DM_update /= 0 ) cycle

       io  = Elecs(iEl)%idx_o
       io  = io - orb_offset(io)
       idx = io + TotUsedOrbs(Elecs(iEl)) - 1

       do ia = 1 , N_tri_part
          if ( has_full_part(N_tri_part,tri_parts,ia,io,idx) ) then
             ! The entire part 'ia' is encapsulated in 
             ! this electrode.
             calc_parts(ia) = .false.
          end if
       end do

    end do

    ! Save the work-space
    ! Now the programmer should "keep a straight tongue"
    ! The zwork points to the array in the zwork_tri
    ! tri-diagonal array. This means that there are two
    ! arrays that point to the same.
    ! Generally the zwork need only to retain the value in
    ! one call!
    zwork => val(zwork_tri)

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
    ! Notice that we DO need it to be the SIESTA size.
#ifdef MPI
    call newDistribution(no_u,MPI_Comm_Self,fdist,name='TS-fake dist')
#else
    call newDistribution(no_u,-1           ,fdist,name='TS-fake dist')
#endif
    
    ! The Hamiltonian and overlap matrices (in Gamma calculations
    ! we will not have any phases, hence, it makes no sense to
    ! have the arrays in complex)
    ! TODO move into a Data2D (could reduce overhead of COMM)
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

    if ( IsVolt ) then
       ! we need only allocate one work-array for
       ! Gf.G.Gf^\dagger
       call re_alloc(GFGGF_work,1,GFGGF_size,routine='transiesta')
    end if

    ! start the itterators
    call itt_init  (SpKp,end1=nspin,end2=ts_nkpnt)
    ! point to the index iterators
    call itt_attach(SpKp,cur1=ispin,cur2=ikpt)
    
    do while ( .not. itt_step(SpKp) )
       
       if ( itt_stepped(SpKp,1) ) then
          call init_DM(sp_dist,sparse_pattern, &
               n_nzs, DM(:,ispin), EDM(:,ispin), &
               tsup_sp_uc, Calc_Forces)
       end if

       ! Include spin factor and 1/(2\pi)
       kpt(:) = ts_kpoint(:,ikpt)
       ! create the k-point in reciprocal space
       call kpoint_convert(ucell,kpt,bkpt,1)
       kw = 0.5_dp / Pi * ts_kweight(ikpt)
       if ( nspin == 1 ) kw = kw * 2._dp
       
#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',1)
#endif

       ! Work-arrays are for MPI distribution...
       call create_HS(sp_dist,sparse_pattern, &
            Ef, &
            N_Elec, Elecs, no_u, n_s, & ! electrodes, SIESTA size
            n_nzs, Hs(:,ispin), Ss, sc_off, &
            spH, spS, kpt, &
            nzwork, zwork)

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',2)
#endif


#ifdef TRANSIESTA_TIMING
       call timer('TS_EQ',1)
#endif

       ! ***************
       ! * EQUILIBRIUM *
       ! ***************
       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
       iE = Nodes - Node
       cE = Eq_E(iE,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(ispin, ikpt, bkpt, &
               cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, .false., forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .false. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(cE, no_u_TS, zwork_tri, &
               N_Elec, Elecs, &
               spH=spH , spS=spS)

          ! *******************
          ! * calc GF         *
          ! *******************
          if ( .not. cE%fake ) then
             call invert_TriMat(zwork_tri,GF_tri,calc_parts)
          end if
          
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
                  GF_tri, &
                  N_Elec, Elecs, &
                  DMidx=mus(imu)%ID)
          end do

          ! step energy-point
          iE = iE + Nodes
          cE = Eq_E(iE,step=Nodes) ! we read them backwards
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
               DM(:,ispin), spuDM, Ef, &
               EDM(:,ispin), spuEDM, kpt, n_s, sc_off)

          ! The remaining code segment only deals with 
          ! bias integration... So we skip instantly

          cycle

       end if

       ! *****************
       ! * only things with non-Equilibrium contour...
       ! *****************

       ! initialize to zero
       ! local sparsity update patterns
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

       ! transfer data to local sparsity arrays
       call add_k_DM(spDM, spuDM, N_mu, &
            spEDM, spuEDM, N_mu, &
            n_s, sc_off, kpt, non_Eq = .false. )

#ifdef TRANSIESTA_TIMING
       call timer('TS_NEQ',1)
#endif

       ! *******************
       ! * NON-EQUILIBRIUM *
       ! *******************
       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
       iE = Nodes - Node
       cE = nEq_E(iE,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(ispin, ikpt, bkpt, &
               cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, .false., forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .true. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(cE, no_u_TS, zwork_tri, &
               N_Elec, Elecs, &
               spH =spH , spS =spS)
          
          ! *******************
          ! * prep GF         *
          ! *******************
          if ( .not. cE%fake ) then
             do iEl = 1 , N_Elec
                has_El(iEl) = has_cE(cE,iEl=iEl)
             end do
             call invert_BiasTriMat_prep(zwork_tri,GF_tri, &
                  N_Elec, Elecs, has_El)
          end if

          ! ** At this point we have calculated the needed
          ! ** information to create the Green's function column
          ! ** for all the electrodes

#ifdef TS_DEV
          if ( .not. cE%fake ) then
          io = 1000 + iE + Node
          open(io,form='unformatted')
          write(io) cE%e
          end if
#endif
          
          ! ****************
          ! * save GF      *
          ! ****************
          do iEl = 1 , N_Elec
             if ( cE%fake ) cycle ! in case we implement an MPI communication solution

             if ( .not. has_El(iEl) ) cycle

             ! ******************
             ! * calc GF-column *
             ! ******************
             call invert_BiasTriMat_col(GF_tri,zwork_tri, &
                  Elecs(iEl), calc_parts)

#ifdef TS_DEV
             ! offset and number of orbitals
             no = TotUsedOrbs(Elecs(iEl))

             idx = 0
             do iid = 1 , N_tri_part
                write(io) tri_parts(iid),no
                write(io) zwork(idx+1:idx+tri_parts(iid)*no)
                idx = idx + tri_parts(iid)*no
             end do
#endif
             
             call GF_Gamma_GF(zwork_tri, Elecs(iEl), calc_parts, &
                  GFGGF_size, GFGGF_work)
#ifdef TRANSIESTA_WEIGHT_DEBUG
             print '(a7,tr1,i3,2(tr1,f10.5),tr5,2(tr1,f10.5))', &
                  trim(Elecs(iEl)%name),iE,zwork(index(zwork_tri,28,28)),cE%e
#endif 

             do iID = 1 , N_nEq_ID
                
                if ( .not. has_cE(cE,iEl=iEl,ineq=iID) ) cycle
                
                call c2weight_neq(cE,kT,iEl,iID, kw,W,imu,ZW)
#ifdef TRANSIESTA_WEIGHT_DEBUG
                print '(a20,2(tr1,i3),2(tr1,e12.5))', &
                     trim(Elecs(iEl)%name),iID,imu,W
#endif 

                call add_DM( spuDM, W, spuEDM, ZW, &
                     zwork_tri, &
                     N_Elec, Elecs, &
                     DMidx=iID, EDMidx=imu )
             end do
          end do

#ifdef TS_DEV
          if ( .not. cE%fake ) &
               close(io)
#endif

          ! step energy-point
          iE = iE + Nodes
          cE = nEq_E(iE,step=Nodes) ! we read them backwards
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
            n_s, sc_off, kpt, non_Eq = .true. )

       if ( TS_W_K_METHOD == TS_W_K_UNCORRELATED ) then
          call weight_DM( N_Elec, Elecs, N_mu, na_u, lasto, &
               spDM, spDMneq, spEDM, n_s, sc_off)
          
#ifdef TRANSIESTA_WEIGHT_DEBUG
          call die('')
#endif
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM(:,ispin), spDM, Ef=Ef, &
               EDM=EDM(:,ispin), spEDM=spEDM, ipnt=ltsup_sc_pnt)
       else if ( itt_last(SpKp,2) ) then ! TS_W_K_METHOD == TS_W_K_CORRELATED
          call weight_DM( N_Elec, Elecs, N_mu, na_u, lasto, &
               spDM, spDMneq, spEDM, n_s, sc_off)
          
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM(:,ispin), spDM, Ef=Ef, &
               EDM=EDM(:,ispin), spEDM=spEDM, ipnt=ltsup_sc_pnt)          
       end if

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',2)
#endif

       ! We don't need to do anything here..

    end do ! spin, k-point

    call itt_destroy(SpKp)

#ifdef TRANSIESTA_DEBUG
    write(*,*) 'Completed TRANSIESTA SCF'
#endif

!***********************
! CLEAN UP
!***********************

    call de_alloc(calc_parts)

    call delete(zwork_tri)
    call delete(GF_tri)

    call delete(spH)
    call delete(spS)

    call delete(spuDM)
    call delete(spuEDM)

    call delete(spDM)
    call delete(spDMneq)
    call delete(spEDM)

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    call clear_TriMat_inversion()
    if ( IsVolt ) then
       call clear_BiasTriMat_inversion()
       call de_alloc(GFGGF_work, routine='transiesta')
    end if
    call clear_mat_inversion()

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS transiesta mem' )
#endif

  end subroutine ts_trik

  subroutine ts_trik_Fermi(N_Elec,Elecs, &
       nq, uGF, ucell, nspin, na_u, lasto, &
       sp_dist, sparse_pattern, &
       no_u, n_nzs, &
       Hs, Ss, DM, Ef, kT, Qtot, converged)

    use units, only : Pi
    use parallel, only : Node, Nodes, IONode

#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D
    use class_dSpData2D
    use class_zSpData2D
    use class_zTriMat

    use m_ts_electype
    ! Self-energy read
    use m_ts_gf
    ! Self-energy expansion
    use m_ts_elec_se

    use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight

    use m_ts_sparse, only : ts_sp_uc, tsup_sp_uc
    use m_ts_sparse, only : ltsup_sp_sc, sc_off

    use m_ts_charge

    use m_ts_cctype
    use m_ts_contour_eq, only : Eq_Eta

    use m_iterator
    use m_mat_invert

    use m_trimat_invert

    use m_ts_tri_scat, only : has_full_part

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

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
    real(dp), intent(in) :: Hs(n_nzs,nspin), Ss(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin)
    real(dp), intent(in) :: kT, Qtot
    real(dp), intent(inout) :: Ef
    logical, intent(out) :: converged

! ******************* Computational arrays *******************
    integer :: nzwork, n_s
    complex(dp), pointer :: zwork(:)
    type(zTriMat) :: zwork_tri, GF_tri
    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D) :: spH, spS
    ! local sparsity pattern in local SC pattern
    type(dSpData2D) :: spDM, spEDM ! EDM for dummy argument
    ! The different sparse matrices that will surmount to the integral
    ! These two lines are in global update sparsity pattern (UC)
    type(zSpData2D) ::  spuDM, spuEDM ! EDM for dummy argument
    ! To figure out which parts of the tri-diagonal blocks we need
    ! to calculate
    logical, pointer :: calc_parts(:) => null()
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    real(dp) :: kw, kpt(3), bkpt(3)
    complex(dp), pointer :: zDM(:,:)
    complex(dp) :: W
! ************************************************************

! ******************** Loop variables ************************
    type(itt2) :: SpKp
    integer, pointer :: ispin, ikpt
    integer :: ia, iEl, io, idx
    integer :: no, no_u_TS
#ifdef MPI
    integer :: MPIerror
#endif
! ************************************************************

    ! Number of supercells
    n_s = size(sc_off,dim=2)
    no_u_TS = no_u - no_Buf

    call newzTriMat(zwork_tri,N_tri_part,tri_parts,'GFinv')
    nzwork = elements(zwork_tri)

    call init_TriMat_inversion(zwork_tri)

    call newzTriMat(GF_tri,N_tri_part,tri_parts,'GF')

    call init_mat_inversion(maxval(tri_parts))

    call re_alloc(calc_parts,1,N_tri_part)
    calc_parts(:) = .true.

    no = 0
    zwork => val(GF_tri)
    do iEl = 1 , N_Elec

       io = TotUsedOrbs(Elecs(iEl))
       Elecs(iEl)%Sigma => zwork(no+1:no+io**2)
       no = no + io ** 2

       if ( Elecs(iEl)%DM_update /= 0 ) cycle

       io  = Elecs(iEl)%idx_o
       io  = io - orb_offset(io)
       idx = io + TotUsedOrbs(Elecs(iEl)) - 1

       do ia = 1 , N_tri_part
          if ( has_full_part(N_tri_part,tri_parts,ia,io,idx) ) then
             calc_parts(ia) = .false.
          end if
       end do
    end do

    zwork => val(zwork_tri)

#ifdef MPI
    call newDistribution(no_u,MPI_Comm_Self,fdist,name='TS-fake dist')
#else
    call newDistribution(no_u,-1           ,fdist,name='TS-fake dist')
#endif
    
    call newzSpData1D(ts_sp_uc,fdist,spH,name='TS spH')
    call newzSpData1D(ts_sp_uc,fdist,spS,name='TS spS')

    call newzSpData2D(tsup_sp_uc,1,fdist, spuDM, name='TS spuDM')
    zDM => val(spuDM)
    call newdSpData2D(ltsup_sp_sc,1, sp_dist,spDM   ,name='TS spDM')
    
    call itt_init  (SpKp,end1=nspin,end2=ts_nkpnt)
    call itt_attach(SpKp,cur1=ispin,cur2=ikpt)
    
    call init_val(spDM)
    do while ( .not. itt_step(SpKp) )

       kpt(:) = ts_kpoint(:,ikpt)
       call kpoint_convert(ucell,kpt,bkpt,1)
       kw = 0.5_dp / Pi * ts_kweight(ikpt)
       if ( nspin == 1 ) kw = kw * 2._dp
       
#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',1)
#endif

       call create_HS(sp_dist,sparse_pattern, &
            Ef, &
            N_Elec, Elecs, no_u, n_s, & ! electrodes, SIESTA size
            n_nzs, Hs(:,ispin), Ss, sc_off, &
            spH, spS, kpt, &
            nzwork, zwork)

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',2)
#endif

       call init_val(spuDM)
       cE%exist = .true.
       cE%fake  = .true.
       cE%e = dcmplx(0._dp, Eq_eta)
       cE%idx = 1
       if ( Node == Nodes - 1 ) cE%fake = .false.

       call read_next_GS(ispin, ikpt, bkpt, &
            cE, N_Elec, uGF, Elecs, &
            nzwork, zwork, .false., forward = .false. )
       do iEl = 1 , N_Elec
          call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
               non_Eq = .false. )
       end do
       
       call prepare_invGF(cE, no_u_TS, zwork_tri, &
            N_Elec, Elecs, &
            spH=spH , spS=spS)
       
       if ( .not. cE%fake ) then
          call invert_TriMat(zwork_tri,GF_tri,calc_parts)
       end if
          
       W  = dcmplx(kw,0._dp)
       call add_DM( spuDM, W, spuEDM, W, &
            GF_tri, N_Elec, Elecs, DMidx=1)

#ifdef MPI
       call timer('TS_comm',1)
       io = size(zDM)
       call MPI_Bcast(zDM(1,1),io,MPI_Double_Complex, &
            Nodes - 1, MPI_Comm_World,MPIerror)
       call timer('TS_comm',2)
#endif

       call add_k_DM(spDM, spuDM, 1, &
            spEDM, spuEDM, 1, &
            n_s, sc_off, kpt, non_Eq = .false. )

    end do ! spin, k-point

    call itt_destroy(SpKp)

#ifdef TRANSIESTA_DEBUG
    write(*,*) 'Completed TRANSIESTA - CHARGE'
#endif

    call ts_qc_Fermi(sp_dist, sparse_pattern, &
         nspin, n_nzs, DM, Ss, Qtot, spDM, Ef, converged)
    
    call de_alloc(calc_parts)

    call delete(zwork_tri)
    call delete(GF_tri)

    call delete(spH)
    call delete(spS)

    call delete(spuDM)
    call delete(spDM)

    call delete(fdist)

    call clear_TriMat_inversion()
    call clear_mat_inversion()

  end subroutine ts_trik_Fermi
  
  subroutine add_DM(DM, DMfact, EDM, EDMfact, &
       GF_tri, &
       N_Elec, Elecs, &
       DMidx,EDMidx)

    use class_zSpData2D
    use class_Sparsity
    use class_zTriMat

    use m_ts_electype

    ! The DM and EDM equivalent matrices
    type(zSpData2D), intent(inout) :: DM
    complex(dp), intent(in) :: DMfact
    type(zSpData2D), intent(inout) :: EDM
    complex(dp), intent(in) :: EDMfact

    ! The Green's function
    type(zTriMat), intent(inout) :: GF_tri
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! the index of the partition
    integer, intent(in) :: DMidx
    integer, intent(in), optional :: EDMidx

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: D(:,:), E(:,:)
    complex(dp), pointer :: Gf(:)
    integer :: io, ind, nr, iu, idx, i1, i2
    logical :: hasEDM

    s  => spar(DM)
    call attach(s, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=nr)
    D => val(DM)
    hasEDM = initialized(EDM)
    if ( hasEDM ) E => val(EDM)
    
    i1 = DMidx
    i2 = i1
    if ( present(EDMidx) ) i2 = EDMidx

    Gf => val(Gf_tri)

    if ( hasEDM ) then

       do io = 1 , nr
          if ( l_ncol(io) == 0 ) cycle

          ! The update region equivalent GF part
          iu = io - orb_offset(io)
          
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
             
             idx = index(Gf_tri,iu,l_col(ind) - orb_offset(l_col(ind)))
             
             D(ind,i1) = D(ind,i1) - GF(idx) * DMfact
             E(ind,i2) = E(ind,i2) - GF(idx) * EDMfact
             
          end do
       end do

    else

       do io = 1 , nr
          if ( l_ncol(io) == 0 ) cycle

          ! The update region equivalent GF part
          iu = io - orb_offset(io)
          
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
             
             idx = index(Gf_tri,iu,l_col(ind) - orb_offset(l_col(ind)))
             
             D(ind,i1) = D(ind,i1) - GF(idx) * DMfact
             
          end do
       end do

    end if

  end subroutine add_DM

  ! creation of the GF^{-1}.
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_invGF(cE, no_u,GFinv_tri, &
       N_Elec, Elecs, spH, spS)

    use class_Sparsity
    use class_zSpData1D
    use class_zTriMat
    use m_ts_electype
    use m_ts_tri_scat, only : insert_Self_Energies
    use m_ts_cctype, only : ts_c_idx

    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    integer, intent(in) :: no_u
    type(zTriMat), intent(inout) :: GFinv_tri
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D), intent(inout) :: spH,  spS
        
    ! Local variables
    complex(dp) :: Z
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: H(:), S(:)
    complex(dp), pointer :: Gfinv(:)
    integer :: io, iu, ind, idx, nr

    if ( cE%fake ) return

    Z = cE%e

    sp => spar(spH)
    H  => val (spH)
    S  => val (spS)

    call attach(sp, n_col=l_ncol, list_ptr=l_ptr, list_col=l_col, &
         nrows_g=nr)

    Gfinv  => val(Gfinv_tri)

    ! Initialize
    GFinv(:) = dcmplx(0._dp,0._dp)

    ! We will only loop in the central region
    ! We have constructed the sparse array to only contain
    ! values in this part...
    do io = 1, nr
       if ( l_ncol(io) == 0 ) cycle

       iu = io - orb_offset(io)

       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 

          ! Notice that we transpose back here...
          ! See symmetrize_HS_kpt
          idx = index(Gfinv_tri,l_col(ind)-orb_offset(l_col(ind)),iu)

          GFinv(idx) = Z * S(ind) - H(ind)
       end do
    end do

    do io = 1 , N_Elec
       call insert_Self_Energies(no_u, Gfinv_tri, Elecs(io))
    end do

  end subroutine prepare_invGF
   
end module m_ts_trik
