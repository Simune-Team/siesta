! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.


! Considerations for optimizing the code:
!  1. Instead of initiating send/recv element wise (done a lot of
!     places) we should utilise a continuous link-up.
!     Utilise MPI_Send_Init/MPI_Recv_Init
!     Best to create a type which takes care of this (there
!     might be several different sizes of communication)
!     At least a dozen of single Send/Recv are used throughout
!     the code. This is sub-optimal, however, I do not know 
!     the exact performance penalty.

module m_tbt_trik

  use precision, only : dp
  use m_timer, only : timer_start, timer_stop, timer_get
  use m_pivot_array, only : Npiv, ipiv

  use m_region
  use m_ts_electype
  use m_ts_sparse_helper, only : create_HS
  use m_ts_tri_common, only : nnzs_tri

  use m_verbosity, only : verbosity
  use m_tbt_hs
  use m_tbt_regions, only : sp_uc, sp_dev, r_aDev
#ifdef NOT_WORKING
  use m_tbt_regions, only : r_oEl
#endif
  use m_tbt_regions, only : r_oDev, r_oEl_alone, r_oElpD
  use m_tbt_regions, only : r_aBuf
  use m_tbt_tri_init, only : ElTri, DevTri

  implicit none

  public :: tbt_trik

  private
  
contains

  subroutine tbt_trik(ispin, N_Elec, Elecs, TSHS, nq, uGf)

    use units, only : Pi, eV
    use parallel, only : Node, Nodes, IONode

#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D
    use class_dSpData1D
    use class_dSpData2D
    use class_zSpData2D
    use class_zTriMat

    use dictionary
    use nf_ncdf, only: hNCDF, ncdf_close

    use m_ts_electype
    ! Self-energy read
    use m_tbt_gf
    ! Self-energy expansion
    use m_ts_elec_se

    use m_tbt_kpoint, only : nkpnt, kpoint, kweight

    use m_tbt_options, only : save_DATA, percent_tracker, N_eigen
#ifdef NCDF_4
    use m_tbt_options, only : cdf_fname, cdf_fname_sigma, cdf_fname_proj
    use m_tbt_sigma_save
#endif

    use m_tbt_contour

    use m_ts_cctype

    use m_iterator
    use m_mat_invert

    use m_trimat_invert

    use m_tbt_tri_init

    ! Gf calculation
    use m_ts_trimat_invert

    ! Gf.Gamma.Gf from ts_tri_scat and all remaining routines
    use m_tbt_tri_scat

    use m_tbt_save
#ifdef NCDF_4
    ! The projections
    use m_tbt_proj, only : N_mol, mols, N_proj_T, proj_T
    use m_tbt_proj, only : N_proj_ME, proj_ME, tLvlMolEl
    use m_tbt_proj, only : proj_LME_assoc
    use m_tbt_proj, only : proj_update, proj_cdf_save_S_D
    use m_tbt_proj, only : proj_bMtk, proj_cdf_save_bGammak
    use m_tbt_proj, only : proj_Mt_mix, proj_cdf_save
    use m_tbt_proj, only : proj_cdf_save_J

    use m_tbt_dH, only : read_next_dH, clean_dH
#endif

    use m_tbt_sparse_helper, only : create_region_HS
    use m_tbt_kregions, only : n_k, r_k
    use m_tbt_kregions, only : kregion_k, kregion_step, calc_GS_k

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in) :: ispin
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    type(tTSHS), intent(inout) :: TSHS
    integer, intent(in) :: nq(N_Elec)
    integer, intent(in) :: uGF(N_Elec)

! ******************* Computational arrays *******************
    integer :: nzwork, nGfwork, nmaxwork
    complex(dp), pointer :: zwork(:), Gfwork(:), maxwork(:)
    type(zTriMat) :: zwork_tri, GF_tri
    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    ! Just as in transiesta, these matrices are transposed.
    type(zSpData1D) :: spH, spS
    real(dp), pointer :: H(:,:), S(:)
    ! To figure out which parts of the tri-diagonal blocks we need
    ! to calculate
    logical :: calc_GF_DOS
    logical, allocatable :: prep_El(:)
    logical, allocatable :: A_parts(:)
! ************************************************************

! ********************** Result arrays ***********************
    real(dp), allocatable, target :: allDOS(:,:)
    real(dp), pointer :: DOS(:,:), DOS_El(:,:)
#ifdef NOT_WORKING
    type tLife
       real(dp), allocatable :: life(:)
    end type tLife
    type(tLife) :: life(N_Elec)
#endif

    ! The transmission values between the different electrodes
    real(dp) :: T(N_Elec+1,N_Elec)
    ! The transmission eigenvalues
    real(dp), allocatable :: Teig(:,:,:)

#ifdef NCDF_4
    ! If the user requests projection molecules
    ! Instead of overwriting the original
    ! electrode, we instead prepare ONE electrode
    ! That can accomodate all projections.
    ! We will at each energy re-create the projection
    ! Matrix. This will take a little time, but it should
    ! be rather fast.
    type(Elec) :: El_p
    type(tLvlMolEl), pointer :: p_E
    real(dp), allocatable :: pDOS(:,:,:)
    real(dp), allocatable :: bTk(:,:), bTkeig(:,:,:)
    type(dSpData1D) :: orb_J
#endif
! ************************************************************

! ****************** Electrode variables *********************
    integer :: nGFGGF ! For the triple-product
    complex(dp), pointer :: GFGGF_work(:) => null()
    integer :: ntt_work
    complex(dp), pointer :: tt_work(:), eig(:)
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    type(tRgn) :: pvt
    real(dp) :: kpt(3), bkpt(3), wkpt
    integer :: info
    logical :: T_all
#ifdef TBT_PHONON
    type(ts_c_idx) :: cOmega
    real(dp) :: omega
#endif
    integer :: n_kpt
! ************************************************************

! ******************** Loop variables ************************
    type(itt1) :: Kp
    integer :: N_E
    integer, pointer :: ikpt
    integer :: iEl, jEl
#ifdef NCDF_4
    integer :: ipt, it
#endif
    integer :: iE, iE_N
    integer :: no, io
    integer :: pad_LHS, pad_RHS
    type(tNodeE) :: nE
! ************************************************************

! ******************* DATA file variables ********************
#ifdef NCDF_4
    ! The netcdf handle for TBT
    type(hNCDF) :: TBTcdf, PROJcdf
#else
    ! Units for IO of ASCII files
    integer, allocatable :: iounits(:), iounits_El(:)
#endif
! ************************************************************

! ******************** Timer variables ***********************
    real(dp) :: loop_time, init_time
    integer  :: last_progress_print
! ************************************************************

    T_all = 'T-all' .in. save_DATA

    ! Create the back-pivoting region
    call rgn_init(pvt,nrows_g(TSHS%sp),val=0)
    do io = 1 , r_oDev%n
       pvt%r(r_oDev%r(io)) = io
    end do

    H => val(TSHS%H_2D)
    S => val(TSHS%S_1D)

    allocate(nE%iE(0:Nodes-1))
    allocate(nE%E(0:Nodes-1))

#ifdef NCDF_4
    ! In case we have projections
    ! Then we allocated the maximum projection
    nullify(El_p%Sigma)
    El_p%no_used = 0
    if ( N_proj_ME > 0 ) then
       no = 0
       do it = 1 , N_proj_ME
          no = max(no,proj_ME(it)%mol%orb%n)
       end do
       El_p%no_used = no
       if ( no == 0 ) then
          call die('Error in determining the maximum size of projection &
               &subspace.')
       end if
       allocate(El_p%Sigma(no**2))

       ! We utilize the Sigma array as we then easily can
       ! reduce the size of Gamma by using pointers
       ! This will leverage MANY copies around the code!

       ! Count maximum number of RHS projections
       io = 0
       do it = 1 , N_proj_T
          io = max(io,size(proj_T(it)%R))
       end do
       ! We allocate data segment for retaining all information
       ! We need to save bGammak, T->N_Elec
       ! In this data array we save 3 quantities:
       allocate(bTk(io+1,N_proj_T))
       allocate(pDOS(r_oDev%n,2,N_proj_T))
       if ( N_eigen > 0 ) then
          allocate(bTkeig(N_eigen,io+1,N_proj_T))
       else
          allocate(bTkeig(1,1,1))
       end if

    end if
#endif

    ! The downfolding need not retain the full structure
    ! Here we calculate the maximum size we need to
    ! be able to down-fold the self-energies to the requested
    ! region.
    pad_LHS = 0
    do iEl = 1 , N_Elec
       if ( ElTri(iEl)%n < 2 ) &
            call die('All regions must be at least 2 blocks big. Error')
       pad_LHS = max(pad_LHS,fold_elements(ElTri(iEl)%n,ElTri(iEl)%r))
       ! In case we are doing in-core calculations, 
       ! we need a certain size for the electrode calculation.
       ! Sadly this is "pretty" big.
       if ( .not. Elecs(iEl)%out_of_core ) then
          no = Elecs(iEl)%no_u ** 2
          ! Determine work array size
          ! H00, H01, S00 and S01 in the work array
          io = no * 4
          if ( Elecs(iEl)%no_u /= Elecs(iEl)%no_used ) io = io + no
          io = io + no * 8
          if ( 'DOS-Elecs' .in. save_DATA ) then
             io = io + no
          end if
          pad_LHS = max(pad_LHS,io)
       end if
    end do

    ! Correct the actual padding by subtracting the 
    ! initial size of the tri-diagonal matrix
    pad_LHS = pad_LHS - nnzs_tri(DevTri%n,DevTri%r)

    ! We now have the maximum size of the electrode down-folding
    ! region.
   
    ! Calculate the size required for the Gf.G.Gf product
    ! Note that we take this with respect to the 
    ! down-folded scattering states (Gamma_L -> Gamma_DL)
    no = maxval(Elecs(:)%o_inD%n)
#ifdef NCDF_4
    if ( N_proj_ME > 0 ) then
       do it = 1 , N_proj_ME
          no = max(no,proj_ME(it)%mol%orb%n)
       end do
       ! Fake a too large electrode
       io = Elecs(1)%o_inD%n
       Elecs(1)%o_inD%n = no
    end if
#endif

    ! Calculate the maximum matrix size which corresponds
    ! to the device region. Here we need to take into account the
    ! work-array size of the GFGGF triple product.
    call GFGGF_needed_worksize(DevTri%n,DevTri%r, &
         N_Elec, Elecs, pad_RHS, nGFGGF)
       
#ifdef NCDF_4
    if ( N_proj_ME > 0 ) then
       ! Reset the too large electrode
       Elecs(1)%o_inD%n = io
    end if
#endif

    ! The minimum padding must be the maximum value of the 
    ! 1) padding required to contain a down-folding region
    ! 2) padding required to calculate the entire Gf column 
    pad_LHS = max(pad_LHS,pad_RHS)

    ! In case the user requests eigenchannel calculations
    ! we need to do a work-query from lapack
    if ( N_Eigen > 0 ) then
       ! 'no' is the maximum Gamma size in device region

       ! In tbtrans we need the GFGGF array
       ! to be at least the size of the maximum electrode
       nGFGGF = max(nGFGGF,no**2)

       ! Allocate Teig
       allocate(Teig(N_eigen,N_Elec,N_Elec))
       
       ! pre-allocate arrays for retrieving sizes
       nullify(zwork)
       allocate(zwork(no))
       eig => zwork(:)
       
       ! 'no' is the maximum size of the electrodes
       ! Work-query of eigenvalue calculations
       call zgeev('N','N',no,zwork,no,eig,zwork,1,zwork,1, &
            zwork,-1,Teig(1,1,1),info)
       if ( info /= 0 ) then
          print *,info
          call die('zgeev: could not determine optimal workarray &
               &size.')
       end if
       ! Minimum padding (eigenvalues) + lwork size
       ! The optimal lwork is stored in zwork(1)
       info = zwork(1)
       info = info + no
       deallocate(zwork)
       pad_LHS = max(pad_LHS,info)

    else

       ! Allocate a dummy Teig for succesfull passing it
       ! to a routine when compiling with warnings
       allocate(Teig(1,1,1))

    end if

    
    ! Note that padding is the extra size to be able to calculate
    ! the spectral function in the BTD format

    ! We allocate for the matrix
    call print_memory('LHS',pad_LHS)
    call newzTriMat(zwork_tri,DevTri%n,DevTri%r,'GFinv', &
         padding = pad_LHS )

    ! The RHS will be the array that retains the 
    ! self-energies and nothing more.
    ! Here we pad with the missing elements to contain
    ! all self-energies.
    pad_RHS = nnzs_tri(DevTri%n,DevTri%r)
    io = 0
    do iEl = 1 , N_Elec
       io = io + max(TotUsedOrbs(Elecs(iEl)),Elecs(iEl)%o_inD%n)**2
    end do
    pad_RHS = io - pad_RHS
    pad_RHS = max(pad_RHS,0) ! truncate at 0
    ! Pad so that GFGGF can be filled, we prefer padding over 
    ! allocating a new region, 
    ! this is only preferred in tbtrans as there might be 
    ! padding any-way.
    pad_RHS = max(pad_RHS,nGFGGF)
    
    call print_memory('RHS',pad_RHS)
    call newzTriMat(GF_tri,DevTri%n,DevTri%r,'GF', &
         padding = pad_RHS )

    ! Create the work-array...
    zwork   => val(zwork_tri,all=.true.)
    nzwork  =  size(zwork)
    Gfwork  => val(Gf_tri,all=.true.)
    nGfwork =  size(Gfwork)
    if ( nGfwork > nzwork ) then
       nmaxwork = nGfwork
       maxwork => Gfwork(:)
    else
       nmaxwork = nzwork
       maxwork => zwork(:)
    end if

    ! Point the work-array for eigenvalue calculation
    if ( N_eigen > 0 ) then
       ! 'no' is still maximum Gamma size in the device region
       ! 'info' is lwork size + eigen size

       ! The back of zwork becomes additional work 
       ! array for eigenvalue calculation
       tt_work => zwork(nzwork-info+1:nzwork-no)
       ntt_work = size(tt_work)
       eig => zwork(nzwork-no+1:nzwork)

    end if

    ! We have increased padding so that GFGGF can be
    ! fitted in the back of the array
    GFGGF_work => GFwork(nGfwork-nGFGGF+1:nGfwork)

    ! Initialize the tri-diagonal inversion routine
    call init_TriMat_inversion(zwork_tri)
    call init_BiasTriMat_inversion(zwork_tri)

    ! initialize the matrix inversion tool
    no = maxval(DevTri%r)
    do iEl = 1 , N_Elec
       no = max(no,maxval(ElTri(iEl)%r))
    end do
    call init_mat_inversion(no)

    ! Figure out how to setup the Green function energy
    ! point
    !call prep_Gfinv_algo(zwork_tri,TSHS%sp,r_oDev,loop_dev)

    ! we use the GF as a placement for the self-energies
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
       io = max(TotUsedOrbs(Elecs(iEl)),Elecs(iEl)%o_inD%n)
       Elecs(iEl)%Sigma => Gfwork(no+1:no+io**2)
       no = no + io ** 2

       ! we have already allocated the H,S, Gamma arrays.

#ifdef NOT_WORKING
       ! Calculate size of life-time array
       io = r_oEl(iEl)%n - TotUsedOrbs(Elecs(iEl))
       allocate(life(iEl)%life(io))
#endif

    end do

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
    ! Notice that we DO need it to be the SIESTA size.
    no = nrows_g(sp_uc)
#ifdef MPI
    call newDistribution(no,MPI_Comm_Self,fdist,name='TBT-fake dist')
#else
    call newDistribution(no,-1           ,fdist,name='TBT-fake dist')
#endif
    
    ! The Hamiltonian and overlap matrices (in Gamma calculations
    ! we will not have any phases, hence, it makes no sense to
    ! have the arrays in complex)
    ! TODO move into a Data2D (could reduce overhead of COMM)
    call newzSpData1D(sp_uc,fdist,spH,name='TBT spH')
    call newzSpData1D(sp_uc,fdist,spS,name='TBT spS')

    ! Allocate data-collecting arrays
    io = r_oDev%n
    jEl = 0
    if ( 'DOS-Elecs' .in. save_DATA ) then
       ! Allocate density of states for the electrodes
       do iEl = 1 , N_Elec
          jEl = max(jEl,Elecs(iEl)%no_u)
       end do
    end if
    ! this allows us to re-use the data arrays for
    ! the DOS storages between electrode DOS and device DOS
    allocate(allDOS(max(io,jEl),N_Elec+1))
    ! Note we point to the same data-segments to reduce memory
    ! requirement
    DOS => allDOS(1:r_oDev%n,1:N_Elec+1)
    DOS_El => allDOS(1:jEl,1:N_Elec)

    ! Initialize which parts to calculate
    allocate(prep_El(N_Elec))
    ! Default to all electrode Green function parts
    prep_El = .true.
    if ( (.not. T_all) .and. &
         ('DOS-A-all' .nin. save_DATA) ) then
       
       ! We only need to prepare the Green function
       ! for Gf.G.Gf product for all electrodes
       ! besides the last one.
       ! Note that calc_GF_DOS has precedence for
       ! calculating the entire Gf
       prep_El(N_elec) = .false.
       
    end if
    
    allocate(A_parts(DevTri%n))
    A_parts = .true.
    ! If the user ONLY wants the transmission function then we
    ! should tell the program not to calculate any more
    ! than needed.
    calc_GF_DOS = 'DOS-Gf' .in. save_DATA
    if ( 'DOS-A' .nin. save_DATA ) then

       ! We have a couple of options that requires
       ! special parts of the Green function

       ! The first partition is whether we want
       ! certain quantities from all electrodes
       A_parts(:) = .false.
       if ( T_all .or. &
            ('T-sum-out' .in. save_DATA) ) then

          ! We need _all_ diagonal blocks of the spectral function
          do iEl = 1 , N_Elec
             do io = 1 , Elecs(iEl)%o_inD%n
                jEl = which_part(Gf_tri, &
                     rgn_pivot(r_oDev,Elecs(iEl)%o_inD%r(io)) )
                A_parts(jEl) = .true.
             end do
          end do

       else
          
          ! we only want the spectral function
          ! for all other electrodes than the first one
          do iEl = 2 , N_Elec
             do io = 1 , Elecs(iEl)%o_inD%n
                jEl = which_part(Gf_tri, &
                     rgn_pivot(r_oDev,Elecs(iEl)%o_inD%r(io)) )
                A_parts(jEl) = .true.
             end do
          end do

       end if
       
    end if

    ! Initialize
    if ( 'T-Gf' .in. save_DATA ) then
       
       ! The only thing we are calculating is the
       ! transmission.
       ! In this case we can limit the calculation
       ! space to speed things up
       prep_El(:) = .true.
       if ( .not. T_all ) then
          prep_El(N_Elec) = .false.
       end if

       ! reset all parts (A_parts is not used when
       ! only using the Green function)
       A_parts = .false.

    end if


#ifdef NCDF_4
    if ( ('orb-current'.in.save_DATA) .or. &
         ('proj-orb-current'.in.save_DATA) ) then
       call newdSpData1D(sp_dev,fdist,orb_J,name='TBT orb_J')
       ! Initialize to 0, it will be overwritten for all
       ! values, so there is no need to initialize it
       ! in the orb_current routine
       call init_val(orb_J)
    end if
#endif

    ! Either we have the k-regions, or the regular k-grid
    if ( n_k == 0 ) then
       ! We have the regular k-grid
       n_kpt = nkpnt
    else
       ! Calculate the permutations of the k-point regions
       n_kpt = 1
       do iEl = 0 , n_k
          n_kpt = n_kpt * size(r_k(iEl)%wkpt)
       end do
       r_k(:)%ik = 1
       r_k(n_k)%ik = 0

    end if
       
    ! start the itterators
    call itt_init  (Kp,end=n_kpt)
    ! point to the index iterators
    call itt_attach(Kp,cur=ikpt)

    ! Number of energy points
    N_E = N_TBT_E()

#ifdef NCDF_4
    ! Open the NetCDF handles

    ! *.TBT.nc file
    call open_cdf_save(cdf_fname, TBTcdf)
    if ( N_proj_ME > 0 ) then
       ! *.TBT.Proj.nc file
       call open_cdf_save(cdf_fname_proj, PROJcdf)
    end if

#else
    ! Allocate units for IO ASCII
    allocate(iounits(1+(N_Elec*2+2)*N_Elec)) ! maximum number of units
    iounits = -100
    call init_save(iounits,ispin,TSHS%nspin,N_Elec,Elecs, &
         save_DATA)
    allocate(iounits_El(N_Elec)) ! maximum number of units
    iounits_El = -100
    call init_save_Elec(iounits_El,ispin,TSHS%nspin,N_Elec,Elecs, &
         save_DATA)
#endif

    ! The current progress is 0%
    last_progress_print = 0

    ! In case we do several spin, then the estimated time
    ! is wrong for the second run, hence we start and stop, and
    ! record the initial time to take that into account
    call timer_start('E-loop')
    call timer_stop('E-loop')
    call timer_get('E-loop',totTime=init_time)
    
    do while ( .not. itt_step(Kp) )

       if ( n_k == 0 ) then

          kpt(:) = kpoint(:,ikpt)
          ! create the k-point in reciprocal space
          call kpoint_convert(TSHS%cell,kpt,bkpt,1)
          wkpt = kweight(ikpt)

       else
          
          ! step the k-region
          call kregion_step()

          ! Get the k-point of the full region
          call kregion_k(-1,bkpt,w=wkpt)

          call kpoint_convert(TSHS%cell,bkpt,kpt,-1)

       end if

#ifndef NCDF_4
       ! Step the k-points in the output files
       call step_kpt_save(iounits,n_kpt,bkpt,wkpt)
       call step_kpt_save(iounits_El,n_kpt,bkpt,wkpt)
#endif

       ! Start timer
       call timer_start('E-loop')
       
#ifdef TBTRANS_TIMING
       call timer('setup-HS',1)
#endif

       ! Work-arrays are for MPI distribution... (reductions)
       iE = size(S)
       if ( n_k == 0 ) then
          call create_HS(TSHS%dit, TSHS%sp, 0._dp, &
               N_Elec, Elecs, TSHS%no_u, product(TSHS%nsc), &
               iE,H(:,1), S(:), TSHS%sc_off, &
               spH, spS, kpt, &
               nmaxwork, maxwork) ! not used...
       else
          call create_region_HS(TSHS%dit,TSHS%sp, 0._dp, &
               TSHS%cell, n_k, r_k, TSHS%na_u, TSHS%lasto, &
               N_Elec, Elecs, TSHS%no_u, product(TSHS%nsc), &
               iE, H(:,1), S(:), TSHS%sc_off, spH, spS, &
               nmaxwork, maxwork)
       end if

#ifdef NCDF_4
       if ( N_mol > 0 ) then

          ! Read in the projections for this k-point
          call proj_update(PROJcdf,N_mol,mols,ikpt)

          ! Calculate the |><| S and save it immediately
          ! For all cases where we need the diagonal of the
          ! sparsity pattern we immediately calculate it...
          !call proj_cdf_save_S_D(cdf_fname_proj, N_mol, mols, ikpt, &
          !     spS,nmaxwork,maxwork)

       end if
#endif

#ifdef TBTRANS_TIMING
       call timer('setup-HS',2)
#endif

       ! Start energy-loop
       iE = 0
       cE = tbt_E(iE+Nodes-Node,step=Nodes) ! we read them backwards
       do while ( cE%exist ) 

          ! The actual energy-loop count
          iE_N = iE + Nodes - Node
          if ( cE%fake ) iE_N = 0

#ifdef TBT_PHONON
          ! Retrieve the frequency before converting energy to 
          ! \omega^2 + i \eta
          omega = real(cE%e,dp)
          ! Copy data, and form: \omega**2 + i \eta
          cOmega = cE
          cOmega%e = dcmplx(omega**2,dimag(cE%e))
#endif

          ! B-cast all nodes current energy segment
          call MPI_BcastNode(iE_N, cE%E, nE)

          ! Print out information about current progress.
          ! We print out a progress report every 5 %
          ! Calculate progress
          jEl = (itt_cur_step(Kp) - 1) * N_E + iE
          iEl = itt_steps(Kp) * N_E
          io = int(100._dp*real(jEl,dp)/real(iEl,dp))
          if ( io - last_progress_print >= percent_tracker ) then
             ! We have passed another 'percent_tracker'% of calculation time

             ! save current progress in integer form
             last_progress_print = io

             ! Stop timer
             call timer_stop('E-loop')

             ! Calculate time passed by correcting for the initial time
             call timer_get('E-loop',totTime=loop_time)
             loop_time = loop_time - init_time

             if ( IONode ) then
                loop_time = iEl * loop_time / (100._dp*real(jEl,dp))
                loop_time = loop_time * ( 100._dp - last_progress_print )
                write(*,'(a,f7.3,'' %, ETA in '',f20.3,'' s'')') &
                     'tbt: Calculated ', &
                     100._dp*real(jEl,dp)/real(iEl,dp), loop_time
             end if

             ! Start the timer again
             call timer_start('E-loop')

          end if

          call timer('read-GS',1)

          ! *******************
          ! * prep Sigma      *
          ! *******************
          ! We have reduced the electrode sizes to only one spin-channel
          ! Hence, it will ALWAYS be the first index
          if ( n_k == 0 ) then
             if ( 'DOS-Elecs' .in. save_DATA ) then
                call read_next_GS(1, ikpt, bkpt, &
                     cE, N_Elec, uGF, Elecs, &
                     nzwork, zwork, .false., forward = .false. , &
                     DOS = DOS_El )

                ! Immediately save the DOS
#ifdef NCDF_4
                call state_cdf_save_Elec(TBTcdf, ikpt, nE, N_Elec, Elecs, &
                     DOS_El, save_DATA)
#else
                call state_save_Elec(iounits_El,nE,N_Elec,Elecs,DOS_El, &
                     save_DATA )
#endif

             else
                call read_next_GS(1, ikpt, bkpt, &
                     cE, N_Elec, uGF, Elecs, &
                     nzwork, zwork, .false., forward = .false. )
             end if
          else
             call calc_GS_k(1, cE, N_Elec, Elecs, uGF, &
                  nzwork, zwork)
          end if

          call timer('read-GS',2)

#ifdef NCDF_4
          ! prepare dH
          call read_next_dH(no,bkpt,nE)
#endif

          call timer('SE-dwn',1)

          do iEl = 1 , N_Elec

             ! We do not have it as a non-equilibrium
             ! contour point as the downfolded Self-energy will
             ! create the Gamma.
             ! Hence it is a waste of time if this is done in this
             ! loop....
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .false. ) 


             ! Down-fold immediately :)
#ifdef TBT_PHONON
             call downfold_SE(cOmega,Elecs(iEl), spH, spS, &
                  r_oElpD(iEl), ElTri(iEl)%n, ElTri(iEl)%r, nzwork,zwork)
#else
             call downfold_SE(cE,Elecs(iEl), spH, spS, &
                  r_oElpD(iEl), ElTri(iEl)%n, ElTri(iEl)%r, nzwork,zwork)
#ifdef NOT_WORKING
             call downfold_SE_life(cE,Elecs(iEl), spH, spS, &
                  r_oElpD(iEl), ElTri(iEl)%n, ElTri(iEl)%r, life(iEl)%life, &
                  nzwork,zwork)
             print *,sum(life(iEl)%life)
#endif
#endif

          end do

          call timer('SE-dwn',2)

#ifdef NCDF_4
          call state_Sigma_save(cdf_fname_Sigma, ikpt, nE, &
               N_Elec, Elecs, nzwork, zwork)


          if ( N_proj_ME > 0 ) then

             call timer('Proj-Gam',1)

             do ipt = 1 , N_proj_ME
                
                io = proj_ME(ipt)%El%o_inD%n ** 2
                call proj_bMtk(proj_ME(ipt)%mol, &
                     proj_ME(ipt)%El%o_inD,proj_ME(ipt)%El%Gamma(1:io), &
                     proj_ME(ipt)%bGk,nzwork,zwork)

             end do

             ! Save the energies
             call cdf_save_E(PROJcdf,nE)
             
             ! Save the projected values
             call proj_cdf_save_bGammak(PROJcdf,N_proj_ME,proj_ME, &
                  ikpt,nE)

             call timer('Proj-Gam',2)

          end if
#endif

          ! Only calculate actual transmission if the user
          ! has requested so...
          if ( ('Sigma-only'.nin.save_DATA) ) then

          call timer('Gf-prep',1)

          ! *******************
          ! * prep GF^-1      *
          ! *******************
#ifdef TBT_PHONON
          call prepare_invGF(cOmega, zwork_tri, r_oDev, pvt, &
               N_Elec, Elecs, spH , spS)
#else
          call prepare_invGF(cE, zwork_tri, r_oDev, pvt, &
               N_Elec, Elecs, spH , spS)
#endif

          ! ********************
          ! * prep GF for scat *
          ! ********************
          if ( .not. cE%fake ) then
             ! We should calculate all nn partitions in case
             ! we want the DOS of the central region.
             ! So possibly we should do:
             !   all_nn = GFDOS
             if ( calc_GF_DOS ) then
                call invert_BiasTriMat_prep(zwork_tri,GF_tri, &
                     N_Elec, Elecs, all_nn = .true. )
#ifdef NCDF_4
             else if ( N_proj_T > 0 ) then
                call invert_BiasTriMat_prep(zwork_tri,GF_tri, &
                     N_Elec, Elecs, all_nn = .true. )
#endif
             else
                call invert_BiasTriMat_prep(zwork_tri,GF_tri, &
                     N_Elec, Elecs, has_El = prep_El )
             end if
             
          end if

          call timer('Gf-prep',2)


          ! Only calculate actual transmission if the user
          ! has requested so...
          if ( ('proj-only'.nin.save_DATA) ) then

          call timer('DOS-Gf-A-T',1)

          ! We have now calculated all block diagonal entries
          ! of the Green's function.
          ! This means that all necessary information to calculate
          ! the entire Green's function resides in GF_tri
          if ( calc_GF_DOS ) then
             if ( .not. cE%fake ) then
                call GF_DOS(r_oDev,Gf_tri,spS,DOS(:,1),nzwork,zwork)
#ifdef TBT_PHONON
                DOS(:,1) = omega * DOS(:,1)
#else
                if ( TSHS%nspin == 1 ) DOS(:,1) = 2._dp * DOS(:,1)
#endif
             end if
          end if

          ! ****************
          ! * Column Gf    *
          ! ****************
          if ( 'T-Gf' .in. save_DATA ) then

             ! We are allowed to calculate the transmission
             ! only by using the diagonal
             if ( .not. cE%fake ) then
                do iEl = 1 , N_Elec
                   if ( iEl == N_Elec .and. .not. T_all ) cycle
                   
                   call invert_BiasTriMat_rgn(GF_tri,zwork_tri, &
                        r_oDev, Elecs(iEl)%o_inD,only_diag=.true.)
                   
                   call GF_T(zwork_tri,Elecs(iEl), &
                        T(N_Elec+1,iEl), T(iEl,iEl), &
                        nGFGGF,GFGGF_work)

                end do
                
                ! Retrieve actual transmissions
                call GF_T_solve(N_Elec,T,T_all)

             end if

          else

          ! We loop over all electrodes
          do iEl = 1 , N_Elec
             if ( iEl == N_Elec .and. ( &
                  (.not. T_all) .and. &
                  ('DOS-A-all' .nin. save_DATA) ) ) cycle

             ! ******************
             ! * calc GF-column *
             ! ******************
             if ( .not. cE%fake ) then
                call invert_BiasTriMat_rgn(GF_tri,zwork_tri, &
                     r_oDev, Elecs(iEl)%o_inD)

                if ( 'T-sum-out' .in. save_DATA ) then
                   call Gf_Gamma(zwork_tri,Elecs(iEl),T(N_Elec+1,iEl))
                end if

             end if

             ! This small conversion of data, ensures
             ! that we do not need to create two EXACT
             ! same functions.
             ! Also, the only reason for doing this
             ! is that the down-folded Gamma will NEVER
             ! have any repetition.
             if ( .not. cE%fake ) then
                call GF_Gamma_GF(zwork_tri, Elecs(iEl), Elecs(iEl)%o_inD%n, &
                     A_parts, &
                     nGFGGF, GFGGF_work)
             end if

             if ( ('DOS-A' .in. save_DATA) .and. .not. cE%fake ) then

                ! Calculate the DOS from the spectral function
                call A_DOS(r_oDev,zwork_tri,spS,DOS(:,1+iEl))
#ifdef TBT_PHONON
                DOS(:,1+iEl) = omega * DOS(:,1+iEl)
#else
                if ( TSHS%nspin == 1 ) DOS(:,1+iEl) = 2._dp * DOS(:,1+iEl)
#endif
                
#ifdef NCDF_4
                if ( 'orb-current' .in. save_DATA ) then

                   call orb_current(spH,zwork_tri,r_oDev,orb_J)

                end if
#endif
             end if

#ifdef NCDF_4
             if ( 'orb-current' .in. save_DATA ) then
                
                ! We need to save it immediately, we
                ! do not want to have several arrays in memory
                call state_cdf_save_J(TBTcdf, ikpt, nE, Elecs(iEl), &
                     orb_J, save_DATA)
                
             end if
#endif
             
             do jEl = 1 , N_Elec
                ! Calculating iEl -> jEl is the
                ! same as calculating jEl -> iEl, hence if we
                ! do not wish to assert this is true, we do not
                ! calculate this.
                if ( (.not. T_all) .and. &
                     jEl < iEl ) cycle
                if ( ('T-sum-out' .nin. save_DATA ) .and. &
                     iEl == jEl ) cycle

                ! Notice that the Gf.G1.Gf.G2 can be performed
                ! for all other electrodes as long as we
                ! have the block diagonal that constitutes the 
                ! "right" electrode.
                
                if ( .not. cE%fake ) then
                   if ( N_eigen > 0 ) then
                      io = Elecs(jEl)%inDpvt%n
                      call A_Gamma_Block(zwork_tri,Elecs(jEl),T(jEl,iEl), &
                           nGFGGF, GFGGF_work)
                      call TT_eigen(io,GFGGF_work, ntt_work, tt_work, eig)
                      ! Copy the eigenvalues over
                      do io = 1 , N_eigen
                         Teig(io,jEl,iEl) = dreal(eig(io))
                      end do
                   else
                      call A_Gamma(zwork_tri,Elecs(jEl),T(jEl,iEl))
                   end if
                end if
                
             end do

          end do

          end if


          call timer('DOS-Gf-A-T',2)

          ! Save the current gathered data
#ifdef NCDF_4
          call state_cdf_save(TBTcdf, ikpt, nE, N_Elec, Elecs, &
               DOS, T, N_eigen, Teig, save_DATA)
#else
          call state_save(iounits,nE,N_Elec,Elecs,DOS, T, &
               N_eigen, Teig, &
               save_DATA )
#endif

          end if ! .not. proj-only
          end if ! .not. Sigma-only


#ifdef NCDF_4

          if ( N_proj_T > 0 ) then

          call timer('T-proj',1)

          ! Calculate the projections
          do ipt = 1 , N_proj_T

             ! Associate projection
             call proj_LME_assoc(p_E,proj_T(ipt)%L)

             if ( cE%fake ) then
#ifdef NCDF_4
                ! We need to fake the IO node to call the save routine
                ! this aint pretty, however it relieves a lot of
                ! superfluous checks in the following block
                if ( ('proj-orb-current' .in. save_DATA) .and. p_E%idx > 0 ) then
                   call proj_cdf_save_J(PROJcdf, ikpt, nE, proj_T(ipt)%L, &
                        orb_J, save_DATA)
                end if
#endif
                cycle
             end if
            ! We have now calculated all block diagonal entries
            ! of the Green's function.
            ! This means that all necessary information to calculate
            ! the entire Green's function resides in GF_tri

            !if ( 'DOS-Gf' .in. save_DATA ) then
            !   if ( .not. cE%fake ) then
            !      call GF_DOS(r_oDev,Gf_tri,spS,DOS(:,1),nzwork,zwork)
#ifdef TBT_PHONON
            !       DOS(:,1) = omega * DOS(:,1)
#else
            !      if ( TSHS%nspin == 1 ) DOS(:,1) = 2._dp * DOS(:,1)
#endif
            !   end if
            !end if

            ! ****************
            ! * Column Gf    *
            ! ****************
            if ( p_E%idx > 0 ) then

               no = p_E%ME%mol%orb%n
               ! We have a Left projection
               ! Insert pointer
               El_p%Gamma => El_p%Sigma(:)
               ! Here we re-create the projection matrix that replaces
               ! the scattering state
               call proj_Mt_mix(p_E%ME%mol,p_E%idx,El_p%Gamma, p_E%ME%bGk)

               call invert_BiasTriMat_rgn(GF_tri,zwork_tri, &
                    r_oDev, p_E%ME%mol%orb)

               if ( 'proj-T-sum-out' .in. save_DATA ) then
                  ! Copy over pivoting table and the size
                  El_p%inDpvt%n =  p_E%ME%mol%pvt%n
                  El_p%inDpvt%r => p_E%ME%mol%pvt%r
                  call Gf_Gamma(zwork_tri,El_p, &
                       bTk(size(proj_T(ipt)%R)+1,ipt))
               end if


            else
               
               iEl = -p_E%idx
               no = Elecs(iEl)%o_inD%n
               El_p%Gamma => Elecs(iEl)%Gamma(:)

               call invert_BiasTriMat_rgn(GF_tri,zwork_tri, &
                    r_oDev, Elecs(iEl)%o_inD)

               if ( 'proj-T-sum-out' .in. save_DATA ) then
                  ! This should work, but I currently do not allow it :(
                  call Gf_Gamma(zwork_tri,Elecs(iEl), &
                       bTk(1+size(proj_T(ipt)%R),ipt))
               end if

            end if

            call GF_Gamma_GF(zwork_tri, El_p, no, &
                 (/(.true.,jEl=1,DevTri%n)/), &
                 nGFGGF, GFGGF_work)

            if ( ('proj-DOS-A' .in. save_DATA) .and. p_E%idx > 0 ) then

               ! Calculate the DOS from the spectral function
               call A_DOS(r_oDev,zwork_tri,spS,pDOS(:,2,ipt))
#ifdef TBT_PHONON
               pDOS(:,2,ipt) = omega * pDOS(:,2,ipt)
#else
               if ( TSHS%nspin == 1 ) pDOS(:,2,ipt) = 2._dp * pDOS(:,2,ipt)
#endif
               
#ifdef NCDF_4
               if ( 'proj-orb-current' .in. save_DATA ) then

                  call orb_current(spH,zwork_tri,r_oDev,orb_J)

                  ! We need to save it immediately, we
                  ! do not want to have several arrays in the
                  ! memory
                  call proj_cdf_save_J(PROJcdf, ikpt, nE, proj_T(ipt)%L, &
                       orb_J, save_DATA)

               end if
#endif
               
            end if

            ! Loop on RHS projections
            do jEl = 1 , size(proj_T(ipt)%R)

               call proj_LME_assoc(p_E,proj_T(ipt)%R(jEl))

               ! Re-create the projection electrode
               if ( p_E%idx > 0 ) then

                  El_p%Gamma => El_p%Sigma(:)
                  El_p%inDpvt%n =  p_E%ME%mol%pvt%n
                  El_p%inDpvt%r => p_E%ME%mol%pvt%r
                  call proj_Mt_mix(p_E%ME%mol,p_E%idx, El_p%Gamma, p_E%ME%bGk)

                  if ( N_eigen > 0 ) then
                     io = El_p%inDpvt%n
                     call A_Gamma_Block(zwork_tri,El_p,bTk(jEl,ipt), &
                          nGFGGF, GFGGF_work)
                     call TT_eigen(io,GFGGF_work, ntt_work, tt_work, eig)
                     do io = 1 , N_eigen
                        bTkeig(io,jEl,ipt) = dreal(eig(io))
                     end do
                  else
                     call A_Gamma(zwork_tri,El_p,bTk(jEl,ipt))
                  end if
                  
               else
                  
                  iEl = -p_E%idx
                  
                  if ( N_eigen > 0 ) then
                     io = Elecs(iEl)%inDpvt%n
                     call A_Gamma_Block(zwork_tri,Elecs(iEl),bTk(jEl,ipt), &
                          nGFGGF, GFGGF_work)
                     call TT_eigen(io,GFGGF_work, ntt_work, tt_work, eig)
                     do io = 1 , N_eigen
                        bTkeig(io,jEl,ipt) = dreal(eig(io))
                     end do
                  else
                     call A_Gamma(zwork_tri,Elecs(iEl),bTk(jEl,ipt))
                  end if
                  
               end if

            end do

          end do ! projections loop

          call timer('T-proj',2)

          ! Save the projections
          call proj_cdf_save(PROJcdf,N_Elec,Elecs, &
               ikpt,nE,N_proj_T,proj_T, &
               pDOS, bTk, N_eigen, bTkeig, save_DATA )

          end if
#endif

          ! For the very first iteration
          ! we print out an estimated time of arrival
          if ( iE == 0 .and. ikpt == 1 ) then

             ! Stop timer
             call timer_stop('E-loop')

             ! Calculate time passed
             call timer_get('E-loop',totTime=loop_time)
             loop_time = loop_time - init_time

             if ( IONode ) then
                iEl = itt_steps(Kp) * N_E
                loop_time = iEl * loop_time / real(Nodes,dp) - loop_time
                write(*,'(a,f20.3,'' s'')') 'tbt: Initial ETA in ', &
                     loop_time
             end if

             ! Start the timer again
             call timer_start('E-loop')

          end if

          ! step energy-point
          iE = iE + Nodes
          cE = tbt_E(iE+Nodes-Node,step=Nodes) ! we read them backwards

       end do

       ! Stop timer
       call timer_stop('E-loop')

    end do ! k-point

#ifdef MPI
    ! Force waiting till everything is done...
    call MPI_Barrier(MPI_Comm_World,iE)
#endif

    call itt_destroy(Kp)

#ifdef TBTRANS_DEBUG
    write(*,*) 'Completed TBTRANS SPIN'
#endif

    !***********************
    ! CLEAN UP
    !***********************
    if ( allocated(Teig) ) deallocate(Teig)

    deallocate(nE%iE,nE%E)

    nullify(DOS,DOS_El)
    deallocate(allDOS)

    deallocate(prep_El)
    deallocate(A_parts)

    call rgn_delete(pvt)

    call delete(zwork_tri)

    call delete(spH)
    call delete(spS)

#ifdef NCDF_4

    ! Close the netcdf file
    call ncdf_close(TBTcdf)
    if ( N_proj_ME > 0 ) then
       call ncdf_close(PROJcdf)
    end if

    call clean_dH( )

    if ( N_proj_ME > 0 ) then
       deallocate(El_p%Sigma)
       deallocate(bTk,pDOS)
       if ( allocated(bTkeig) ) deallocate(bTkeig)
    end if

    call delete(orb_J)

    ! Before we delete the Gf tri-diagonal matrix
    ! we need to create the sigma mean if requested.
    ! This is because %Sigma => Gfwork(:)
    call state_Sigma2mean(cdf_fname_sigma,N_Elec,Elecs)
#endif
    call delete(GF_tri)

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    call clear_BiasTriMat_inversion()
    call clear_TriMat_inversion()
    call clear_mat_inversion()

#ifdef NCDF_4
    ! Once we have cleaned up we can easily do the
    ! conversion of the TBT.nc file to the regular txt files
    ! We should have plenty of memory to do this.
    if ( ('proj-only'.nin.save_DATA).and.('Sigma-only'.nin.save_DATA) ) then

       ! We will guesstimate the current using the weights
       ! First we need to copy them over, we use S
       iE_N = N_tbt_E()
       nullify(S)
       allocate(S(iE_N))
       do iE = 1 , iE_N
          cE = tbt_E(iE)
          call c2weight(cE,S(iE))
       end do
       call state_cdf2ascii(cdf_fname,TSHS%nspin,ispin,N_Elec,Elecs, &
            iE_N,S,save_DATA)
       deallocate(S)
    end if
#else
    call end_save(iounits)
    call end_save(iounits_El)
    deallocate(iounits)
#endif

  contains

    subroutine print_memory(name,padding)
      use m_verbosity, only : verbosity
      use precision, only : i8b
      use m_ts_tri_common, only : nnzs_tri_i8b
      character(len=*), intent(in) :: name
      integer, intent(in) :: padding
      integer(i8b) :: nsize
      real(dp) :: mem
      if ( verbosity > 4 .and. IONode ) then
         nsize = nnzs_tri_i8b(DevTri%n,DevTri%r)
         nsize = nsize + padding
         mem = real(nsize,dp) * 16._dp / 1024._dp ** 2
         if ( mem > 600._dp ) then
            write(*,'(3a,i0,a,f8.3,a)') 'tbtrans: ',name,' Green function size / memory: ', &
                 nsize,' / ',mem / 1024._dp,' GB'
         else
            write(*,'(3a,i0,a,f8.2,a)') 'tbtrans: ',name,' Green function size / memory: ', &
                 nsize,' / ',mem,' MB'
         end if
      end if
    end subroutine print_memory

  end subroutine tbt_trik
  
  ! creation of the GF^{-1} for a certain region
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_invGF(cE, GFinv_tri, r, pvt, &
       N_Elec, Elecs, spH, spS)

    use class_Sparsity
    use class_zSpData1D
    use class_zTriMat
    use m_ts_cctype, only : ts_c_idx
    use m_tbt_tri_scat, only : insert_Self_Energy_Dev
    use intrinsic_missing, only : SFIND
#ifdef NCDF_4
    use m_tbt_dH, only : dH, add_zdH_TriMat
#endif

    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    type(zTriMat), intent(inout) :: GFinv_tri
    type(tRgn), intent(in) :: r, pvt
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D), intent(inout) :: spH,  spS

    ! Local variables
    complex(dp) :: Z
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: H(:), S(:)
    complex(dp), pointer :: Gfinv(:)
    integer :: io, iu, ind, idx, ju

    if ( cE%fake ) return

    Z = cE%e

    sp => spar(spH)
    H  => val (spH)
    S  => val (spS)

    call attach(sp, n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    Gfinv => val(Gfinv_tri)

!$OMP parallel default(shared), private(iu,io,ind,ju,idx)

    ! Initialize
!$OMP workshare
    GFinv(:) = dcmplx(0._dp,0._dp)
!$OMP end workshare

    ! We will only loop in the central region
    ! We have constructed the sparse array to only contain
    ! values in this part...
!$OMP do 
    do iu = 1, r%n
       io = r%r(iu) ! get the orbital in the big sparsity pattern
       if ( l_ncol(io) /= 0 ) then
          
       ! Loop over non-zero entries here
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          ju = pvt%r(l_col(ind))
          ! If it is zero, then *must* be electrode
          ! or fold down region.
          if ( ju == 0 ) cycle
          
          ! Notice that we transpose back here...
          ! See symmetrize_HS_kpt
          idx = index(Gfinv_tri,ju,iu)
          
          GFinv(idx) = Z * S(ind) - H(ind)
          
       end do
             
       end if
    end do
!$OMP end do nowait

!$OMP end parallel

    do io = 1 , N_Elec
       call insert_self_energy_dev(Gfinv_tri,Gfinv,r,Elecs(io))
    end do

#ifdef NCDF_4
    if ( dH%lvl > 0 ) then
       ! Add dH
       call add_zdH_TriMat(dH%dH, Gfinv_tri, r)
    end if
#endif

  end subroutine prepare_invGF

#ifdef NOT_USED
  ! There are two variants of populating the
  ! Green function for each energy point.
  ! If you have a very small device region it is
  ! beneficial to search the sparsity pattern for the
  ! indices.
  ! On the other hand for a very large device region
  ! it is better to loop the sparsity pattern and search
  ! the device region.
  ! This routine decides which is better by doing a
  ! one-time loop for comparison
  ! * NOT USED
  ! This was used in the old setup.
  ! When one wants to reduce memory requirement
  ! this can be reinstantiated.
  ! However, then we are talking about more
  ! than 2,000,000 elements...
  subroutine prep_Gfinv_algo(dev_tri,sp,r,loop_dev)

    use class_Sparsity
    use class_zTriMat
    use intrinsic_missing, only : SFIND

    type(zTriMat), intent(inout) :: dev_tri
    type(Sparsity), intent(inout) :: sp
    type(tRgn), intent(in) :: r
    logical, intent(out) :: loop_dev

    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: io, iu, ind
    integer :: itot
#ifdef _OPENMP
    real(dp) :: t_d, t_s
#else
    real :: t_d, t_s, tt
#endif

    call attach(sp, n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    ! search some middle orbital
    io = r%r(r%n/2)

    ! Time: loop_dev
    itot = 0
#ifdef _OPENMP
    t_d = omp_get_wtime( )
#else
    call cpu_time( tt )
    t_d = tt
#endif
    do iu = 1 , r%n
          
       ind = SFIND(l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io)),r%r(iu))
       if ( ind == 0 ) cycle
       ind = l_ptr(io) + ind

       ! This is only to dis-allow the optimizer from removing
       ! this entire loop
       itot = itot + ind

    end do
#ifdef _OPENMP
    t_d = omp_get_wtime( ) - t_d
#else
    call cpu_time( tt )
    t_d = tt - t_d
#endif
    
    itot = 0
    ! Time: loop_sparsity
#ifdef _OPENMP
    t_s = omp_get_wtime( )
#else
    call cpu_time( tt )
    t_s = tt
#endif
    do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
       
       iu = rgn_pivot(r,l_col(ind))
       ! Important to check this as we cannot ensure
       ! the entire sparsity pattern in this
       ! tri-diagonal matrix
       if ( iu == 0 ) cycle
       itot = itot + iu
       
    end do
#ifdef _OPENMP
    t_s = omp_get_wtime( ) - t_s
#else
    call cpu_time( tt )
    t_s = tt - t_s
#endif

    ! if the timing for the searching of
    ! sparsity is slower than looping the
    ! device, we prefer to loop the device
    loop_dev = t_s > t_d

  end subroutine prep_Gfinv_algo
#endif

  subroutine downfold_SE(cE, El, spH, spS,r,np,p,nwork,work)

    use class_Sparsity
    use class_zSpData1D
    use m_ts_electype
    use m_ts_cctype, only : ts_c_idx

    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    ! Electrode, this *REQUIRES* that the down-folding region
    ! only contains the self-energies of this electrode... :(
    type(Elec), intent(inout) :: El
    ! The hamiltonian and overlap
    type(zSpData1D), intent(in) :: spH, spS
    ! The region of downfolding... (+ the region connecting to the device..)
    type(tRgn), intent(in) :: r
    ! number of parts that constitute the tri-diagonal region
    integer, intent(in) :: np
    ! parts associated
    integer, intent(in) :: p(np)
    ! Work-arrays...
    integer, intent(in) :: nwork
    complex(dp), intent(inout), target :: work(nwork)

    ! All our work-arrays...
    complex(dp), pointer :: A(:), B(:), C(:), Y(:)

    integer :: no, off, i, ii, j, ierr
    integer :: ip, itmp

    if ( cE%fake ) return

    ip = maxval(p)

    ! Copy down the downfolded size...
    no = El%o_inD%n

    if ( p(np) /= no ) then
       call die('Something went wrong... The last segment MUST be &
            &equivalent to the down-folded region. No more, no less.')
    end if

    ! Check that there is space enough in the work array.
    do ip = 1 , np - 1
       itmp = p(ip) * p(ip+1) * 2 ! B,C
       itmp = itmp  + p(ip) ** 2   ! A
       itmp = itmp  + p(ip+1) ** 2 ! Y
       if ( itmp > nwork ) then
          call die('Work array is too small... A, B, C, Y')
       end if
    end do

    ! check that the last segment also holds...
    itmp = p(np-1) ** 2 ! A
    itmp = itmp + p(np) ** 2 ! Y
    if ( itmp > nwork ) then ! we do not need y here...
       call die('Work array is too small..., A, Y-next')
    end if
    ! We start by pointing the Y array to the far back of the
    ! work-array. In that way we ensure that no elements overlap
    ! and hence we can re-use that array directly
    
    ! loop and convert..
    off = 0
    do ip = 2 , np 

       ! Set up pointers
       i = 1
       A => work(i:i-1+p(ip-1)**2)
       i = i + p(ip-1)**2
       B => work(i:i-1+p(ip)*p(ip-1))
       i = i + p(ip)*p(ip-1)
       C => work(i:i-1+p(ip)*p(ip-1))

       call prep_HS(cE%E,El,spH,spS,r,off,p(ip-1),off,p(ip-1),A)

       if ( ip > 2 ) then
!$OMP parallel workshare default(shared)
         A(:) = A(:) - Y(:)
!$OMP end parallel workshare
       end if

       ! Ensures that Y is not overwritten
       i = off + p(ip-1)
       call prep_HS(cE%E,El,spH,spS,r,i,p(ip),off,p(ip-1),B)
       call prep_HS(cE%E,El,spH,spS,r,off,p(ip-1),i,p(ip),C)

       ! increment offset
       off = off + p(ip-1)

       ! re-point Y
       if ( ip == np ) then
          ! Sigma should have been emptied by the previous loops :)
          Y => El%Sigma(:)
          if ( no /= p(np) ) call die('must be enforced')
       else
          Y => work(nwork-p(ip)**2+1:nwork)
       end if

       ! Calculate: [An-1 - Yn-1] ^-1 Cn
       call zgesv(p(ip-1),p(ip),A,p(ip-1),ipiv,C,p(ip-1),ierr)
       if ( ierr /= 0 ) then
          write(*,'(a,i0)') 'Inversion of down-projection failed: ',ierr
       end if
       
       ! Calculate: Bn-1 [An-1 - Yn-1] ^-1 Cn
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',p(ip),p(ip),p(ip-1),dcmplx(1._dp,0._dp), &
            B,p(ip),C,p(ip-1),dcmplx(0._dp,0._dp),Y,p(ip))
       
    end do

    ! At this point we should be left with the last segment
    ! which is the self-energy projected into the device region...
    if ( r%n /= off + p(np) ) then
       print *,r%n,off+p(np)
       call die('Error in regional size, should not be encountered')
    end if
        
    ! Create Gamma...
    ! \Gamma ^ T = (\Sigma - \Sigma^\dagger)^T
!$OMP parallel do default(shared), private(j,i,ii,ip)
    do j = 1 , no
       ii = no * ( j - 1 )
       ip = j - no
       do i = 1 , j - 1
          ii = ii + 1
          ip = ip + no
          El%Gamma(ii) = El%Sigma(ip) - dconjg( El%Sigma(ii) )
          El%Gamma(ip) = El%Sigma(ii) - dconjg( El%Sigma(ip) )
       end do
       ii = no*(j-1) + j
       El%Gamma(ii) = El%Sigma(ii) - dconjg( El%Sigma(ii) )
    end do
!$OMP end parallel do

    ! aaaaannnnnD DONE!

  end subroutine downfold_SE

#ifdef NOT_WORKING

  subroutine downfold_SE_life(cE, El, spH, spS,r,np,p,life,nwork,work)

    use class_Sparsity
    use class_zSpData1D
    use m_ts_electype
    use m_ts_cctype, only : ts_c_idx

    use m_mat_invert

    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    ! Electrode, this *REQUIRES* that the down-folding region
    ! only contains the self-energies of this electrode... :(
    type(Elec), intent(inout) :: El
    ! The hamiltonian and overlap
    type(zSpData1D), intent(in) :: spH, spS
    ! The region of downfolding... (+ the region connecting to the device..)
    type(tRgn), intent(in) :: r
    ! number of parts that constitute the tri-diagonal region
    integer, intent(in) :: np
    ! parts associated
    integer, intent(in) :: p(np)
    ! The lifetime on the respective orbitals, all the way down
    real(dp), intent(out) :: life(r%n-El%o_inD%n)
    ! Work-arrays...
    integer, intent(in) :: nwork
    complex(dp), intent(inout), target :: work(nwork)

    ! All our work-arrays...
    complex(dp), pointer :: A(:), B(:), C(:), Y(:)

    integer :: no, off, i, ii, j, jj, ierr, o_life
    integer :: ip, itmp

    complex(dp), external :: zdotu

    if ( cE%fake ) return

    ip = maxval(p)

    ! Copy down the downfolded size...
    no = El%o_inD%n

    if ( p(np) /= no ) then
       call die('Something went wrong... The last segment MUST be &
            &equivalent to the down-folded region. No more, no less.')
    end if

    ! Check that there is space enough in the work array.
    do ip = 1 , np - 1
       itmp = p(ip) * p(ip+1) * 2 ! B,C
       itmp = itmp  + p(ip) ** 2   ! A
       itmp = itmp  + p(ip+1) ** 2 ! Y
       if ( itmp > nwork ) then
          call die('Work array is too small... A, B, C, Y')
       end if
    end do

    ! check that the last segment also holds...
    itmp = p(np-1) ** 2 ! A
    itmp = itmp + p(np) ** 2 ! Y
    if ( itmp > nwork ) then ! we do not need y here...
       call die('Work array is too small..., A, Y-next')
    end if
    ! We start by pointing the Y array to the far back of the
    ! work-array. In that way we ensure that no elements overlap
    ! and hence we can re-use that array directly
    
    ! loop and convert..
    off = 0
    do ip = 2 , np 

       ! Set up pointers
       i = 1
       A => work(i:i-1+p(ip-1)**2)

       call prep_HS(cE%E,El,spH,spS,r,off,p(ip-1),off,p(ip-1),A)

       ! Now we can calculate the life-time of the electrons
       ! The life-time is the g_ii (local Green function)
       !   Tr[Gf_ii^0\Gamma Gf_ii^0^\dagger \Gamma]
       ! First calculate Gf, but we retain A
       
       ierr = p(ip-1)
       if ( size(El%Sigma) >= ierr**2 ) then
          C => El%Sigma
       else
          call die('Sigma array too small')
       end if
       
       i = 1 + ierr**2
       B => work(i:i-1+ierr**2)

       if ( ip > 2 ) then
!$OMP parallel workshare default(shared)
          A(:) = A(:) - Y(:)
!$OMP end parallel workshare
       else
          Y => El%Sigma
       end if
       
       ! Calculate Gf^0
       call mat_invert(A,B,ierr,MI_IN_PLACE_LAPACK)
       
       ! Calculate the Gamma function
       ! \Gamma ^ T = (\Sigma - \Sigma^\dagger)^T
!$OMP parallel do default(shared), private(j,i,ii,jj)
       do j = 1 , ierr
          ii = ierr * ( j - 1 )
          jj = j - ierr
          do i = 1 , j - 1
             ii = ii + 1
             jj = jj + ierr
             C(ii) = Y(jj) - dconjg( Y(ii) )
             C(jj) = Y(ii) - dconjg( Y(jj) )
          end do
          ii = ierr*(j-1) + j
          C(ii) = Y(ii) - dconjg( Y(ii) )
       end do
!$OMP end parallel do

       ! Calculate matrix product
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','T',ierr,ierr,ierr,dcmplx(1._dp,0._dp), &
            A,ierr,C,ierr,dcmplx(0._dp,0._dp),B,ierr)
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','C',ierr,ierr,ierr,dcmplx(1._dp,0._dp), &
            B,ierr,A,ierr,dcmplx(0._dp,0._dp),El%Sigma,ierr)

       ! Calculate trace
       do i = 1 , ierr
          ! Calculate offset
          j = (i-1) * ierr + 1
          life(off+i) = - real( zdotu(ierr,El%Sigma(j),1,C(j),1) ,dp)
       end do

       ! Reassign pointers
       i = 1 + p(ip-1)**2
       B => work(i:i-1+p(ip)*p(ip-1))
       i = i + p(ip)*p(ip-1)
       C => work(i:i-1+p(ip)*p(ip-1))

       ! Ensures that Y is not overwritten
       i = off + p(ip-1)
       call prep_HS(cE%E,El,spH,spS,r,off,p(ip-1),i,p(ip),B)

       ! increment offset
       off = off + p(ip-1)

       ! re-point Y
       if ( ip == np ) then
          ! Sigma should have been emptied by the previous loops :)
          Y => El%Sigma(:)
          if ( no /= p(np) ) call die('must be enforced')
       else
          Y => work(nwork-p(ip)**2+1:nwork)
       end if

       ! Calculate: [An-1 - Yn-1] ^-1 Cn
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',p(ip-1),p(ip),p(ip-1),dcmplx(1._dp,0._dp), &
            A,p(ip-1),B,p(ip-1),dcmplx(0._dp,0._dp),C,p(ip))

       call prep_HS(cE%E,El,spH,spS,r,i,p(ip),off,p(ip-1),B)

       ! Calculate: Bn-1 [An-1 - Yn-1] ^-1 Cn
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',p(ip),p(ip),p(ip-1),dcmplx(1._dp,0._dp), &
            B,p(ip),C,p(ip-1),dcmplx(0._dp,0._dp),Y,p(ip))

    end do

    ! At this point we should be left with the last segment
    ! which is the self-energy projected into the device region...
    if ( r%n /= off + p(np) ) then
       print *,r%n,off+p(np)
       call die('Error in regional size, should not be encountered')
    end if
        
    ! Create Gamma...
    ! \Gamma ^ T = (\Sigma - \Sigma^\dagger)^T
!$OMP parallel do default(shared), private(j,i,ii,ip)
    do j = 1 , no
       ii = no * ( j - 1 )
       ip = j - no
       do i = 1 , j - 1
          ii = ii + 1
          ip = ip + no
          El%Gamma(ii) = El%Sigma(ip) - dconjg( El%Sigma(ii) )
          El%Gamma(ip) = El%Sigma(ii) - dconjg( El%Sigma(ip) )
       end do
       ii = no*(j-1) + j
       El%Gamma(ii) = El%Sigma(ii) - dconjg( El%Sigma(ii) )
    end do
!$OMP end parallel do

    ! aaaaannnnnD DONE!

  end subroutine downfold_SE_life

#endif

  
  ! creation of the GF^{-1} for a certain region
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prep_HS(Z, El, spH, spS, &
       r,off1,n1,off2,n2,M)

    use class_Sparsity
    use class_zSpData1D
    use m_tbt_tri_scat, only : insert_Self_Energy
    use intrinsic_missing, only : SFIND

#ifdef NCDF_4
    use m_tbt_dH, only : dH, add_zdH_Mat
#endif

    ! the current energy point
    complex(dp), intent(in) :: Z
    ! Electrodes...
    type(Elec), intent(inout) :: El
    ! The Hamiltonian and overlap sparse matrices
    type(zSpData1D), intent(in) :: spH,  spS
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    ! The sizes and offsets of the matrix
    integer, intent(in) :: off1, n1, off2, n2
    complex(dp), intent(out) :: M(n1,n2)

    ! Local variables
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: H(:), S(:)
    integer :: io, iu, ind, ju

    sp => spar(spH)
    H  => val (spH)
    S  => val (spS)

    call attach(sp, n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    ! We will only loop in the region
!$OMP parallel do default(shared), private(iu,io,ind,ju)
    do iu = 1 , n2
       io = r%r(off2+iu) ! get the orbital in the sparsity pattern

       ! Initialize to zero
       M(:,iu) = dcmplx(0._dp,0._dp)

       if ( l_ncol(io) /= 0 ) then

       ! Loop on entries here...
       do ju = 1 , n1
 
          ! Check if the orbital exists in the region
          ! We are dealing with a UC sparsity pattern.
          ind = SFIND(l_col(l_ptr(io)+1:l_ptr(io) + l_ncol(io)),r%r(off1+ju))
          if ( ind == 0 ) cycle
          ind = l_ptr(io) + ind

          ! Notice that we transpose back here...
          ! See symmetrize_HS_kpt
          M(ju,iu) = Z * S(ind) - H(ind)

       end do
       end if
    end do
!$OMP end parallel do

    call insert_Self_Energy(n1,n2,M,r,El,off1,off2)

#ifdef NCDF_4
    if ( dH%lvl > 0 ) then
       ! Add dH
       call add_zdH_Mat(dH%dH, r, off1,n1,off2,n2, M)
    end if
#endif

  end subroutine prep_HS

end module m_tbt_trik
