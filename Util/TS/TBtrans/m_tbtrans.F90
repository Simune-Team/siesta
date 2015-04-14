! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_tbtrans

  use precision, only : dp

  use m_tbt_hs

  use m_tbt_tri_init, only : tbt_tri_init, tbt_tri_print_opti

  use m_tbt_trik

  implicit none

  public :: tbt
  private

contains

  subroutine tbt(TSHS, kT)

#ifdef MPI
    use mpi_siesta, only : MPI_Barrier, MPI_Comm_World
#endif
    use units, only : eV
    use alloc, only : re_alloc, de_alloc

    use fdf, only: fdf_get

    use parallel, only : IONode

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use m_tbt_contour

    use m_ts_electype

    use m_tbt_options, only : N_Elec, Elecs
#ifdef NCDF_4
    use m_tbt_kpoint, only : nkpnt, kpoint, kweight
    use m_tbt_options, only : save_DATA
    use m_tbt_options, only : cdf_fname, cdf_fname_sigma, cdf_fname_proj
    use m_tbt_regions, only : r_aDev, r_aBuf, r_oDev, sp_dev

    use m_tbt_save
    use m_tbt_proj, only : N_mol, mols, init_proj_save
    use m_tbt_sigma_save, only : init_Sigma_save
#else
    use m_tbt_kpoint, only : nkpnt
#endif
    use m_tbt_kregions, only : n_k, r_k, kregion_step, kregion_k

    use m_ts_gf, only : read_Green

    use dictionary

! ********************
! * INPUT variables  *
! ********************
    type(tTSHS), intent(inout) :: TSHS
    real(dp), intent(in) :: kT

! ******************** IO descriptors ************************
    integer, allocatable :: uGF(:)
! ************************************************************

! ****************** Electrode variables *********************
    integer, allocatable :: nq(:)
! ************************************************************
    ! Temporary variables
    integer :: nkpt
    real(dp), pointer :: kpt(:,:), wkpt(:)
    real(dp) :: k(3)

! * local variables
    integer :: iEl, NEn, no_used, no_used2, ispin, ils, i

    ! Total number of energy-points...
    NEn = N_TBT_E()

#ifdef NCDF_4
    ! Initialize the tri-diagonal matrices!
    if ( N_mol > 0 ) then
       call tbt_tri_init( TSHS%dit, TSHS%sp, mols(:)%orb )
    else
       call tbt_tri_init( TSHS%dit, TSHS%sp )
    end if
#else
    call tbt_tri_init( TSHS%dit, TSHS%sp )
#endif

    ! Suggest to the user an optimal device region for
    ! fastest calculation
    call tbt_tri_print_opti(TSHS%na_u,TSHS%lasto,r_oDev)

    if ( fdf_get('TBT.Analyze',.false.) ) then
#ifdef MPI
       call MPI_Barrier(MPI_Comm_World,i)
#endif
       call die('Stopping TBtrans on purpose after analyzation step...')
    end if

    ! Open GF files...
    ! Read-in header of Green's functions
    ! Prepare for the calculation
    ! We read in the k-points that the electrode was generated with.
    ! Furthermore we read in the expansion q-points
    ! They are communicated in the routine

    ! We need to initialize tbtrans
    if ( IONode ) then
       write(*,*) ! new-line
    end if

    call timer('TBT',1)

    ! in case the file-descriptor is negative it basically 
    ! means "out-of-core" calculation.
    allocate(uGF(N_Elec),nq(N_Elec))
    uGF(:) = -1
    do iEl = 1 , N_Elec

       nq(iEl) = product(Elecs(iEl)%Rep)

       ! Allocate the electrode quantities
       nullify(Elecs(iEl)%HA,Elecs(iEl)%SA,Elecs(iEl)%Gamma)

       ! We allocate for once as much space as needed,

       ! Allocate the non-repeated hamiltonian and overlaps...
       no_used = Elecs(iEl)%no_used
       if ( Elecs(iEl)%pre_expand > 1 ) then ! > 1 also expand H, S before writing
          no_used = TotUsedOrbs(Elecs(iEl))
          nq(iEl) = 1
       end if

       ! If we using bulk electrodes, we need not the Hamiltonian, 
       ! nor the overlap...
       if ( .not. Elecs(iEl)%Bulk ) then
          call re_alloc(Elecs(iEl)%HA,1,no_used,1,no_used,1,nq(iEl),routine='tbtrans')
          call re_alloc(Elecs(iEl)%SA,1,no_used,1,no_used,1,nq(iEl),routine='tbtrans')
       end if

       ! We need to retain the maximum size of Gamma
       ! I.e. take into account that the down-projected region
       ! could be larger than the read in self-energy
       no_used = TotUsedOrbs(Elecs(iEl))
       no_used2 = no_used
       if ( Elecs(iEl)%pre_expand == 0 ) then
          no_used2 = Elecs(iEl)%no_used
       end if
       ! We need to assure that the entire down-folded Gamma
       ! can be saved
       if ( no_used * no_used2 < Elecs(iEl)%o_inD%n ** 2 ) then
          no_used2 = Elecs(iEl)%o_inD%n ** 2 / no_used + 1
       end if

       call re_alloc(Elecs(iEl)%Gamma,1,no_used*no_used2,routine='tbtrans')

       ! This seems stupid, however, we never use the expansion array and
       ! GammaT at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called:
       ! first the GAA is "emptied" of information and then
       ! Gamma is filled.
       no_used2 = no_used
       if ( Elecs(iEl)%pre_expand == 0 ) no_used2 = Elecs(iEl)%no_used
       Elecs(iEl)%GA => Elecs(iEl)%Gamma(1:no_used*no_used2)

    end do

    call open_GF(N_Elec,Elecs,uGF,nkpnt,NEn,spin_idx)
    if ( n_k > 0 ) then
       ! We initialize the region that belongs to
       ! each electrode
       do iEl = 1 , N_Elec
          uGF(iEl) = -1
          do ils = 0 , n_k
             if ( in_rgn(r_k(ils)%atm,Elecs(iEl)%idx_a) ) then
                uGF(iEl) = ils
                exit
             end if
          end do
          if ( uGF(iEl) == -1 ) then
             call die('Error in setting up different k-regions, &
                  &all electrodes are not fully contained.')
          end if
       end do
    end if
    
    do ils = 1 , TSHS%nspin

       ! Determine the actual spin-index
       if ( spin_idx == 0 ) then
          ispin = ils
       else
          ispin = spin_idx
       end if

       ! The initial spin has already been 
       ! setup for the first spin, hence we
       ! only need to re-read them for 
       ! following spin calculations
       if ( ils > 1 ) then
          call prep_next_HS(ispin,Volt)
          do iEl = 1 , N_Elec
             ! Re-read in the electrode 
             ! Hamiltonian and build the H/S
             ! for this spin.
             if ( Elecs(iEl)%out_of_core ) cycle
             call init_Electrode_HS(Elecs(iEl),ispin)
          end do
       end if

#ifdef NCDF_4
       if ( n_k == 0 ) then
          nkpt =  nkpnt
          kpt  => kpoint(:,:)
          wkpt => kweight(:)
       else
          nullify(kpt,wkpt)
          nkpt = 1
          do i = 0 , n_k
             nkpt = nkpt * size(r_k(i)%wkpt)
          end do
          allocate(kpt(3,nkpt),wkpt(nkpt))
          r_k(:)%ik   = 1
          r_k(n_k)%ik = 0
          do i = 1 , nkpt
             call kregion_step( )
             call kregion_k(-1, k, w = wkpt(i) )
             call kpoint_convert(TSHS%cell,k,kpt(:,i),1)
          end do
       end if

       ! If the user has requested only to calculate
       ! the projections, we do not initialize the cdf_fname
       if ( ('proj-only'.nin.save_DATA).and.('Sigma-only'.nin.save_DATA) ) then

          ! Initialize data files
          call name_save( ispin, TSHS%nspin,cdf_fname, end = 'nc')
          call init_cdf_save(cdf_fname,TSHS,r_oDev,ispin,N_Elec, Elecs, &
               nkpt, kpt, wkpt, NEn, r_aDev, r_aBuf, sp_dev, save_DATA )
       end if
       
       call name_save( ispin, TSHS%nspin,cdf_fname_sigma, end = 'Sigma.nc')
       call init_Sigma_save(cdf_fname_sigma,TSHS,r_oDev,ispin,N_Elec, Elecs, &
            nkpt, kpt, wkpt, NEn, r_aDev, r_aBuf )

       if ( ('Sigma-only'.nin.save_DATA) ) then
       
          call name_save( ispin, TSHS%nspin, cdf_fname_proj, end = 'Proj.nc' )
          call init_Proj_save( cdf_fname_proj, TSHS , r_oDev, ispin, N_Elec, Elecs, &
               nkpt, kpt, wkpt, NEn , r_aDev, r_aBuf, sp_dev, save_DATA )
       end if

       if ( n_k /= 0 ) then
          deallocate(kpt,wkpt)
          nullify(kpt,wkpt)
       end if
#endif

       call tbt_trik(ispin,N_Elec, Elecs, TSHS, nq, uGF)

       ! the spin-index is zero for all,
       ! and one of the allowed spin indices if
       ! a specific one is requested...
       if ( ispin == spin_idx ) exit

    end do

    ! Close files
    do iEl = 1 , N_Elec
       if ( IONode .and. Elecs(iEl)%out_of_core ) then
          call io_close(uGF(iEl))
       end if
    end do
       
    !***********************
    !       Clean up
    !***********************
    do iEl = 1 , N_Elec
       if ( .not. Elecs(iEl)%out_of_core ) then
          call delete(Elecs(iEl))
       end if
    end do

    !***********************
    !  Clean up electrodes
    !***********************
    do iEl = 1 , N_Elec
       if ( associated(Elecs(iEl)%HA) ) then
          call de_alloc(Elecs(iEl)%HA,routine='transiesta')
          call de_alloc(Elecs(iEl)%SA,routine='transiesta')
       end if
       call de_alloc(Elecs(iEl)%Gamma,routine='transiesta')
    end do

    deallocate(uGF,nq)

    call timer('TBT',2)

  contains

    subroutine init_Electrode_HS(El,spin_idx)
      use class_Sparsity
      use class_dSpData1D
      use class_dSpData2D
      use alloc, only : re_alloc
      type(Elec), intent(inout) :: El
      integer, intent(in) :: spin_idx

      ! If already initialized, return immediately
      if ( initialized(El%sp) ) return

      ! Read-in and create the corresponding transfer-matrices
      call delete(El) ! ensure clean electrode
      call read_Elec(El,Bcast=.true., IO = .false., ispin = spin_idx )
      
      if ( .not. associated(El%isc_off) ) then
         call die('An electrode file needs to be a non-Gamma calculation. &
              &Ensure at least two k-points in the T-direction.')
      end if
      
      call create_sp2sp01(El, IO = .false.)

      ! Clean-up, we will not need these!
      ! we should not be very memory hungry now, but just in case...
      call delete(El%H)
      call delete(El%S)
      
      ! We do not accept onlyS files
      if ( .not. initialized(El%H00) ) then
         call die('An electrode file must contain the Hamiltonian')
      end if

      call delete(El%sp)

    end subroutine init_Electrode_HS

    subroutine open_GF(N_Elec,Elecs,uGF,nkpnt,NEn,spin_idx)
      integer, intent(in) :: N_Elec
      type(Elec), intent(inout) :: Elecs(N_Elec)
      integer, intent(out) :: uGF(N_Elec)
      integer, intent(in) :: nkpnt, NEn, spin_idx

      ! Local variables
      integer :: iEl
      
      do iEl = 1 , N_Elec

         ! Initialize k-points
         Elecs(iEl)%bkpt_cur(:) = huge(1._dp)

         if ( Elecs(iEl)%out_of_core ) then
            
            if ( IONode ) then
               call io_assign(uGF(iEl))
               open(file=Elecs(iEl)%GFfile,unit=uGF(iEl),form='unformatted')
            end if
            
            call read_Green(uGF(iEl),Elecs(iEl), nkpnt, NEn )

!            if ( IONode .and. spin_idx > 1 ) then
!               ! In case the user has requested to only use
!               ! one of the spin-channels we step forward to that one
!               do is = 1 , spin_idx - 1
!                  ! Skip all H and S arrays
!                  do i = 1 , nkpnt
!                     read(uGF(iEl)) ! H
!                     read(uGF(iEl)) ! S
!                  end do
!                  ! Skip all header lines AND GS lines
!                  do i = 1 , nkpnt * NEn
!                     read(uGF(iEl)) ! Header line
!                     read(uGF(iEl)) ! GS
!                  end do
!               end do
!            end if

         else
            
            ! prepare the electrode to create the surface self-energy
            call init_Electrode_HS(Elecs(iEl),max(1,spin_idx))
            
         end if

      end do
      
    end subroutine open_GF

  end subroutine tbt

end module m_tbtrans

