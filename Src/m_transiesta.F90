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

module m_transiesta

  use precision, only : dp

  use files, only : slabel

  use m_ts_sparse, only : ts_sparse_init
  use m_ts_method, only : ts_method, TS_SPARSITY, TS_SPARSITY_TRI

  use m_ts_tri_init, only : ts_tri_init

  use m_ts_fullg
  use m_ts_fullk

  use m_ts_trig
  use m_ts_trik

  implicit none

  public :: transiesta
  private

contains

  subroutine transiesta(TSiscf,nspin, &
       sp_dist, sparse_pattern, &
       Gamma, ucell, no_u, na_u, lasto, xa, n_nzs, &
       xij, H, S, DM, EDM, Ef, kT, &
       Qtot)

    use alloc, only : re_alloc, de_alloc

    use parallel, only : IONode

    use class_OrbitalDistribution
    use class_Sparsity

    use m_ts_kpoints, only : ts_nkpnt

    use m_ts_electype

    use m_ts_options, only : N_Elec, Elecs
    use m_ts_options, only : IsVolt, Calc_Forces
    use m_ts_options, only : no_BufL, no_BufR

    use m_ts_contour_eq , only : N_Eq_E
    use m_ts_contour_neq, only : N_nEq_E

    use m_ts_charge
    
    use m_ts_gf, only : read_Green

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in)  :: TSiscf
    integer, intent(in)  :: nspin
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    logical, intent(in)  :: Gamma
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in)  :: no_u, na_u
    integer, intent(in)  :: lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: xij(3,n_nzs)
    real(dp), intent(in) :: H(n_nzs,nspin), S(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Ef, kT, Qtot

! ******************** IO descriptors ************************
    integer, allocatable :: uGF(:)
! ************************************************************

! ****************** Electrode variables *********************
    integer, allocatable :: nq(:)
! ************************************************************

! * local variables
    integer :: iEl, no, NEn, no_used, no_used2

    ! Open GF files...
    ! Read-in header of Green's functions
    ! Prepare for the calculation
    ! We read in the k-points that the electrode was generated with.
    ! Furthermore we read in the expansion q-points
    ! They are communicated in the routine

    if ( TSiscf == 1 ) then
       ! We need to initialize TRANSIESTA

       call timer('TS_init',1)

       call ts_sparse_init(slabel,Gamma, sp_dist, sparse_pattern, &
            na_u, lasto)

       if ( ts_method == TS_SPARSITY_TRI ) then
          ! initialize the tri-diagonal partition
          call ts_tri_init()
       end if

       ! print out estimated memory usage...
       call ts_print_memory(Gamma)

       call ts_print_charges(Elecs, sp_dist, sparse_pattern, &
            nspin, n_nzs, DM, S)

       if ( .not. Calc_Forces .and. IONode ) then
          write(*,'(a)') 'transiesta: Notice that the forces are NOT updated'
       end if

       call timer('TS_init',2)

    end if


    call timer('TS',1)

    ! Total number of energy-points...
    NEn = N_Eq_E() + N_nEq_E()
    
    allocate(uGF(N_Elec))
    allocate(nq(N_Elec))
    do iEl = 1 , N_Elec
       if ( IONode ) then
          call io_assign(uGF(iEl))
          open(file=GFFile(Elecs(iEl)),unit=uGF(iEl),form='unformatted')
       end if
       call read_Green(uGF(iEl),Elecs(iEl), ts_nkpnt, NEn, .false. )
       nq(iEl) = Rep(Elecs(iEl))


       ! Allocate the electrode quantities
       nullify(Elecs(iEl)%HA,Elecs(iEl)%SA,Elecs(iEl)%Gamma)

       ! We allocate for once as much space as needed,

       ! Allocate the non-repeated hamiltonian and overlaps...
       no_used = UsedOrbs(Elecs(iEl))
       call re_alloc(Elecs(iEl)%HA,1,no_used,1,no_used,1,nq(iEl),routine='transiesta')
       call re_alloc(Elecs(iEl)%SA,1,no_used,1,no_used,1,nq(iEl),routine='transiesta')

       no_used = TotUsedOrbs(Elecs(iEl))
       if ( IsVolt ) then
          ! We need Gamma's with voltages (now they are both GAA and GammaT)
          no_used2 = no_used
       else
          ! This is only for having space for GA
          no_used2 = UsedOrbs(Elecs(iEl))
       end if
       call re_alloc(Elecs(iEl)%Gamma,1,no_used,1,no_used2,routine='transiesta')

       ! This seems stupid, however, we never use the expansion array and
       ! GammaT at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called:
       ! first the GAA is "emptied" of information and then
       ! Gamma is filled.
       no_used2 = UsedOrbs(Elecs(iEl))
       Elecs(iEl)%GA => Elecs(iEl)%Gamma(1:no_used,1:no_used2)

    end do


    if ( ts_method == TS_SPARSITY ) then
       if ( Gamma ) then ! we can't even do this in TS_Gamma it would result in errorneous connections
          call ts_fullg(N_Elec,Elecs, &
               nq,uGF, &
               nspin, &
               sp_dist, sparse_pattern, &
               no_u, n_nzs, &
               H, S, DM, EDM, Ef, kT)
       else
          call ts_fullk(N_Elec,Elecs, &
               nq,uGF, &
               nspin, &
               sp_dist, sparse_pattern, &
               no_u, n_nzs, &
               H, S, xij, DM, EDM, Ef, kT)
       end if
    else if ( ts_method == TS_SPARSITY_TRI ) then
       if ( Gamma ) then
          call ts_trig(N_Elec,Elecs, &
               nq, uGF, nspin, &
               sp_dist, sparse_pattern, &
               no_u, n_nzs, &
               H, S, DM, EDM, Ef, kT)
       else
          call ts_trik(N_Elec,Elecs, &
               nq,uGF, &
               nspin, &
               sp_dist, sparse_pattern, &
               no_u, n_nzs, &
               H, S, xij, DM, EDM, Ef, kT)
       end if
    else

       call die('Error in code')
    end if

!***********************
!     Close Files
!***********************
    if ( IONode ) then
       do iEl = 1 , N_Elec
          call io_close(uGF(iEl))
       end do
    end if
    
!***********************
!     Clean up electrodes
!***********************
    do iEl = 1 , N_Elec
       call de_alloc(Elecs(iEl)%HA,routine='transiesta')
       call de_alloc(Elecs(iEl)%SA,routine='transiesta')
       call de_alloc(Elecs(iEl)%Gamma,routine='transiesta')
    end do

    deallocate(uGF,nq)

    ! We do the charge correction of the transiesta
    ! computation here (notice that the routine will automatically
    ! return if no charge-correction is requested)
    call ts_charge_correct(no_BufL+no_BufR, Elecs, sp_dist, &
         sparse_pattern, nspin, n_nzs, DM, EDM, S, Qtot, &
         TS_RHOCORR_METHOD)

    call ts_print_charges(Elecs, sp_dist, sparse_pattern, &
         nspin, n_nzs, DM, S, method = TS_INFO_SCF)

    call timer('TS',2)

#ifdef TS_DEV
    call die('to not disturb the TSDE')
#endif

  end subroutine transiesta

  subroutine ts_print_memory(Gamma)
    
    use parallel, only : IONode

#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Max
    use mpi_siesta, only : MPI_Double_Precision
#endif 

    use class_Sparsity
    use m_ts_options, only : IsVolt, Calc_Forces
    use m_ts_options, only : N_mu, N_Elec, Elecs
    use m_ts_options, only : no_BufL, no_BufR
    use m_ts_contour_neq, only : N_nEq_id
    use m_ts_tri_init, only : tri_parts, N_tri_part
    use m_ts_sparse, only : ts_sp_uc, tsup_sp_uc, ltsup_sp_sc
    use m_ts_electype

    logical, intent(in) :: Gamma ! SIESTA Gamma
    integer :: i, no_E
    real(dp) :: mem, tmp_mem
#ifdef MPI
    integer :: MPIerror
#endif

    ! estimate the amount of memory used...
    ! H and S
    i = nnzs(ts_sp_uc)
    mem = i * 2

    ! global sparsity update
    i = nnzs(tsup_sp_uc)
    mem = mem + i * max(N_mu,N_nEq_id)
    if ( Calc_Forces ) mem = mem + i * N_mu
    if ( Gamma ) then
       mem = mem * 8._dp
    else
       mem = mem * 16._dp
    end if

    ! local sparsity update
    if ( IsVolt ) then
       i = nnzs(ltsup_sp_sc)
       if ( Calc_Forces ) mem = mem + i * N_mu * 8._dp ! always in double 
       mem = mem + i * ( N_mu + N_nEq_id ) * 8._dp ! always in double
    end if

    ! Add electrode sizes
    tmp_mem = 0._dp
    do i = 1 , N_Elec
       if ( IsVolt ) then
          tmp_mem = tmp_mem + TotUsedOrbs(Elecs(i)) * UsedOrbs(Elecs(i)) * 2
          tmp_mem = tmp_mem + TotUsedOrbs(Elecs(i)) ** 2
       else
          tmp_mem = tmp_mem + TotUsedOrbs(Elecs(i)) * UsedOrbs(Elecs(i)) * 3
       end if
    end do
    mem = mem + tmp_mem * 16._dp

#ifdef MPI
    call MPI_Reduce(mem,tmp_mem,1, MPI_Double_Precision, MPI_Max, 0, &
         MPI_Comm_World,MPIerror)
    mem = tmp_mem
#endif

    mem = mem / 1024._dp ** 2
    if ( IONode ) then
       write(*,'(/,a,f10.2,a)') &
            'transiesta: Memory usage of sparse arrays and electrodes (static): ', &
            mem,'MB'
    end if

    if ( ts_method == TS_SPARSITY_TRI ) then
       ! Calculate size of the tri-diagonal matrix
       mem = tri_parts(N_tri_part)**2
       do i = 1 , N_tri_part - 1
          mem = mem + tri_parts(i)*( tri_parts(i) + 2 * tri_parts(i+1) )
       end do
       mem = mem * 16._dp * 2 / 1024._dp ** 2
       if ( IONode ) &
            write(*,'(a,f10.2,a)') &
            'transiesta: Memory usage of tri-diagonal matrices: ', &
            mem,'MB'
    else
       ! Calculate size of the full matrices
       no_E = sum(TotUsedOrbs(Elecs),.not. Elecs(:)%DM_CrossTerms)
       i = nrows_g(ts_sp_uc) - no_BufL - no_BufR
       ! LHS
       mem = i ** 2
       ! RHS
       if ( IsVolt ) then
          mem = mem + i * max(i-no_E,sum(TotUsedOrbs(Elecs)))
       else
          mem = mem + i * (i-no_E)
       end if
       mem = mem * 16._dp / 1024._dp ** 2
       if ( IONode ) &
            write(*,'(a,f10.2,a)') &
            'transiesta: Memory usage of full matrices: ', &
            mem,'MB'
    end if

  end subroutine ts_print_memory

end module m_transiesta

