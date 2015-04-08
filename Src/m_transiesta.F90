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
  use m_ts_method, only : ts_method
  use m_ts_method, only : TS_FULL, TS_BTD
#ifdef MUMPS
  use m_ts_method, only : TS_MUMPS
#endif

  use m_ts_tri_init, only : ts_tri_init

  use m_ts_fullg
  use m_ts_fullk

  use m_ts_trig
  use m_ts_trik

#ifdef MUMPS
  use m_ts_mumpsg
  use m_ts_mumpsk
#endif

  implicit none

  public :: transiesta
  private

contains

  subroutine transiesta(TSiscf,nspin, &
       sp_dist, sparse_pattern, &
       Gamma, ucell, nsc, isc_off, no_u, na_u, lasto, xa, n_nzs, &
       H, S, DM, EDM, Ef, kT, &
       Qtot, Fermi_correct)

    use units, only : eV
    use alloc, only : re_alloc, de_alloc

    use parallel, only : IONode

    use class_OrbitalDistribution
    use class_Sparsity

    use m_ts_kpoints, only : ts_nkpnt, ts_Gamma

    use m_ts_electype

    use m_ts_options, only : N_Elec, Elecs
    use m_ts_options, only : IsVolt, Calc_Forces

    use m_ts_options, only : opt_TriMat_method

    use m_ts_contour_eq , only : N_Eq_E
    use m_ts_contour_neq, only : N_nEq_E

    use m_ts_charge
    
    use m_ts_gf, only : read_Green
    use m_interpolate

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in)  :: TSiscf
    integer, intent(in)  :: nspin
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    logical, intent(in)  :: Gamma
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in)  :: nsc(3), no_u, na_u
    integer, intent(in) :: isc_off(3,product(nsc))
    integer, intent(in)  :: lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: H(n_nzs,nspin), S(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: kT, Qtot
    real(dp), intent(inout) :: Ef
    logical, intent(in) :: Fermi_correct

! ******************** IO descriptors ************************
    integer, allocatable :: uGF(:)
! ************************************************************

! ****************** Electrode variables *********************
    integer, allocatable :: nq(:)
! ************************************************************

! * local variables
    integer :: iEl, NEn, no_used, no_used2
    logical :: converged
    ! In case of Fermi-correction, we save the previous steps
    ! and do a spline interpolation... :)
    integer :: N_F, i_F, ioerr
    real(dp), pointer :: Q_Ef(:,:) => null()

    ! Open GF files...
    ! Read-in header of Green functions
    ! Prepare for the calculation
    ! We read in the k-points that the electrode was generated with.
    ! Furthermore we read in the expansion q-points
    ! They are communicated in the routine

    if ( TSiscf == 1 ) then
       ! We need to initialize TRANSIESTA

       call timer('TS_init',1)

       ! For the fermi-correction, we need the 
       ! local sparsity pattern...
       converged = IsVolt .or. TS_RHOCORR_METHOD == TS_RHOCORR_FERMI
       call ts_sparse_init(slabel,converged, N_Elec, Elecs, &
            ucell, nsc, na_u, xa, lasto, sp_dist, sparse_pattern, Gamma, &
            isc_off)

       if ( ts_method == TS_BTD ) then
          ! initialize the tri-diagonal partition
          call ts_tri_init( sp_dist, sparse_pattern , N_Elec, &
               Elecs, IsVolt, ucell, na_u, lasto ,nsc, isc_off, &
               opt_TriMat_method )
       end if

       ! print out estimated memory usage...
       call ts_print_memory(ts_Gamma)

       call ts_print_charges(N_Elec,Elecs, sp_dist, sparse_pattern, &
            nspin, n_nzs, DM, S)

       if ( .not. Calc_Forces .and. IONode ) then
          write(*,'(a)') 'transiesta: *** Notice that the forces are NOT updated ***'
          write(*,'(a)') 'transiesta: *** Will set the forces to zero ***'
       end if
       if ( .not. Calc_Forces ) then
          if ( IONode ) then
             write(*,'(a)') 'transiesta: *** The forces are NOT updated ***'
             write(*,'(a)') 'transiesta: ***  Will set the forces to 0  ***'
          end if
!$OMP parallel workshare default(shared)
          EDM(:,:) = 0._dp
!$OMP end parallel workshare
       end if

       call timer('TS_init',2)

    end if


    call timer('TS',1)

    ! Total number of energy-points...
    NEn = N_Eq_E() + N_nEq_E()

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

       if ( IsVolt .or. .not. Elecs(iEl)%Bulk ) then
          ! If we using bulk electrodes, we need not the Hamiltonian, 
          ! nor the overlap...
          call re_alloc(Elecs(iEl)%HA,1,no_used,1,no_used,1,nq(iEl),routine='transiesta')
          call re_alloc(Elecs(iEl)%SA,1,no_used,1,no_used,1,nq(iEl),routine='transiesta')
       end if

       no_used = TotUsedOrbs(Elecs(iEl))
       if ( IsVolt ) then
          ! We need Gamma's with voltages (now they are both GAA and GammaT)
          no_used2 = no_used
       else 
          ! This is only for having space for GA
          if ( Elecs(iEl)%pre_expand > 0 ) then
             no_used2 = no_used
          else
             no_used2 = Elecs(iEl)%no_used
          end if
       end if
       call re_alloc(Elecs(iEl)%Gamma,1,no_used*no_used2,routine='transiesta')

       ! This seems stupid, however, we never use the expansion array and
       ! GammaT at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called:
       ! first the GAA is "emptied" of information and then
       ! Gamma is filled.
       if ( Elecs(iEl)%pre_expand == 0 ) no_used2 = Elecs(iEl)%no_used
       Elecs(iEl)%GA => Elecs(iEl)%Gamma(1:no_used*no_used2)

    end do

    ! start calculation
    converged = .false.
    if ( Fermi_correct ) then

       ! we will utilize the old Fermi-level to correct the 
       ! EDM matrix (just in case the 
       ! electrode region elements are not taken care of)

       ! Allocate for interpolation
       N_F = 10
       i_F = 0
       call re_alloc(Q_Ef,1,N_F,1,2)

    end if

    do while ( .not. converged ) 

       call open_GF(N_Elec,Elecs,uGF,NEn,.false.)
       
       if ( ts_method == TS_FULL ) then
          if ( ts_Gamma ) then
             call ts_fullg(N_Elec,Elecs, &
                  nq, uGF, nspin, na_u, lasto, &
                  sp_dist, sparse_pattern, &
                  no_u, n_nzs, &
                  H, S, DM, EDM, Ef, kT)
          else
             call ts_fullk(N_Elec,Elecs, &
                  nq, uGF, &
                  ucell, nspin, na_u, lasto, &
                  sp_dist, sparse_pattern, &
                  no_u, n_nzs, &
                  H, S, DM, EDM, Ef, kT)
          end if
       else if ( ts_method == TS_BTD ) then
          if ( ts_Gamma ) then
             call ts_trig(N_Elec,Elecs, &
                  nq, uGF, nspin, na_u, lasto, &
                  sp_dist, sparse_pattern, &
                  no_u, n_nzs, &
                  H, S, DM, EDM, Ef, kT)
          else
             call ts_trik(N_Elec,Elecs, &
                  nq, uGF, &
                  ucell, nspin, na_u, lasto, &
                  sp_dist, sparse_pattern, &
                  no_u, n_nzs, &
                  H, S, DM, EDM, Ef, kT)
          end if
#ifdef MUMPS
       else if ( ts_method == TS_MUMPS ) then
          if ( ts_Gamma ) then
             call ts_mumpsg(N_Elec,Elecs, &
                  nq, uGF, nspin, na_u, lasto, &
                  sp_dist, sparse_pattern, &
                  no_u, n_nzs, &
                  H, S, DM, EDM, Ef, kT)
          else
             call ts_mumpsk(N_Elec,Elecs, &
                  nq, uGF, &
                  ucell, nspin, na_u, lasto, &
                  sp_dist, sparse_pattern, &
                  no_u, n_nzs, &
                  H, S, DM, EDM, Ef, kT)
          end if
#endif
       else

          call die('Error in code')
       end if

       ! Close files
       do iEl = 1 , N_Elec
          if ( IONode .and. Elecs(iEl)%out_of_core ) then
             call io_close(uGF(iEl))
          end if
       end do
       
       if ( Fermi_correct ) then

          i_F = i_F + 1
          if ( N_F < i_F ) then
             N_F = N_F + 10
             call re_alloc(Q_Ef,1,N_F,1,2,copy=.true.)
          end if

          ! Save current fermi level and charge
          call ts_get_charges(N_Elec, sp_dist, sparse_pattern, &
               nspin, n_nzs, DM, S, Qtot = Q_Ef(i_F,1) )
          Q_Ef(i_F,2) = Ef

          if ( i_F < 2 ) then

          call open_GF(N_Elec,Elecs,uGF,1,.true.)
          
          if ( ts_method == TS_FULL ) then
             if ( Q_Ef(i_F,1) > Qtot ) then
                Ef = Ef - 0.01_dp * eV
             else
                Ef = Ef + 0.01_dp * eV
             end if
          else if ( ts_method == TS_BTD ) then
             if ( ts_Gamma ) then
                call ts_trig_Fermi(N_Elec,Elecs, &
                     nq, uGF, nspin, na_u, lasto, &
                     sp_dist, sparse_pattern, &
                     no_u, n_nzs, &
                     H, S, DM, Ef, kT, Qtot, converged)
             else
                call ts_trik_Fermi(N_Elec,Elecs, &
                     nq, uGF, &
                     ucell, nspin, na_u, lasto, &
                     sp_dist, sparse_pattern, &
                     no_u, n_nzs, &
                     H, S, DM, Ef, kT, Qtot, converged)
             end if
#ifdef MUMPS
          else if ( ts_method == TS_MUMPS ) then
             if ( Q_Ef(i_F,1) > Qtot ) then
                Ef = Ef - 0.01_dp * eV
             else
                Ef = Ef + 0.01_dp * eV
             end if
#endif
          else
             
             call die('Error in code')
          end if

          ! Close files
          do iEl = 1 , N_Elec
             if ( IONode .and. Elecs(iEl)%out_of_core ) then
                call io_close(uGF(iEl))
             end if
          end do

          else

             ! In case we have accumulated 2 or more points
             call interp_spline(i_F,Q_Ef(1:i_F,1),Q_Ef(1:i_F,2),Qtot,Ef)

             ! Truncate to the maximum allowed change in Fermi-level
             converged = ts_qc_Fermi_truncate(Q_Ef(i_F,2), &
                  TS_RHOCORR_FERMI_MAX, Ef)

             if ( IONode ) then
                write(*,'(a,e11.4,a)') 'transiesta: cubic spline. dEf = ', &
                     (Ef-Q_Ef(i_F,2))/eV, ' eV'
             end if

             ! Even if we have converged we allow the interpolation
             ! to do a final step. If dQ is very small it should be very
             ! close to the found value.
             ! If the truncation already is reached we stop as that
             ! *MUST* be the maximal change.
             if ( .not. converged ) &
                  converged = abs(Q_Ef(i_F,1) - Qtot) < &
                  TS_RHOCORR_FERMI_TOLERANCE

          end if

       else
          
          ! If no Fermi-correction, we are converged
          converged = .true.

       end if

    end do

    if ( IONode .and. Fermi_correct ) then

       ! After converge we write out the convergence
       call io_assign(iEl)
       inquire(file='TS_FERMI', exist=converged)
       if ( converged ) then
          open(unit=iEl,file='TS_FERMI',position='append',form='formatted', &
               status='old',iostat=ioerr)
          write(iEl,'(/,a,i0)') '# TSiscf = ',TSiscf
       else
          open(unit=iEl,file='TS_FERMI',form='formatted', &
               status='new')
          write(iEl,'(a,i0)') '# TSiscf = ',TSiscf
       end if
       N_F = i_F
       write(iEl,'(a,i0)')'# ',N_F ! Number of iterations
       do i_F = 1 , N_F
          write(iEl,'(2(tr1,e15.6))') Q_Ef(i_F,2)/eV,Q_Ef(i_F,1) - Qtot
       end do

       call io_close(iEl)

    end if
    if ( Fermi_correct ) then

       ! Guess-stimate the actual Fermi-shift
       ! typically will the above be "too" little
       ! So we interpolate between all previous 
       ! estimations for this geometry...
       call ts_qc_Fermi_file(Ef)

       ! We have now calculated the new Ef
       ! We shift it EDM to the correct level
       Q_Ef(1,2) = Ef - Q_Ef(1,2)
       call daxpy(n_nzs*nspin,Q_Ef(1,2),DM(1,1),1,EDM(1,1),1)

       call de_alloc(Q_Ef)

    end if

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

    ! We do the charge correction of the transiesta
    ! computation here (notice that the routine will automatically
    ! return if no charge-correction is requested)
    call ts_qc(N_Elec,Elecs, sp_dist, &
         sparse_pattern, nspin, n_nzs, DM, EDM, S, Qtot, &
         TS_RHOCORR_METHOD)

    call ts_print_charges(N_Elec,Elecs, sp_dist, sparse_pattern, &
         nspin, n_nzs, DM, S, method = TS_INFO_SCF)

    call timer('TS',2)

#ifdef TS_DEV
    call die('to not disturb the TSDE')
#endif

  contains

    subroutine init_Electrode_HS(El)
      use class_Sparsity
      use class_dSpData1D
      use class_dSpData2D
      use alloc, only : re_alloc
      type(Elec), intent(inout) :: El
      
      ! If already initialized, return immediately
      if ( initialized(El%sp) ) return

      ! Read-in and create the corresponding transfer-matrices
      call delete(El) ! ensure clean electrode
      call read_Elec(El,Bcast=.true., IO = .false.)
      
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

    subroutine open_GF(N_Elec,Elecs,uGF,NEn,Fermi_correct)
      integer, intent(in) :: N_Elec
      type(Elec), intent(inout) :: Elecs(N_Elec)
      integer, intent(out) :: uGF(N_Elec)
      integer, intent(in) :: NEn
      logical, intent(in) :: Fermi_correct

      ! Local variables
      integer :: iEl
      
      do iEl = 1 , N_Elec

         ! Initialize k-points
         Elecs(iEl)%bkpt_cur(:) = huge(1._dp)

         if ( .not. Fermi_correct ) then
            if ( Elecs(iEl)%out_of_core ) then
               
               if ( IONode ) then
                  call io_assign(uGF(iEl))
                  open(file=Elecs(iEl)%GFfile,unit=uGF(iEl),form='unformatted')
               end if
               
            else

               ! prepare the electrode to create the surface self-energy
               call init_Electrode_HS(Elecs(iEl))
               
            end if
         else

            if ( Elecs(iEl)%out_of_core ) then
               if ( IONode ) then
                  call io_assign(uGF(iEl))
                  open(file=trim(Elecs(iEl)%GFfile)//'-Fermi', &
                       unit=uGF(iEl),form='unformatted')
               end if
            end if
            
         end if

         if ( Elecs(iEl)%out_of_core ) then
            call read_Green(uGF(iEl),Elecs(iEl), ts_nkpnt, NEn )
         end if
         
      end do
      
    end subroutine open_GF

  end subroutine transiesta

  subroutine ts_print_memory(ts_Gamma)
    
    use parallel, only : IONode

#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Max
    use mpi_siesta, only : MPI_Double_Precision
#endif 

    use class_Sparsity
    use m_ts_options, only : IsVolt, Calc_Forces
    use m_ts_options, only : N_mu, N_Elec, Elecs
    use m_ts_contour_neq, only : N_nEq_id
    use m_ts_sparse, only : ts_sp_uc, tsup_sp_uc, ltsup_sp_sc
    use m_ts_electype

    use m_ts_tri_init, only : c_Tri
    use m_ts_tri_common, only : GFGGF_needed_worksize, nnzs_tri
    use m_ts_method, only : no_Buf

    logical, intent(in) :: ts_Gamma ! transiesta Gamma
    integer :: i, f, no_E
    integer :: padding, worksize
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
    if ( ts_Gamma ) then
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
       f = 1
       if ( IsVolt ) then
          if ( Elecs(i)%pre_expand > 1 ) f = product(Elecs(i)%Rep)
          tmp_mem = tmp_mem + f * TotUsedOrbs(Elecs(i)) * Elecs(i)%no_used * 2 !H,S
          tmp_mem = tmp_mem + TotUsedOrbs(Elecs(i)) ** 2 ! GS/Gamma
       else
          if ( .not. Elecs(i)%Bulk ) then
             if ( Elecs(i)%pre_expand > 1 ) f = product(Elecs(i)%Rep) ! H,S
             tmp_mem = tmp_mem + f * TotUsedOrbs(Elecs(i)) * Elecs(i)%no_used * 2
          end if
          if ( Elecs(i)%pre_expand > 0 ) f = product(Elecs(i)%Rep) ! GS
          tmp_mem = tmp_mem + f * TotUsedOrbs(Elecs(i)) * Elecs(i)%no_used
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

    if ( ts_method == TS_BTD ) then

       ! Calculate size of the tri-diagonal matrix
       if ( IsVolt ) then
          call GFGGF_needed_worksize(c_Tri%n,c_Tri%r, &
               N_Elec, Elecs, padding, worksize)
       else
          padding = 0
          worksize = 0
       end if

       mem = nnzs_tri(c_Tri%n,c_Tri%r)
       mem = (mem * 2 + padding + worksize ) * 16._dp / 1024._dp ** 2
       if ( IONode ) &
            write(*,'(a,f10.2,a)') &
            'transiesta: Memory usage of tri-diagonal matrices: ', &
            mem,'MB'
    else if ( ts_method == TS_FULL ) then
       ! Calculate size of the full matrices
       ! Here we calculate number of electrodes not needed to update the cross-terms
       no_E = sum(TotUsedOrbs(Elecs),Elecs(:)%DM_update==0)
       i = nrows_g(ts_sp_uc) - no_Buf
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
#ifdef MUMPS
    else if ( ts_method == TS_MUMPS ) then
       if ( IONode ) then
          write(*,'(a)')'transiesta: Memory usage is determined by MUMPS.'
          write(*,'(a)')'transiesta: Search in TS_MUMPS_<Node>.dat for: ### Minimum memory.'
       end if
#endif
    end if

  end subroutine ts_print_memory

end module m_transiesta

