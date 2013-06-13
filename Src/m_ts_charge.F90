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

! Module for correcting the density matrix for retaining a constant charge density
! The idea is to introduce several different schemes of charge corrections.

module m_ts_charge

  use precision, only: dp
  
  implicit none

  public 

  ! Info parameters for obtaining charge calculations (mulliken charges in regions)
  integer, parameter :: TS_INFO_FULL = 0
  integer, parameter :: TS_INFO_SCF = 1


  ! Method parameters for the charge-correction
  integer, save :: TS_RHOCORR_METHOD = 0
  integer, parameter :: TS_RHOCORR_BUFFER = 1
  integer, parameter :: TS_RHOCORR_UPDATE = 2
  real(dp), save :: TS_RHOCORR_FACTOR = 0.75_dp

  private :: dp

contains

  ! Retrive the mulliken charges in each region of the transiesta setup
  subroutine ts_get_charges(ElLeft, ElRight, no_BufL, no_BufR, &
       dit, sp, nspin, n_nzs, DM, S, Q)

    use m_ts_electype
    ! left stuff
    use m_ts_method, only : get_scat_region
    use parallel, only : IONode, Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB

! **********************
! * INPUT variables    *
! **********************
    type(OrbitalDistribution), intent(inout) :: dit
    ! The electrodes
    type(Elec), intent(in) :: ElLeft, ElRight
    ! The buffer regions
    integer, intent(in) :: no_BufL, no_BufR
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and overlap
    real(dp), intent(in) :: DM(n_nzs,nspin), S(n_nzs)
    ! The charge in the regions
    real(dp), intent(out) :: Q(0:9, nspin)

! **********************
! * LOCAL variables    *
! **********************
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, jo, ispin, r
    integer :: no_u_TS
    integer :: no_L, no_R
    integer :: lmethod
#ifdef MPI
    real(dp) :: tmp(0:9, nspin)
    integer :: MPIerror
#endif

    ! Calculate the buffer region and electrode orbitals
    no_L = TotUsedOrbs(ElLeft)
    no_R = TotUsedOrbs(ElRight)

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_lo,nrows_g=no_u)
    
    no_u_TS = no_u - no_BufL - no_BufR

    ! Initialize charges
#ifdef MPI
    tmp(:,:) = 0._dp
#else
    Q(:,:)   = 0._dp
#endif

    do ispin = 1 , nspin
       do lio = 1 , no_lo

          ! obtain the global index of the orbital.
          io = index_local_to_global(dit,lio,Node) - no_BufL

          ! Loop number of entries in the row... (index frame)
          do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             ! as the local sparsity pattern is a super-cell pattern,
             ! we need to check the unit-cell orbital
             ! The unit-cell column index
             jo = UCORB(l_col(ind),no_u) - no_BufL
             
             r = get_scat_region(io,no_L,jo,no_R,no_u_TS)
#ifdef MPI
             tmp(r,ispin) = tmp(r,ispin) + &
                  DM(ind,ispin) * S(ind)
#else
             Q(r,ispin) = Q(r,ispin) + &
                  DM(ind,ispin) * S(ind)
#endif
          end do
       end do
    end do

#ifdef MPI
    call MPI_AllReduce(tmp(0,1),Q(0,1),10*nspin,MPI_Double_Precision,MPI_SUM, &
         MPI_Comm_World,MPIerror)
#endif

  end subroutine ts_get_charges

    

  ! A subroutine for printing out the charge distribution in the cell
  ! it will currently only handle the full charge distribution, and
  ! not per k-point.
  subroutine ts_print_charges(ElLeft, ElRight, no_BufL, no_BufR, &
       dit, sp, &
       nspin, n_nzs, DM, S, &
       method)
    use m_ts_electype
    ! left stuff
    use m_ts_method, only : get_scat_region
    use parallel, only : IONode, Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB

! **********************
! * INPUT variables    *
! **********************
    type(Elec), intent(in) :: ElLeft, ElRight
    integer, intent(in) :: no_BufL, no_BufR
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and overlap
    real(dp), intent(in) :: DM(n_nzs,nspin), S(n_nzs)
    ! The method by which it should be printed out...
    integer, intent(in), optional :: method

! **********************
! * LOCAL variables    *
! **********************
    real(dp) :: Q(0:9,nspin)
    integer :: ispin, lmethod

    lmethod = TS_INFO_FULL
    if ( present(method) ) lmethod = method

    call ts_get_charges(ElLeft, ElRight, no_BufL, no_BufR, &
         dit, sp, nspin, n_nzs, DM, S, Q)
         
    ! it will only be the IONode which will write out...
    if ( .not. IONode ) return

    if ( lmethod == TS_INFO_FULL ) then

       write(*,'(/,a)') 'transiesta: Charge distribution:'
       if ( nspin > 1 ) then
          write(*,'(a,2(f12.5,tr1))') &
               'Total charge                  [Q]    :', &
               sum(Q(:,1)),sum(Q(:,2))
          if ( no_BufL > 0 ) write(*,'(a,2(f12.5,tr1),/,a,2(f12.5,tr1))') &
               'Left buffer                   [LB]   :',Q(1,1), Q(1,2), &
               'Left buffer/left electrode    [LB-L] :',Q(2,1), Q(2,2)
          write(*,'(a,2(f12.5,tr1),4(/,a,2(f12.5,tr1)))') &
               'Left electrode                [L]    :',Q(3,1), Q(3,2), &
               'Left electrode/device         [L-C]  :',Q(4,1), Q(4,2), &
               'Device                        [C]    :',Q(5,1), Q(5,2), &
               'Device/right electrode        [C-R]  :',Q(6,1), Q(6,2), &
               'Right electrode               [R]    :',Q(7,1), Q(7,2)
          if ( no_BufR > 0 ) write(*,'(a,2(f12.5,tr1),/,a,2(f12.5,tr1))') &
               'Right electrode/right buffer  [R-RB] :',Q(8,1), Q(8,2), &
               'Right buffer                  [RB]   :',Q(9,1), Q(9,2)
          write(*,'(a,2(f12.5,tr1),/)') &
               'Other                         [O]    :',Q(0,1), Q(0,2)

       else
          write(*,'(a,f12.5)') &
               'Total charge                  [Q]    :',sum(Q(:,1))
          if ( no_BufL > 0 ) write(*,'(a,f12.5,/,a,f12.5)') &
               'Left buffer                   [LB]   :',Q(1,1), &
               'Left buffer/left electrode    [LB-L] :',Q(2,1)
          write(*,'(a,f12.5,4(/,a,f12.5))') &
               'Left electrode                [L]    :',Q(3,1), &
               'Left electrode/device         [L-C]  :',Q(4,1), &
               'Device                        [C]    :',Q(5,1), &
               'Device/right electrode        [C-R]  :',Q(6,1), &
               'Right electrode               [R]    :',Q(7,1)
          if ( no_BufR > 0 ) write(*,'(a,f12.5,/,a,f12.5)') &
               'Right electrode/right buffer  [R-RB] :',Q(8,1), &
               'Right buffer                  [RB]   :',Q(9,1)
          write(*,'(a,f12.5,/)') &
               'Other                         [O]    :',Q(0,1)
       end if

    else if ( lmethod == TS_INFO_SCF ) then

       ! We write out the information from the SCF cycle...

       write(*,'(a,7(1x,a9))') 'ts-charge:','O','L','L-C','C', &
            'C-R','R','Qt'
       do ispin = 1 , nspin
          write(*,'(a,7(1x,f9.3))') 'ts-charge:', &
               Q(0,ispin), &
               Q(3,ispin),Q(4,ispin), &
               Q(5,ispin),Q(6,ispin), &
               Q(7,ispin),sum(Q(:,ispin))
       end do

    end if
    
  end subroutine ts_print_charges

  subroutine ts_charge_correct(ElLeft, ElRight,&
       no_BufL, no_BufR, &
       dit, sp, nspin, n_nzs, DM, EDM, S, Qtot, &
       method)

    use m_ts_electype
    use parallel, only : IONode, Node
    use class_OrbitalDistribution
    use class_Sparsity

! **********************
! * INPUT variables    *
! **********************
    ! The electrodes
    type(Elec), intent(in) :: ElLeft, ElRight
    ! The buffer regions
    integer, intent(in) :: no_BufL, no_BufR
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrices and overlap
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: S(n_nzs)
    real(dp), intent(in) :: Qtot
    integer, intent(in) :: method

    if ( method == TS_RHOCORR_BUFFER ) then
       call ts_charge_correct_buffer(ElLeft, ElRight,&
            no_BufL, no_BufR, &
            dit, sp, nspin, n_nzs, DM, EDM, S, Qtot)
    end if

  end subroutine ts_charge_correct


  subroutine ts_charge_correct_buffer(ElLeft, ElRight,&
       no_BufL, no_BufR, &
       dit, sp, nspin, n_nzs, DM, EDM, S, Qtot)

    use m_ts_electype
    ! left stuff
    use m_ts_method, only : get_scat_region
    use parallel, only : IONode, Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB

! **********************
! * INPUT variables    *
! **********************
    ! The electrodes
    type(Elec), intent(in) :: ElLeft, ElRight
    ! The buffer regions
    integer, intent(in) :: no_BufL, no_BufR
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and energy density matrix
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    ! The overlap
    real(dp), intent(in) :: S(n_nzs)
    ! Total charge of the system
    real(dp), intent(in) :: Qtot

! **********************
! * LOCAL variables    *
! **********************
    ! The charge in the regions
    real(dp) :: Q(0:9, nspin)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, jo, ispin, r
    integer :: no_u_TS
    integer :: no_L, no_R
    real(dp) :: reD(2), addQ(2,nspin)
#ifdef MPI
    real(dp) :: tmp(2, nspin)
    integer :: MPIerror
#endif

    call ts_get_charges(ElLeft, ElRight, no_BufL, no_BufR, &
       dit, sp, nspin, n_nzs, DM, S, Q)

    ! Calculate the buffer region and electrode orbitals
    no_L = TotUsedOrbs(ElLeft)
    no_R = TotUsedOrbs(ElRight)

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_lo,nrows_g=no_u)
    
    no_u_TS = no_u - no_BufL - no_BufR

    ! Calculate the density factor for obtaining the correct charge
    ! For the left buffer region
    reD(1) = real(no_BufL+no_BufR,dp) / real(no_BufL,dp)
    reD(1) = (Qtot-sum(Q(:,:))) / (reD(1)*sum(Q(1,:)))

    ! Calculate the density factor for obtaining the correct charge
    ! For the right buffer region
    reD(2) = real(no_BufL+no_BufR,dp) / real(no_BufR,dp)
    reD(2) = (Qtot-sum(Q(:,:))) / (reD(2) * sum(Q(9,:)))

    ! Reset to zero if not existing
    if ( no_BufL == 0 ) reD(1) = 0.0_dp
    if ( no_BufR == 0 ) reD(2) = 0.0_dp

    ! Apply charge-correction factor 
    ! This will reduce "heavy" charge fluctuations and
    ! should guard against this.
    reD(:) = reD(:) * TS_RHOCORR_FACTOR

    addQ(:,:) = 0.0_dp
    
    do ispin = 1 , nspin
       do lio = 1 , no_lo
          
          ! obtain the global index of the orbital.
          io = index_local_to_global(dit,lio,Node) - no_BufL

          ! Loop number of entries in the row... (index frame)
          do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             r = get_scat_region(io,no_L,jo,no_R,no_u_TS)

             if ( r == 1 .and. no_BufL > 0 ) then

                DM(ind,ispin) = DM(ind,ispin)*reD(1)+ &
                     DM(ind,ispin)
! We are not sure to do with the energy matrix, hence we don't do anything
!                EDM(ind,ispin) = EDM(ind,ispin)*reD(1)+ &
!                     EDM(ind,ispin)
                addQ(1,ispin) = addQ(1,ispin) + &
                     DM(ind,ispin)*reD(1)
             else if ( r == 9 .and. no_BufR > 0 ) then
                DM(ind,ispin) = DM(ind,ispin)*reD(2) + &
                     DM(ind,ispin)
! We are not sure to do with the energy matrix, hence we don't do anything
!                EDM(ind,ispin) = EDM(ind,ispin)*reD(2) + &
!                     EDM(ind,ispin)
                addQ(2,ispin) = addQ(2,ispin) + &
                     DM(ind,ispin)*reD(2)
             end if
          end do
       end do
       
    end do

  end subroutine ts_charge_correct_buffer

#ifdef TS_NOT_IMPLEMENTED_YET
  subroutine ts_charge_correct_update(ElLeft, ElRight,&
       no_BufL, no_BufR, &
       dit, sp, nspin, n_nzs, DM, EDM, S, Qtot)

    use m_ts_electype
    ! left stuff
    use m_ts_method, only : get_scat_region
    use parallel, only : IONode, Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB

! **********************
! * INPUT variables    *
! **********************
    type(OrbitalDistribution), intent(inout) :: dit
    ! The electrodes
    type(Elec), intent(in) :: ElLeft, ElRight
    ! The buffer regions
    integer, intent(in) :: no_BufL, no_BufR
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and overlap
    real(dp), intent(in) :: DM(n_nzs,nspin), S(n_nzs)
    ! The charge in the regions
    real(dp), intent(out) :: Q(0:9, nspin)

! **********************
! * LOCAL variables    *
! **********************
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, jo, ispin, r
    integer :: no_u_TS
    integer :: no_L, no_R
    real(dp) :: reD(2), addQ(2,nspin)
#ifdef MPI
    real(dp) :: tmp(2, nspin)
    integer :: MPIerror
#endif

    call ts_get_charges(ElLeft, ElRight, no_BufL, no_BufR, &
       dit, sp, nspin, n_nzs, DM, S, Q)

    ! Calculate the buffer region and electrode orbitals
    no_L = TotUsedOrbs(ElLeft)
    no_R = TotUsedOrbs(ElRight)

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_lo,nrows_g=no_u)
    
    no_u_TS = no_u - no_BufL - no_BufR

    ! Calculate the density factor for obtaining the correct charge
    ! For the left buffer region
    reD(1) = real(no_BufL+no_BufR,dp) / real(no_BufL,dp)
    reD(1) = (Qtot-sum(Q(:,:))) / (reD(1)*sum(Q(1,:)))

    ! Calculate the density factor for obtaining the correct charge
    ! For the right buffer region
    reD(2) = real(no_BufL+no_BufR,dp) / real(no_BufR,dp)
    reD(2) = (Qtot-sum(Q(:,:))) / (reD(2) * sum(Q(9,:)))

    ! Reset to zero if not existing
    if ( no_BufL == 0 ) reD(1) = 0.0_dp
    if ( no_BufR == 0 ) reD(2) = 0.0_dp

    ! Apply charge-correction factor 
    ! This will reduce "heavy" charge fluctuations and
    ! should guard against this.
    reD(:) = reD(:) * ChargeCorr_factor

    addQ(:,:) = 0.0_dp
    
    do ispin = 1 , nspin
       do lio = 1 , no_lo
          
          ! obtain the global index of the orbital.
          io = index_local_to_global(dit,lio,Node) - no_BufL

          ! Loop number of entries in the row... (index frame)
          do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             r = get_scat_region(io,no_L,jo,no_R,no_u_TS)

             if ( r == 1 .and. no_BufL > 0 ) then

                DM(ind,ispin) = DM(ind,ispin)*reD(1)+ &
                     DM(ind,ispin)
! We are not sure to do with the energy matrix, hence we don't do anything
!                EDM(ind,ispin) = EDM(ind,ispin)*reD(1)+ &
!                     EDM(ind,ispin)
                addQ(1,ispin) = addQ(1,ispin) + &
                     DM(ind,ispin)*reD(1)
             else if ( r == 9 .and. no_BufR > 0 ) then
                DM(ind,ispin) = DM(ind,ispin)*reD(2) + &
                     DM(ind,ispin)
! We are not sure to do with the energy matrix, hence we don't do anything
!                EDM(ind,ispin) = EDM(ind,ispin)*reD(2) + &
!                     EDM(ind,ispin)
                addQ(2,ispin) = addQ(2,ispin) + &
                     DM(ind,ispin)*reD(2)
             end if
          end do
       end do
       
    end do

  end subroutine ts_charge_correct_update
#endif

end module m_ts_charge
