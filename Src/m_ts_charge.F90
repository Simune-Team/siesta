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
  subroutine ts_get_charges(Nelecs,dit, sp, nspin, n_nzs, DM, S, Q)

    use m_ts_method
    use parallel, only : IONode, Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB
    use m_ts_electype

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: Nelecs
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and overlap
    real(dp), intent(in) :: DM(n_nzs,nspin), S(n_nzs)
    ! The charge in the regions
    real(dp), intent(out) :: Q(0:1+1+Nelecs*2, nspin)

! **********************
! * LOCAL variables    *
! **********************
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, jo, ir, jr, r
#ifdef MPI
    real(dp) :: tmp(0:1+1+Nelecs*2, nspin)
    integer :: MPIerror
#endif

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_lo,nrows_g=no_u)
    
    ! Initialize charges
#ifdef MPI
    tmp(:,:) = 0._dp
#else
    Q(:,:)   = 0._dp
#endif

    do lio = 1 , no_lo

       ! obtain the global index of the orbital.
       io = index_local_to_global(dit,lio,Node)
       ir = get_orb_type(io)

       ! Loop number of entries in the row... (index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)
          jr = get_orb_type(jo)

          if      ( all((/ir,jr/) == TYP_BUFFER) ) then
             r = 1 ! buffer
          else if ( any((/ir,jr/) == TYP_BUFFER) ) then
             r = 0 ! other
          else if ( all((/ir,jr/) == TYP_DEVICE) ) then
             r = 2 ! device
          else if ( any((/ir,jr/) == TYP_DEVICE) ) then
             r = 4+(ir+jr-1)*2 ! device/electrode
          else if ( ir == jr ) then
             r = 3+(ir-1)*2 ! electrode/electrode
          else
             r = 0 ! other
          end if
#ifdef MPI
          tmp(r,:) = tmp(r,:) + DM(ind,:) * S(ind)
#else
          Q(r,:) = Q(r,:) + DM(ind,:) * S(ind)
#endif
       end do
    end do

#ifdef MPI
    call MPI_AllReduce(tmp(0,1),Q(0,1),size(tmp), &
         MPI_Double_Precision,MPI_SUM, MPI_Comm_World,MPIerror)
#endif

  end subroutine ts_get_charges

    

  ! A subroutine for printing out the charge distribution in the cell
  ! it will currently only handle the full charge distribution, and
  ! not per k-point.
  subroutine ts_print_charges(Elecs,dit, sp, &
       nspin, n_nzs, DM, S, &
       method)
    use parallel, only : IONode, Node
    use m_ts_electype
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB

! **********************
! * INPUT variables    *
! **********************
    type(Elec), intent(in) :: Elecs(:)
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
    integer :: Nelecs, i
    real(dp), allocatable :: Q(:,:)
    integer :: ispin, lmethod

    lmethod = TS_INFO_FULL
    if ( present(method) ) lmethod = method

    Nelecs = size(Elecs)
    allocate(Q(0:2+Nelecs*2,nspin))

    call ts_get_charges(Nelecs,dit, sp, nspin, n_nzs, DM, S, Q)
         
    ! it will only be the IONode which will write out...
    if ( .not. IONode ) return

    if ( lmethod == TS_INFO_FULL ) then

       write(*,'(/,a)') 'transiesta: Charge distribution:'
       if ( nspin > 1 ) then
          write(*,'(a,2(f12.5,tr1))') &
               'Total charge                  [Q]  :', &
               sum(Q(:,1)),sum(Q(:,2))
          if ( Q(1,1) > 0._dp .or. Q(1,2) > 0._dp ) then
             write(*,'(a,2(f12.5,tr1))') &
               'Buffer                        [B]  :',Q(1,1), Q(1,2)
          end if
          write(*,'(a,2(f12.5,tr1))') &
               'Device                        [C]  :',Q(2,1), Q(2,2)
          do i = 1 , Nelecs 
             write(*,'(a,t31,a,i0,a,f12.5)') &
                  trim(name(Elecs(i))),'[E',i,'] :', &
                  Q(3+(i-1)*2,1), Q(3+(i-1)*2,2)
             write(*,'(a,t22,a,i0,a,2(f12.5,tr1))') &
                  trim(name(Elecs(i))),'/ device [C',i,'] :', &
                  Q(4+(i-1)*2,1), Q(4+(i-1)*2,2)
          end do
          write(*,'(a,2(f12.5,tr1),/)') &
               'Other                         [O]  :',Q(0,1), Q(0,2)

       else
          write(*,'(a,f12.5)') &
               'Total charge                  [Q]  :', sum(Q(:,1))
          if ( Q(1,1) > 0._dp ) then
             write(*,'(a,f12.5)') &
               'Buffer                        [B]  :',Q(1,1)
          end if
          write(*,'(a,f12.5)') &
               'Device                        [C]  :',Q(2,1)
          do i = 1 , Nelecs 
             write(*,'(a,t31,a,i0,a,f12.5)') &
               trim(name(Elecs(i)))         ,'[E',i,'] :',Q(3+(i-1)*2,1)
             write(*,'(a,t22,a,i0,a,f12.5)') &
               trim(name(Elecs(i))),'/ device [C',i,'] :',Q(4+(i-1)*2,1)
          end do
          write(*,'(a,f12.5,/)') &
               'Other                         [O]  :',Q(0,1)
       end if

    else if ( lmethod == TS_INFO_SCF ) then

       ! We write out the information from the SCF cycle...
       write(*,'(a,1x,a9)',advance='no') 'ts-q:','D'
       do i = 1 , Nelecs
          write(*,'(1x,a8,i0,1x,a8,i0)',advance='no') 'E',i,'C',i
       end do
       write(*,'(1x,a9)') 'Q'
       do ispin = 1 , nspin
          write(*,'(a,1x,f9.3)',advance='no') 'ts-q:', Q(2,ispin)
          do i = 1 , Nelecs
             write(*,'(2(1x,f9.3))',advance='no') Q(3+(i-1)*2,ispin),Q(4+(i-1)*2,ispin)
          end do
          write(*,'(1x,f9.3)') sum(Q(:,ispin))
       end do
       
    end if

    deallocate(Q)
    
  end subroutine ts_print_charges

  subroutine ts_charge_correct(no_Buf,Elecs, &
       dit, sp, nspin, n_nzs, DM, EDM, S, Qtot, &
       method)

    use m_ts_electype
    use class_OrbitalDistribution
    use class_Sparsity

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: no_Buf
    ! The electrodes
    type(Elec), intent(in) :: Elecs(:)
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
       call ts_charge_correct_buffer(no_Buf,Elecs, &
            dit, sp, nspin, n_nzs, DM, EDM, S, Qtot)
    end if

  end subroutine ts_charge_correct


  subroutine ts_charge_correct_buffer(no_Buf,Elecs, &
       dit, sp, nspin, n_nzs, DM, EDM, S, Qtot)

    use m_ts_electype
    use m_ts_method
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
    integer, intent(in) :: no_Buf
    ! The electrodes
    type(Elec), intent(in) :: Elecs(:)
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
    integer :: Nelecs
    real(dp), allocatable :: Q(:,:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, ispin
    real(dp) :: reD, addQ(nspin)

    Nelecs = size(Elecs)
    allocate(Q(0:2+Nelecs*2,nspin))
    
    call ts_get_charges(Nelecs, &
       dit, sp, nspin, n_nzs, DM, S, Q)

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_lo,nrows_g=no_u)
    
    ! Calculate the density factor for obtaining the correct charge
    ! For the left buffer region
    reD = (Qtot-sum(Q(:,:))) / sum(Q(1,:))
    
    ! immediately deallocate charge
    deallocate(Q)

    addQ(:) = 0.0_dp

    ! Reset to zero if not existing
    if ( no_Buf == 0 ) return

    ! Apply charge-correction factor 
    ! This will reduce "heavy" charge fluctuations and
    ! should guard against this.
    reD = reD * TS_RHOCORR_FACTOR
    
    do ispin = 1 , nspin
       do lio = 1 , no_lo
          
          ! obtain the global index of the orbital.
          io = index_local_to_global(dit,lio,Node)

          if ( get_orb_type(io) /= TYP_BUFFER ) cycle

          ! Loop number of entries in the row... (index frame)
          do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             if ( get_orb_type(l_col(ind)) /= TYP_BUFFER ) cycle
             
             DM(ind,ispin) = DM(ind,ispin) * reD + &
                  DM(ind,ispin)
             ! We are not sure to do with the energy matrix, hence we don't do anything
             ! As the energy density matrix is an integral over the 
             ! density times energy I expect this is the "correct" way to introduce 
             ! this...
             EDM(ind,ispin) = EDM(ind,ispin) * reD + &
                  EDM(ind,ispin)
             addQ(ispin) = addQ(ispin) + &
                  DM(ind,ispin) * reD

          end do
       end do
       
    end do

  end subroutine ts_charge_correct_buffer

end module m_ts_charge
