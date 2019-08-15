!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2019, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! Module for calculating energy contributions for the transiesta
! calculations.
! Effectively many of the energy contributions stem from Hamiltonian
! and/or DM elements.
! However, in TranSiesta we only update/use a subset of the full
! matrices. As such the energies should only be calculated
! on the used elements rather than the full.
!
! There are many more elements that need to be split, while currently
! we only calculate a subset of what is actually needed.
! The main problem in calculating NEGF energies is:
!  1. the open-boundary problem, and
!  2. how to handle real-space energies such as Exc from dhscf.
module ts_energies_m

  implicit none

contains

  subroutine ts_compute_energies()

    use precision, only: dp
    use m_ts_options, only: N_Elec, Elecs
    use m_spin, only: spin

    use sparse_matrices, only: dit => block_dist
    use sparse_matrices, only: n_nzs => maxnh
    use sparse_matrices, only: Dold, Dscf, H
    use sparse_matrices, only: H_kin_1D, H_vkb_1D, H_so_2D
    use sparse_matrices, only: sp => sparse_pattern
    
    use m_ts_method
    use parallel, only : Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use geom_helper, only : UCORB
    use m_ts_electype, only: Elec
    use m_energies, only: NEGF_Ebs, NEGF_Ekin, NEGF_Eso, NEGF_Etot, NEGF_Enl
    use m_energies, only: NEGF_DEharr

! **********************
! * LOCAL variables    *
! **********************
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, jo, ir, jr, r, ispin
    real(dp), pointer :: H_vkb(:), H_kin(:), H_so(:,:)
    real(dp) :: Etmp(5,0:1+1+N_Elec*2)
#ifdef MPI
    real(dp) :: tmp(5,0:1+1+N_Elec*2)
    integer :: MPIerror
#endif

    H_vkb => val(H_vkb_1D)
    H_kin => val(H_kin_1D)
    if ( spin%SO ) then
      H_so => val(H_so_2D)
    end if

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
        n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
        nrows=no_lo,nrows_g=no_u)
    
    ! Initialize energies
    Etmp(:,:) = 0._dp
    
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ir,ind,jo,jr,r,ispin), &
!$OMP&reduction(+:Etmp)
    do lio = 1 , no_lo

      ! obtain the global index of the orbital.
      io = index_local_to_global(dit,lio,Node)
      ir = orb_type(io)
      
      ! Loop number of entries in the row... (index frame)
      do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
        
        ! as the local sparsity pattern is a super-cell pattern,
        ! we need to check the unit-cell orbital
        ! The unit-cell column index
        jo = UCORB(l_col(ind),no_u)
        jr = orb_type(jo)
        
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

        ! Ebs
        if ( spin%SO ) then
          Etmp(1,r) = Etmp(1,r) + H(ind,1) * Dscf(ind,1) &
              + H(ind,2) * Dscf(ind,2) &
              + H(ind,3) * Dscf(ind,7) &
              + H(ind,4) * Dscf(ind,8) &
              - H(ind,5) * Dscf(ind,5) &
              - H(ind,6) * Dscf(ind,6) &
              + H(ind,7) * Dscf(ind,3) &
              + H(ind,8) * Dscf(ind,4)
        else if ( spin%NCol ) then
          Etmp(1,r) = Etmp(1,r) + H(ind,1) * Dscf(ind,1) &
              + H(ind,2) * Dscf(ind,2) &
              + (H(ind,3) * Dscf(ind,3) &
              + H(ind,4) * Dscf(ind,4) ) * 2
        else
          Etmp(1,r) = Etmp(1,r) + sum(H(ind,:) * Dscf(ind,:) )
        end if

        do ispin = 1, spin%spinor
          ! Ekin
          Etmp(2,r) = Etmp(2,r) + H_kin(ind) * Dscf(ind,ispin)
          ! Enl
          Etmp(3,r) = Etmp(3,r) + H_vkb(ind) * Dscf(ind,ispin)
        end do

        ! Eso
        if ( spin%SO ) then
          Etmp(4,r) = Etmp(4,r) + H_so(ind,1)*Dscf(ind,7) + H_so(ind,2)*Dscf(ind,8) &
              + H_so(ind,5)*Dscf(ind,3) + H_so(ind,6)*Dscf(ind,4) &
              - H_so(ind,3)*Dscf(ind,5) - H_so(ind,4)*Dscf(ind,6)

        end if

        ! DEharr
        if ( spin%SO ) then
          Etmp(5,r) = Etmp(5,r) + H(ind,1) * ( Dscf(ind,1) - Dold(ind,1) )  &
                          + H(ind,2) * ( Dscf(ind,2) - Dold(ind,2) )  &
                          + H(ind,3) * ( Dscf(ind,7) - Dold(ind,7) )  &
                          + H(ind,4) * ( Dscf(ind,8) - Dold(ind,8) )  &
                          - H(ind,5) * ( Dscf(ind,5) - Dold(ind,5) )  &
                          - H(ind,6) * ( Dscf(ind,6) - Dold(ind,6) )  &
                          + H(ind,7) * ( Dscf(ind,3) - Dold(ind,3) )  &
                          + H(ind,8) * ( Dscf(ind,4) - Dold(ind,4) )
        else if ( spin%NCol ) then
          Etmp(5,r) = Etmp(5,r) + H(ind,1) * ( Dscf(ind,1) - Dold(ind,1) )  &
              + H(ind,2) * ( Dscf(ind,2) - Dold(ind,2) )  &
              + 2.0_dp * H(ind,3) * ( Dscf(ind,3) - Dold(ind,3) )  &
              + 2.0_dp * H(ind,4) * ( Dscf(ind,4) - Dold(ind,4) )
        else
          do ispin = 1, spin%spinor
            Etmp(5,r) = Etmp(5,r) + H(ind,ispin) * (Dscf(ind,ispin) - Dold(ind,ispin))
          end do
        end if

      end do
    end do
!$OMP end parallel do

    ! Now add the *other* contributions
    do ir = 1, N_Elec
      if ( Elecs(ir)%DM_update >= 1 ) then
        ! add cross-terms
        r = 4 + (TYP_DEVICE+ir - 1) * 2
        Etmp(:,2) = Etmp(:,2) + Etmp(:,r)
      end if
      if ( Elecs(ir)%DM_update == 2 ) then
        r = 3 + (ir - 1) * 2
        ! add diagonal electrode terms
        Etmp(:,2) = Etmp(:,2) + Etmp(:,r)
      end if
    end do

#ifdef MPI
    call MPI_AllReduce(Etmp(1,2),tmp(1,2),size(Etmp, 1), &
        MPI_Double_Precision,MPI_SUM, MPI_Comm_World,MPIerror)
    Etmp(:,2) = tmp(:,2)
#endif

    ! Copy data over
    NEGF_Ebs = Etmp(1,2)
    NEGF_Ekin = Etmp(2,2)
    NEGF_Enl = Etmp(3,2)
    NEGF_Eso = Etmp(4,2)
    NEGF_DEharr = Etmp(5,2)

  end subroutine ts_compute_energies

end module ts_energies_m

  
