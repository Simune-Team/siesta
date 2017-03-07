! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_compute_energies
  implicit none
  public :: compute_energies
CONTAINS
  subroutine  compute_energies(iscf)

    ! This routine computes the Harris energy from E_KS(DM_in):
    !
    !    E_Harris = E_KS(DM_in) + Tr[H_in*(DM_out-DM_in)]
    !
    ! and possibly E_KS(DM_out) as
    !
    !    E_KS(DM_out) = Tr[H_in*DM_out] + E_HXC(DM_out)
    !
    ! Note that E_KS(DM_in), as computed in setup_hamiltonian, is not
    ! variational if mixing the DM, as the kinetic energy term is
    ! using a DM (DM_in) which is not derived from wave-functions. It
    ! is also not consistent if we are using the charge density as
    ! primary variable for mixing. In this case the Enl and Ekin parts
    ! are computed with the previous step's DM, whereas the "scf" part
    ! comes from the (mixed) rho.
    !
    ! E_KS(DM_out) computed in setup_hamiltonian with DM_out is
    ! variational. Rather than calling again setup_hamiltonian, it is
    ! enough to compute the "scf" part of the energy by calling dhscf
    ! with DM_out.

    ! The routine checks the mixing type and proceeds accordingly

      use precision,       only: dp
      use fdf,             only: fdf_get
      use siesta_options,  only: g2cut, Temp
      use siesta_options,  only: mix_charge, mixH
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H
      use sparse_matrices, only: Dscf, Dold
      use m_dhscf,         only: dhscf
      use m_energies
      use atomlist,        only: no_u, iaorb, iphkb, qtot, indxuo, datm,   &
                                 lastkb, no_s, rmaxv, indxua, iphorb, lasto, &
                                 rmaxo, no_l
      use m_ntm,           only: ntm

! CC RC  Added for the offSpOrb
      use m_spin,          only: spin 
      use sparse_matrices, only: H0_SO
      use sparse_matrices, only: listht
      use parallel,        only: IONode 
!      use m_spin,          only: NoMagn, SPpol, NonCol, SpOrb 
!      use m_spin,          only: nspin, h_spin_dim, spinor_dim
! CC RC  Added for the offSpOrb

      use m_dipol,         only: dipol
      use siesta_geom,     only: na_u, na_s, xa, isa
      use m_rhog,          only: rhog
#ifdef MPI
      use m_mpi_utils,     only: globalize_sum
#endif

      integer, intent(in)   :: iscf

      integer               :: ihmat, ifa, istr, ispin, io
#ifdef MPI
      real(dp) :: buffer1
#endif

      real(dp) :: const
      real(dp) :: dummy_stress(3,3), dummy_fa(1,1)
      real(dp) :: dummy_E, g2max, dummy_H(1,1)
      logical  :: mixDM

! CC RC  Added for the offSpOrb
      real(dp)    :: Ebs_tmp(4), DEharr_tmp(4)
!      integer     :: iot 

      complex(dp):: Hc, Dc, Ebs_Daux(2,2), Ebs_Haux(2,2)
      complex(dp):: DEharr_Daux(2,2), DEharr_Haux(2,2), DEharr_Daux_old(2,2)
! CC RC  Added for the offSpOrb

      mixDM = (.not. (mixH .or. mix_charge))

!     Compute the band-structure energy

      call compute_EBS()

!      stop'Stopping after compute_EBS...'
    
      ! These energies were calculated in the latest call to
      ! setup_hamiltonian, using as ingredient D_in 

      ! Ecorrec comes from O(N)...
      ! DUext is the energy of the charge from DM_in in a field.
      ! Emad, Emm, Emeta are extra terms that are added for
      ! consistency of the total energy.

      ! E0 = Ena + Ekin + Enl + Eso - Eions

      call update_DEna()
      call update_Etot()


! Harris energy
! It does not make sense if mixing the charge density, as there is no DM_in
! If mixing the Hamiltonian its usefulness is questionable also.

      if (mix_charge) then   ! possibly add mixH here
         EHarrs = 0.0_dp
      else
         call compute_DEharr()
         Eharrs = Etot + DEharr
!         write(6,'(2(a,f14.8))') ' DEharr = ', DEharr
      endif

! Possible correction to Etot if mixing the DM. This is purely
! cosmetic, to show uniformly converged values during the scf
! cycle. The final energy will be completely correct if the DM is not
! mixed after scf convergence.

! The correction is mandatory if mixing the charge density. In this
! case the call to dhscf is needed to generate rho(G)_out

! If mixing H the KS energy is already variational, as it is
! computed with DM_out

      if (mixDM) then
         if (fdf_get("SCF.Want.Variational.EKS",.false.)) then
            call compute_correct_EKS()
         endif
      else if (mixH) then
         ! not needed
      else if (mix_charge) then
         call compute_correct_EKS()
      endif

      call update_FreeE(Temp)

 CONTAINS

   subroutine compute_EBS()

      integer :: i, j, ind

      Ebs = 0.0_dp
      Ebs_tmp = 0.0_dp  ! CC RC  Added for the offSpOrb
!       const factor takes into account that there are two nondiagonal
!       elements in non-collinear spin density matrix, stored as
!       ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
!       ispin=4 => Imag(D12)

      if ( spin%SO_off ) then

        Ebs_Daux = dcmplx(0.0_dp, 0.0_dp)
        Ebs_Haux = dcmplx(0.0_dp, 0.0_dp)

        do io = 1, maxnh

!         iot = listht(io)

         Ebs_Haux(1,1) = dcmplx(H(io,1), H(io,5))  ! Huu_ij
         Ebs_Haux(2,2) = dcmplx(H(io,2), H(io,6))  ! Hdd_ij
         Ebs_Haux(1,2) = dcmplx(H(io,3), H(io,4))  ! Hud_ij
         Ebs_Haux(2,1) = dcmplx(H(io,3),-H(io,4))  ! Hdu_ij

         Ebs_Daux(1,1) = dcmplx(Dscf(io,1),-Dscf(io,5))  ! Duu_ji
         Ebs_Daux(2,2) = dcmplx(Dscf(io,2),-Dscf(io,6))  ! Ddd_ji
         Ebs_Daux(1,2) = dcmplx(Dscf(io,7),-Dscf(io,8))  ! Dud_ji
         Ebs_Daux(2,1) = dcmplx(Dscf(io,3),-Dscf(io,4))  ! Ddu_ji
!         Ebs_Daux(1,1) = dcmplx(Dscf(iot,1), Dscf(iot,5))  ! Duu_ji
!         Ebs_Daux(2,2) = dcmplx(Dscf(iot,2), Dscf(iot,6))  ! Ddd_ji
!         Ebs_Daux(1,2) = dcmplx(Dscf(iot,3), Dscf(iot,4))  ! Dud_ji
!         Ebs_Daux(2,1) = dcmplx(Dscf(iot,7), Dscf(iot,8))  ! Ddu_ji

         Ebs_tmp(1) = Ebs_tmp(1) + real( Ebs_Haux(1,1)*Ebs_Daux(1,1))
         Ebs_tmp(2) = Ebs_tmp(2) + real( Ebs_Haux(2,2)*Ebs_Daux(2,2))
         Ebs_tmp(3) = Ebs_tmp(3) + real( Ebs_Haux(1,2)*Ebs_Daux(2,1))
         Ebs_tmp(4) = Ebs_tmp(4) + real( Ebs_Haux(2,1)*Ebs_Daux(1,2))

        enddo

        Ebs = sum ( Ebs_tmp )

        if ( IONode .and. spin%deb_offSO .or. spin%deb_P ) then
         write(spin%iout_SO,'(a,f16.10)') ' compute_Ebs: Ebs    = ', Ebs
        endif

      elseif ( spin%SO .and. .not. spin%SO_off ) then
        do io = 1,maxnh
          Ebs = Ebs      + H(io,1) * ( Dscf(io,1) )   &
                         + H(io,2) * ( Dscf(io,2) )   &
                         + H(io,3) * ( Dscf(io,7) )   &
                         + H(io,4) * ( Dscf(io,8) )   &
                         - H(io,5) * ( Dscf(io,5) )   &
                         - H(io,6) * ( Dscf(io,6) )   &
                         + H(io,7) * ( Dscf(io,3) )   &
                         + H(io,8) * ( Dscf(io,4) )
        enddo
      else if ( spin%NCol ) then
        do io = 1,maxnh
          Ebs    = Ebs    + H(io,1) * ( Dscf(io,1)  )   &
                          + H(io,2) * ( Dscf(io,2)  )   &
                 + 2.0_dp * H(io,3) * ( Dscf(io,3)  )   &
                 + 2.0_dp * H(io,4) * ( Dscf(io,3)  )
        enddo
      else if ( spin%Col )  then
        do io = 1,maxnh
          Ebs    = Ebs    + H(io,1) * Dscf(io,1)  &
                          + H(io,2) * Dscf(io,2)
        enddo
      else if ( spin%none ) then
        do io = 1,maxnh
          Ebs    = Ebs + H(io,1) *  Dscf(io,1)
        enddo
      endif

#ifdef MPI
!     Global reduction 
      call globalize_sum(Ebs,buffer1)
      Ebs = buffer1
#endif
    end subroutine compute_EBS

    subroutine compute_DEharr()

      DEharr = 0.0_dp
      DEharr_tmp = 0.0_dp ! CC RC  Added for the offSpOrb
             !       const factor takes into account that there are two nondiagonal
             !       elements in non-collinear spin density matrix, stored as
             !       ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
             !       ispin=4 => Imag(D12)

      if ( spin%SO_off ) then

        DEharr_Daux     = dcmplx(0.0_dp, 0.0_dp)
        DEharr_Daux_old = dcmplx(0.0_dp, 0.0_dp)
        DEharr_Haux     = dcmplx(0.0_dp, 0.0_dp)

        do io = 1, maxnh

!          iot = listht(io)

          DEharr_Haux(1,1) = dcmplx(H(io,1), H(io,5))
          DEharr_Haux(2,2) = dcmplx(H(io,2), H(io,6))
          DEharr_Haux(1,2) = dcmplx(H(io,3), H(io,4))
          DEharr_Haux(2,1) = dcmplx(H(io,3),-H(io,4))

!          DEharr_Daux(1,1) = dcmplx(Dscf(iot,1), Dscf(iot,5))
!          DEharr_Daux(2,2) = dcmplx(Dscf(iot,2), Dscf(iot,6))
!          DEharr_Daux(1,2) = dcmplx(Dscf(iot,3), Dscf(iot,4))
!          DEharr_Daux(2,1) = dcmplx(Dscf(iot,7), Dscf(iot,8))
!
!          DEharr_Daux_old(1,1) = dcmplx(Dold(iot,1), Dold(iot,5))
!          DEharr_Daux_old(2,2) = dcmplx(Dold(iot,2), Dold(iot,6))
!          DEharr_Daux_old(1,2) = dcmplx(Dold(iot,3), Dold(iot,4))
!          DEharr_Daux_old(2,1) = dcmplx(Dold(iot,7), Dold(iot,8))

          DEharr_Daux(1,1) = dcmplx(Dscf(io,1),-Dscf(io,5))
          DEharr_Daux(2,2) = dcmplx(Dscf(io,2),-Dscf(io,6))
          DEharr_Daux(1,2) = dcmplx(Dscf(io,7),-Dscf(io,8))
          DEharr_Daux(2,1) = dcmplx(Dscf(io,3),-Dscf(io,4))

          DEharr_Daux_old(1,1) = dcmplx(Dold(io,1),-Dold(io,5))
          DEharr_Daux_old(2,2) = dcmplx(Dold(io,2),-Dold(io,6))
          DEharr_Daux_old(1,2) = dcmplx(Dold(io,7),-Dold(io,8))
          DEharr_Daux_old(2,1) = dcmplx(Dold(io,3),-Dold(io,4))

          DEharr_tmp(1) = DEharr_tmp(1) &
             + real( DEharr_Haux(1,1)*(DEharr_Daux(1,1)-DEharr_Daux_old(1,1)))
          DEharr_tmp(2) = DEharr_tmp(2) &
             + real( DEharr_Haux(2,2)*(DEharr_Daux(2,2)-DEharr_Daux_old(2,2)))
          DEharr_tmp(3) = DEharr_tmp(3) &
             + real( DEharr_Haux(1,2)*(DEharr_Daux(2,1)-DEharr_Daux_old(2,1)))
          DEharr_tmp(4) = DEharr_tmp(4) &
             + real( DEharr_Haux(2,1)*(DEharr_Daux(1,2)-DEharr_Daux_old(1,2)))

         enddo

         DEharr = sum ( DEharr_tmp )

        if ( IONode .and. spin%deb_offSO .or. spin%deb_P ) then
         write(spin%iout_SO,'(a,f16.10)') & 
            ' compute_DEharr: DEharr = ', DEharr
        endif

      elseif ( spin%SO .and. .not.spin%SO_off ) then
        do io = 1,maxnh
          DEharr = DEharr + H(io,1) * ( Dscf(io,1) - Dold(io,1) )  &
                          + H(io,2) * ( Dscf(io,2) - Dold(io,2) )  &
                          + H(io,3) * ( Dscf(io,7) - Dold(io,7) )  &
                          + H(io,4) * ( Dscf(io,8) - Dold(io,8) )  &
                          - H(io,5) * ( Dscf(io,5) - Dold(io,5) )  &
                          - H(io,6) * ( Dscf(io,6) - Dold(io,6) )  &
                          + H(io,7) * ( Dscf(io,3) - Dold(io,3) )  &
                          + H(io,8) * ( Dscf(io,4) - Dold(io,4) )
        enddo
      else if ( spin%NCol ) then
        do io = 1,maxnh
          DEharr = DEharr + H(io,1) * ( Dscf(io,1) - Dold(io,1) )  &
                          + H(io,2) * ( Dscf(io,2) - Dold(io,2) )  &
                 + 2.0_dp * H(io,3) * ( Dscf(io,3) - Dold(io,3) )  &
                 + 2.0_dp * H(io,4) * ( Dscf(io,3) - Dold(io,3) )
        enddo
      elseif ( spin%Col )  then
        do io = 1,maxnh
          DEharr = DEharr + H(io,1) * ( Dscf(io,1) - Dold(io,1) )  &
                          + H(io,2) * ( Dscf(io,2) - Dold(io,2) )
        enddo
      elseif ( spin%none ) then
        do io = 1,maxnh
          DEharr = DEharr + H(io,1) * ( Dscf(io,1) - Dold(io,1) )
        enddo
      endif

#ifdef MPI
      !     Global reduction of DEharr
      call globalize_sum( DEharr, buffer1 )
      DEharr = buffer1
#endif
    end subroutine compute_DEharr
      
    subroutine compute_correct_EKS()

      use files, only : filesOut_t    ! derived type for output file names
      use class_dSpData1D, only : val
      use class_dSpData2D, only : val
      use sparse_matrices, only: H_kin_1D, H_vkb_1D, H_so_2D

      type(filesOut_t)  :: filesOut  ! blank output file names
      real(dp), pointer :: H_vkb(:), H_kin(:), H_so(:,:)

      ! Compute E_KS(DM_out)

      g2max = g2cut
      ifa  = 0
      istr = 0
      ihmat = 0

      ! Pass DM_out to compute E_HXC(out)

      ! Remove unwanted arguments...

      call dhscf( spin%Grid, no_s, iaorb, iphorb, no_l,      &
                  no_u, na_u, na_s, isa, xa, indxua,                        &
                  ntm, ifa, istr, ihmat, filesOut,                          &
                  maxnh, numh, listhptr, listh, Dscf, Datm,                 &
                  maxnh, dummy_H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, &
                  Exc, Dxc, dipol, dummy_stress, dummy_fa, dummy_stress)
      ! (when mix_charge is true, rhog will contain the output rho(G))

      call update_DEna()

!     Compute Tr[H_0*DM_out] = Ekin + Enl + Eso with DM_out

      H_vkb => val(H_vkb_1D)
      H_kin => val(H_kin_1D)
      
      Ekin = 0.0_dp
      Enl  = 0.0_dp
      do ispin = 1, spin%spinor
         do io = 1,maxnh
            Ekin = Ekin + H_kin(io) * Dscf(io,ispin)
            Enl  = Enl  + H_vkb(io) * Dscf(io,ispin)
         enddo
      enddo

#ifdef MPI
      ! Global reduction of Ekin, Enl
      call globalize_sum( Ekin, buffer1 )
      Ekin = buffer1
      call globalize_sum( Enl, buffer1 )
      Enl = buffer1
#endif

      Eso = 0._dp
      Enl_SO = 0.0_dp
      if ( spin%SO_off ) then

        do io = 1, maxnh

!          iot = listht(io)

!-------- Enl_SO(u,u)
!          Dc = cmplx(Dscf(iot,1), Dscf(iot,5))
          Dc = cmplx(Dscf(io,1), -Dscf(io,5))
          Hc = H0_SO(io,1)
          Enl_SO = Enl_SO + real( Hc*Dc )
!-------- Enl_SO(d,d)
!          Dc = cmplx(Dscf(iot,2), Dscf(iot,6))
          Dc = cmplx(Dscf(io,2), -Dscf(io,6))
          Hc = H0_SO(io,2)
          Enl_SO = Enl_SO + real( Hc*Dc )
!-------- Enl_SO(u,d)
!          Dc = cmplx(Dscf(iot,3), Dscf(iot,4))
          Dc = cmplx(Dscf(io,7), -Dscf(io,8))
          Hc = H0_SO(io,4)
          Enl_SO = Enl_SO + real( Hc*Dc )
!-------- Enl_SO(d,u)
!          Dc = cmplx(Dscf(iot,7), Dscf(iot,8))
          Dc = cmplx(Dscf(io,3), -Dscf(io,4))
          Hc = H0_SO(io,3)
          Enl_SO = Enl_SO + real( Hc*Dc )

        enddo

      elseif ( spin%SO .and. .not.spin%SO_off ) then
         ! Sadly some compilers (g95), does
         ! not allow bounds for pointer assignments :(
         H_so => val(H_so_2D)
         do io = 1,maxnh
            Eso = Eso + H_so(io,1)*Dscf(io,7) + H_so(io,2)*Dscf(io,8) &
                 + H_so(io,5)*Dscf(io,3) + H_so(io,6)*Dscf(io,4) &
                 - H_so(io,3)*Dscf(io,5) - H_so(io,4)*Dscf(io,6)
         end do
#ifdef MPI
         ! Global reduction of Eso
         call globalize_sum( Eso, buffer1 )
         Eso = buffer1
#endif
      end if
      
      ! E0 = Ena + Ekin + Enl + Eso - Eions

      ! Clarify: Ecorrec (from O(N))
      call update_Etot()
      
    end subroutine compute_correct_EKS

  end subroutine compute_energies

end module m_compute_energies
