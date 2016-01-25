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
      use sparse_matrices, only: H_kin, H_vkb
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H
      use sparse_matrices, only: Dscf, Dold
      use m_dhscf,         only: dhscf
      use m_energies
      use atomlist,        only: no_u, iaorb, iphkb, qtot, indxuo, datm,   &
                                 lastkb, no_s, rmaxv, indxua, iphorb, lasto, &
                                 rmaxo, no_l
      use m_ntm,           only: ntm
      use m_spin,          only: nspin
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

      real(dp) :: const, Escf_out
      real(dp) :: dummy_stress(3,3), dummy_fa(1,1)
      real(dp) :: dummy_E, g2max, dummy_H(1,1)
      logical  :: mixDM

      mixDM = (.not. (mixH .or. mix_charge))

!     Compute the band-structure energy

      call compute_EBS()


      ! These energies were calculated in the latest call to
      ! setup_hamiltonian, using as ingredient D_in 

      ! Ecorrec comes from O(N)...
      ! DUext is the energy of the charge from DM_in in a field.
      ! Emad, Emm, Emeta are extra terms that are added for
      ! consistency of the total energy.

      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec + Emad + Emm + Emeta


! Harris energy
! It does not make sense if mixing the charge density, as there is no DM_in
! If mixing the Hamiltonian its usefulness is questionable also.

      if (mix_charge) then   ! possibly add mixH here
         EHarrs = 0.0_dp
      else
         call compute_DEharr()
         Eharrs = Etot + DEharr
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
         
      FreeE  = Etot - Temp * Entropy

 CONTAINS

   subroutine compute_EBS()

      Ebs = 0.0_dp
      do ispin = 1,nspin
!       const factor takes into account that there are two nondiagonal
!       elements in non-collinear spin density matrix, stored as
!       ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
!       ispin=4 => Imag(D12)
        const = 1._dp
        if (ispin .gt. 2) const = 2._dp
        do io = 1,maxnh
          Ebs = Ebs + H(io,ispin) * const * Dscf(io,ispin)
        enddo
      enddo
#ifdef MPI
!     Global reduction 
      call globalize_sum(Ebs,buffer1)
      Ebs = buffer1
#endif
    end subroutine compute_EBS

    subroutine compute_DEharr()

          DEharr = 0.0_dp
          do ispin = 1,nspin
             !       const factor takes into account that there are two nondiagonal
             !       elements in non-collinear spin density matrix, stored as
             !       ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
             !       ispin=4 => Imag(D12)
             const = 1._dp
             if (ispin .gt. 2) const = 2._dp
             do io = 1,maxnh
                DEharr = DEharr +  &
                         H(io,ispin) * const * ( Dscf(io,ispin) - Dold(io,ispin) )
             enddo
          enddo
#ifdef MPI
          !     Global reduction of DEharr
          call globalize_sum( DEharr, buffer1 )
          DEharr = buffer1
#endif
    end subroutine compute_DEharr
      
    subroutine compute_correct_EKS()

      use files, only : filesOut_t    ! derived type for output file names

      type(filesOut_t)    :: filesOut  ! blank output file names

      ! Compute E_KS(DM_out)

      g2max = g2cut
      ifa  = 0
      istr = 0
      ihmat = 0

      ! Pass DM_out to compute E_HXC(out)

      ! Remove unwanted arguments...

      call dhscf( nspin, no_s, iaorb, iphorb, no_l,                         &
                  no_u, na_u, na_s, isa, xa, indxua,                        &
                  ntm, ifa, istr, ihmat, filesOut,                          &
                  maxnh, numh, listhptr, listh, Dscf, Datm,                 &
                  maxnh, dummy_H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, &
                  Exc, Dxc, dipol, dummy_stress, dummy_fa, dummy_stress)
      ! (when mix_charge is true, rhog will contain the output rho(G))

      DEna     = Enascf - Enaatm

      ! Clarify: 
      !          DUext (external electric field) -- should it be in or out?

      Escf_out = DEna + DUscf + DUext + Exc 

!     Compute Tr[H_0*DM_out] = Ekin + Enl with DM_out

      Ekin = 0.0_dp
      Enl  = 0.0_dp
      do ispin = 1,min(nspin,2)
        do io = 1,maxnh
          Ekin = Ekin + H_kin(io) * Dscf(io,ispin)
          Enl  = Enl  + H_vkb(io)* Dscf(io,ispin)
        enddo
      enddo
#ifdef MPI
!     Global reduction of Ekin, Enl
      call globalize_sum( Ekin, buffer1 )
      Ekin = buffer1
      call globalize_sum( Enl, buffer1 )
      Enl = buffer1
#endif

      ! E0 = Ekin + Enl - Eions + Ena

      ! Clarify: Ecorrec (from O(N))
      !          
      Etot = Ekin + Enl - Eions + Ena + Escf_out + Ecorrec + Emad + Emm + Emeta
    end subroutine compute_correct_EKS

  end subroutine compute_energies

end module m_compute_energies
