! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_setup_H0
      private
      public :: setup_H0
      CONTAINS

      subroutine setup_H0(G2max)
      
C     Computes non-self-consistent part of the Hamiltonian
C     and initializes data structures on the grid.
      
      USE siesta_options,  only: g2cut
      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: H_so_2D
      use sparse_matrices, only: Dscf

! CC RC  Added for the offSpOrb
      use sparse_matrices, only: H0_offsiteSO, listhptr, numh ! , listht
      use m_nlefsm,        only: nlefsm_offsiteSO
      use m_energies,      only: Enl_offsiteSO, Enl, Ekin
      use m_spin,          only: spin
      use parallel,        only: IONode
      use atomlist,        only: no_u
! CC RC  Added for the offSpOrb

      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, indxuo, datm, 
     &                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     &                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm
      use m_kinefsm,    only: kinefsm
      use m_naefs,      only: naefs
      use m_dnaefs,     only: dnaefs
      use m_dhscf,      only: dhscf_init
      use m_energies,   only: Eions, Ena, DEna, Emm, Emeta, Eso
      use m_ntm
      use m_spin,       only: spin
      use spinorbit,    only: spinorb
      use alloc, only: re_alloc, de_alloc
      use class_dSpData1D, only: val
      use class_dSpData2D, only: val

#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      implicit none
      real(dp), intent(inout) :: g2max
      
      real(dp) :: dummy_stress(3,3), dummy_fa(1,1), dummy_dm(1,1)
      real(dp) :: dummy_E
      integer  :: ia, is

! CC RC  Added for the offSpOrb
      real(dp)    :: dummy_E_offsiteSO
      real(dp) :: eV = 1._dp / 13.60580_dp
      character(len=1) :: cindent = ' '
      integer     :: io, ispin, i, j !, ind, iot
!      integer     :: ij, ji, js
      complex(dp) :: Hc, Dc
#ifdef MPI
      real(dp)            :: buffer1
#endif
! CC RC  Added for the offSpOrb

      real(dp), pointer :: H_val(:), H_so(:,:)

      real(dp) :: err_H0, err_H0_offsiteSO(6)

#ifdef DEBUG
      call write_debug( '    PRE setup_H0' )
#endif

!----------------------------------------------------------------------BEGIN
      call timer('Setup_H0',1)

C     Self-energy of isolated ions
      Eions = 0.0_dp
      do ia = 1,na_u
        is = isa(ia)
        Eions = Eions + uion(is)
      enddo

!     In these routines, add a flag to tell them NOT to compute
!     forces and stresses in this first pass, only energies.

!     Neutral-atom: energy

      if( IONode .and. spin%deb_offSO ) write(spin%iout_offsiteSO,'(a)')
     &   '    setup_H0: Calling naefs & dnaefs...'
      call naefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &           isa, Ena, dummy_fa, dummy_stress,
     &           forces_and_stress=.false.)
      call dnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, DEna, dummy_fa, dummy_stress,
     &            forces_and_stress=.false.) 
      Ena = Ena + DEna
 
C     Metadynamics energy
      if (lMetaForce) then
        call meta(xa,na_u,ucell,Emeta,dummy_fa,dummy_stress,
     $           .false.,.false.)
      endif

C     Add on force field contribution to energy
      call twobody( na_u,xa,isa,ucell,Emm,
     &              ifa=0,fa=dummy_fa,istr=0,stress=dummy_stress)

!
!     Now we compute matrix elements of the Kinetic and Non-local
!     parts of H

!     Kinetic: matrix elements only
      H_val => val(H_kin_1D)
!$OMP parallel workshare default(shared)
      H_val(:) = 0.0_dp
!$OMP end parallel workshare
      if( IONode .and. spin%deb_offSO ) write(spin%iout_offsiteSO,'(a)')
     &  '    setup_H0: Calling kinefsm...'
      call kinefsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &             maxnh, maxnh, lasto, iphorb, isa, 
     &             numh, listhptr, listh, numh, listhptr, listh, 
     &             1,
     &             dummy_dm, dummy_E, dummy_fa, dummy_stress,
     &             H_val,
     &             matrix_elements_only=.true.) 

!     Non-local-pseudop:  matrix elements only
      H_val => val(H_vkb_1D)
!$OMP parallel workshare default(shared)
      H_val(:) = 0.0_dp
!$OMP end parallel workshare


      if ( .not.spin%SO_offsite ) then
       if(IONode .and. spin%deb_offSO) write(spin%iout_offsiteSO,'(a)')
     &  '    setup_H0: Calling nlefsm...'
       call nlefsm(scell, na_u, na_s, isa, xa, indxua, 
     &             maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     &             numh, listhptr, listh, numh, listhptr, listh, 
     &             1,
     &             dummy_dm, dummy_E, dummy_fa, dummy_stress,
     &             H_val,
     &             matrix_elements_only=.true.) 
      elseif ( spin%SO_offsite ) then
        H0_offsiteSO=dcmplx(0.0d0,0.0d0)
        if(IONode .and. spin%deb_offSO) write(spin%iout_offsiteSO,'(a)')
     &   '    setup_H0: Calling nlefsm_offsiteSO...'
        call nlefsm_offsiteSO(scell, na_u, na_s, isa, xa, indxua,
     &                 maxnh, maxnh, lasto, lastkb, iphorb, iphKB,
     &                 numh, listhptr, listh, numh, listhptr, listh,
     &                 spin%Grid,
     &                 dummy_E, dummy_E_offsiteSO, dummy_fa,
     &                 dummy_stress, H_val,
     &                 matrix_elements_only=.true.)


        Enl_offsiteSO = 0.0d0
        do io = 1, maxnh

!-------- Enl_offsiteSO(u,u)
          Dc = cmplx(Dscf(io,1),Dscf(io,5))
          Hc = H0_offsiteSO(io,1)
          Enl_offsiteSO = Enl_offsiteSO + real( Hc*Dc )
!-------- Enl_offsiteSO(d,d)
          Dc = cmplx(Dscf(io,2),Dscf(io,6))
          Hc = H0_offsiteSO(io,2)
          Enl_offsiteSO = Enl_offsiteSO + real( Hc*Dc )
!-------- Enl_offsiteSO(u,d)
          Dc = cmplx(Dscf(io,3),Dscf(io,4))
          Hc = H0_offsiteSO(io,4)
          Enl_offsiteSO = Enl_offsiteSO + real( Hc*Dc )
!-------- Enl_offsiteSO(d,u)
          Dc = cmplx(Dscf(io,7), -Dscf(io,8))
          Hc = H0_offsiteSO(io,3)
          Enl_offsiteSO = Enl_offsiteSO + real( Hc*Dc )

        enddo

#ifdef MPI
! Global reduction of Ekin, Enl and Eso 
      call globalize_sum(Enl_offsiteSO,buffer1)
      Enl_offsiteSO = buffer1
#endif

        if ( IONode .and. spin%deb_offSO ) 
     &   write(spin%iout_offsiteSO,'(a,f14.10)')
     &   '    setup_H0: Enl_offsiteSO = ', Enl_offsiteSO

      endif

!CC RC
c---- check symmetries for H0
!      err_H0 = 0.0d0; err_H0_offsiteSO(1:6) = 0.0d0
!      do i = 1, no_u
!       do js = 1, numh(i)
!        ij = listhptr(i)+js
!        j = listh(ij)
!        ji = listhptr(j)+i
!        if ( ji.lt.ij ) then
!         err_H0 = max( err_H0, abs(H_val(ij)-H_val(ji)) )
!         H_val(ij) = 0.5d0*( H_val(ij) + H_val(ji) )
!         H_val(ji) = H_val(ij)
!        endif
!        if ( spin%SO_offsite ) then
!
!c------- hermiticity
!         if ( ji.lt.ij ) then
!          err_H0_offsiteSO(1) = max( err_H0_offsiteSO(1),
!     &                   cdabs(H0_offsiteSO(ij,1)-conjg(H0_offsiteSO(ji,1))) )
!          err_H0_offsiteSO(2) = max( err_H0_offsiteSO(2),
!     &                   cdabs(H0_offsiteSO(ij,2)-conjg(H0_offsiteSO(ji,2))) )
!         endif
!         err_H0_offsiteSO(3) = max( err_H0_offsiteSO(3),
!     &                  cdabs(H0_offsiteSO(ij,3)-conjg(H0_offsiteSO(ji,4))) )
!         err_H0_offsiteSO(4) = max( err_H0_offsiteSO(4),
!     &                   cdabs(H0_offsiteSO(ij,4)-conjg(H0_offsiteSO(ji,3))) )
!
!c------- spin-box anti-hermitian 
!         err_H0_offsiteSO(5) = max( err_H0_offsiteSO(5),
!     &                  cdabs(H0_offsiteSO(ij,3)+conjg(H0_offsiteSO(ij,4))) )
!         err_H0_offsiteSO(6) = max( err_H0_offsiteSO(6),
!     &                  cdabs(H0_offsiteSO(ij,4)+conjg(H0_offsiteSO(ij,3))) )
!        endif
!       enddo
!      enddo
!
!c---- Print errors
!      write(6,'(a,f8.5)') trim(cindent)//
!     & '    : H0 sym error [eV]=',err_H0/eV
!      if ( spin%SO_offsite ) then
!       write(6,'(4x,a,7f8.5)')trim(cindent)//
!     & ': H0_offsiteSO max/sym error  =',maxval(cdabs(H0_offsiteSO))/eV,err_H0_offsiteSO/eV
!       if ( (err_H0+maxval(err_H0_offsiteSO,1))/eV .gt. 1.0d-4 ) then
!        write(6,'(a)') 'siesta_H0: ERROR, H0 seems not hermitian!'
!        stop
!       endif
!      endif
!CC RC

! ..................

! If in the future the spin-orbit routine is able to compute
! forces and stresses, then "last" will be needed. If we are not
! computing forces and stresses, calling it in the first iteration
! should be enough
!
      if ( spin%SO .and. .not.spin%SO_offsite ) then
         H_so => val(H_so_2D)
!$OMP parallel workshare default(shared)
         H_so = 0._dp
!$OMP end parallel workshare
         call spinorb(no_u,no_l,iaorb,iphorb,isa,indxuo,
     &        maxnh,numh,listhptr,listh,Dscf,H_so,Eso)
      else
         Eso = 0._dp
      end if

      if ( IONode .and. spin%deb_offSO ) 
     &  write(spin%iout_offsiteSO,'(a)') ' setup_H0: Before dhscf_init'

C     This will take care of possible changes to the mesh and atomic-related
C     mesh structures for geometry changes
      g2max = g2cut
      call dhscf_init( spin%Grid, no_s, iaorb, iphorb,
     &                 no_l, no_u, na_u, na_s,
     &                 isa, xa, indxua, ucell,
     &                 mscell, G2max, ntm,
     &                 maxnh, numh, listhptr, listh, datm,
     &                 dummy_fa, dummy_stress)

      if ( IONode .and. spin%deb_offSO ) 
     &   write(spin%iout_offsiteSO,'(a)') ' setup_H0: After dhscf_init'
  
      call timer('Setup_H0',2)

#ifdef DEBUG
      call write_debug( '    POS setup_H0' )
#endif

!---------------------------------------------------------------------- END
      END subroutine setup_H0
      END module m_setup_H0
