      module m_final_H_f_stress
      private
      public :: final_H_f_stress
      CONTAINS

      subroutine final_H_f_stress( iscf )

      USE siesta_options, only: g2cut, savehs, temp
      USE siesta_options, only: recompute_H_after_scf
      use sparse_matrices, only: numh, listh, listhptr
      use sparse_matrices, only: S, Dscf, Escf, maxnh, xijo
      use sparse_matrices, only: H, H0, H_kin, H_vkb ! CC RC  Added 

      use siesta_geom
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     &                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     &                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm
      use m_overfsm,    only: overfsm
      use m_kinefsm,    only: kinefsm
      use m_naefs,      only: naefs
      use m_dnaefs,     only: dnaefs
      use m_grdsam,     only: grdsam
      use m_stress
      use m_energies
      use m_steps, only: istp
      use m_ntm
      use m_spin,         only: nspin, SpOrb ! CC RC  Added for osSO
      use m_spin,         only: h_spin_dim, spinor_dim
      use spinorbit,      only: spinorb
      use m_dipol
      use m_forces,          only: fa
      use alloc, only: re_alloc, de_alloc
      use m_gamma
      use m_hsx, only: write_hsx
      use sys, only: die
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
      use parallel, only: IOnode
#ifdef TRANSIESTA
! TSS Begin
      use files,        only : slabel, label_length
      use siesta_options, only : idyn, ia1
      use m_ts_global_vars, only : ts_istep
      use m_ts_options, only : savetshs
      use m_ts_kpoints, only : ts_gamma_scf, ts_kscell, ts_kdispl
      use m_ts_io,      only : ts_iohs
! TSS End
#endif /* TRANSIESTA */
      implicit none

      integer,   intent(in) :: iscf
      real(dp)              :: stressl(3,3)
      real(dp), pointer     :: fal(:,:)   ! Local-node part of atomic F
      integer               :: io, is
      integer               :: ifa     ! Calculate forces? 0=>no, 1=>yes
      integer               :: istr    ! Calculate stress? 0=>no, 1=>yes
      real(dp)              :: g2max
      ! To avoid overwriting the current Hamiltonian
      real(dp), pointer     :: H_tmp(:,:) => null()  
#ifdef MPI
      real(dp)              :: stresstmp(3,3)
      real(dp), pointer     :: fatmp(:,:)
      real(dp)              :: buffer1
#endif
#ifdef TRANSIESTA
!     TSS Begin
      character(len=label_length+13) :: fname
      character(8)                   :: numstring
      character(8)                   :: iotask
!     SIESTA External
      character(len=label_length+5)  :: PASTE
!     TSS End
#endif /* TRANSIESTA */
!--------------------------------------------------------------------- BEGIN

!     Initialize Hamiltonian ........................................
      call re_alloc(H_tmp, 1,maxnh, 1,nspin,
     $                       'H_tmp','final_H_f_stress')
      H_tmp = 0.0_dp

!     Initialize forces and stress ...................
      nullify(fal)
      call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'final_H_f_stress' )

      fa(1:3,1:na_u)   = 0.0_dp
      fal(1:3,1:na_u)  = 0.0_dp
      stress(1:3,1:3)  = 0.0_dp
      stressl(1:3,1:3) = 0.0_dp  ! fal already initialized

! Neutral-atom: energy, forces and stress ............................
!      write(6,'(a)') '----final_H_f_stress: Calling naefs & dnaefs '
      call naefs( na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, Ena, fa, stress, forces_and_stress=.true. )
      call dnaefs( na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, DEna, fa, stress, forces_and_stress=.true. )
      Ena = Ena + DEna

!      write(6,'(a)') '----final_H_f_stress: Calling kinefsm '
!     Kinetic: energy, forces, stress and matrix elements .................
      call kinefsm( na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &              maxnh, maxnh, lasto, iphorb, isa, 
     &              numh, listhptr, listh, numh, listhptr, listh, 
     &              spinor_dim, Dscf, Ekin, fal, stressl, H_tmp,
     &              matrix_elements_only=.false. ) 
#ifdef MPI
! Global reduction of energy terms
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
#endif
! ..................

!      write(6,'(a)') '----final_H_f_stress: Calling nlefsm '
! Non-local-pseudop: energy, forces, stress and matrix elements .......
      call nlefsm( scell, na_u, na_s, isa, xa, indxua, 
     &             maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     &             numh, listhptr, listh, numh, listhptr, listh, 
     &             spinor_dim, Dscf, Enl, fal, stressl, H_tmp,
     &             matrix_elements_only=.false. ) 

#ifdef MPI
! Global reduction of energy terms
      call globalize_sum(Enl,buffer1)
      Enl = buffer1
#endif
! ..................
! CC RC  Calculation of the osSO contribution to the last H call

! If in the future the spin-orbit routine is able to compute
! forces and stresses, then "last" will be needed. If we are not
! computing forces and stresses, calling it in the first iteration
! should be enough
!
!      write(6,'(a)') '----final_H_f_stress: Calling spinorbit '
       if (SpOrb) then
          call spinorb(no_u,no_l,iaorb,iphorb,isa,indxuo,
     .                 maxnh,numh,listhptr,listh,Dscf,H,Eso)
#ifdef MPI
! Global reduction of energy terms
          call globalize_sum(Eso,buffer1)
          Eso = buffer1
#endif
        endif

! Save or get partial Hamiltonian (non-SCF part) ......................
          if (SpOrb) then
              H0(:,:) = H_kin(:,:) + H_vkb(:,:) + H(:,:)
           else
              H0(:,1) = H_kin(:,1) + H_vkb(:,1) 
           endif
! ..................


!     Non-SCF part of total energy
      E0 = -Eions + Ena + Ekin + Enl + Eso ! CC RC Added Eso

!     Non-local-pseudop: energy, forces, stress and matrix elements
!     Add SCF contribution to energy and matrix elements 
      g2max = g2cut
!     Last call to dhscf and grid-cell sampling if requested
      ifa  = 1
      istr = 1

      ! This will call dhscf with the final DM coming out of
      ! the scf cycle

      call grdsam( nspin, h_spin_dim, SpOrb, no_s, iaorb, iphorb, 
     &             no_l, no_u, na_u, na_s, isa, xa, indxua,
     &             ucell, ntm, ifa, istr, maxnh,
     &             maxnh, numh, listhptr, listh, Dscf, Datm, H_tmp,
     &             Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     &             Exc, Dxc, dipol, fa, stress, fal, stressl )

!     Orthonormalization forces ...........................................
!     **Note that this will re-calculate S also, but should get 
!     the same result
      call overfsm( na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &              maxnh, maxnh, lasto, iphorb, isa, 
     &              numh, listhptr, listh, numh, listhptr, listh, 
     &              spinor_dim, Escf, fal, stressl, S ) 

!     Metadynamics forces
      if (lMetaForce) then
        call meta(xa,na_u,ucell,Emeta,fa,stress,.true.,.true.)
      endif

!     Add on force field contribution if required
      call twobody( na_u, xa, isa, ucell, Emm,
     &              ifa=1, fa=fa, istr=1, stress=stress )

#ifdef MPI
!     Global reduction of forces and stresses
      nullify(fatmp)
      call re_alloc( fatmp, 1, 3, 1, na_u, 'fatmp', 'final_H_f_stress' )
      call globalize_sum(stressl(1:3,1:3),stresstmp(1:3,1:3))
      call globalize_sum( fal(1:3,1:na_u), fatmp(1:3,1:na_u) )
      stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
      call de_alloc( fatmp, 'fatmp', 'final_H_f_stress' )
#else
      stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif

      ! If backward compatibility is requested, recover
      ! the old behavior with respect to H

      if (recompute_H_after_scf) then
         if (ionode) then
            write(6,"(a)") ":!: Updating H after scf cycle" //
     $                     " as requested by compat. option"
         endif
         H = H_tmp
      endif

      call de_alloc( fal, 'fal', 'final_H_f_stress' )
      call de_alloc( H_tmp, 'H_tmp', 'final_H_f_stress' )

! Save Hamiltonian and overlap matrices ............................
! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file
!
! Note that we save the Hamiltonian coming out of the scf cycle,
! not the one computed (from DM_out) in this routine.
! This call could be moved to a more appropriate place
!
      if (savehs) then
         call write_hsx( gamma, no_u, no_s, h_spin_dim, indxuo,
     &        maxnh, numh, listhptr, listh, H, S, qtot,
     &        temp, xijo)
      endif
#ifdef TRANSIESTA
! TSS Begin
      if (savetshs) then
         if (idyn.eq.6) then ! Force constants run
            write(numstring,'(2i4)') ia1, ts_istep
            fname = paste( slabel, numstring )
            fname = paste( fname, '.TSHS' )
            iotask='write'
            call ts_iohs( iotask, gamma, .false., no_u, no_s, nspin,
     &                    indxuo, maxnh, numh, listhptr, listh, H, S, 
     &                    qtot, temp, xijo, label_length+5, fname, na_u,
     &                    lasto, isa, ef, ucell, ts_kscell, ts_kdispl,
     &                    ts_gamma_scf, xa, ts_istep, ia1)
         else
            fname  = paste( slabel, '.TSHS' )
            iotask ='write'
            call ts_iohs( iotask, gamma, .false., no_u, no_s, nspin,
     &                    indxuo, maxnh, numh, listhptr, listh, H, S, 
     &                    qtot, temp, xijo, label_length+5, fname, na_u,
     &                    lasto, isa, ef, ucell,ts_kscell, ts_kdispl,
     &                    ts_gamma_scf, xa, 0, 0 )
         endif
      endif
! TSS End
#endif /* TRANSIESTA */
!----------------------------------------------------------------------- END
      END subroutine final_H_f_stress
      END module m_final_H_f_stress
