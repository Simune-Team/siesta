      module m_final_H_f_stress
      private
      public :: final_H_f_stress
      CONTAINS

      subroutine final_H_f_stress( iscf )
      
      USE m_normalize_dm, only: normalize_dm
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     .                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     .                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm
      use m_overfsm,    only: overfsm
      use m_kinefsm,    only: kinefsm
      use m_naefs,      only: naefs
      use m_dnaefs,     only: dnaefs
      use siesta_cml
      use m_stress
      use m_energies
      use m_steps, only: istp
      use m_ntm
      use m_spin,         only: nspin
      use m_dipol
      use m_forces,          only: fa
      use alloc, only: re_alloc, de_alloc
      use m_gamma
      use sys, only: die
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      implicit none

      integer, intent(in) :: iscf

      real(dp)              :: stressl(3,3)
      real(dp), pointer     :: fal(:,:)   ! Local-node part of atomic forces

#ifdef MPI
      real(dp)              :: stresstmp(3,3)
      real(dp), pointer     :: fatmp(:,:)
      real(dp):: buffer1   ! Temporary for the result of a globalize operation
#endif

      integer :: ia, io, is, ispin

      integer:: ifa     ! Calculate forces? 0=>no, 1=>yes
      integer:: istr    ! Calculate stress? 0=>no, 1=>yes
      integer:: ihmat   ! Calculate hamiltonian? 0=>no, 1=>yes

      real(dp) :: g2max



!------------------------------------------------------------------------- BEGIN
       ! Do we need this here??  *****
   
      if (cml_p) then
         call cmlStartStep(xf=mainXML, type='SCF', index=iscf)
      endif

      ! Normalize density matrix to exact charge
      call normalize_dm(first=.false.)

      ! Initialize Hamiltonian ........................................
      H = 0.0_dp

! Initialize forces and stress ...................

      nullify(fal)
      call re_alloc(fal,1,3,1,na_u,name="fal",
     $              routine="setup_hamiltonian")

      fa(1:3,1:na_u) = 0.0_dp
      stress(1:3,1:3) = 0.0_dp
      stressl(1:3,1:3) = 0.0_dp  ! fal already initialized
! ..................


! Neutral-atom: energy, forces and stress ............................
          call naefs(na_u, na_s, scell, xa, indxua, rmaxv,
     .               isa, Ena, fa, stress, forces_and_stress=.true.)
          call dnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     .                isa, DEna, fa, stress, forces_and_stress=.true.) 
          Ena = Ena + DEna
! ..................

! Kinetic: energy, forces, stress and matrix elements .................
          call kinefsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo, no_l,
     .                 maxnh, maxnh, lasto, iphorb, isa, 
     .                 numh, listhptr, listh, numh, listhptr, listh, 
     .                 min(nspin,2), Dscf, Ekin, fal, stressl, H,
     .                 forces_and_stress=.true.) 
#ifdef MPI
! Global reduction of energy terms
          call globalize_sum(Ekin,buffer1)
          Ekin = buffer1
#endif
! ..................

! Non-local-pseudop: energy, forces, stress and matrix elements .......
          call nlefsm(scell, na_u, na_s, isa, xa, indxua, 
     .                maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                min(nspin,2), Dscf, Enl, fal, stressl, H,
     .                forces_and_stress=.true.)
        
#ifdef MPI
! Global reduction of energy terms
          call globalize_sum(Enl,buffer1)
          Enl = buffer1
#endif
! ..................

! Non-SCF part of total energy .......................................
          E0 = -Eions + Ena + Ekin + Enl

! Non-local-pseudop: energy, forces, stress and matrix elements .......
! Add SCF contribution to energy and matrix elements ..................
        g2max = g2cut
! Last call to dhscf and grid-cell sampling if requested
          ifa  = 1
          istr = 1
          call grdsam( nspin, no_s, iaorb, iphorb, 
     .                 no_l, no_u, na_u, na_s, isa, xa, indxua,
     .                 ucell, mscell, g2max, ntm, ifa, istr, maxnh,
     .                 maxnh, numh, listhptr, listh, Dscf, Datm, H,
     .                 Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                 Exc, Dxc, dipol, fa, stress, fal, stressl)



! Orthonormalization forces ...........................................
           call overfsm(na_u, na_s, no_s, scell, xa, indxua,rmaxo,no_l,
     .                 maxnh, maxnh, lasto, iphorb, isa, 
     .                 numh, listhptr, listh, numh, listhptr, listh, 
     .                 min(nspin,2), Escf, fal, stressl, S ) 
! ..................

        ! Metadynamics forces
        if (lMetaForce) then
          call meta(xa,na_u,ucell,Emeta,fa,stress,.true.,.true.)
        endif

        ! Add on force field contribution if required
           call twobody(na_u,xa,isa,ucell,Emm,
     $                  ifa=1,fa=fa,istr=1,stress=stress)

#ifdef MPI
           !Global reduction of forces and stresses
           nullify(fatmp)
           call re_alloc(fatmp,1,3,1,na_u,name="fatmp",
     $              routine="setup_hamiltonian")
           call globalize_sum(stressl(1:3,1:3),stresstmp(1:3,1:3))
           call globalize_sum(fal(1:3,1:na_u),fatmp(1:3,1:na_u))
           stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
           fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
           call de_alloc(fatmp,routine="setup_hamiltonian")
#else
           stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
           fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif

        call de_alloc(fal,routine="setup_hamiltonian")

! Find entropy ........................................................
! Entropy is returned from the call to diagon. To add to the energy
! the entropy computed from the input charge, here it is assigned to the one
! of the former SCF step

        Entropy = 0.0_dp
        if (isolve .eq. 0) then
          if (istp.gt.1 .or. iscf.gt.1) then
            Entropy = Entrop
          endif
        endif
        
!!!! Possible BUG: Dold should not be overwritten in the last iteration
! Save present density matrix ........................................
        do is = 1,nspin
          do io = 1,nh
            Dold(io,is) = Dscf(io,is)
            Eold(io,is) = Escf(io,is)
          enddo
        enddo        

! Save Hamiltonian and overlap matrices ............................
        if (savehs) then
          call iohs( 'write', gamma, no_u, no_s, nspin, indxuo,
     $               maxnh, numh, listhptr, listh, H, S, qtot, temp,
     $               xijo )
        endif
!--------------------------------------------------------------------------- END
      END subroutine final_H_f_stress
      END module m_final_H_f_stress
