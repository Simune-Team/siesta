      subroutine setup_hamiltonian( first, last, iscf )
      USE SIESTA_TODO
      implicit none
      logical :: first, last
      integer :: iscf

      external :: automatic_cell,
     .  bonds, cgvc, cgvc_zmatrix, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, iozm,
     .  kgrid, kgridinit, kinefsm, mulliken, naefs, neighb,
     .  pulayx, 
     .  reinit, shaper, spnvec, 
     .  timer, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion, eggbox


!------------------------------------------------------------------------- BEGIN
            if (cml_p) then
              call cmlStartStep(xf=mainXML, type='SCF', index=iscf)
            endif

            ! Normalize density matrix to exact charge
            qsol = 0.0_dp
            do ispin = 1,min(nspin,2)
              do io = 1,nh
                qsol = qsol + Dscf(io,ispin) * s(io)
              enddo
            enddo
#ifdef MPI
            call globalize_sum(qsol,buffer1)
            qsol = buffer1
#endif
        if (IOnode) then
          if (.not.first .and.
     .       abs(qsol/qtot-1._dp).gt.1.d-2) write(6,'(a,2f15.6)')
     .      'siesta: WARNING: Qtot, Tr[D*S] =', qtot, qsol
        endif
        do ispin = 1,nspin
          do io = 1,nh
            Dscf(io,ispin) = Dscf(io,ispin) * qtot/qsol
            Escf(io,ispin) = Escf(io,ispin) * qtot/qsol
          enddo
        enddo
! ..................

! Initialize Hamiltonian ........................................
        H = 0.0_dp

! Initialize forces and stress ...................
        if (first.or.last) then
          fa(1:3,1:na_u) = 0.0_dp
          fal(1:3,1:na_u) = 0.0_dp
          stress(1:3,1:3) = 0.0_dp
          stressl(1:3,1:3) = 0.0_dp
        endif
! ..................

! Self-energy of isolated ions ........................................
        if (first) then
          Eions = 0.0_dp
          do ia = 1,na_u
            is = isa(ia)
            Eions = Eions + uion(is)
          enddo
        endif
! ..................

! Neutral-atom: energy, forces and stress ............................
! First time for energy, last time for forces
        if (first.or.last) then
          call naefs(na_u, na_s, scell, xa, indxua, rmaxv,
     .               maxna, isa, jna, xij, r2ij,
     .               Ena, fa, stress)
          call dnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     .               maxna, isa, jna, xij, r2ij,
     .               DEna, fa, stress) 
          Ena = Ena + DEna
        endif
! ..................

! Kinetic: energy, forces, stress and matrix elements .................
        if (first.or.last) then
          call kinefsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo, no_l,
     .                 maxna, maxnh, maxnh, lasto, iphorb, isa, 
     .                 numh, listhptr, listh, numh, listhptr, listh, 
     .                 min(nspin,2), Dscf, jna, xij, r2ij,
     .                 Ekin, fal, stressl, H ) 
#ifdef MPI
! Global reduction of energy terms
          call globalize_sum(Ekin,buffer1)
          Ekin = buffer1
#endif
        endif
! ..................

! Non-local-pseudop: energy, forces, stress and matrix elements .......
        if (first.or.last) then
          call nlefsm(scell, na_u, na_s, isa, xa, indxua, maxna,
     .                maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                min(nspin,2), Dscf, Enl, fal, stressl, H)
        
! Check whether maxna has been increased during nlefsm & resize if needed
          if (maxna.gt.size(jna)) then
            call memory('D','I',size(jna),'siesta')
            deallocate(jna)
            call memory('D','D',size(r2ij),'siesta')
            deallocate(r2ij)
            call memory('D','D',size(xij),'siesta')
            deallocate(xij)
            allocate(jna(maxna))
            call memory('A','I',maxna,'siesta')
            allocate(r2ij(maxna))
            call memory('A','D',maxna,'siesta')
            allocate(xij(3,maxna))
            call memory('A','D',3*maxna,'siesta')
          endif

#ifdef MPI
! Global reduction of energy terms
          call globalize_sum(Enl,buffer1)
          Enl = buffer1
#endif
        endif
! ..................

! Save or get partial Hamiltonian (non-SCF part) ......................
        if (first.or.last) then
          do io = 1,nh
            H0(io) = H(io,1)
          enddo
        else
          do ispin = 1,nspin
            if (ispin .le. 2) then
              do io = 1,nh
                H(io,ispin) = H0(io)
              enddo
            else
              do io = 1,nh
                H(io,ispin) = 0.0_dp
              enddo
            endif
          enddo          
        endif
! ..................

! Non-SCF part of total energy .......................................
        if (first.or.last) then
          E0 = -Eions + Ena + Ekin + Enl
        else
          E0 = 0.0_dp
          do ispin = 1,min(nspin,2)
            do io = 1,nh
              E0 = E0 + H0(io) * Dscf(io,ispin)
            enddo
          enddo
#ifdef MPI
! Global reduction of E0
          call globalize_sum(E0,buffer1)
          E0 = buffer1
#endif
          E0 = E0 - Eions + Ena
        endif
! ..................

! Non-local-pseudop: energy, forces, stress and matrix elements .......
! Add SCF contribution to energy and matrix elements ..................
        g2max = g2cut
        if (last) then
! Last call to dhscf and grid-cell sampling if requested
          ifa  = 1
          istr = 1
          call grdsam( nspin, no_s, iaorb, iphorb, 
     .                 no_l, no_u, na_u, na_s, isa, xa, indxua,
     .                 ucell, mscell, g2max, ntm, ifa, istr, maxnh,
     .                 maxnh, numh, listhptr, listh, Dscf, Datm, H,
     .                 Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                 Exc, Dxc, dipol, fa, stress, fal, stressl)
        else
          ifa  = 0
          istr = 0
          ihmat = 1
          call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .                no_u, na_u, na_s, isa, xa, indxua, 
     .                ucell, mscell, g2max, ntm,
     .                ifa, istr, ihmat, ' ', ' ', ' ', ' ', ' ', ' ',
     .                maxnh, numh, listhptr, listh, Dscf, Datm,
     .                maxnh, numh, listhptr, listh, H,
     .                Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                Exc, Dxc, dipol, fa, stress, fal, stressl)
        endif
            
! Output memory use after first call to dhscf
        if (istp.eq.1 .and. iscf.eq.1) call printmemory( 6, 0 )

*       if (istp.eq.1 .and. iscf.eq.1) write(6,'(/,a,f10.3,a)')
*    .    'siesta: dhscf mesh cutoff =', g2max, ' Ry'

! ..................

! Orthonormalization forces ...........................................
        if (last) then
           call overfsm(na_u, na_s, no_s, scell, xa, indxua,rmaxo,no_l,
     .                 maxna, maxnh, maxnh, lasto, iphorb, isa, 
     .                 numh, listhptr, listh, numh, listhptr, listh, 
     .                 min(nspin,2), Escf, jna, xij, r2ij,
     .                 fal, stressl, S ) 
        endif
! ..................

! Metadynamics forces
        if (lMetaForce.and.(first.or.last)) then
          call meta(xa,na_u,ucell,Emeta,fa,stress,last,last)
        endif

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
      END subroutine setup_hamiltonian
