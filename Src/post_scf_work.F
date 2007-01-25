      subroutine post_scf_work( first, last, iscf )
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
          ! If converged, make last iteration to find forces
          call setup_hamiltonian( first, last, iscf )

          ! Add on force field contribution if required
          call twobody(na_u,xa,isa,ucell,Emm,ifa,fa,istr,stress)

          ! Print energies
          DEna = Enascf - Enaatm
          Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec + Emad+Emm+
     .         Emeta
          Eharrs = Etot + DEharr
          FreeE  = Etot - Temp * Entropy

          ! If last iteration, exit SCF loop 
          do ispin = 1,nspin
            do io = 1,nh
              Dscf(io,ispin) = Dold(io,ispin)
              Escf(io,ispin) = Eold(io,ispin)
            enddo
          enddo
          if (dumpcharge) then
            call plcharge( no_s, na_s, no_u, maxnh, maxna, nspin,
     .                      isa, iphorb, indxuo, lasto,
     .                      scell, nsc, xa, rmaxo, datm )
          endif
          if (cml_p) call cmlEndStep(mainXML)
!--------------------------------------------------------------------------- END
      END subroutine post_scf_work
