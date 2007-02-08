      subroutine scfconvergence_test( first, iscf, last )
      USE SIESTA_TODO
      USE siesta_options
      implicit none

      logical :: first
      integer :: iscf
      logical :: last

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
          ! convergence test
          DEna = Enascf - Enaatm
          Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec+Emad+Emm+
     &           Emeta
          Eharrs = Etot + DEharr
          FreeE  = Etot - Temp * Entropy

          ! Recalculating the energy in the last iter (for gridcellsampling)
          ! but preserving the value of Eharrs1
          Eharrs1 = Eharrs

          dEmax = abs(Etot - Elast)
          Elast = Etot
          if (require_energy_convergence) then
            if (dDmax.lt.dDtol.and.dEmax.lt.dEtol) last = .true.
          else
            if (dDmax.lt.dDtol) last = .true.
          endif


          ! Print energies
          if (IOnode) then
            call siesta_write_energies( first, iscf )

            if (harrisfun) then
              write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ',Eharrs/eV
              if (cml_p) then
                call cmlStartPropertyList(mainXML, title='SCF Cycle')
                call cmlAddProperty(xf=mainXML, property=Eharrs/eV,
     .           units="siestaUnits:eV", dictRef="siesta:Eharrs", 
     .           fmt="(f14.7)")
                call cmlEndPropertyList(mainXML)
              endif
            endif
            if (cml_p) call cmlEndStep(mainXML)
            ! flush stdout
            call pxfflush(6)
            call wallclock("-------------- end of scf step")
          endif

!--------------------------------------------------------------------------- END
      END subroutine scfconvergence_test
