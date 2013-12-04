      module m_scfconvergence_test
      private
      public :: scfconvergence_test
      CONTAINS

      subroutine scfconvergence_test( first, iscf, dDmax,
     $                                conv_harris, conv_etot,
     $                                converged)
      USE siesta_options
      use siesta_cml
      use m_wallclock, only : wallclock
      use parallel, only: IOnode
      use write_subs
      use m_energies
      use m_convergence, only: converger_t
      use m_convergence, only: add_value, is_converged

      implicit none

      integer :: iscf
      logical :: first

      real(dp), intent(in) :: dDmax     ! Max. change in dens. matrix
      type(converger_t), intent(inout)  :: conv_harris, conv_etot
      logical, intent(out)              :: converged

!------------------------------------------------------------------- BEGIN

      call timer( 'SCFconv', 1 )
      ! convergence test

      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec+Emad+Emm+
     &       Emeta
      Eharrs = Etot + DEharr

      call add_value(conv_harris, eharrs)
      call add_value(conv_etot, etot)

      FreeE  = Etot - Temp * Entropy

      converged = .false.

      Elast = Etot
      if (require_energy_convergence) then
         if (   dDmax.lt.dDtol
     &        .and. is_converged(conv_etot)) then
            converged = .true.
            if (IOnode) then
               write(6,"(a)") "SCF Convergence by Etot+DM criteria"
            endif
         endif
      else if (require_harris_convergence) then
         if (is_converged(conv_harris)) then
            converged = .true.
            if (IOnode) then
               write(6,"(a)") "SCF Convergence by Harris criterion"
            endif
         endif
      else
        if (dDmax.lt.dDtol) converged = .true.
      endif

      ! Print energies
      if (IOnode) then
        call siesta_write_energies( first, iscf, dDmax )

        if (harrisfun) then
          write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ',Eharrs/eV
          if (cml_p) then
            call cmlStartPropertyList(mainXML, title='SCF Cycle')
            call cmlAddProperty(xf=mainXML, value=Eharrs/eV,
     .       units="siestaUnits:eV", dictRef="siesta:Eharrs", 
     .       fmt="r7")
            call cmlEndPropertyList(mainXML)
          endif
        endif
       
        ! flush stdout
        call pxfflush(6)
        call wallclock("-------------- end of scf step")
      endif

      call timer( 'SCFconv', 2 )

!---------------------------------------------------------------- END
      END subroutine scfconvergence_test
      end module m_scfconvergence_test
