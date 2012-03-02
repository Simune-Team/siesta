      module m_check_scf_convergence
      private
      public :: check_scf_convergence
      CONTAINS

      subroutine check_scf_convergence( iscf, last, E_Harris, E_KS_good,
     $                                  dDmax, dHmax,
     $                                  conv_harris, conv_etot)

      use precision,             only: dp
      use siesta_cml,            only: cml_p, mainXML
      use siesta_cml,            only: cmlAddProperty
      use siesta_cml,            only: cmlStartPropertyList
      use siesta_cml,            only: cmlEndPropertyList
      use siesta_options,        only: Temp, harrisfun, dDtol
      use siesta_options,        only: require_harris_convergence
      use siesta_options,        only: require_energy_convergence
      use siesta_options,        only: energy_tolerance
      use siesta_options,        only: harris_tolerance
      use m_wallclock,           only: wallclock
      use parallel,              only: IOnode
      use units,                 only: eV
      use write_subs,            only: siesta_write_energies
      use m_energies,            only: Etot, EHarrs, Eharrs1
      use m_energies,            only: FreeE, Entropy
      use m_convergence,         only: converger_t
      use m_convergence,         only: add_value, is_converged

      implicit none

      integer, intent(in)  :: iscf
      logical, intent(out) :: last
      real(dp), intent(in) :: E_Harris, E_KS_Good
      real(dp), intent(in) :: dDmax     ! Max. change in density matrix
      real(dp), intent(in) :: dHmax     ! Max. change in  H
      type(converger_t), intent(inout)  :: conv_harris, conv_etot

!----------------------------------------------------------------- BEGIN

      call timer( 'SCFconv', 1 )
      ! convergence test

      ! Monitor Etot and Eharrs for optional convergence
      call add_value(conv_harris, E_Harris)
      call add_value(conv_etot, E_KS_Good)

      Etot = E_KS_Good
      Eharrs = E_Harris

      ! The electronic entropy is computed in compute_dm.
      ! It is used here to make the energy variational if
      ! occupation numbers change. Note that this is the
      ! entropy of the band structure computed from H_in,
      ! which might be claimed to be an "out" quantity.
      ! In previous versions the entropy of the previous
      ! SCF step was used instead (with no clear reason).

      FreeE  = Etot - Temp * Entropy

      ! Recalculating the energy in the last iter (for
      ! gridcellsampling) but preserving the value of Eharrs1
      Eharrs1 = Eharrs

      if (require_energy_convergence) then
         if (   dDmax.lt.dDtol
     &        .and. is_converged(conv_etot)) then
            last = .true.
            if (IOnode) then
               write(6,"(a)") "SCF Convergence by Etot+DM criteria"
            endif
         endif
      else if (require_harris_convergence) then
         if (is_converged(conv_harris)) then
            last = .true.
            if (IOnode) then
               write(6,"(a)") "SCF Convergence by Harris criterion"
            endif
         endif
      else
        if (dDmax.lt.dDtol) last = .true.
      endif

      ! Print energies
      if (IOnode) then
        call siesta_write_energies( iscf, dDmax )
        write(6,"(a,f12.7)") "Max(|H_out - H_in|): ", dHmax

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

!--------------------------------------------------------------------------- END
      END subroutine check_scf_convergence
      end module m_check_scf_convergence
