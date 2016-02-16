      module m_scfconvergence_test
      private
      public :: scfconvergence_test
      CONTAINS

      subroutine scfconvergence_test( first, iscf, dDmax, dHmax,
     $                                conv_harris, conv_freeE,
     $                                converged)
      USE siesta_options
      use siesta_cml
      use m_wallclock, only : wallclock
      use parallel, only: IOnode
      use write_subs
      use m_energies
      use units, only: eV
      use m_convergence, only: converger_t, tolerance
      use m_convergence, only: add_value, is_converged

      implicit none

      integer :: iscf
      logical :: first

      real(dp), intent(in) :: dDmax     ! Max. change in dens. matrix
      real(dp), intent(in) :: dHmax     ! Max. change in H
      type(converger_t), intent(inout)  :: conv_harris, conv_freeE
      logical, intent(out)              :: converged


      logical :: default_DM_convergence
      logical :: default_NoEnergy_convergence
      logical :: is_Harris_conv, is_FreeEnergy_conv
      logical :: is_H_conv, is_DM_conv, is_Delta_conv

!------------------------------------------------------------------- BEGIN

      ! convergence test

      call add_value(conv_harris, eharrs)
      call add_value(conv_freeE, freeE)

      ! Print energies
      if (IOnode) then
        call siesta_write_energies( iscf, dDmax, dHmax )

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

      is_DM_conv = (dDmax.lt.dDtol)
      is_H_conv =  (dHmax.lt.dHtol)
      is_Harris_conv = is_converged(conv_harris)
      is_FreeEnergy_conv = is_converged(conv_freeE)

      default_DM_convergence =
     $      (.not. require_hamiltonian_convergence)
      default_NoEnergy_convergence =
     $      (.not. require_energy_convergence)

      is_Delta_conv = ((require_hamiltonian_convergence
     $                  .and. is_H_conv) .or.
     $                 (default_DM_convergence
     $                  .and. is_DM_conv) ) 

      converged = ( (require_harris_convergence
     $               .and. is_Harris_conv) .or.
     $              (require_energy_convergence
     $               .and. is_FreeEnergy_conv 
     $               .and. is_Delta_conv)  .or.
     $              (default_NoEnergy_convergence
     $               .and. is_Delta_conv)      )

      if (converged .and. IOnode) then

         if (require_harris_convergence) then
            write(6,"(/,a)") "SCF Convergence by Harris criterion"
            write(6,"(a,f14.8)") "|EH_i-EH_(i-1)| (eV) < ",
     $           tolerance(conv_harris)/eV
         endif
      
         if (require_hamiltonian_convergence) then
            if (require_energy_convergence) then
              write(6,"(/,a)") "SCF Convergence by FreeE+dHmax criteria"
              write(6,"(a,f14.8)") "|FreeE_i-FreeE_(i-1)| (eV) < ",
     $              tolerance(conv_freeE)/eV
            else
              write(6,"(/,a)") "SCF Convergence by dHmax criterion"
            endif
            write(6,"(a,f14.8)") "max |H_out - H_in| (eV): ", dHmax/eV
         endif
            
         if (default_DM_convergence) then
            if (require_energy_convergence) then
              write(6,"(/,a)") "SCF Convergence by FreeE+dDmax criteria"
              write(6,"(a,f14.8)") "|FreeE_i-FreeE_(i-1)| (eV) < ",
     $                             tolerance(conv_freeE)/eV
            else 
              write(6,"(/,a)") "SCF Convergence by dDmax criteria"
            endif

            if (mix_charge) then
               write(6,"(a,f14.8)") "max |DM_i - DM_(i-1)|: ", dDmax
            else
               write(6,"(a,f14.8)") "max |DM_out - DM_in|: ", dDmax
            endif
         endif
      
         write(6,"(a,i0,a)") "SCF cycle converged after ",
     $           iscf," iterations"

         call pxfflush(6)

      endif

!---------------------------------------------------------------- END
      END subroutine scfconvergence_test
      end module m_scfconvergence_test
