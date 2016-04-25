      module m_scfconvergence_test
      private
      public :: scfconvergence_test
      CONTAINS

      subroutine scfconvergence_test( first, iscf, dDmax,
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
      use ldau_specs,    only: switch_ldau
      use m_ldau,        only: conv_ldaupop

      implicit none

      integer :: iscf
      logical :: first

      real(dp), intent(in) :: dDmax     ! Max. change in dens. matrix or H
      type(converger_t), intent(inout)  :: conv_harris, conv_freeE
      logical, intent(out)              :: converged

!------------------------------------------------------------------- BEGIN

      ! convergence test

      call add_value(conv_harris, eharrs)
      call add_value(conv_freeE, freeE)

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

      converged = .false.

      if (require_energy_convergence) then
         if( switch_ldau ) then
           if (   dDmax.lt.dDtol
     $            .and. is_converged(conv_freeE)
     $            .and. conv_ldaupop) then
             converged = .true.
             if (IOnode) then
               write(6,"(/,a)") "SCF Convergence by FreeE+dmax criteria"
               if (mixH) then
                  write(6,"(a,f14.8)") "max |H_out - H_in| (Ry): ",
     $                                  dDmax/eV
               else if (mix_charge) then
                  write(6,"(a,f14.8)") "max |DM_i - DM_(i-1)|: ", dDmax
               else
                  write(6,"(a,f14.8)") "max |DM_out - DM_in|: ", dDmax
               endif
               write(6,"(a,f14.8)") "|FreeE_i-FreeE_(i-1)| (eV) < ",
     $                              tolerance(conv_freeE)/eV
               call pxfflush(6)
               write(6,"(/,a)") 
     $           "Population in the LDA+U projector convergence "
               write(6,"(a)") 
     $           "by dtol_ldaupop criterion"
             endif
           endif 
         endif
         if (   dDmax.lt.dDtol
     &        .and. is_converged(conv_freeE)) then
            converged = .true.
            if (IOnode) then
              write(6,"(/,a)") "SCF Convergence by FreeE+dmax criteria"
              if (mixH) then
                 write(6,"(a,f14.8)") "max |H_out - H_in| (Ry): ",
     $                                 dDmax/eV
              else if (mix_charge) then
                 write(6,"(a,f14.8)") "max |DM_i - DM_(i-1)|: ", dDmax
              else
                 write(6,"(a,f14.8)") "max |DM_out - DM_in|: ", dDmax
              endif
              write(6,"(a,f14.8)") "|FreeE_i-FreeE_(i-1)| (eV) < ",
     $                             tolerance(conv_freeE)/eV
              call pxfflush(6)
            endif
         endif
      else if (require_harris_convergence) then
         if (is_converged(conv_harris)) then
            converged = .true.
            if (IOnode) then
               write(6,"(/,a)") "SCF Convergence by Harris criterion"
               write(6,"(a,f14.8)") "|EH_i-EH_(i-1)| (eV) < ",
     $                             tolerance(conv_harris)/eV
               call pxfflush(6)
            endif
         endif
      else
        if ( switch_ldau ) then 
          if (dDmax.lt.dDtol.and.conv_ldaupop) then
            converged = .true.
            if (IOnode) then
              write(6,"(/,a)") "SCF Convergence by dMax criterion"
              if (mixH) then
                 write(6,"(a,f14.8)") "max |H_out - H_in| (eV): ",
     $                                 dDmax/eV
              else if (mix_charge) then
                 write(6,"(a,f14.8)") "max |DM_i - DM_(i-1)|: ", dDmax
              else
                 write(6,"(a,f14.8)") "max |DM_out - DM_in|: ", dDmax
              endif
              write(6,"(/,a)") 
     $          "Population in the LDA+U projector convergence "
              write(6,"(a)") 
     $          "by dtol_ldaupop criterion"
              call pxfflush(6)
            endif
          endif 
        else 
          if (dDmax.lt.dDtol) then
            converged = .true.
            if (IOnode) then
               write(6,"(/,a)") "SCF Convergence by dMax criterion"
               if (mixH) then
                  write(6,"(a,f14.8)") "max |H_out - H_in| (eV): ",
     $                                  dDmax/eV
               else if (mix_charge) then
                  write(6,"(a,f14.8)") "max |DM_i - DM_(i-1)|: ", dDmax
               else
                  write(6,"(a,f14.8)") "max |DM_out - DM_in|: ", dDmax
               endif
               call pxfflush(6)
            endif
          endif 
        endif
      endif
      if (converged) then
         if (IOnode) then
            write(6,"(a,i4,a)") "SCF cycle converged after ",
     $           iscf," iterations"
         endif
      endif


!---------------------------------------------------------------- END
      END subroutine scfconvergence_test
      end module m_scfconvergence_test
