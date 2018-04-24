! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
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

      use m_spin,      only: nspin
      use sparse_matrices, only:  Dscf, S, numh, listhptr, listh, maxnh
      use siesta_geom, only: na_u, isa
      use atomlist, only:    no_u, lasto, iaorb, iphorb

      implicit none

      integer :: iscf
      logical :: first

      real(dp), intent(in) :: dDmax     ! Max. change in dens. matrix or H
      type(converger_t), intent(inout)  :: conv_harris, conv_freeE
      logical, intent(out)              :: converged

      real(dp) :: dummy_qspin(8)
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

!     Print spin polarization at each SCF step if requested before mixing
      if (spndeb) then
         call print_spin(dummy_qspin)
      endif
!     Print populations at each SCF step if requested before mixing ......
      if (muldeb) then
        if (ionode) write (6,"(/a)")
     &    'siesta: Mulliken populations (DM_out)'
        call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     &                 numh, listhptr, listh, S, Dscf, isa,
     &                 lasto, iaorb, iphorb )
      endif

      converged = .false.

      if (require_energy_convergence) then
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
        if (dDmax.lt.dDtol) then
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
              call pxfflush(6)
           endif
        endif
      endif
      if (converged) then
         if (IOnode) then
            write(6,"(a,i4,a)") "SCF cycle converged after ",
     $           iscf," iterations"
         endif
      endif
      if (harrisfun) then
         converged = .true.
         if (IOnode) then
            write(6,"(a)") "Harris-functional calculation " //
     $                     "considered 'converged'"
         endif
      endif

!---------------------------------------------------------------- END
      END subroutine scfconvergence_test
      end module m_scfconvergence_test
