      module m_scfconvergence_test
      private
      public :: scfconvergence_test
      CONTAINS

      subroutine scfconvergence_test( iscf, last,
     $                                conv_harris, conv_etot)

      use precision,             only: dp
      use siesta_cml,            only: cml_p, mainXML
      use siesta_cml,            only: cmlAddProperty
      use siesta_cml,            only: cmlStartPropertyList
      use siesta_cml,            only: cmlEndPropertyList
      use siesta_options,        only: mixH, Temp, harrisfun, dDtol
      use siesta_options,        only: require_harris_convergence
      use siesta_options,        only: require_energy_convergence
      use siesta_options,        only: energy_tolerance
      use siesta_options,        only: harris_tolerance
      use m_wallclock,           only: wallclock
      use parallel,              only: IOnode
      use units,                 only: eV
      use write_subs,            only: siesta_write_energies
      use m_energies            
      use atomlist,              only: no_l
      use m_spin,                only: nspin
      use sparse_matrices,       only: numh, listhptr
      use sparse_matrices,       only: H, Hold, Dscf, Dold
#ifdef MPI
      use m_mpi_utils,           only: globalize_max
#endif
      use m_convergence,         only: converger_t
      use m_convergence,         only: add_value, is_converged
#ifdef TRANSIESTA
      use m_ts_global_vars,      only: TSinit, TSrun, TSiscf
      use m_ts_options,          only: TSmode
      use m_pulay,               only: init_pulay_arrays
#endif /* TRANSIESTA */
      implicit none

      integer, intent(in)  :: iscf
      logical, intent(out) :: last
      type(converger_t), intent(inout)  :: conv_harris, conv_etot

      real(dp) :: dDmax     ! Max. change in dens. matrix or H
      real(dp), dimension(:,:), pointer :: Xin => null()
      real(dp), dimension(:,:), pointer :: Xout => null()  

      integer               :: ind, in, is, i
#ifdef MPI
      real(dp)              :: buffer1
#endif

!----------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '  PRE scfconvergence_test' )
#endif
      call timer( 'SCFconv', 1 )
      ! convergence test

      if (mixH) then
         ! Mixing Hamiltonian
         Xin => Hold
         Xout => H
      else
         ! Mix density matrix
         Xin => Dold
         Xout => Dscf
      endif

    !  Compute current maximum deviation ...........
      dDmax = 0.0_dp
      do is = 1,nspin
         do i = 1,no_l
            do in = 1,numh(i)
               ind = listhptr(i) + in
               dDmax = max(dDmax, abs(Xout(ind,is) - Xin(ind,is)))
            enddo
         enddo
      enddo
    ! .......
#ifdef MPI
!     Ensure that dDmax is the same on all nodes for convergence test/output
      call globalize_max(dDmax,buffer1)
      dDmax = buffer1
#endif


      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec+Emad+Emm+
     &       Emeta
      Eharrs = Etot + DEharr

      call add_value(conv_harris, eharrs)
      call add_value(conv_etot, etot)

      FreeE  = Etot - Temp * Entropy

      ! Recalculating the energy in the last iter (for
      ! gridcellsampling) but preserving the value of Eharrs1
      Eharrs1 = Eharrs

      Elast = Etot
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

#ifdef TRANSIESTA
!     TSS Begin
      if (TSmode .and. last .and. TSinit) then
         last   = .false.     ! we're only done with the initial diagon run
         TSinit = .false.     ! now we start the TRANSIESTA run
         TSrun = .true.
           !
           ! From now on, a new mixing cycle starts,
           ! so we have to reset the Pulay record keeping
           !
           call init_pulay_arrays()  
           !
         if(IONode) then
           write(6,'(/a)')
     &'                          ************************       '
           write(6,'(a)')
     &'                          *   TRANSIESTA BEGIN   *       '
           write(6,'(a)')
     &'                          ************************       '
        end if
      end if
      if (TSiscf .eq. 1) call timer( 'IterSCF', 3 ) !TSS
      if (TSiscf .eq. 1) call timer( 'TS', 3 ) !TSS
      if (last) TSiscf = 0

!     TSS End
#endif /* TRANSIESTA */

      call timer( 'SCFconv', 2 )
#ifdef DEBUG
      call write_debug( '  POS scfconvergence_test' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine scfconvergence_test
      end module m_scfconvergence_test
