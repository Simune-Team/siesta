      module m_save_density_matrix
      public :: save_density_matrix

      CONTAINS

      subroutine save_density_matrix(when,file)

      ! Stores DM and EDM on files
      ! This version uses module variables.
      ! Eventually it can be a cleaner routine

      use precision,            only:  dp
      use sparse_matrices,      only:  maxnh, numh
      use sparse_matrices,      only:  listh, listhptr
      use m_iodm,               only:  write_dm
      use m_spin,               only:  nspin
      use atomlist,             only:  no_l
      use siesta_options,       only:  writedm, writedm_cdf
      use siesta_options,       only:  writedm_cdf_history
      use m_steps,              only:  istp
      use sparse_matrices,      only:  Dscf
#ifdef CDF
      use iodm_netcdf,          only:  write_dm_netcdf
#endif
#ifdef TRANSIESTA
      use siesta_options,       only:  idyn, writedm
      use sparse_matrices,      only:  Escf
      use m_steps,              only:  istp
      use m_ts_iodm,            only:  write_ts_dm
      use m_ts_global_vars,     only:  TSrun
      use m_energies,           only:  ef
#endif

      logical, intent(in), optional  :: when
      character(len=*), intent(in), optional   :: file

      logical :: do_it

      ! Control the condition to write the DM externally
      ! This gives more flexibility
      ! We retain the "writedm" default condition for compatibility

      if (present(when)) then
         do_it = when
      else
         do_it = writedm
      endif

#ifdef CDF
!     Save density matrix on disk, after mixing, to
!     be used for re-starting the SCF cycle.
      if (writedm_cdf_history) then
         call write_dm_netcdf( no_l, maxnh, nspin, Dscf,
     &                         overwrite=.false. )
      else if (writedm_cdf) then
         call write_dm_netcdf( no_l, maxnh, nspin, Dscf,
     &                         overwrite=.true. )
      endif
#endif
    
#ifndef TRANSIESTA
      if (do_it) then
        if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
          if (istp.eq.1)
     &    call write_dm (maxnh, no_l, nspin,
     &               numh, listhptr, listh, Dscf,
     $               userfile=file)
        else
          call write_dm (maxnh, no_l, nspin,
     &               numh, listhptr, listh, Dscf,
     $               userfile=file)
        endif
      endif ! writedm
#else /* TRANSIESTA */
! TSS Begin
      if (do_it) then
        if (.not.TSrun) then   !TSS save Dscf
          if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
            if (istp.eq.1)
     &        call write_dm( maxnh, no_l, nspin, numh,
     &                       listhptr, listh, Dscf,
     $                       userfile=file)
          else
            call write_dm( maxnh, no_l, nspin, numh,
     &                     listhptr, listh, Dscf,
     $                     userfile=file)
          endif
        else
!         TSS write DscfE
             if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
               if (istp .eq.1)
     .           call write_ts_dm (maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf, Escf,ef)
             else
                call write_ts_dm (maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf, Escf,ef)
             endif            ! idyn.eq.6
        endif ! TSrun
      endif ! writedm
#endif /* TRANSIESTA */
      end subroutine save_density_matrix

      end module m_save_density_matrix
