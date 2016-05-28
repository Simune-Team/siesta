! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_save_density_matrix
      public :: save_density_matrix

      CONTAINS

      subroutine save_density_matrix(id)

      ! Stores DM and EDM on files
      ! This version uses module variables.
      ! Eventually it can be a cleaner routine

      use precision,            only:  dp
      use sparse_matrices,      only:  maxnh, numh
      use sparse_matrices,      only:  listh, listhptr
      use m_iodm,               only:  write_dm
      use m_spin,               only:  nspin
      use atomlist,             only:  no_l
      use siesta_options,       only:  idyn, writedm, writedm_cdf
      use siesta_options,       only:  writedm_cdf_history
      use m_steps,              only:  istp
      use sparse_matrices,      only:  Dscf
#ifdef CDF
      use iodm_netcdf,          only:  write_dm_netcdf
#endif
#ifdef TRANSIESTA
      use sparse_matrices,      only:  Escf
      use m_steps,              only:  istp
      use m_ts_iodm,            only:  write_ts_dm
      use m_ts_global_vars,     only:  TSrun
      use m_energies,           only:  ef
#endif

      implicit none

      ! This is a just a comment in the caller side, unused for now
      character(len=*), intent(in)   :: id

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
      if (writedm) then
        if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
          if (istp.eq.1)
     &    call write_dm (maxnh, no_l, nspin,
     &               numh, listhptr, listh, Dscf)
        else
          call write_dm (maxnh, no_l, nspin,
     &               numh, listhptr, listh, Dscf)
        endif
      endif ! writedm
#else /* TRANSIESTA */
! TSS Begin
      if (writedm) then
        if (.not.TSrun) then   !TSS save Dscf
          if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
            if (istp.eq.1)
     &        call write_dm( maxnh, no_l, nspin, numh,
     &                       listhptr, listh, Dscf )
          else
            call write_dm( maxnh, no_l, nspin, numh,
     &                     listhptr, listh, Dscf)
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
