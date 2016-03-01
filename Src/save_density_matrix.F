! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_save_density_matrix

      implicit none
      public :: save_density_matrix

      CONTAINS

      subroutine save_density_matrix(id)

      ! Stores DM and EDM on files
      ! This version uses module variables.
      ! Eventually it can be a cleaner routine

      use precision,            only:  dp
      use sparse_matrices,      only:  maxnh, numh
      use sparse_matrices,      only:  listh, listhptr
      use files, only : slabel
      use m_iodm,               only:  write_dm
      use m_spin,               only:  nspin
      use atomlist,             only:  no_l
      use siesta_options,       only:  writedm, writedm_cdf
      use siesta_options,       only:  writedm_cdf_history
      use siesta_options,       only:  idyn
      use m_steps,              only:  istp
      use sparse_matrices,      only:  Dscf, DM_2D
#ifdef CDF
      use iodm_netcdf,          only:  write_dm_netcdf
#ifdef NCDF_4
      use dictionary
      use siesta_options, only:  write_cdf
      use m_ncdf_siesta, only : cdf_save_state
#endif
#endif
#ifdef TRANSIESTA
      use sparse_matrices,      only:  EDM_2D
      use m_ts_iodm,            only:  write_ts_dm
      use m_ts_global_vars,     only:  TSrun
      use m_ts_options,         only:  TS_DE_save
      use m_energies,           only:  Ef
#endif

      implicit none

      ! This is a just a comment in the caller side, unused for now
      character(len=*), intent(in)   :: id
#ifdef CDF
#ifdef NCDF_4
      type(dict) :: dic_save
#endif
#endif
#ifdef TIMING_IO
      integer :: i
#endif

#ifdef DEBUG
      call write_debug( '    PRE save_density_matrix' )
#endif

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
#ifdef NCDF_4
      if ( write_cdf ) then

         dic_save = ('Ef'.kv.1)
         dic_save = dic_save//('DM'.kv.1)//('EDM'.kv.1)
         call cdf_save_state('SIESTA.nc',dic_save)

         call delete(dic_save)

      end if
#endif
#endif
    
#ifndef TRANSIESTA
      if (writedm) then
#ifdef TIMING_IO
         call timer('IO-W-DM',1)
         do i = 1 , 100
#endif
        if ((idyn .eq. 6) .or. (idyn .eq. 7).or.(idyn.eq.9)) then
          if (istp.eq.1)
     &          call write_dm(trim(slabel)//'.DM', nspin, DM_2D)
        else
           call write_dm(trim(slabel)//'.DM', nspin, DM_2D)
        end if
#ifdef TIMING_IO
      end do
      call timer('IO-W-DM',2)
      call timer('IO-W-DM',3)
#endif
      endif ! writedm
#else /* TRANSIESTA */
! TSS Begin
      if (writedm) then
        if (.not.TSrun) then   !TSS save Dscf
#ifdef TIMING_IO
         call timer('IO-W-DM',1)
         do i = 1 , 100
#endif
          if ((idyn .eq. 6) .or. (idyn .eq. 7).or.(idyn==9)) then
            if (istp.eq.1)
     &            call write_dm(trim(slabel)//'.DM', nspin, DM_2D)
          else
             call write_dm(trim(slabel)//'.DM', nspin, DM_2D)
          end if
#ifdef TIMING_IO
      end do
      call timer('IO-W-DM',2)
      call timer('IO-W-DM',3)
#endif
        else
!         TSS write DscfE
#ifdef TIMING_IO
         call timer('IO-W-TS-DE',1)
         do i = 1 , 100
#endif
             if ((idyn .eq. 6) .or. (idyn .eq. 7).or.(idyn==9)) then
               if (istp .eq.1)
     &               call write_ts_dm(trim(slabel)//'.TSDE', nspin,
     &               DM_2D, EDM_2D, Ef)
             else
                call write_ts_dm(trim(slabel)//'.TSDE', nspin,
     &               DM_2D, EDM_2D, Ef)
             end if 
#ifdef TIMING_IO
      end do
      call timer('IO-W-TS-DE',2)
      call timer('IO-W-TS-DE',3)
#endif
        end if ! TSrun
        if ( TS_DE_save ) then
           call write_ts_dm(trim(slabel)//'.TSDE', nspin,
     &          DM_2D, EDM_2D, Ef)
        end if
      end if ! writedm
#endif /* TRANSIESTA */

#ifdef DEBUG
      call write_debug( '    POS save_density_matrix' )
#endif
      
      end subroutine save_density_matrix

      end module m_save_density_matrix
