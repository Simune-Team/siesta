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

      subroutine save_density_matrix(when,file)

      ! Stores DM and EDM on files
      ! This version uses module variables.
      ! Eventually it can be a cleaner routine

      use precision,            only:  dp
      use sparse_matrices,      only:  maxnh, numh
      use sparse_matrices,      only:  listh, listhptr
      use files, only : slabel
      use m_iodm,               only:  write_dm
      use m_matio,              only:  write_mat
      use m_spin,               only:  nspin, h_spin_dim
      use atomlist,             only:  no_l
      use siesta_options,       only:  writedm, writedm_cdf
      use siesta_options,       only:  writedm_cdf_history
      use siesta_options,       only:  idyn
      use m_steps,              only:  istp
      use sparse_matrices,      only:  Dscf, DM_2D
      use fdf,                  only:  fdf_get
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

      logical, intent(in), optional  :: when
      character(len=*), intent(in), optional   :: file

      logical :: do_it

#ifdef CDF
#ifdef NCDF_4
      type(dict) :: dic_save
#endif
#endif
#ifdef TIMING_IO
      integer :: i
#endif

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
         call write_dm_netcdf( no_l, maxnh, h_spin_dim, Dscf,
     &                         overwrite=.false. )
      else if (writedm_cdf) then
         call write_dm_netcdf( no_l, maxnh, h_spin_dim, Dscf,
     &                         overwrite=.true. )
      endif
#ifdef NCDF_4
      if ( write_cdf ) then

         dic_save = ('Ef'.kv.1)
         dic_save = dic_save//('DM'.kv.1)//('EDM'.kv.1)
         call cdf_save_state(trim(slabel)//'.nc',dic_save)

         call delete(dic_save)

      end if
#endif
#endif
    
#ifndef TRANSIESTA
      if (do_it) then
#ifdef TIMING_IO
         call timer('IO-W-DM',1)
         do i = 1 , 100
#endif
        if (((idyn .eq. 6) .or. (idyn .eq. 7).or.(idyn.eq.9))
     $        .and. (istp .ne. 1)                      ) then
            ! do not write
        else
          if (fdf_get("Use.Blocked.WriteMat",.false.)) then
             call write_mat (maxnh, no_l, h_spin_dim,
     &               numh, listhptr, listh, Dscf,
     $               userfile=trim(file)//'.blocked',compatible=.false.)
          else
             ! The filename falls back to the standard one     
             !
             call write_dm(trim(slabel)//'.DM', DM_2D)
	  endif
        end if
#ifdef TIMING_IO
      end do
      call timer('IO-W-DM',2)
      call timer('IO-W-DM',3)
#endif
      endif ! writedm
#else /* TRANSIESTA */
! TSS Begin
      if (do_it) then
        if (.not.TSrun) then   !TSS save Dscf
#ifdef TIMING_IO
         call timer('IO-W-DM',1)
         do i = 1 , 100
#endif
          if ((idyn .eq. 6) .or. (idyn .eq. 7).or.(idyn==9)) then
            if (istp.eq.1)
     &            call write_dm(trim(slabel)//'.DM', DM_2D)
          else
             call write_dm(trim(slabel)//'.DM', DM_2D)
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
     &               call write_ts_dm(trim(slabel)//'.TSDE', 
     &               DM_2D, EDM_2D, Ef)
             else
                call write_ts_dm(trim(slabel)//'.TSDE',
     &               DM_2D, EDM_2D, Ef)
             end if 
#ifdef TIMING_IO
      end do
      call timer('IO-W-TS-DE',2)
      call timer('IO-W-TS-DE',3)
#endif
        end if ! TSrun
        if ( TS_DE_save ) then
           call write_ts_dm(trim(slabel)//'.TSDE',
     &          DM_2D, EDM_2D, Ef)
        end if
      end if ! writedm
#endif /* TRANSIESTA */

      end subroutine save_density_matrix

      end module m_save_density_matrix
