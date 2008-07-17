      MODULE m_new_dm
      private
      public :: new_dm

      CONTAINS

      subroutine new_dm( auxchanged )
      use sparse_matrices
      USE siesta_options
      use siesta_geom
      use atomlist, only: datm, no_s, iaorb, lasto, no_u, no_l
      use m_steps,  only: istp
      use m_spin,   only: nspin
      use m_extrapol, only: extrapol

      implicit none

      logical:: auxchanged ! Has auxiliary supercell changed?

      logical:: dminit     ! Initialize density matrix?
      integer:: iord       ! Order of DM extrapol. between geom. steps
      integer:: istpsave   ! Temporary to save and restore istp
      logical:: itest      ! Check if DM have just been read
      logical :: found

      ! Initialize density matrix
      ! Set density matrix for first step
      found = .false.
      dminit = .false.
      if (istp .eq. 1) dminit = .true.
      if (istp .ne. 1 .and. harrisfun) dminit = .true.
      if (istp .ne. 1 .and. (idyn .eq. 6)
     $         .and. usesavedm .and. writedm)  dminit = .true.
      if (istp .ne. 1 .and. (idyn .eq. 7)
     $         .and. usesavedm)  dminit = .true.

      ! If auxiliary cell has changed, optionally reset density matrix
      ! and set usesavedata to false to avoid reading back saved copy
      if (initdmaux.and.auxchanged) then
        dminit = .true.
        usesavedmloc = .false.
      else
        usesavedmloc = usesavedm
      endif

      if (dminit)
     .   call initdm(Datm, Dscf, Dold, lasto, na_s,
     .               maxnh, no_s, no_l, nspin, na_u, no_l, nspin,
     .               numh, numhold, listhptr, listhptrold,
     .               listh, listhold, iaorb, found, inspn, 
     .               usesavedmloc )


      ! Initialize energy-density matrix to zero
      Escf(1:maxnh,1:nspin) = 0.0_dp

      ! Extrapolate density matrix between steps
      itest = .false.
      istpsave = 0
      iord = 1
      if (idyn .eq. 0) iord = 0
      if (idyn .eq. 6) iord = 0
      if (idyn .eq. 7) iord = 0
      !  If DM has just been read from disk, 
      !  call extrapol with istep = 2 and iord = 0
      !  to make it update the structure of DM, if needed
      if (found .and. ((istp .eq. 1) .or. (idyn .eq. 6)
     .                               .or. (idyn .eq. 7))) then
        istpsave = istp
        istp = 2
        iord = 0
        itest = .true.
      endif
      if (.not.harrisfun)
     . call extrapol(istp, iord, nspin, no_s, no_l, maxnh, 
     .              numh, listhptr, listh, numhold, listhptrold,
     .              listhold, Dscfsave, Dscf)
      !  If DM have just been read, restore istp
      if (itest) istp = istpsave
      itest = .false.
      END subroutine new_dm

      END MODULE m_new_dm




