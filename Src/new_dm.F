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

#ifdef TRANSIESTA
      !TSS Begin
      ! SIESTA Modules
      use parallel, only: IOnode 
      use m_energies, only: ef
      use sys, only : die

      ! TS Modules
      use m_ts_options, only : VIn, VOut, TSinit, TSrun, foundts
      use m_ts_io, only : TSiodm 

      ! TSS End 

#endif /* TRANSIESTA */
      implicit none

      logical:: auxchanged ! Has auxiliary supercell changed?

      logical:: dminit     ! Initialize density matrix?
      integer:: iord       ! Order of DM extrapolation between geom. steps
      integer:: istpsave   ! Temporary to save and restore istp when reading DM from file
      logical:: itest      ! Temporary used to check if DM have just been read
      logical :: found
#ifdef TRANSIESTA
      ! TSS Begin
      ! Internal Variables
      integer :: ts_alloc_stat
      ! TSS End
#endif /* TRANSIESTA */

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

#ifdef TRANSIESTA
      ! TSS Begin: Mixing arrays
      if(.not.allocated(VIn)) then
       allocate(VIn(maxnh,nspin), STAT=ts_alloc_stat)
       if(ts_alloc_stat /= 0) 
     .   call die('Failled to Allocate VIn in new_dm !!')
       allocate(VOut(maxnh,nspin))
       if(ts_alloc_stat /= 0) 
     .   call die('Failled to Allocate VOut in new_dm !!')
      endif 
      ! TSS End

#endif /* TRANSIESTA */
      if (dminit)
     .   call initdm(Datm, Dscf, Dold, lasto, na_s,
     .               maxnh, no_s, no_l, nspin, na_u, no_l, nspin,
     .               numh, numhold, listhptr, listhptrold,
     .               listh, listhold, iaorb, found, inspn, 
     .               usesavedmloc, no_u)

#ifdef TRANSIESTA
      ! TSS Begin
      ! TSS check of TSDM stuff

      if (istp .eq. 1 .and. isolve .eq. 2) then ! try to read TSDM
         call TSiodm('read', maxnh, no_l, nspin,
     .        numh, listhptr, listh, Dscf, Escf, ef, foundts)

         if(.not. foundts )  then ! not a TS continuation run
            TSinit = .true.     ! start to converge a diagon run
            TSrun = .false.

            if(IONode) then
               write(6,'(a)')
     .              'TRANSIESTA: No TS-DensityMatrix file found'
               write(6,'(a)')
     .              'TRANSIESTA: Initialization runs using diagon'
            endif
         else
            TSinit = .false.
            TSrun = .true.

            if(IONode) then
               write(6,'(a)') 'TRANSIESTA: Continuation run'
               write(6,'(/a)')
     .              '                     ************************  '
               write(6,'(a)')
     .              '                     *   TRANSIESTA BEGIN   *  '
               write(6,'(a)')
     .              '                     ************************  '
            endif


      !     copy Dsfc -> Dold

            Dold(1:nh,1:nspin) = Dscf(1:nh,1:nspin)
            Eold(1:nh,1:nspin) = Escf(1:nh,1:nspin)

         end if                 !found TSDM-file
      end if                    !istp == 1  and isolve == 2

      ! TSS End

#endif /* TRANSIESTA */

      ! Initialize energy-density matrix to zero for first call to overfsm
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




