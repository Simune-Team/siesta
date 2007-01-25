      subroutine state_init( istep )
      USE SIESTA_TODO
      implicit none

      integer :: istep

      external :: automatic_cell,
     .  bonds, cgvc, cgvc_zmatrix, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, iozm,
     .  kgrid, kgridinit, kinefsm, mulliken, naefs, neighb,
     .  pulayx, 
     .  reinit, shaper, spnvec, 
     .  timer, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion, eggbox


!------------------------------------------------------------------------- BEGIN
      call timer( 'IterMD', 1 )
      istp = istp + 1
      if (IOnode) then
        write(6,'(/2a)') 'siesta:                 ',
     .                    '=============================='
        select case (idyn)
        case (0)
          write(6,'(28(" "),a,i6)') 'Begin CG move = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='CG', index=istp)
        case (1:5)
          write(6,'(28(" "),a,i6)') 'Begin MD step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='MD', index=istep)
        case (6)
          write(6,'(28(" "),a,i6)') 'Begin FC step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FC', index=istep)
          if (istep .eq. 0) then
            write(6,'(28(" "),a)') 'Undisplaced coordinates'
          else
            iadispl = (istep-mod(istep-1,6))/6+ia1
            write(6,'(28(" "),a,i6)') 'displace atom   ',
     .        iadispl
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(28(" "),a,i6)') 'in direction    ',
     .        ixdispl
            dx=-dx
            write(6,'(28(" "),a,f8.4,a)') 'by       ',
     .                      dx, ' Bohr'
! Displace atom by dx...
            xa(ixdispl,iadispl) = xa(ixdispl,iadispl) + dx
          endif
        case (7)
          call phonon_set_coords(istep,xa,ucell)
        case (8)
          write(6,'(28(" "),a,i6)') 'Begin Server step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FS', index=istep)

        end select
        write(6,'(2a)') '                        ',
     .                    '=============================='
      else ! not IOnode
        select case (idyn)
! We don't need to do anything for 0<idyn<6
        case(6)
          if (istep .ne. 0) then
            iadispl = (istep-mod(istep-1,6))/6 + ia1
            ix = mod(istep-1,6) + 1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            dx = - dx
! Displace atom by dx...
            xa(ixdispl,iadispl) = xa(ixdispl,iadispl) + dx
          endif
        case(7)
          call phonon_set_coords(istep,xa,ucell)
        end select
      endif

! Get coordinates from driver program trough pipe
      if (idyn.eq.8 .and. istep.ne.inicoor) then
        call coordsFromPipe( na_u, xa, ucell )
        if (volcel(ucell) < 1.0e-8_dp) then
           call automatic_cell(ucell,scell,na_u,xa,isa,charnet)
        endif
      end if

      if (IOnode) then
! Print Z-matrix coordinates
        if (lUseZmatrix) then
          call write_Zmatrix
        endif
! Print atomic coordinates ............................................
        call outcoor( ucell, xa, na_u, ' ', writec )
        call siesta_write_positions()
      endif
! ...................

! Actualize things if variable cell ...................................
      auxchanged = .false.
      if ( (varcel .or. idyn.eq.8)
     $     .and. (istep.ne.inicoor) .and. (.not.gamma) ) then

! k-grid 
        nkpnt = maxk
        call kgridinit( ucell, kscell, kdispl, kcutof, nkpnt )

! If number of k points has changed - re-size arrays
        if (nkpnt .gt. maxk) then
          call memory('D','D',size(kpoint),'siesta')
          deallocate(kpoint)
          call memory('D','D',size(kweight),'siesta')
          deallocate(kweight)
          allocate(kpoint(3,nkpnt))
          call memory('A','D',3*nkpnt,'siesta')
          allocate(kweight(nkpnt))
          call memory('A','D',nkpnt,'siesta')
          kpoint(1:3,1)=0.0_dp
          kweight(1)=1.0_dp
          call kgrid( ucell, kscell, kdispl,
     .      nkpnt, kpoint, kweight )
          call memory('D','D',size(eo),'siesta')
          deallocate(eo)
          call memory('D','D',size(qo),'siesta')
          deallocate(qo)
          allocate(eo(no_u,nspin,nkpnt))
          call memory('A','D',no_u*nspin*nkpnt,'siesta')
          allocate(qo(no_u,nspin,nkpnt))
          call memory('A','D',no_u*nspin*nkpnt,'siesta')
          maxk = nkpnt
       else
          call kgrid( ucell, kscell, kdispl,
     .                nkpnt, kpoint, kweight )

       endif
 
! Print k-points 
       if (IOnode) call siesta_write_k_points()
 
! Find required supercell

        auxchanged = .false.

        if (gamma) then
           nsc(1:3) = 1
        else if (fixauxcell) then
           nsc(i) = nscold(i)
        else
           do i=1,3
              veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
              nsc(i) = ceiling( 2 * rmaxh / veclen )
           enddo
           if (.not. naiveauxcell)
     $         call check_sc_factors(ucell,nsc,2*rmaxh)
        endif

        mscell = 0.0_dp
        do i = 1, 3
           mscell(i,i) = nsc(i)
           if (nsc(i).ne.nscold(i)) auxchanged = .true.
           nscold(i) = nsc(i)
        enddo

! Madelung correction for charged systems .............................
        if (charnet .ne. 0.0_dp) then
          call madelung(ucell, shape, charnet, Emad)
        endif

      endif
! End variable cell actualization

! Auxiliary supercell
      call superc(ucell, scell, nsc)

! Print unit cell and find its volume
      if (IOnode) call outcell(ucell)
      volume = volcel( ucell )
! ...................

! Initialize neighb subroutine ........................................
  144 ia = 0
      isel = 0
      rmax = max( 2._dp*rmaxv, 2._dp*rmaxo, rmaxo+rmaxkb )
      nnia = maxna
      if (allocated(jna)) then
        call memory('D','I',size(jna),'siesta')
        deallocate(jna)
      endif
      if (allocated(r2ij)) then
        call memory('D','D',size(r2ij),'siesta')
        deallocate(r2ij)
      endif
      if (allocated(xij)) then
        call memory('D','D',size(xij),'siesta')
        deallocate(xij)
      endif
      allocate(jna(maxna))
      call memory('A','I',maxna,'siesta')
      allocate(r2ij(maxna))
      call memory('A','D',maxna,'siesta')
      allocate(xij(3,maxna))
      call memory('A','D',3*maxna,'siesta')
      call neighb( scell, rmax, na_s, xa, ia, isel,
     .             nnia, jna, xij, r2ij )
      nnamax = 0
      do ia = 1,na_s
        nnia = 0
        call neighb( scell, rmax, na_s, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        nnamax = max( nnamax, nnia )
      enddo
      if (nnamax .gt. maxna) then
! Increase maxna with safety margin when atoms move
        maxna = nnamax + 0.10 * nnamax + 10
        overflow = .true.
      else
        overflow = .false.
      endif
      if (overflow) goto 144
! ..................

! Check if any two atoms are unreasonably close .......................
      do ia = 1,na_s
        r2min = huge(1._dp)
        jamin = 0
        nnia = maxna
        call neighb( scell, rmax, na_s, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        do j = 1,nnia
          ja = jna(j)
          if ( r2ij(j).lt.r2min .and. ja.ge.ia ) then
!           Check that it is not the same atom
            if ( ja.ne.ia .or. r2ij(j).gt.1.d-12 ) then
              r2min = r2ij(j)
              jamin = ja
            endif
          endif
        enddo
        rmin = sqrt( r2min )
        if (IOnode) then
          if ( rmin .lt. rijmin ) write(6,'(a,2i6,a,f12.6,a)')
     .      'siesta: WARNING: Atoms', ia, jamin, ' too close: rij =',
     .       rmin/Ang, ' Ang'
        endif
      enddo
! ..................

! List of nonzero Hamiltonian matrix elements .........................
      overflow=.true.
      overflowed=.false.
      do while (overflow)
         nh = maxnh
         call re_alloc(listh,1,maxnh,name='listh',routine='siesta',
     .                 copy=.false.)
         call hsparse( negl, scell, nsc, na_s, isa, xa, lasto, lastkb, 
     .                 iphorb, iphKB, nh, numh, listhptr, listh )
         if (nh .gt. maxnh) then
      ! Increase maxnh with safety margin for when atoms move
            maxnh = 1.05 * nh + 40
            overflowed=.true.
         else
            overflow=.false.
         endif
      enddo
      ! In first step, allocate anyway (to catch corner case
      ! where one node has nh=0, and doesn't overflow)
      if (istp==1) then
        call re_alloc(Dscf,1,maxnh,1,nspin,name='Dscf',
     .                routine='siesta',copy=.false.)
        call re_alloc(Dscfsave,1,maxnh,1,nspin,name='Dscfsave',
     .                routine='siesta',copy=.false.)
        call re_alloc(listhold,1,maxnh,name='listhold',
     .                routine='siesta',copy=.false.)
        ! Initialise Dscfsave to avoid problems in extrapol
        Dscfsave(1:maxnh,1:nspin) = 0.0_dp
        listhold(1:maxnh) = 0
      elseif (overflowed) then
        ! We need to preserve the contents of these arrays
        ! when reallocating.
        call re_alloc(Dscf,1,maxnh,1,nspin,name='Dscf',
     .                routine='siesta',copy=.true.)
        call re_alloc(Dscfsave,1,maxnh,1,nspin,name='Dscfsave',
     .                routine='siesta',copy=.true.)
        call re_alloc(listhold,1,maxnh,name='listhold',
     .                routine='siesta',copy=.true.)
      endif
      if (istp==1.or.overflowed) then
        call re_alloc(Dold,1,maxnh,1,nspin,name='Dold',
     .                routine='siesta',copy=.false.)
        call re_alloc(Eold,1,maxnh,1,nspin,name='Eold',
     .                routine='siesta',copy=.false.)
        call re_alloc(Escf,1,maxnh,1,nspin,name='Escf',
     .                routine='siesta',copy=.false.)
      endif

! Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      call re_alloc(H,1,maxnh,1,nspin,name='H',
     .              routine='siesta',shrink=.false.,copy=.false.)
      call re_alloc(H0,1,maxnh,name='H0',routine='siesta',
     .              shrink=.false.,copy=.false.)
      call re_alloc(S,1,maxnh,name='S',routine='siesta',
     .              shrink=.false.,copy=.false.)

! ..................

! Some printout for debugging ........................................
*     if (IOnode) then
*       write(6,'(/,a)') 'siesta: connected orbitals'
*       do io = 1,no_u
*         call GlobalToLocalOrb(io,Node,Nodes,iio)
*         if (iio.gt.0) then
*           write(6,'(i6,4x,15i4)') 
*    .        io, (listh(listhptr(iio)+j),j=1,numh(iio))
*         endif
*#ifdef MPI
*         call MPI_Barrier(MPI_Comm_World,MPIerror)
*#endif
*       enddo
*       write(6,*) ' '
*     endif
! ..................

! Find vectors between orbital centers ................................
      if (allocated(xijo)) then
        call memory('D','D',size(xijo),'siesta')
        deallocate(xijo)
      endif
      if (.not.gamma) then
        nxij = maxnh
        allocate(xijo(3,nxij))
        call memory('A','D',3*nxij,'siesta')
        call xijorb( negl, scell, na_u, na_s, xa,
     .               lasto, lastkb, rco, rckb,
     .               maxnh, numh, listhptr, listh, xijo )
      else
        nxij = 1
        allocate(xijo(3,1))
        call memory('A','D',3,'siesta')
      endif
! ..................

! Initialize density matrix ...........................................
! set density matrix for first step
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
     .               usesavedmloc, no_u)


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
! ..................

! Check for Pulay auxiliary matrices sizes ...................................
      if (pulfile .or. maxsav .le. 0) then
        nauxpul = 1
        if (.not.allocated(auxpul)) then
          allocate(auxpul(nauxpul,2))
          call memory('A','D',2*nauxpul,'siesta')
        endif
      else
        nauxpul = 0
        do io = 1,no_l
          nauxpul = nauxpul + numh(io)
        enddo
        nauxpul = nauxpul * nspin * maxsav
#ifdef MPI
        call globalize_max(nauxpul,ntmp)
        nauxpul = ntmp
#endif
! Increase nauxpul with safety margin when atoms move
        nauxpul = 1.1 * nauxpul + 10
        if (allocated(auxpul)) then
          if (size(auxpul,1).ne.nauxpul) then
            call memory('D','D',size(auxpul),'siesta')
            deallocate(auxpul)
            allocate(auxpul(nauxpul,2))
            call memory('A','D',2*nauxpul,'siesta')
          endif
        else
          allocate(auxpul(nauxpul,2))
          call memory('A','D',2*nauxpul,'siesta')
        endif
      endif
! ....................

! Find overlap matrix ...............................................
      call overfsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo, no_l,
     .             maxna, maxnh, maxnh, lasto, iphorb, isa, 
     .             numh, listhptr, listh, numh, listhptr, listh, 
     .             min(nspin,2), Escf, jna, xij, r2ij,
     .             fal, stress, S )
! ..................

!--------------------------------------------------------------------------- END
      END subroutine state_init
