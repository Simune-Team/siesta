      MODULE m_state_init

      private
      public :: state_init

      CONTAINS

      subroutine state_init( istep )
      USE SIESTA_TODO
      USE Kpoint_grid, only: setup_Kpoint_grid
      USE m_new_dm, only: new_dm
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: rmaxv, iphorb, iphkb, rco, rckb, indxua,no_s,
     .                    rmaxo, no_u, lastkb, lasto, rmaxkb, superc
      use atomlist, only: no_l
      use alloc, only: re_alloc, alloc_report
      use m_hsparse,    only: hsparse
      use m_overlap,    only: overlap
      use m_check_supercell, only: check_sc_factors
      use siesta_cml
      use siesta_cmlsubs, only: siesta_cml_init
      use zmatrix,      only: lUseZmatrix, write_zmatrix
      use Kpoint_grid
      use m_energies, only: Emad
      use write_subs
      use m_steps
      use parallel, only: IOnode
      use alloc, only: re_alloc
      use m_numbs_neighb
      use m_spin, only: nspin
      use m_rmaxh
      use m_auxpul
      use m_eo

      implicit none

      integer :: istep

      real(dp):: veclen     ! Length of a unit-cell vector
      real(dp):: rmax
      logical :: cell_can_change

      integer :: i, ia, ix, io, iadispl, ixdispl, isel, nnamax
      integer :: nnia      ! Number of neighbors of a given atom

      external :: automatic_cell,
     .  bonds, cgvc, cgvc_zmatrix, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, iozm,
     .  kinefsm, mulliken, naefs, neighb, !kgrid, kgridinit,
     .  pulayx, 
     .  reinit, shaper, spnvec, 
     .  timer, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion, eggbox

      real(dp), external :: volcel
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
            write(6,'(28(" "),a,i6)') 'displace atom   ', iadispl
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(28(" "),a,i6)') 'in direction    ', ixdispl
            write(6,'(28(" "),a,f8.4,a)') 'by       ', dx, ' Bohr'
          endif

        case (8)
          write(6,'(28(" "),a,i6)') 'Begin Server step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FS', index=istep)

        end select
        write(6,'(2a)') '                        ',
     .                    '=============================='

        ! Print Z-matrix coordinates
        if (lUseZmatrix) then
          call write_Zmatrix()
        endif

        ! Print atomic coordinates
        call outcoor( ucell, xa, na_u, ' ', writec )
        call siesta_write_positions()
      endif


! Actualize things if variable cell ...................................
      ! This might be updated with PR MD options
      cell_can_change = (varcel .or. idyn.eq.8)

      auxchanged = .false.
      if ( cell_can_change
     $     .and. (istep.ne.inicoor) .and. (.not.gamma) ) then

        ! Will print k-points also
        call setup_Kpoint_grid( ucell )

        call re_alloc(eo,1,no_u,1,nspin,1,maxk,name="eo",
     $                routine="state_init")
        call re_alloc(qo,1,no_u,1,nspin,1,maxk,name="qo",
     $                routine="state_init")

        ! Find required supercell
        auxchanged = .false.

        if (gamma) then
           nsc(1:3) = 1
        else if (fixauxcell) then
           nsc(1:3) = nscold(1:3)
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
      ! Do not move from here, as the coordinates might have changed
      ! even if not the unit cell
      call superc(ucell, scell, nsc)

! Print unit cell and compute cell volume
! Possible BUG: 
! Note that this volume is later used in write_subs and the md output
! routines, even if the cell later changes.
      if (IOnode) call outcell(ucell)
      volume_of_some_cell = volcel(ucell)
! ...................

! Initialize neighb subroutine ........................................
  144 continue
      isel = 0  
      ! AG ** POSIBLE BUG 
      ! The last argument should be twice as big, and there
      ! should be a check for negl (neglect of KB overlap)
      ! Compare with rmaxh in siesta_init...
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
      ! Initialize neighb subroutine  (ia=0)
      call neighb( scell, rmax, na_s, xa, 0, isel,
     .             nnia, jna, xij, r2ij )
      nnamax = 0
      do ia = 1,na_s
        nnia = 0  ! First pass without filling arrays
                  ! just to get maximum size needed
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


      ! Check if any two atoms are unreasonably close
      call proximity_check(rmax)

      ! List of nonzero Hamiltonian matrix elements
      if (istp==1) then 
        ! Nullify arrays
        nullify(Dscf,Dscfsave,Dold,Eold,Escf)
        nullify(listhold,listh)
        nullify(H,H0,S)
        nullify(xijo)
      endif

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


      if (overflowed.or.istp==1) then
        ! In first step, allocate anyway (to catch corner case
        ! where one node has nh=0, and doesn't overflow)
        call re_alloc(Dscf,1,maxnh,1,nspin,name='Dscf',
     .                routine='siesta',copy=.true.)
        call re_alloc(Dscfsave,1,maxnh,1,nspin,name='Dscfsave',
     .                routine='siesta',copy=.true.)
        call re_alloc(listhold,1,maxnh,name='listhold',
     .                routine='siesta',copy=.true.)

        if (istp==1) then
          ! Initialise Dscfsave to avoid problems in extrapol
          Dscfsave(1:maxnh,1:nspin) = 0.0_dp
          listhold(1:maxnh) = 0
        endif

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

      if (.not.gamma) then
        call re_alloc(xijo,1,3,1,maxnh,name="xijo",
     $        routine="state_init",shrink=.false.,copy=.false.)
        call xijorb( negl, scell, na_u, na_s, xa,
     .               lasto, lastkb, rco, rckb,
     .               maxnh, numh, listhptr, listh, xijo )
      else
        call re_alloc(xijo,1,3,1,1,name="xijo",
     $        routine="state_init",shrink=.false.,copy=.false.)
      endif
! ..................

!     Initialize density matrix
      call new_dm()

! Check for Pulay auxiliary matrices sizes ...................................
      if (pulfile .or. maxsav .le. 0) then
        nauxpul = 1
        if (.not. allocated(auxpul)) then
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
      call overlap(na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     .             maxna, maxnh, lasto, iphorb, isa, 
     .             numh, listhptr, listh, 
     .             jna, xij, r2ij, S )
! ..................

!--------------------------------------------------------------------------- END
      END subroutine state_init
      END module m_state_init





