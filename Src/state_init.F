      MODULE m_state_init

      private
      public :: state_init

      CONTAINS

      subroutine state_init( istep )
      USE Kpoint_grid, only: setup_Kpoint_grid
      USE m_new_dm, only: new_dm
      USE m_proximity_check, only: proximity_check
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: iphorb, iphkb, rco, rckb, indxua,no_s, rmaxkb,
     .                    rmaxo, no_u, lastkb, lasto, superc, rmaxv
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
      use m_spin, only: nspin
      use m_rmaxh
      use m_pulay, only: init_pulay_arrays
      use m_eo
      use m_gamma
#ifdef MPI
      use m_mpi_utils, only : globalize_max
#endif

      implicit none

      integer :: istep

      real(dp):: veclen     ! Length of a unit-cell vector
      real(dp):: rmax
      logical :: cell_can_change

      integer :: i, ix, iadispl, ixdispl

      logical:: auxchanged ! Has auxiliary supercell changed?

      external ::  madelung, timer

#ifdef MPI
      integer :: ntmp      ! Temporary to call an MPI routine
#endif
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


      ! Actualize things if variable cell 
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

        ! Madelung correction for charged systems 
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

      ! Use largest possible range in program, except hsparse...
      ! 2 * rmaxv: Vna overlap
      ! 2 * rmaxo: orbital overlap
      ! rmaxo + rmaxkb: Non-local KB action
      ! 2.0_dp * (rmaxo+rmaxkb) : Orbital interaction through KB projectors

      rmax = max( 2._dp*rmaxv, 2._dp*rmaxo, rmaxo+rmaxkb)

      ! AG: Why not use the largest range and be done with it?
      if (.not. negl) then
         rmax = max(rmax, 2.0_dp * (rmaxo+rmaxkb) )
      endif

      ! Check if any two atoms are unreasonably close
      call proximity_check(rmax)

      if (istp==1) then 
        nullify(Dscf,Dscfsave,Dold,Eold,Escf)
        nullify(listhold,listh)
        nullify(H,H0,S)
        nullify(xijo)
      endif

      ! List of nonzero Hamiltonian matrix elements
      ! and, if applicable,  vectors between orbital centers

      ! Listh (and xijo if .not. gamma) are allocated inside hsparse
      ! nh is set by the routine also
      nh = maxnh
      call hsparse( negl, scell, nsc, na_s, isa, xa, lasto, lastkb, 
     .                 iphorb, iphKB, nh, numh, listhptr, listh,
     $                 gamma, xijo )

      maxnh = nh
      !OBSOLETE: 
      !  safety margin for when atoms move 
      !!!!!!!!   if (nh > maxnh)  maxnh = 1.05 * nh + 40 

      ! Use re_alloc regardless of whether nh is greater than 
      ! the old maxnh or not, BUT DO NOT SHRINK, as some of
      ! these arrays are later used in extrapol.

      call re_alloc(Dscf,1,maxnh,1,nspin,name='Dscf',
     .     routine='siesta',copy=.true.,SHRINK=.false.)
      call re_alloc(Dscfsave,1,maxnh,1,nspin,name='Dscfsave',
     .     routine='siesta',copy=.true.,SHRINK=.false.)
      call re_alloc(listhold,1,maxnh,name='listhold',
     .     routine='siesta',copy=.true.,SHRINK=.false.)
!==============================================================

      if (istp==1) then
         ! Initialise Dscfsave to avoid problems in extrapol
         ! This might be superfluous, as alloc fills array with zeroes
         Dscfsave(1:maxnh,1:nspin) = 0.0_dp
         listhold(1:maxnh) = 0
      endif

      call re_alloc(Dold,1,maxnh,1,nspin,name='Dold',
     .     routine='siesta',copy=.false.)
      call re_alloc(Eold,1,maxnh,1,nspin,name='Eold',
     .     routine='siesta',copy=.false.)
      call re_alloc(Escf,1,maxnh,1,nspin,name='Escf',
     .     routine='siesta',copy=.false.)

      ! Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      call re_alloc(H,1,maxnh,1,nspin,name='H',
     .              routine='siesta',shrink=.false.,copy=.false.)
      call re_alloc(H0,1,maxnh,name='H0',routine='siesta',
     .              shrink=.false.,copy=.false.)
      call re_alloc(S,1,maxnh,name='S',routine='siesta',
     .              shrink=.false.,copy=.false.)


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

      ! Initialize density matrix
      call new_dm( auxchanged )

      ! Check for size of Pulay auxiliary matrices
      call init_pulay_arrays()

      ! Find overlap matrix 
      call overlap(na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     .             maxnh, lasto, iphorb, isa, 
     .             numh, listhptr, listh, S)

!--------------------------------------------------------------------------- END
      END subroutine state_init
      END module m_state_init





