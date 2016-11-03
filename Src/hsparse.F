! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_hsparse
!
!     We need an explicit interface since assumed-shape arrays are used
!
      CONTAINS

C *********************************************************************
C Routine to find nonzero hamiltonian matrix elements.
C Writen by J.Soler and P.Ordejon, June 1997
C Optimized by A. Garcia, August 2009
C **************************** INPUT **********************************
C logical negl         : Working option: Neglect interactions
C                        between non-overlaping orbitals
C real*8  cell(3,3)    : Supercell vectors CELL(IXYZ,IVECT)
C integer nsc(3)       : Num. of unit cells in each supercell direction
C integer na           : Number of atoms in supercell
C real*8  xa(3,na)     : Atomic positions in cartesian coordinates
C integer lasto(0:na)  : Last orbital of each atom in array iphorb
C integer lastkb(0:na) : Last KB proj. of each atom in array iphkb
C integer iphorb(no)   : Orbital index of each orbital in its atom,
C                        where no=lasto(na)
C integer iphkb(nkb)   : Index of each KB projector in its atom
C                        (negative). Here nkb=lastkb(na)
C integer isa(na)      : Species index of each atom
C set_xijo             : logical
C **************************** OUTPUT *********************************
C xijo                  : vectors between orbital centers
C integer nlhmax        : Maximum number
C                         of basis orbitals interacting, either directly
C                         or through a KB projector, with any orbital
C **************************** BEHAVIOR *******************************
C Equivalent pairs of atoms are assigned the same sparse index, i.e.
C if atoms i1 and j1 are equivalent, and so are i2 and j2, and the
C vector from i1 to i2 is the same as that from j1 to j2, and if
C listh(i,i1)=i2, then listh(i,j1)=j2. It is thus also implied that
C the matrix elements must be equal, i.e. H(i,i1)=H(i,i2).
C *********************************************************************
C
C  Modules
C
      subroutine hsparse( negl, cell, nsc, na, isa, xa,
     &                    lasto, lastkb, iphorb, iphkb,     
     &                    nlhmax, gamma,
     $                    set_xijo, folding)

      use precision
      use parallel,        only : Node, Nodes
      use parallelsubs,    only : GetNodeOrbs, GlobalToLocalOrb
      use atmfuncs,        only : rcut, nofis, nkbfis
      use atm_types,       only : nspecies, species_info, species
      use listsc_module,   only : listsc_init
      use sorting
      use neighbour,       only : jna=>jan, xij, r2ij, maxna=>maxnna
      use neighbour,       only : mneighb, reset_neighbour_arrays
      use sys,             only : die
      use alloc,           only : re_alloc, de_alloc
      use atomlist,        only : no_l
      use atomlist,        only : in_kb_orb_u_range
      use sparse_matrices, only : listhptr, numh, listh, xijo
      use radial,          only : rad_func
      use ldau_specs, only  : switch_ldau       ! Switch that determines whether
                                                !   and LDA+U simulation is 
                                                !   required or not

      implicit none

      integer,           intent(in) :: na
      integer,           intent(in) :: iphkb(:), iphorb(:)
      integer,           intent(in) :: isa(na), lastkb(0:na),
     $                                 lasto(0:na)
      integer,           intent(in) :: nsc(3)
      real(dp),          intent(in) :: cell(3,3), xa(3,na)
      logical,           intent(in) :: negl
      integer,        intent(out)   :: nlhmax
      logical, intent(in)           :: gamma
      logical, intent(in)           :: set_xijo
      logical, intent(out)    :: folding

      external               timer


      real(dp), parameter  :: tol = 1.0d-8   ! tolerance

      real(dp), allocatable     :: rkbmax(:) ! maximum KB radius of each species
      real(dp), allocatable     :: rorbmax(:) ! maximum ORB radius of each species
      real(dp), allocatable     :: rldaumax(:)! maximum LDA+U radius of each species
      integer      :: maxnkb  = 500          ! max no. of atoms with
                                             ! KB projectors which 
                                             ! overlap another
                                             ! atom's orbitals.
      integer
     &  ia, iio, ikb, ind, inkb, io, ioa, is, isel, 
     &  j, ja, jnat, jo, joa, js, 
     &  ka, kna, ko, koa, ks, 
     &  ncells, nlh, nna, nnkb, no, nua, nuo, nuotot


      real(dp)
     &  rci, rcj, rck, rij, rik, rjk,
     &  rmax, rmaxkb, rmaxo, rmaxldau

      logical, save :: warn1 = .false.

      real(dp), dimension(:), pointer :: rckb
      logical, dimension(:),  pointer :: conect
      integer, dimension(:),  pointer :: index
      integer, dimension(:),  pointer :: knakb
      integer, dimension(:),  pointer :: listhtmp

      logical :: connected
      type(species_info), pointer :: spp
      type(rad_func),     pointer :: pp

C -------------------------------------
#ifdef DEBUG
      call write_debug( '    PRE hsparse' )
#endif
C     Start time counter
      call timer( 'hsparse', 1 )

      ! Reduce check to only IONode
      if ( .not. warn1 ) warn1 = Node /= 0

C     Check size of internal arrays
      ncells = nsc(1) * nsc(2) * nsc(3)
      nua = na / ncells
      nuotot = lasto(nua)
      no = nuotot * ncells

!     Number of orbitals in local node
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

C     Allocate local arrays
      nullify(conect,listhtmp)
      call re_alloc( conect, 1, no, 'conect', 'hsparse' )
      call re_alloc( listhtmp, 1, no, 'listhtmp', 'hsparse' )

      ! We force nnkb == 0
      ! for the case where no LDAU and negl == .true.
      nnkb = 0
      
C     Allocate local arrays that depend on parameters
      nullify( knakb, rckb, index )
      call re_alloc( knakb, 1, maxnkb, 'knakb', 'hsparse' )
      call re_alloc( rckb, 1, maxnkb, 'rckb', 'hsparse' )
      call re_alloc( index, 1, maxna, 'index', 'hsparse' )
!
!     Find maximum radius of the KB and LDA+U projectors of each species (and total)
!     Also for orbitals (for possible future extra optimizations)
!     
      allocate(rkbmax(nspecies),rorbmax(nspecies),rldaumax(nspecies))
      rmaxo    = 0.0_dp
      rmaxkb   = 0.0_dp
      rmaxldau = 0.0_dp
      do is = 1, nspecies

         ! Species orbital range
         rorbmax(is) = 0.0_dp
         do io = 1, nofis(is)
            rorbmax(is) = max(rorbmax(is),rcut(is,io))
         enddo
         rmaxo = max(rmaxo, rorbmax(is))
         
         ! Species KB range
         rkbmax(is) = 0.0_dp
         do ikb = 1, nkbfis(is)
            rkbmax(is) = max(rkbmax(is),rcut(is,-ikb))
         enddo
         rmaxkb = max(rmaxkb, rkbmax(is))
         
         ! Species LDAU range
         rldaumax(is) = 0.0_dp
         if ( switch_ldau ) then
           spp => species(is)
           do io = 1, spp%n_pjldaunl
             pp => spp%pjldau(io)
             rldaumax(is) = max( rldaumax(is), pp%cutoff )
           enddo
         endif
         rmaxldau = max(rmaxldau, rldaumax(is))

      enddo

      
      ! Now calculate the actual extend of the basis orbitals
      if ( negl ) then
         rmax = rmaxo + rmaxldau
      else
         rmax = rmaxo + max(rmaxkb, rmaxldau)
      end if
      rmax = rmax * 2._dp

      
      isel = 0
C     Initialize internal data structures in neighb
      call mneighb( cell, rmax, na, xa, 0, isel, nna )

C     Initialize number of neighbour orbitals
      do io = 1,nuo 
        numh(io) = 0
      enddo 

C     Initialize vector switch only once
      do io = 1,no
        conect(io) = .false.
        listhtmp(io) = 0
      enddo
      do ia = 1 , na
         in_kb_orb_u_range(ia) = .false.
      end do
      
      folding = .false.

C------------------------------------C
C     Find number of non-zeros in H  C
C------------------------------------C
C     Loop on atoms in unit cell
      do ia = 1,nua

C       Find neighbour atoms within maximum range
        call mneighb( cell, rmax, na, xa, ia, isel, nna )
C       In case neighbor arrays have expanded
        call re_alloc( index, 1, maxna, 'index', 'hsparse' )

        
C       Order neighbours in a well defined way
        call ordvec( tol, 3, nna, xij, index )
        call iorder( jna, 1, nna, index )
        call order ( r2ij, 1, nna, index )

C       Loop on orbitals of atom ia
        do io = lasto(ia-1)+1,lasto(ia)
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then
            is = isa(ia)
            ioa = iphorb(io)
            rci = rcut(is,ioa)

            ! To be used in nlefsm
            ! Mark atoms whose KB projs overlap this process's
            ! base orbitals
            do kna = 1,nna
               ka = jna(kna)
               rik = sqrt( r2ij(kna) )
               ks = isa(ka)
               if ( rik < (rkbmax(ks) + rci)) then
                  in_kb_orb_u_range(ka) = .true.
               endif
            enddo

            if ( negl .and. switch_ldau ) then
              ! We ONLY add the LDAU projectors
              nnkb = 0
              do kna = 1,nna
                ka = jna(kna)
                rik = sqrt( r2ij(kna) )
                ks = isa(ka)
                rck = rldaumax(ks)
                if ( rci + rck > rik ) then
                  call extend_projector(nnkb, kna, rck)
                end if
              end do
            else if ( .not. negl ) then
              ! Find overlaping KB and LDA+U projectors
              nnkb = 0
              do kna = 1,nna
                ka = jna(kna)
                rik = sqrt( r2ij(kna) )
                ks = isa(ka)
                ! It is only necessary to check with
                ! the *largest* KB projector
                rck = rkbmax(ks)
                if ( rci + rck > rik ) then
                  call extend_projector(nnkb, kna, rck)
                endif
                if( switch_ldau ) then
                  rck = rldaumax(ks)
                  if ( rci + rck > rik ) then
                    call extend_projector(nnkb, kna, rck)
                  end if
                end if
              end do
            end if

C           Find orbitals connected by direct overlap or
C           through a KB or LDA+U projector
            do jnat = 1,nna
              ja = jna(jnat)
              js = isa(ja)
              rij = sqrt( r2ij(jnat) )
              do jo = lasto(ja-1)+1,lasto(ja)

                !If not yet connected (we only allow one connection, and
                !reserve space for it)
                if (.not.conect(jo)) then
                  joa = iphorb(jo)
                  rcj = rcut(js,joa)
C                 Find if there is direct overlap
                  if ( rci+rcj > rij ) then
                    conect(jo) = .true.
                  else
C                   Find if jo overlaps with a KB/LDAU projector
                    do inkb = 1,nnkb
                      rck = rckb(inkb)
                      kna = knakb(inkb)
                      rjk = sqrt( (xij(1,kna)-xij(1,jnat))**2 +
     &                            (xij(2,kna)-xij(2,jnat))**2 +
     &                            (xij(3,kna)-xij(3,jnat))**2 )
                      if ( rcj + rck > rjk ) then
                        conect(jo) = .true.
                        exit
                      endif
                    end do
                  endif
                  if ( conect(jo) ) then
                    numh(iio) = numh(iio) + 1
                    listhtmp(numh(iio)) = jo
                  endif
                endif
              enddo
            enddo

C           Restore conect array for next orbital io
            do j = 1,numh(iio)
              jo = listhtmp(j)
              conect(jo) = .false.
            enddo
          endif
        enddo
      enddo

C     Find optimum value for nlhmax
      if (nuo .gt. 0) then
        nlh = 0 
        do io = 1,nuo
          nlh = nlh + numh(io)
        enddo
      else    
        nlh = 0 
      endif     

      nlhmax = nlh
      ! I suspect that having nlhmax == 1 will result in errors
      ! It would be much better to not allow this "oversubscription".
      if ( nlhmax < 1 ) then
         call die('Sparse pattern is oversubscribed with nodes, 
     &please reduce number of nodes.')
      end if

      call re_alloc( listh, 1, max(1,nlhmax), 'listh',
     &     'sparseMat', SHRINK=.true. )
      if (set_xijo) then
         call re_alloc( xijo, 1, 3, 1, max(1,nlhmax),
     &   'xijo', 'sparseMat' )
      else
         call re_alloc( xijo, 1, 3, 1, 1, 'xijo', 'sparseMat' )
      endif

C     Now fill in the arrays, including xijo if using k-points
C     Set up listhptr
      listhptr(1) = 0
      do io = 2,nuo
        listhptr(io) = listhptr(io-1) + numh(io-1)
      enddo

C------------------------------------C
C     Find full H sparsity pattern   C
C------------------------------------C
C     Loop on atoms in unit cell
      do ia = 1,nua
C       Find neighbour atoms within maximum range
        call mneighb( cell, rmax, na, xa, ia, isel, nna )
C       In case neighbor arrays have expanded
        call re_alloc( index, 1, maxna, 'index', 'hsparse' )

C       Order neighbours in a well defined way
        call ordvec( tol, 3, nna, xij, index )
        call iorder( jna, 1, nna, index )
        call order(  r2ij, 1, nna, index )

C       Loop on orbitals of atom ia
        do io = lasto(ia-1)+1,lasto(ia)
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then
            numh(iio) = 0

            is = isa(ia)
            ioa = iphorb(io)
            rci = rcut(is,ioa)

C           Find overlaping KB and LDA+U projectors
            if ( negl .and. switch_ldau ) then
              ! We ONLY add the LDAU projectors
              nnkb = 0
              do kna = 1,nna
                ka = jna(kna)
                rik = sqrt( r2ij(kna) )
                ks = isa(ka)
                rck = rldaumax(ks)
                if ( rci + rck > rik ) then
                  call add_projector(nnkb, kna, rck)
                end if
              end do
            else if ( .not. negl ) then
              nnkb = 0
              do kna = 1,nna
                ka = jna(kna)
                rik = sqrt( r2ij(kna) )
                ks = isa(ka)
                ! It is only necessary to check with
                ! the *largest* KB/LDAU projector
                rck = rkbmax(ks)
                if ( rci + rck > rik ) then
                   call add_projector(nnkb, kna, rck)
                end if
                if( switch_ldau ) then
                   rck = rldaumax(ks)
                   if ( rci + rck > rik ) then
                      call add_projector(nnkb, kna, rck)
                   endif
                endif
              enddo
            endif
        
C           Find orbitals connected by direct overlap or
C           through a KB or LDA+U projector
            do jnat = 1,nna
              ja = jna(jnat)
              js = isa(ja)
              rij = sqrt( r2ij(jnat) )
              do jo = lasto(ja-1)+1,lasto(ja)

                joa = iphorb(jo)
                rcj = rcut(js,joa)

                ! Check connectivity
                connected = rci + rcj > rij
                
                if ( .not. connected ) then

! Find if jo overlaps with a KB projector
! nnkb is zero if no KB/LDAU projector exists
                   do inkb = 1,nnkb
                    rck = rckb(inkb)
                    kna = knakb(inkb)
                    rjk = sqrt( (xij(1,kna)-xij(1,jnat))**2 +
     .                   (xij(2,kna)-xij(2,jnat))**2 +
     .                   (xij(3,kna)-xij(3,jnat))**2 )
                    if ( rcj+rck > rjk ) then
                       connected = .true.
                       exit
                    end if
                  end do

                end if

                if (connected) then
                  if (conect(jo)) then
                    folding = .true.
                       
                    ! This test is now deferred to be able
                    ! to catch multiple images while avoiding
                    ! false positives (i.e., we test first
                    ! whether there is indeed a connection).

                    ! If already connected and using supercell, 
                    ! the latter might not be big enough...
                    ! We warn the user and keep the first instance
                    ! of xij (same behavior as the old xijorb, as
                    ! earlier jnats are closer)
                    ! Warn also if Gamma-point calculation, just
                    ! in case
                    if (.not.warn1) then
                      if ( gamma ) then
                        call check_cohp(io,jo)
                      else
                        write(6,'(/,a,2i6,a,/,a)')
     .                        'WARNING: orbital pair ',io,jo,
     .                        ' is multiply connected'
                      endif
                      warn1 = .true.
                    endif
                       
                  else
                     
                    conect(jo) = .true.
                    numh(iio) = numh(iio) + 1
                    ind = listhptr(iio)+numh(iio)
                    listh(ind) = jo
                    if (set_xijo) then
                      xijo(1:3,ind) = xij(1:3,jnat)
                    endif
                  endif
                endif
              enddo
            enddo
          
            ! Restore conect array for next orbital io
            do j = 1,numh(iio)
              jo = listh(listhptr(iio)+j)
              conect(jo) = .false.
            end do
           
          end if                 ! iio > 0
        end do                   ! io
      end do                     ! ia

!!        print "(a5,i3,a40,3i8)",
!!     $         "Node: ", Node, "in hsparse nuo, nuotot, nlhmax: ",
!!     $         nuo, nuotot, nlhmax

C     Initialize listsc
      call LISTSC_INIT( nsc, nuotot )


C     Deallocate local arrays
      call reset_neighbour_arrays( )

      call de_alloc( listhtmp, 'listhtmp', 'hsparse' )
      call de_alloc( conect, 'conect', 'hsparse' )
      call de_alloc( index, 'index', 'hsparse' )
      call de_alloc( knakb, 'knakb', 'hsparse' )
      call de_alloc( rckb, 'rckb', 'hsparse' )

      call timer( 'hsparse', 2 )
#ifdef DEBUG
      call write_debug( '    POS hsparse' )
#endif

      contains

      ! Possibly increase the arrays used for creating
      ! the data arrays
      ! ALSO add data to the arrays
      subroutine extend_projector(nnkb, kna, rck)
      integer, intent(inout) :: nnkb
      integer, intent(in)    :: kna
      real(dp), intent(in)   :: rck
      
      if ( nnkb ==  maxnkb ) then
         maxnkb = maxnkb + 10
         call re_alloc( knakb, 1, maxnkb, 'knakb',
     &        'hsparse', .true. )
         call re_alloc( rckb, 1, maxnkb, 'rckb',
     &        'hsparse', .true. )
      end if
      call add_projector(nnkb, kna, rck)
      end subroutine extend_projector

      subroutine add_projector(nnkb, kna, rck)
      integer, intent(inout) :: nnkb
      integer, intent(in)    :: kna
      real(dp), intent(in)   :: rck
      nnkb        = nnkb + 1
      knakb(nnkb) = kna
      rckb(nnkb)  = rck
      end subroutine add_projector
      
      end subroutine hsparse
      
      subroutine check_cohp(io,jo)
      use siesta_options, only: write_coop

      integer, intent(in) :: io, jo
 
      if (write_coop) then
         write(6,'(/,a,2i6,a,/,a)')
     .     'NOTE: orbital pair ',io,jo,
     .     ' (at least) is multiply connected.',
     .     'NOTE: Your COOP/COHP analysis might ' //
     $     'be affected by folding.'
         write(0,'(/,a,2i6,a,/,a)')
     .     'NOTE: orbital pair ',io,jo,
     .     ' (at least) is multiply connected.',
     .     'NOTE: Your COOP/COHP analysis might ' //
     $     'be affected by folding.'
      endif
      end subroutine check_cohp

      end module m_hsparse
