      module atomlist

      use precision
      use alloc
      use ionew, only: IOnode

      implicit none

!
!     Instead of "generic" na, no, and nokb, we use:
!
! For "supercell" (intended for k-point calcs)
      integer, save                  :: na_s         ! Number of atoms
      integer, save                  :: no_s         ! Number of orbitals
      integer, save                  :: nokb_s       ! Number of KB projs

! Same for "unit", or "real" cell:
      integer, save                  ::  na_u, no_u, nokb_u

! Here 'na' is a generic number. It could be na_u or na_s, depending
! on whether we need a supercell or not.

C integer isa(na)           : Species index of each atom
C integer lasto(0:na)       : Position of last orbital of each atom
C integer lastkb(0:na)      : Position of last KB proj. of each atom
C integer iza(na)           : Atomic number of each atom
C real*8 amass(na)          : Atomic mass of each atom
C real*8 qa(na)             : Neutral atom charge of each atom

      integer, pointer, save  :: isa(:) ! 
      integer, pointer, save  :: iza(:) ! 
      integer, pointer, save  :: lasto(:) ! 
      integer, pointer, save  :: lastkb(:)
      real(dp), pointer, save  :: amass(:), qa(:)

      real(dp), pointer, save  :: xa(:,:)
!        Atomic coordinates
      real(dp), pointer, save  :: xalast(:,:)
!        Atomic coordinates (it doesn't really belong here)

      integer, pointer, save           :: indxua(:)
!        Index of equivalent atom in "u" cell

      real(dp), save         :: rmaxv  ! Max cutoff for local pot Vna
      real(dp), save         :: rmaxo  ! Max cuoff for atomic orbitals
      real(dp), save         :: rmaxkb ! Max cuoff for KB projectors

      real(dp), save         :: qtot ! Total number of electrons


      integer, pointer, save  :: iaorb(:)  ! Atomic index of each orbital
      integer, pointer, save  :: iphorb(:) ! Orbital index of each 
                                               ! orbital in its atom
      real(dp), pointer, save :: Datm(:) ! Neutral atom charge 
                                             ! of each orbital
      real(dp), pointer, save :: rco(:)  ! Cutoff radius of each orbital

      integer, pointer, save           :: indxuo(:)
!        Index of equivalent orbital in "u" cell

      integer, pointer, save     :: iakb(:)
!         Atomic index of each KB projector
      integer, pointer, save     :: iphKB(:)
!         Index of each KB projector in its atom (negative)
      real(dp), pointer, save   :: rckb(:)
!         Cutoff radius of each KB projector
!
      CONTAINS

      subroutine initatomlists

C Routine to initialize the atomic lists.
C
      use  atmfuncs, only: nofis, nkbfis, izofis, massfis,
     $                     rcut, atmpopfio

      integer  ia, io, is, nkba, noa, nol, nokbl, ioa, ikb
!
      na_s = na_u
      no_s = no_u
      nokb_s = nokb_u

      nullify(iaorb, indxuo, iphorb, Datm, rco)
      call realloc(iaorb, 1, no_u, routine='initatomlists')
      call realloc(indxuo, 1, no_u, routine='initatomlists')
      call realloc(iphorb, 1, no_u, routine='initatomlists')
      call realloc(Datm, 1, no_u, routine='initatomlists')
      call realloc(rco, 1, no_u, routine='initatomlists')
!
!
      nullify(iaKB, iphKB, rckb)
      call realloc(iaKB, 1, nokb_u, routine='initatomlists')
      call realloc(iphKB, 1, nokb_u, routine='initatomlists')
      call realloc(rckb, 1, nokb_u, routine='initatomlists')

c Initialize atomic lists
      nol = 0
      nokbl = 0
      qtot = 0._dp
      rmaxv  = 0._dp
      rmaxo  = 0._dp
      rmaxkb = 0._dp
      lasto(0) = 0
      lastkb(0) = 0
      do ia = 1,na_u
        is = isa(ia)
        noa  = nofis(is)
        nkba = nkbfis(is)
        lasto(ia)  = lasto(ia-1)  + noa
        lastkb(ia) = lastkb(ia-1) + nkba
        rmaxv = max( rmaxv, rcut(is,0) )
        iza(ia) = izofis(is)
        amass(ia) = massfis(is)
        qa(ia) = 0.0_dp
        do io = 1,noa
          nol = nol + 1
          rmaxo = max( rmaxo, rcut(is,io) )
          iaorb(nol) = ia
          iphorb(nol) = io
          Datm(nol) = atmpopfio(is,io)
          qa(ia) = qa(ia) + Datm(nol)
          qtot = qtot + Datm(nol)
        enddo
        do io = 1,nkba
          nokbl = nokbl + 1
          rmaxkb = max( rmaxkb, rcut(is,-io) ) 
          iaKB(nokbl) = ia
          iphKB(nokbl) = -io
        enddo
      enddo

! Find rco and rckb .............................

      do ia = 1,na_u
        is = isa(ia)
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rco(io) = rcut(is,ioa)
        enddo
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rckb(ikb) = rcut(is,ioa)
        enddo
      enddo
      
      if (IOnode)
     $   write(6,'(a,3(1x,i4))')
     $   'initatomlists: Number of atoms, orbitals, and projectors: ',
     $     na_u, no_u, nokb_u

      end subroutine initatomlists


      subroutine superc( ucell, scell, nsc)

C Finds the supercell required to avoid multiple image overlaps,
C and expands arrays from unit cell to supercell
C Written by J.M.Soler. August 1998.
! Rewritten Alberto Garcia, May 2000.

      implicit none
      integer, intent(in)  :: nsc(3)      ! Diagonal elements of mscell
      real(dp), intent(in) :: ucell(3,3)  ! Unit cell vectors
      real(dp), intent(out) :: scell(3,3) ! Supercell vectors

      external superx

C Internal variables
      integer           ia, io, iua, iuo, ja, ncells,
     $                  na, no, nokb

!
!      Find number of cells, atoms and orbitals in supercell
      ncells = nsc(1) * nsc(2) * nsc(3)
      na    = na_u   * ncells
      no    = no_u   * ncells
      nokb  = nokb_u * ncells

!
!     Reallocate arrays if needed
!
      if (na.gt.na_s) then
        call realloc(indxua, 1, na, routine='superc',copy=.true.)
        call realloc(isa, 1, na, routine='superc',copy=.true.)
        call realloc(iza, 1, na, routine='superc',copy=.true.)
        call realloc(lastkb, 0, na, routine='superc',copy=.true.)
        call realloc(lasto, 0, na, routine='superc',copy=.true.)
        call realloc(qa, 1, na, routine='superc',copy=.true.)
        call realloc(xa, 1,3, 1,na, routine='superc',copy=.true.)
        call realloc(xalast, 1,3, 1,na, routine='superc',copy=.true.)
      endif

      na_s  = na

C Find supercell vectors and atomic coordinates in supercell 
      call superx( ucell, nsc, na_u, na_s, xa, scell )

C Find indxua and expand isa, iza, lasto and lastkb to supercell 
      do ia = 1,na_s
        ja = mod(ia-1,na_u) + 1
        indxua(ia) = ja
        isa(ia)    = isa(ja)
        iza(ia)    = iza(ja)
        lasto(ia)  = lasto(ia-1)  + lasto(ja)  - lasto(ja-1)
        lastkb(ia) = lastkb(ia-1) + lastkb(ja) - lastkb(ja-1)
      enddo

! Reallocate orbital arrays

      if (no.gt.no_s) then
        call realloc(iaorb, 1,no, routine='superc',copy=.true.)
        call realloc(indxuo, 1,no, routine='superc',copy=.true.)
        call realloc(iphorb, 1,no, routine='superc',copy=.true.)
        call realloc(Datm, 1,no, routine='superc',copy=.true.)
        call realloc(rco, 1,no, routine='superc',copy=.true.)
      endif

      no_s = no

C Find indxuo and expand iaorb, iphorb, and rco 
      do io = 1,no_s
        indxuo(io) = mod(io-1,no_u) + 1
      enddo
      do ia = 1,na_s
        do io = lasto(ia-1)+1,lasto(ia)
          iuo = indxuo(io)
          iaorb(io)  = ia
          iphorb(io) = iphorb(iuo)
          rco(io)    = rco(iuo)
        enddo
      enddo

! Reallocate projector arrays

      if (nokb .gt. nokb_s) then
        call realloc(iaKB, 1,nokb, routine='superc',copy=.true.)
        call realloc(iphKB, 1,nokb, routine='superc',copy=.true.)
        call realloc(rckb, 1,nokb, routine='superc',copy=.true.)
      endif

      nokb_s = nokb

C Expand iakb and iphKB and rckb

      do ia = 1,na_s
        iua = indxua(ia)
        iuo = lastkb(iua-1)
        do io = lastkb(ia-1)+1,lastkb(ia)
          iuo = iuo + 1
          iakb(io)  = ia
          iphKB(io) = iphKB(iuo)
          rckb(io)  = rckb(iuo)
        enddo
      enddo

      if (IOnode .and. ncells.gt.1) then
         write(6,'(/,a,i6,a,i6,a,i6,a,i8)')
     .    'superc: Internal auxiliary supercell:',
     .     nsc(1), ' x', nsc(2), ' x', nsc(3), '  =', ncells

         write(6,'(a,3(1x,i4))')
     $     'superc: Number of atoms, orbitals, and projectors: ',
     $     na_s, no_s, nokb_s
           
      endif

      end subroutine superc


      end module atomlist



