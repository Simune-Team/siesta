! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_hsparse
!
!     We need an explicit interface since assumed-shape arrays are used
!
      CONTAINS

C *********************************************************************
C Routine to find nonzero hamiltonian matrix elements.
C Writen by J.Soler and P.Ordejon, June 1997
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
C **************************** INPUT and OUTPUT ***********************
C integer nlhmax: Input : Dimension of listh
C                 Output: Actual value needed for nlhmax (Maximum number
C                         of basis orbitals interacting, either directly
C                         or through a KB projector, with any orbital)
C **************************** OUTPUT *********************************
C integer numh(*)       : Number of nonzero elements of each row of the
C                         Hamiltonian matrix between atomic orbitals
C integer listhptr(*)   : Pointer to where each row of listh starts - 1
C                         The reason for pointing to the element before
C                         the first one is so that when looping over the
C                         elements of a row there is no need to shift by
C                         minus one.
C integer listh(nlhmax) : Nonzero Hamiltonian-matrix element 
C                         column indexes for each matrix row
C                         For parallel execution, listh contains the
C                         elements for rows that involve
C                         any locally stored orbitals. In the case
C                         where parallelisation is over K points then
C                         the full listh matrix is needed on every
C                         Node.
C gamma                 : logical
C xijo                  : vectors between orbital centers
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
     .                    lasto, lastkb, iphorb, iphkb,
     .                    nlhmax, numh, listhptr, listh, gamma, xijo)
      use precision
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GetNodeOrbs, GlobalToLocalOrb
      use atmfuncs,      only : rcut
      use listsc_module, only : listsc_init
      use sorting
      use m_neighbors,   only : maxna, jna, xij, r2ij
      use sys,           only : die
      use alloc,         only : re_alloc, de_alloc

      implicit none

      integer,  intent(in)    :: na
      integer,  intent(in)    :: iphkb(:), iphorb(:)
      integer,  intent(in)    :: isa(na), lastkb(0:na),
     $                           lasto(0:na)
      integer,  intent(in)    :: nsc(3)
      real(dp), intent(in)    :: cell(3,3), xa(3,na)
      logical,  intent(in)    :: negl

      integer,  intent(inout) :: nlhmax
      integer,  intent(out)   :: listhptr(:)
      integer,  intent(out)   :: numh(:)

      integer,  pointer       :: listh(:)
      logical, intent(in)     :: gamma
      real(dp),  pointer      :: xijo(:,:)

      external                   neighb, timer, memory


      real(dp), parameter        :: tol = 1.0d-8   ! tolerance for comparing vector-coordinates


      integer      :: maxnkb  = 1000         ! max no. of neighboring KB projectors

      integer
     .  ia, iio, ikb, ind, inkb, io, ioa, is, isel, 
     .  j, ja, jnat, jo, joa, js, 
     .  ka, kna, ko, koa, ks, 
     .  ncells, nlh, nna, nnkb, no, nua, nuo, nuotot


      real(dp)
     .  rci, rcj, rck, rij, rik, rjk,
     .  rmax, rmaxkb, rmaxo

      logical, save :: warn1 = .false.

      real(dp), dimension(:), pointer :: rckb
      logical, dimension(:),  pointer :: conect
      integer, dimension(:),  pointer :: index, knakb, listhtmp


C -------------------------------------

C Start time counter
      call timer( 'hsparse', 1 )

C Check size of internal arrays
      ncells = nsc(1) * nsc(2) * nsc(3)
      nua = na / ncells
      nuotot = lasto(nua)
      no = nuotot * ncells
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)  !Numero de orbitales del proceso MPI

C Allocate local arrays
      nullify(conect,listhtmp)
      call re_alloc(conect,1,no,name="conect",routine="hsparse")
      call re_alloc(listhtmp,1,no,name="listhtmp",routine="hsparse")

C Find maximum range of basis orbitals and KB projectors
      rmaxo = 0.0d0
      rmaxkb = 0.0d0
      do ia = 1,na
        is = isa(ia)
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphkb(ikb)
          rmaxkb = max( rmaxkb, rcut(is,ioa) )
        enddo
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rmaxo = max( rmaxo, rcut(is,ioa) )
        enddo
      enddo

      if (negl) then
        rmax = 2.0d0 * rmaxo
      else
        rmax = 2.0d0 * (rmaxo+rmaxkb)
      endif

      ! Allocate local arrays that depend on parameters
      call re_alloc(knakb,1,maxnkb, name="knakb")
      call re_alloc(rckb,1,maxnkb, name="rckb")
      call re_alloc(index,1,maxna, name="index")

      isel = 0
      ! Initialize internal data structures in neighb
      ! Note that this is necessary if, for example,
      ! we switch from considering unit-cell numbering
      ! to supercell numbering...
      nna = maxna
      call neighb( cell, rmax, na, xa, 0, isel,
     .             nna, jna, xij, r2ij )

C Initialize number of neighbour orbitals
      do io = 1,nuo 
        numh(io) = 0
      enddo 

C Initialize vector switch only once
      do io = 1,no
        conect(io) = .false.
        listhtmp(io) = 0
      enddo

C----------------------------------
C  Find number of non-zeros in H  _
C----------------------------------
C Loop on atoms in unit cell

      do ia = 1,nua

C Find neighbour atoms within maximum range
        nna = maxna
        call neighb( cell, rmax, na, xa, ia, isel,
     .               nna, jna, xij, r2ij )

C Order neighbours in a well defined way
          call ordvec( tol, 3, nna, xij, index )
          call iorder( jna, 1, nna, index )
          call order ( r2ij, 1, nna, index )

C Loop on orbitals of atom ia
          do io = lasto(ia-1)+1,lasto(ia)

            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then

              is = isa(ia)
              ioa = iphorb(io)
              rci = rcut(is,ioa)

C Find overlaping KB projectors
              if (.not.negl) then
                nnkb = 0
                do kna = 1,nna
                  ka = jna(kna)
                  rik = sqrt( r2ij(kna) )
                  do ko = lastkb(ka-1)+1,lastkb(ka)
                    ks = isa(ka)
                    koa = iphkb(ko)
                    rck = rcut(ks,koa)
                    if (rci+rck .gt. rik) then
C Check maxnkb - if too small then increase array sizes
                      if (nnkb.eq.maxnkb) then
                         maxnkb = maxnkb + 100
                         call re_alloc(knakb,1,maxnkb,
     $                        name="knakb",copy=.true.)
                         call re_alloc(rckb,1,maxnkb,
     $                        name="rckb",copy=.true.)
                      endif
                      nnkb = nnkb + 1
                      knakb(nnkb) = kna
                      rckb(nnkb) = rck
                    endif
                  enddo
                enddo
              endif

C Find orbitals connected by direct overlap or
C through a KB projector
              do jnat = 1,nna
                ja = jna(jnat)
                rij = sqrt( r2ij(jnat) )
                do jo = lasto(ja-1)+1,lasto(ja)

C If not yet connected 
                  if (.not.conect(jo)) then
                    js = isa(ja)
                    joa = iphorb(jo)
                    rcj = rcut(js,joa)
C Find if there is direct overlap
                    if (rci+rcj .gt. rij) then
                      conect(jo) = .true.
                    elseif (.not.negl) then
C Find if jo overlaps with a KB projector
                      do inkb = 1,nnkb
                        rck = rckb(inkb)
                        kna = knakb(inkb)
                        rjk = sqrt( (xij(1,kna)-xij(1,jnat))**2 +
     .                              (xij(2,kna)-xij(2,jnat))**2 +
     .                              (xij(3,kna)-xij(3,jnat))**2 )
                        if (rcj+rck .gt. rjk) then
                          conect(jo) = .true.
                          goto 50
                        endif
                      enddo
   50                 continue
                    endif
                    if (conect(jo)) then
                      numh(iio) = numh(iio) + 1
                      listhtmp(numh(iio)) = jo
                    endif
                  endif
                enddo
              enddo

C Restore conect array for next orbital io
              do j = 1,numh(iio)
                jo = listhtmp(j)
                conect(jo) = .false.
              enddo
            endif
          enddo
      enddo


C Find optimum value for nlhmax
      if (nuo .gt. 0) then
        nlh = 0 
        do io = 1,nuo
          nlh = nlh + numh(io)
        enddo
      else    
        nlh = 0 
      endif     

      nlhmax = nlh
      call re_alloc(listh,1,nlhmax,name='listh',routine='hsparse')
      if (.not. gamma) then
         call re_alloc(xijo,1,3,1,nlhmax,name='xijo',routine='hsparse')
      else
         call re_alloc(xijo,1,3,1,1,name='xijo',routine='hsparse')
      endif

!
!     Now fill in the arrays, including xijo if using k-points
!
C Set up listhptr
        do io = 1,nuo
          if (io.gt.1) then
            listhptr(io) = listhptr(io-1) + numh(io-1)
          else
            listhptr(io) = 0
          endif
        enddo
C---------------------------------
C  Find full H sparsity pattern  -
C---------------------------------
C Loop on atoms in unit cell
        do ia = 1,nua

C Find neighbour atoms within maximum range
          nna = maxna
          call neighb( cell, rmax, na, xa, ia, isel,
     .                 nna, jna, xij, r2ij )

C Order neighbours in a well defined way
          call ordvec( tol, 3, nna, xij, index )
          call iorder( jna, 1, nna, index )
          call order(  r2ij, 1, nna, index )

C Loop on orbitals of atom ia
          do io = lasto(ia-1)+1,lasto(ia)

            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then
              numh(iio) = 0

              is = isa(ia)
              ioa = iphorb(io)
              rci = rcut(is,ioa)

C Find overlaping KB projectors
              if (.not.negl) then
                nnkb = 0 
                do kna = 1,nna
                  ka = jna(kna)
                  rik = sqrt( r2ij(kna) )
                  do ko = lastkb(ka-1)+1,lastkb(ka)
                    ks = isa(ka)
                    koa = iphkb(ko)
                    rck = rcut(ks,koa)
                    if (rci+rck .gt. rik) then
                      nnkb = nnkb + 1
                      knakb(nnkb) = kna
                      rckb(nnkb) = rck
                    endif
                  enddo
                enddo
              endif
          
C Find orbitals connected by direct overlap or
C through a KB projector
              do jnat = 1,nna
                ja = jna(jnat)
                rij = sqrt( r2ij(jnat) )
                do jo = lasto(ja-1)+1,lasto(ja)

                   if (conect(jo)) then

                      if (.not. gamma) then
                         ! If already connected and using supercell, 
                         ! the latter might not be big enough...
                         ! We warn the user and keep the first instance
                         ! of xij (same behavior as the old xijorb, as
                         ! earlier jnats are closer)
                         if (.not.warn1) then
                            if (Node.eq.0) then
                               write(6,'(/,a,2i6,a,/)')
     .                              'WARNING: orbital pair ',io,jo,
     .                              ' is multiply connected'
                            endif
                            warn1 = .true.
                         endif
                      endif

                   else         ! not conect(jo)

                    js = isa(ja)
                    joa = iphorb(jo)
                    rcj = rcut(js,joa)
C Find if there is direct overlap
                    if (rci+rcj .gt. rij) then
                      conect(jo) = .true.
                    elseif (.not.negl) then
C Find if jo overlaps with a KB projector
                      do inkb = 1,nnkb
                        rck = rckb(inkb)
                        kna = knakb(inkb)
                        rjk = sqrt( (xij(1,kna)-xij(1,jnat))**2 +
     .                              (xij(2,kna)-xij(2,jnat))**2 +
     .                              (xij(3,kna)-xij(3,jnat))**2 )
                        if (rcj+rck .gt. rjk) then
                          conect(jo) = .true.
                          goto 55
                        endif
                      enddo
   55                 continue
                    endif
                    if (conect(jo)) then
                      numh(iio) = numh(iio) + 1
                      ind = listhptr(iio)+numh(iio)
                      listh(ind) = jo
                      if (.not. gamma) then
                         xijo(1:3,ind) = xij(1:3,jnat)
                      endif
                    endif
                   endif
                enddo
              enddo

C Restore conect array for next orbital io

              do j = 1,numh(iio)
                 jo = listh(listhptr(iio)+j)
                 conect(jo) = .false.
              enddo   
     
            endif     ! iio > 0
          enddo       ! io
        enddo         ! ia
   
C Initialize listsc
      call listsc_init( nsc, nuotot )

C Deallocate local arrays

      call de_alloc(listhtmp,name="listhtmp")
      call de_alloc(conect,name="conect")
      call de_alloc(index,name="index")
      call de_alloc(knakb,name="knakb")
      call de_alloc(rckb,name="rckb")

      call timer( 'hsparse', 2 )

      end subroutine hsparse

      end module m_hsparse





