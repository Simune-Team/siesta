      subroutine hsparse( negl, cell, na, isa, xa,
     .                    lasto, lastkb, iphorb, iphkb,
     .                    nomax, numh, listh )
C *********************************************************************
C Routine to find nonzero hamiltonian matrix elements.
C Writen by J.Soler and P.Ordejon, June 1997
C **************************** INPUT **********************************
C logical negl         : Working option: Neglect interactions
C                        between non-overlaping orbitals
C integer na           : Total number of atoms
C real*8  cell(3,3)    : Unit cell vectors CELL(IXYZ,IVECT)
C real*8  xa(3,na)     : Atomic positions in cartesian coordinates
C integer lasto(0:na)  : Last orbital of each atom in array iphorb
C integer lastkb(0:na) : Last KB proj. of each atom in array iphkb
C integer iphorb(no)   : Orbital index of each orbital in its atom,
C                        where no=lasto(na)
C integer iphkb(nkb)   : Index of each KB projector in its atom
C                        (negative). Here nkb=lastkb(na)
C integer isa(na)      : Species index of each atom
C **************************** INPUT and OUTPUT ***********************
C integer nomax : Input : First dimension of listh
C                 Output: Actual value needed for nomax (Maximum number
C                         of basis orbitals interacting, either directly
C                         or through a KB projector, with any orbital)
C **************************** OUTPUT *********************************
C integer numh(no)        : Number of nonzero elements of each row of the
C                           Hamiltonian matrix between atomic orbitals
C integer listh(nomax,no) : Nonzero Hamiltonian-matrix element 
C                           column indexes for each matrix row
C                           (only if input_nomax.ge.output_nomax)
C **************************** BEHAVIOR *******************************
C Equivalent pairs of atoms are assigned the same sparse index, i.e.
C if atoms i1 and j1 are equivalent, and so are i2 and j2, and the
C vector from i1 to i2 is the same as that from i2 to j2, and if
C listh(i,i1)=i2, then listh(i,j1)=j2. It is thus also implied that
C the matrix elements must be equal, i.e. H(i,i1)=H(i,i2).
C *********************************************************************
      implicit          none
      integer           na, nomax
      integer           iphkb(*), iphorb(*), isa(na),
     .                  lastkb(0:na), lasto(0:na), listh(nomax,*),
     .                  numh(*)
      double precision  cell(3,3), rcut, xa(3,na)
      logical           negl
      external          chkdim, neighb, ordvec, rcut, timer

C Internal variables -----------------
C maxna  = maximum number of neighbour atoms of any atom
C maxnkb = maximum number of neighbour KB projectors of any orbital
C maxo   = maximum number of basis orbitals
C tol    = tolerance for comparing vector-coordinates
      integer maxna, maxnkb, maxo
      double precision tol
      parameter ( maxna  =   500 )
      parameter ( maxnkb =  1000 )
      parameter ( maxo   = 10000 )
      parameter ( tol    = 1.d-8 )

      integer
     .  ia, ikb, index(maxna), inkb, io, ioa, is, isel, 
     .  j, ja, jna, jana(maxna), jo, joa, js, 
     .  ka, kna, knakb(maxnkb), ko, koa, ks,
     .  nna, nnkb, no

      double precision
     .  r2ij(maxna), rci, rcj, rck, rckb(maxnkb), rij, rik, rjk,
     .  rmax, rmaxkb, rmaxo, xij(3,maxna)

      logical
     .  conect(maxo)
C -------------------------------------

C     Start time counter
*     call timer( 'hsparse', 1 )

C     Check size of internal arrays
      no = lasto(na)
      call chkdim( 'hsparse', 'maxo', maxo, no, 1 )

C     Find maximum range of basis orbitals and KB projectors
      rmaxo = 0.d0
      rmaxkb = 0.d0
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
        rmax = 2.d0 * rmaxo
      else
        rmax = 2.d0 * (rmaxo+rmaxkb)
      endif

c     Initialize neighb subroutine
      isel = 0
      ia = 0
      nna = maxna
      call neighb( cell, rmax, na, xa, ia, isel,
     .             nna, jana, xij, r2ij )

c     Initialize number of neighbour orbitals
      do io= 1,no 
        numh(io)=0
      enddo 

c     Initialize vector switch only once
      do io = 1,no
        conect(io) = .false.
      enddo

c     Loop on atoms
      do ia = 1,na

c       Find neighbour atoms within maximum range
        nna = maxna
        call neighb( cell, rmax, na, xa, ia, isel,
     .               nna, jana, xij, r2ij )
        call chkdim( 'hsparse', 'maxna', maxna, nna, 1 )

c       Order neighbours in a well defined way
        call ordvec( tol, 3, nna, xij, index )
        call iorder( jana, 1, nna, index )
        call order(  r2ij, 1, nna, index )

c       Loop on orbitals of atom ia
        do io = lasto(ia-1)+1,lasto(ia)
          is = isa(ia)
          ioa = iphorb(io)
          rci = rcut(is,ioa)

c         Find overlaping KB projectors
          if (.not.negl) then
            nnkb = 0
            do kna = 1,nna
              ka = jana(kna)
              rik = sqrt( r2ij(kna) )
              do ko = lastkb(ka-1)+1,lastkb(ka)
                ks = isa(ka)
                koa = iphkb(ko)
                rck = rcut(ks,koa)
                if (rci+rck .gt. rik) then
                  nnkb = nnkb + 1
                  call chkdim( 'hsparse', 'maxnkb', maxnkb, nnkb, 1 )
                  knakb(nnkb) = kna
                  rckb(nnkb) = rck
                endif
              enddo
            enddo
          endif

c         Find orbitals connected by direct overlap or
c         through a KB projector
          do jna = 1,nna
            ja = jana(jna)
            rij = sqrt( r2ij(jna) )
            do jo = lasto(ja-1)+1,lasto(ja)
c             If not yet connected
              if (.not.conect(jo)) then
                js = isa(ja)
                joa = iphorb(jo)
                rcj = rcut(js,joa)
c               Find if there is direct overlap
                if (rci+rcj .gt. rij) then
                  conect(jo) = .true.
                elseif (.not.negl) then
c                 Find if jo overlaps with a KB projector
                  do inkb = 1,nnkb
                    rck = rckb(inkb)
                    kna = knakb(inkb)
                    rjk = sqrt( (xij(1,kna)-xij(1,jna))**2 +
     .                          (xij(2,kna)-xij(2,jna))**2 +
     .                          (xij(3,kna)-xij(3,jna))**2 )
                    if (rcj+rck .gt. rjk) then
                      conect(jo) = .true.
                      goto 50
                    endif
                  enddo
   50             continue
                endif
                if (conect(jo)) then
                  numh(io) = numh(io) + 1
                  if (numh(io).le.nomax) listh(numh(io),io) = jo
                endif
              endif
            enddo
          enddo

c         Restore conect array for next orbital io
          if (numh(io).gt.nomax) then
            do jna = 1,nna
              ja = jana(jna)
              do jo = lasto(ja-1)+1,lasto(ja)
                conect(jo) = .false.
              enddo
            enddo
          else
            do j = 1,numh(io)
              jo = listh(j,io)
              conect(jo) = .false.
            enddo
          endif

        enddo
      enddo

C     Find optimum value for nomax
      nomax = 0
      do io = 1,no
        nomax = max( nomax, numh(io) )
      enddo

*     call timer( 'hsparse', 2 )
      end
