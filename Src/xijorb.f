C $Id: xijorb.f,v 1.7 1999/05/05 17:25:34 emilio Exp $

      subroutine xijorb( negl, scell, nua, na, xa,
     .                   lasto, lastkb, rco, rckb,
     .                   nomax, numh, listh, xijo )
C *********************************************************************
C Finds vectors between orbital centers.
C Writen by J.Soler. July 1997
C **************************** INPUT **********************************
C logical negl         : Working option: Neglect interactions
C                        between non-overlaping orbitals
C real*8  scell(3,3)   : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua          : Number of atoms in unit cell (first of list)
C integer na           : Total number of atoms
C real*8  xa(3,na)     : Atomic positions in cartesian coordinates
C integer lasto(0:na)  : Last orbital of each atom in array rco
C integer lastkb(0:na) : Last KB proj. of each atom in array rckb
C real*8  rco(no)      : Cutoff radius of basis orbitals 
C                        where no=lasto(na)
C real*8  rckb(nkb)    : Cutoff radius of KB projectors
C                        where nkb=lastkb(na)
C integer nomax        : First/second dimension of listh/xijo
C integer numh(no)        : Number of nonzero elements of each row of the
C                           Hamiltonian matrix between atomic orbitals
C integer listh(nomax,no) : Nonzero Hamiltonian-matrix element 
C                           column indexes for each matrix row
C                           (only if input_nomax.ge.output_nomax)
C **************************** OUTPUT *********************************
C real*8  xijo(3,nomax,no) : Vectors between orbital centers
C *********************************************************************
      implicit          none
      integer           na, nomax
      integer           lastkb(0:na), lasto(0:na), listh(nomax,*),
     .                  nua, numh(*)
      double precision  scell(3,3), rckb(*), rco(*), xa(3,na),
     .                  xijo(3,nomax,*)
      logical           negl
      external          chkdim, neighb, timer

C Internal variables -----------------
C maxna  = maximum number of neighbour atoms of any atom
C maxnkb = maximum number of neighbour KB projectors of any orbital
C maxo   = maximum number of basis orbitals
      integer maxna, maxnkb, maxo
      parameter ( maxna  =  1000 )
      parameter ( maxnkb =  2000 )
      parameter ( maxo   = 20000 )

      integer
     .  ia, ikb, inkb, io, isel, 
     .  j, ja, jna, jana(maxna), jnao(maxo), jo,
     .  ka, kna, knakb(maxnkb), ko,
     .  nkb, nna, nnkb, no

      double precision
     .  r2ij(maxna), rci, rcj, rck, rcnkb(maxnkb), rij, rik, rjk,
     .  rmax, rmaxkb, rmaxo, xija(3,maxna)

      logical
     .  conect, warn1, warn2
     
      save warn1, warn2
      data warn1, warn2 /2*.false./
C -------------------------------------

C     Start time counter
*     call timer( 'xijorb', 1 )

C     Check size of internal arrays
      nkb = lastkb(na)
      no = lasto(na)
      call chkdim( 'xijorb', 'maxo', maxo, no, 1 )

C     Find maximum range of basis orbitals and KB projectors
      rmaxo = 0.d0
      rmaxkb = 0.d0
      do io = 1,no
        rmaxo = max( rmaxo, rco(io) )
      enddo
      do ikb = 1,nkb
        rmaxkb = max( rmaxkb, rckb(ikb) )
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
      call neighb( scell, rmax, na, xa, ia, isel,
     .             nna, jana, xija, r2ij )

c     Initialize vector jnao only once
      do io = 1,no
        jnao(io) = 0
      enddo

c     Loop on atoms (only within unit cell)
      do ia = 1,nua

c       Find neighbour atoms within maximum range
        nna = maxna
        call neighb( scell, rmax, na, xa, ia, isel,
     .               nna, jana, xija, r2ij )
        call chkdim( 'xijorb', 'maxna', maxna, nna, 1 )

c       Loop on orbitals of atom ia
        do io = lasto(ia-1)+1,lasto(ia)
          rci = rco(io)

c         Find overlaping KB projectors
          if (.not.negl) then
            nnkb = 0
            do kna = 1,nna
              ka = jana(kna)
              rik = sqrt( r2ij(kna) )
              do ko = lastkb(ka-1)+1,lastkb(ka)
                rck = rckb(ko)
                if (rci+rck .gt. rik) then
                  nnkb = nnkb + 1
                  call chkdim( 'hsparse', 'maxnkb', maxnkb, nnkb, 1 )
                  knakb(nnkb) = kna
                  rcnkb(nnkb) = rck
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
              rcj = rco(jo)

c             Find if orbitals io and jo are connected
              conect = .false.
c             Find if there is direct overlap
              if (rci+rcj .gt. rij) then
                conect = .true.
              elseif (.not.negl) then
c               Find if jo overlaps with a KB projector
                do inkb = 1,nnkb
                  rck = rcnkb(inkb)
                  kna = knakb(inkb)
                  rjk = sqrt( (xija(1,kna)-xija(1,jna))**2 +
     .                        (xija(2,kna)-xija(2,jna))**2 +
     .                        (xija(3,kna)-xija(3,jna))**2 )
                  if (rcj+rck .gt. rjk) then
                    conect = .true.
                    goto 50
                  endif
                enddo
   50           continue
              endif

c             Add to list of connected orbitals
              if (conect) then
                if (jnao(jo) .eq. 0) then
                  jnao(jo) = jna
                else
                  if (.not.warn1) then
                    write(6,'(/,a,2i6,a,/)')
     .                'xijorb: WARNING: orbital pair ', io, jo,
     .                ' is multiply connected'
                    warn1 = .true.
                  endif
                  kna = jnao(jo)
                  if (r2ij(jna) .lt. r2ij(kna)) jnao(jo) = jna
                endif
              endif
            enddo
          enddo

c         Copy interatomic vectors into xijo
          do j = 1,numh(io)
            jo = listh(j,io)
            jna = jnao(jo)
            if (jna .eq. 0) then
              if (.not.warn2) then
                write(6,'(/,a,2i6,a,/)')
     .            'xijorb: WARNING: orbital pair ', io, jo,
     .            ' is in listh but not really connected'
                warn2 = .true.
              endif
              xijo(1,j,io) = 0.d0
              xijo(2,j,io) = 0.d0
              xijo(3,j,io) = 0.d0
            else
              xijo(1,j,io) = xija(1,jna)
              xijo(2,j,io) = xija(2,jna)
              xijo(3,j,io) = xija(3,jna)
            endif
          enddo

c         Restore jnao array for next orbital io
          do jna = 1,nna
            ja = jana(jna)
            do jo = lasto(ja-1)+1,lasto(ja)
              jnao(jo) = 0
            enddo
          enddo

        enddo
      enddo

*     call timer( 'xijorb', 2 )
      end
