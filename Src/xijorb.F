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
C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi
#endif

      implicit          none
      integer           na, nomax
      integer           lastkb(0:na), lasto(0:na), listh(nomax,*),
     .                  nua, numh(*)
      double precision  scell(3,3), rckb(*), rco(*), xa(3,na),
     .                  xijo(3,nomax,*)
      logical           negl
      external          neighb, timer

C Internal variables -----------------
C maxna  = maximum number of neighbour atoms of any atom
C maxnkb = maximum number of neighbour KB projectors of any orbital
      integer, save ::
     .  maxna, maxnkb

      integer
     .  ia, ikb, inkb, io, isel, iio,
     .  j, ja, jna, jo,
     .  ka, kna, ko, maxnain,
     .  nkb, nna, nnkb, no, Node, Nodes

      integer, dimension(:), allocatable ::
     .  jana, jnao, knakb, ibuffer

#ifdef MPI
      integer
     .  MPIerror
#endif

      double precision
     .  rci, rcj, rck, rij, rik, rjk,
     .  rmax, rmaxkb, rmaxo

      double precision, dimension(:), allocatable ::
     .  r2ij, rcnkb, dpbuffer

      double precision, dimension(:,:), allocatable ::
     .  xija

      logical
     .  conect, warn1, warn2, overflow
     
      external
     .  memory

      save warn1, warn2
      data warn1, warn2 /2*.false./
      data maxna  / 1000 /
      data maxnkb / 2000 /
C -------------------------------------

C     Start time counter
*     call timer( 'xijorb', 1 )

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

C Set local variables
      nkb = lastkb(na)
      maxnkb = max(nkb,maxnkb)
      no = lasto(na)

C Find maximum range of basis orbitals and KB projectors
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

C Allocate local arrays that depend on input values
      allocate(jnao(no))
      call memory('A','I',no,'xijorb')

C Allocate local arrays that depend on parameters
      allocate(knakb(maxnkb))
      call memory('A','I',maxnkb,'xijorb')
      allocate(rcnkb(maxnkb))
      call memory('A','D',maxnkb,'xijorb')
  100 if (.not.allocated(jana)) then
        allocate(jana(maxna))
        call memory('A','I',maxna,'xijorb')
      endif
      if (.not.allocated(r2ij)) then
        allocate(r2ij(maxna))
        call memory('A','D',maxna,'xijorb')
      endif
      if (.not.allocated(xija)) then
        allocate(xija(3,maxna))
        call memory('A','D',3*maxna,'xijorb')
      endif

C Initialize neighb subroutine
      isel = 0
      ia = 0
      nna = maxna
      call neighb( scell, rmax, na, xa, ia, isel,
     .             nna, jana, xija, r2ij )
      overflow = (nna.gt.maxna)
      if (overflow) then
        call memory('D','I',size(jana),'xijorb')
        deallocate(jana)
        call memory('D','D',size(r2ij),'xijorb')
        deallocate(r2ij)
        call memory('D','D',size(xija),'xijorb')
        deallocate(xija)
        maxna = nna + nint(0.1*nna)
        goto 100
      endif

C Initialize vector jnao only once
      do io = 1,no
        jnao(io) = 0
      enddo

C Loop on atoms (only within unit cell)
      overflow = .false.
      maxnain = maxna
      do ia = 1,nua

C Find neighbour atoms within maximum range
        nna = maxnain
        call neighb( scell, rmax, na, xa, ia, isel,
     .               nna, jana, xija, r2ij )
        if (.not.overflow) overflow = (nna.gt.maxna)
        if (overflow) then
          maxna = max(maxna,nna)
        endif

C Don't do the actual work if the neighbour arrays are too small
        if (.not.overflow) then

C Loop on orbitals of atom ia
          do io = lasto(ia-1)+1,lasto(ia)

C Is this orbital on this node?
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then

              rci = rco(io)

C Find overlaping KB projectors
              if (.not.negl) then
                nnkb = 0
                do kna = 1,nna
                  ka = jana(kna)
                  rik = sqrt( r2ij(kna) )
                  do ko = lastkb(ka-1)+1,lastkb(ka)
                    rck = rckb(ko)
                    if (rci+rck .gt. rik) then
C Check maxnkb - if too small then increase array sizes
                      if (nnkb.eq.maxnkb) then
                        allocate(ibuffer(maxnkb))
                        call memory('A','I',maxnkb,'xijorb')
                        ibuffer(1:maxnkb) = knakb(1:maxnkb)
                        call memory('D','I',size(knakb),'xijorb')
                        deallocate(knakb)
                        allocate(knakb(maxnkb+100))
                        call memory('A','I',maxnkb+100,'xijorb')
                        knakb(1:maxnkb) = ibuffer(1:maxnkb)
                        call memory('D','I',size(ibuffer),'xijorb')
                        deallocate(ibuffer)
                        allocate(dpbuffer(maxnkb))
                        call memory('A','D',maxnkb,'xijorb')
                        dpbuffer(1:maxnkb) = rcnkb(1:maxnkb)
                        call memory('D','D',size(rcnkb),'xijorb')
                        deallocate(rcnkb)
                        allocate(rcnkb(maxnkb+100))
                        call memory('A','D',maxnkb+100,'xijorb')
                        rcnkb(1:maxnkb) = dpbuffer(1:maxnkb)
                        call memory('D','D',size(dpbuffer),'xijorb')
                        deallocate(dpbuffer)
                        maxnkb = maxnkb + 100
                      endif
                      nnkb = nnkb + 1
                      knakb(nnkb) = kna
                      rcnkb(nnkb) = rck
                    endif
                  enddo
                enddo
              endif

C Find orbitals connected by direct overlap or
C through a KB projector
              do jna = 1,nna
                ja = jana(jna)
                rij = sqrt( r2ij(jna) )
                do jo = lasto(ja-1)+1,lasto(ja)
                  rcj = rco(jo)

C Find if orbitals io and jo are connected
                  conect = .false.
C Find if there is direct overlap
                  if (rci+rcj .gt. rij) then
                    conect = .true.
                  elseif (.not.negl) then
C Find if jo overlaps with a KB projector
                    do inkb = 1,nnkb
                      rck = rcnkb(inkb)
                      kna = knakb(inkb)
                      rjk = sqrt( (xija(1,kna)-xija(1,jna))**2 +
     .                            (xija(2,kna)-xija(2,jna))**2 +
     .                            (xija(3,kna)-xija(3,jna))**2 )
                      if (rcj+rck .gt. rjk) then
                        conect = .true.
                        goto 50
                      endif
                    enddo
   50               continue
                  endif

C Add to list of connected orbitals
                  if (conect) then
                    if (jnao(jo) .eq. 0) then
                      jnao(jo) = jna
                    else
                      if (.not.warn1) then
                        if (Node.eq.0) then
                          write(6,'(/,a,2i6,a,/)')
     .                      'xijorb: WARNING: orbital pair ',io,jo,
     .                      ' is multiply connected'
                        endif
                        warn1 = .true.
                      endif
                      kna = jnao(jo)
                      if (r2ij(jna) .lt. r2ij(kna)) jnao(jo) = jna
                    endif
                  endif
                enddo
              enddo

C Copy interatomic vectors into xijo
              do j = 1,numh(io)
                jo = listh(j,io)
                jna = jnao(jo)
                if (jna .eq. 0) then
                  if (.not.warn2) then
                    if (Node.eq.0) then
                      write(6,'(/,a,2i6,a,/)')
     .                  'xijorb: WARNING: orbital pair ', io, jo,
     .                  ' is in listh but not really connected'
                    endif
                    warn2 = .true.
                  endif
                  xijo(1,j,iio) = 0.d0
                  xijo(2,j,iio) = 0.d0
                  xijo(3,j,iio) = 0.d0
                else
                  xijo(1,j,iio) = xija(1,jna)
                  xijo(2,j,iio) = xija(2,jna)
                  xijo(3,j,iio) = xija(3,jna)
                endif
              enddo

C Restore jnao array for next orbital io
              do jna = 1,nna
                ja = jana(jna)
                do jo = lasto(ja-1)+1,lasto(ja)
                  jnao(jo) = 0
                enddo
              enddo

            endif

          enddo

        endif

      enddo

C If maxna dimension was no good then reset arrays and start again
      if (overflow) then
        call memory('D','I',size(jana),'xijorb')
        deallocate(jana)
        call memory('D','D',size(r2ij),'xijorb')
        deallocate(r2ij)
        call memory('D','D',size(xija),'xijorb')
        deallocate(xija)
        maxna = maxna + nint(0.1*maxna)
        goto 100
      endif

C Deallocate local arrays
      call memory('D','I',size(jnao),'xijorb')
      deallocate(jnao)
      call memory('D','I',size(jana),'xijorb')
      deallocate(jana)
      call memory('D','I',size(knakb),'xijorb')
      deallocate(knakb)
      call memory('D','D',size(rcnkb),'xijorb')
      deallocate(rcnkb)
      call memory('D','D',size(r2ij),'xijorb')
      deallocate(r2ij)
      call memory('D','D',size(xija),'xijorb')
      deallocate(xija)

*     call timer( 'xijorb', 2 )
      end
