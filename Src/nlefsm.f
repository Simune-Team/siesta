      subroutine nlefsm( scell, nua, na, isa, xa, indxua, maxna,
     .                   maxnh, maxnd, lasto, lastkb, iphorb, 
     .                   iphKB, numd, listdptr, listd, numh, 
     .                   listhptr, listh, nspin, Dscf, Enl, 
     .                   fa, stress, H )
C *********************************************************************
C Calculates non-local (NL) pseudopotential contribution to total 
C energy, atomic forces, stress and hamiltonian matrix elements.
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, June 1997.
C **************************** INPUT **********************************
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer isa(na)          : Species index of each atom
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C integer indxua(na)       : Index of equivalent atom in unit cell
C integer maxna            : Maximum number of atoms for neighb
C integer maxnh            : First dimension of H and listh
C integer maxnd            : Maximum number of elements of the
C                            density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer lastkb(0:na)     : Position of last KB projector of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
C integer iphKB(nokb)      : Index of each KB projector in its atom,
C                            where nokb=lastkb(na)
C integer numd(nuo)        : Number of nonzero elements of each row of the
C                            density matrix
C integer listdptr(nuo)    : Pointer to the start of each row (-1) of the
C                            density matrix
C integer listd(maxnd)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer numh(nuo)        : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listhptr(nuo)    : Pointer to the start of each row (-1) of the
C                            hamiltonian matrix
C integer listh(maxnh)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Number of spin components
C real*8  Dscf(maxnd,nspin): Density matrix
C ******************* INPUT and OUTPUT *********************************
C real*8 fa(3,na)          : NL forces (added to input fa)
C real*8 stress(3,3)       : NL stress (added to input stress)
C real*8 H(maxnh,nspin)    : NL Hamiltonian (added to input H)
C **************************** OUTPUT *********************************
C real*8 Enl               : NL energy
C *********************************************************************
C
C  Modules
C
      use precision
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut, epskb

      implicit none

      integer
     .  maxna, maxnh, na, maxnd, nspin, nua

      integer
     .  indxua(na), iphKB(*), iphorb(*), isa(na),  
     .  lasto(0:na), lastkb(0:na), listd(maxnd), listh(maxnh),
     .  numd(*), numh(*), listdptr(*), listhptr(*)

      real(dp)
     .  scell(3,3), Dscf(maxnd,nspin), Enl,
     .  fa(3,nua), H(maxnh,nspin),
     .  stress(3,3), xa(3,na), volcel

      external
     .  neighb, timer, volcel, memory

C Internal variables ................................................
C maxno  = maximum number of basis orbitals overlapping a KB projector
      integer, save ::
     .  maxno
  
      integer
     .  ia, ikb, ina, ind, ino,
     .  io, iio, ioa, is, ispin, ix, 
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nkb, nna, nno, no, nuo, nuotot, maxnain, maxkba

      integer, dimension(:), allocatable, save ::
     .  iano, iana, iono, ibuffer

      real(dp)
     .  Cijk, epsk, fik, rki, rmax, rmaxkb, rmaxo, 
     .  Sik, Sjk, volume

      real(dp), dimension(:), allocatable, save ::
     .  Di, Vi, r2ki

      real(dp), dimension(:,:), allocatable, save ::
     .  Ski, xki, xno, dpbuffer2

      real(dp), dimension(:,:,:), allocatable, save ::
     .  grSki, dpbuffer3

      logical
     .  within, overflow
     
      logical, dimension(:), allocatable, save ::
     .  listed, listedall

      data maxno  / 2000 /
C ......................

C Start time counter
      call timer( 'nlefsm', 1 )

C Find unit cell volume
      volume = volcel( scell ) * nua / na

C Find maximum range
      rmaxo = 0.0d0
      rmaxkb = 0.0d0
      do ia = 1,na
        is = isa(ia)
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rmaxkb = max( rmaxkb, rcut(is,ioa) )
        enddo
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rmaxo = max( rmaxo, rcut(is,ioa) )
        enddo
      enddo
      rmax = rmaxo + rmaxkb

C Initialize arrays Di and Vi only once
      no = lasto(na)
      nuotot = lasto(nua)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

C Allocate local memory
      allocate(Di(no))
      call memory('A','D',no,'nlefsm')
      allocate(Vi(no))
      call memory('A','D',no,'nlefsm')
      allocate(listed(no))
      call memory('A','L',no,'nlefsm')
      allocate(listedall(no))
      call memory('A','L',no,'nlefsm')

      Enl = 0.0d0
      do jo = 1,no
        Di(jo) = 0.0d0
        Vi(jo) = 0.0d0
        listed(jo) = .false.
        listedall(jo) = .false.
      enddo

C Make list of all orbitals needed for this node
      do io = 1,nuo
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        listedall(iio) = .true.
        do j = 1,numh(io)
          jo = listh(listhptr(io)+j)
          listedall(jo) = .true.
        enddo
      enddo

C Find maximum number of KB projectors of one atom = maxkba
      maxkba = 0
      do ka = 1,na
        nkb = lastkb(ka) - lastkb(ka-1)
        maxkba = max(maxkba,nkb)
      enddo

C Allocate local arrays that depend on saved parameters
      allocate(iano(maxno))
      call memory('A','I',maxno,'nlefsm')
      allocate(iono(maxno))
      call memory('A','I',maxno,'nlefsm')
      allocate(xno(3,maxno))
      call memory('A','D',3*maxno,'nlefsm')
      allocate(Ski(maxkba,maxno))
      call memory('A','D',maxkba*maxno,'nlefsm')
      allocate(grSki(3,maxkba,maxno))
      call memory('A','D',3*maxkba*maxno,'nlefsm')
  100 if (.not.allocated(iana)) then
        allocate(iana(maxna))
        call memory('A','I',maxna,'nlefsm')
      endif
      if (.not.allocated(r2ki)) then
        allocate(r2ki(maxna))
        call memory('A','D',maxna,'nlefsm')
      endif
      if (.not.allocated(xki)) then
        allocate(xki(3,maxna))
        call memory('A','D',3*maxna,'nlefsm')
      endif

C Initialize neighb subroutine
      nna = maxna
      call neighb( scell, rmax, na, xa, 0, 0,
     .             nna, iana, xki, r2ki )
      overflow = (nna.gt.maxna)
      if (overflow) then
        call memory('D','I',size(iana),'nlefsm')
        deallocate(iana)
        call memory('D','D',size(r2ki),'nlefsm')
        deallocate(r2ki)
        call memory('D','D',size(xki),'nlefsm')
        deallocate(xki)
        maxna = nna + nint(0.1*nna)
        goto 100
      endif

C Loop on atoms with KB projectors      
      overflow = .false.
      maxnain = maxna
      do ka = 1,na
        kua = indxua(ka)
        ks = isa(ka)
        nkb = lastkb(ka) - lastkb(ka-1)

C Find neighbour atoms
        nna = maxnain
        call neighb( scell, rmax, na, xa, ka, 0,
     .               nna, iana, xki, r2ki )
        if (.not.overflow) overflow = (nna.gt.maxna)
        if (overflow) then
          maxna = max(maxna,nna)
        endif

C Don't do the actual work if the neighbour arrays are too small
        if (.not.overflow) then

C Find neighbour orbitals
          nno = 0
          do ina = 1,nna
            ia = iana(ina)
            is = isa(ia)
            rki = sqrt(r2ki(ina))
            do io = lasto(ia-1)+1,lasto(ia)

C Only calculate if needed locally
              if (listedall(io)) then

                ioa = iphorb(io)
          
C Find if orbital is within range
                within = .false.
                do ko = lastkb(ka-1)+1,lastkb(ka)
                  koa = iphKB(ko)
                  if ( rki .lt. rcut(is,ioa)+rcut(ks,koa) ) 
     .              within = .true.
                enddo

C Find overlap between neighbour orbitals and KB projectors
                if (within) then
C Check maxno - if too small then increase array sizes
                  if (nno.eq.maxno) then
                    allocate(ibuffer(maxno))
                    call memory('A','I',maxno,'nlefsm')
                    ibuffer(1:maxno) = iano(1:maxno)
                    call memory('D','I',size(iano),'nlefsm')
                    deallocate(iano)
                    allocate(iano(maxno+100))
                    call memory('A','I',maxno+100,'nlefsm')
                    iano(1:maxno) = ibuffer(1:maxno)
                    call memory('D','I',size(ibuffer),'nlefsm')
                    ibuffer(1:maxno) = iono(1:maxno)
                    call memory('D','I',size(iono),'nlefsm')
                    deallocate(iono)
                    allocate(iono(maxno+100))
                    call memory('A','I',maxno+100,'nlefsm')
                    iono(1:maxno) = ibuffer(1:maxno)
                    deallocate(ibuffer)
                    allocate(dpbuffer2(3,maxno))
                    call memory('A','D',3*maxno,'nlefsm')
                    dpbuffer2(1:3,1:maxno) = xno(1:3,1:maxno)
                    call memory('D','D',size(xno),'nlefsm')
                    deallocate(xno)
                    allocate(xno(3,maxno+100))
                    call memory('A','D',3*(maxno+100),'nlefsm')
                    xno(1:3,1:maxno) = dpbuffer2(1:3,1:maxno)
                    call memory('D','D',size(dpbuffer2),'nlefsm')
                    deallocate(dpbuffer2)
                    allocate(dpbuffer2(maxkba,maxno))
                    call memory('A','D',maxkba*maxno,'nlefsm')
                    dpbuffer2(1:maxkba,1:maxno) = Ski(1:maxkba,1:maxno)
                    call memory('D','D',size(Ski),'nlefsm')
                    deallocate(Ski)
                    allocate(Ski(maxkba,maxno+100))
                    call memory('A','D',maxkba*(maxno+100),'nlefsm')
                    Ski(1:maxkba,1:maxno) = dpbuffer2(1:maxkba,1:maxno)
                    call memory('D','D',size(dpbuffer2),'nlefsm')
                    deallocate(dpbuffer2)
                    allocate(dpbuffer3(3,maxkba,maxno))
                    call memory('A','D',3*maxkba*maxno,'nlefsm')
                    dpbuffer3(1:3,1:maxkba,1:maxno) = grSki(1:3,
     .                1:maxkba,1:maxno)
                    call memory('D','D',size(grSki),'nlefsm')
                    deallocate(grSki)
                    allocate(grSki(3,maxkba,maxno+100))
                    call memory('A','D',3*maxkba*(maxno+100),'nlefsm')
                    grSki(1:3,1:maxkba,1:maxno) = dpbuffer3(1:3,
     .                1:maxkba,1:maxno)
                    call memory('D','D',size(dpbuffer3),'nlefsm')
                    deallocate(dpbuffer3)
                    maxno = maxno + 100
                  endif
                  nno = nno + 1
                  iono(nno) = io
                  iano(nno) = ia
                  do ix = 1,3
                    xno(ix,nno) = xki(ix,ina)
                  enddo
                  ikb = 0
                  do ko = lastkb(ka-1)+1,lastkb(ka)
                    ikb = ikb + 1
                    ioa = iphorb(io)
                    koa = iphKB(ko)
                    call matel( 'S', ks, is, koa, ioa, xki(1,ina),
     .                    Ski(ikb,nno), grSki(1,ikb,nno) )
                  enddo

                endif

              endif

            enddo

          enddo

C Loop on neighbour orbitals
          do ino = 1,nno
            io = iono(ino)
            ia = iano(ino)

            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then
C Valid orbital

              if (ia .le. nua) then

C Scatter density matrix row of orbital io
                do j = 1,numd(iio)
                  ind = listdptr(iio)+j
                  jo = listd(ind)
                  do ispin = 1,nspin
                    Di(jo) = Di(jo) + Dscf(ind,ispin)
                  enddo
                enddo
          
C Scatter filter of desired matrix elements
                do j = 1,numh(iio)
                  jo = listh(listhptr(iio)+j)
                  listed(jo) = .true.
                enddo

C Find matrix elements with other neighbour orbitals
                do jno = 1,nno
                  jo = iono(jno)
                  if (listed(jo)) then

C Loop on KB projectors
                    ikb = 0
                    do ko = lastkb(ka-1)+1,lastkb(ka)
                      ikb = ikb + 1
                      koa = iphKB(ko)
                      epsk = epskb(ks,koa)
                      Sik = Ski(ikb,ino)
                      Sjk = Ski(ikb,jno)
                      Vi(jo) = Vi(jo) + epsk * Sik * Sjk
                      Cijk = Di(jo) * epsk
                      Enl = Enl + Cijk * Sik * Sjk
                      do ix = 1,3
                        fik = 2.d0 * Cijk * Sjk * grSki(ix,ikb,ino)
                        fa(ix,ia)  = fa(ix,ia)  - fik
                        fa(ix,kua) = fa(ix,kua) + fik
                        do jx = 1,3
                          stress(jx,ix) = stress(jx,ix) +
     .                                xno(jx,ino) * fik / volume
                        enddo
                      enddo
                    enddo

                  endif
                enddo

C Pick up contributions to H and restore Di and Vi
                do j = 1,numh(iio)
                  ind = listhptr(iio)+j
                  jo = listh(ind)
                  do ispin = 1,nspin
                    H(ind,ispin) = H(ind,ispin) + Vi(jo)
                  enddo
                  Vi(jo) = 0.0d0
                  listed(jo) = .false.
                enddo
                do j = 1,numd(iio)
                  jo = listd(listdptr(iio)+j)
                  Di(jo) = 0.0d0
                enddo

              endif

            endif

          enddo

        endif

      enddo

      if (overflow) then
        call memory('D','I',size(iana),'nlefsm')
        deallocate(iana)
        call memory('D','D',size(r2ki),'nlefsm')
        deallocate(r2ki)
        call memory('D','D',size(xki),'nlefsm')
        deallocate(xki)
        maxna = nna + nint(0.1*nna)
        goto 100
      endif

C Deallocate local memory
      call memory('D','D',size(Di),'nlefsm')
      deallocate(Di)
      call memory('D','D',size(Vi),'nlefsm')
      deallocate(Vi)
      call memory('D','L',size(listed),'nlefsm')
      deallocate(listed)
      call memory('D','L',size(listedall),'nlefsm')
      deallocate(listedall)
      call memory('D','I',size(iano),'nlefsm')
      deallocate(iano)
      call memory('D','I',size(iono),'nlefsm')
      deallocate(iono)
      call memory('D','D',size(xno),'nlefsm')
      deallocate(xno)
      call memory('D','D',size(Ski),'nlefsm')
      deallocate(Ski)
      call memory('D','D',size(grSki),'nlefsm')
      deallocate(grSki)
      call memory('D','I',size(iana),'nlefsm')
      deallocate(iana)
      call memory('D','D',size(r2ki),'nlefsm')
      deallocate(r2ki)
      call memory('D','D',size(xki),'nlefsm')
      deallocate(xki)

      call timer( 'nlefsm', 2 )

      end
